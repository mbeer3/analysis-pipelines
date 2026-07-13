# --------------------------------------------------

library(openxlsx)
library(idr)
library(GenomicRanges)
library(Gviz)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)
library(ArchR)

library(BSgenome.Hsapiens.UCSC.hg19)
hg19.info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)


# --------------------------------------------------
SAMPLE <- sapply(c("SAMPLE1", "SAMPLE2"), function(sheet){
  lapply(c("Chrom1", "Chrom3", "Chrom8", "Chrom18"), function(chr){
    df <- read.xlsx(paste0("../../sequencing-results/", chr, ".xlsx"), sheet = sheet,
                    check.names = F, rowNames = F, colNames = T)
    df <- strsplit(df$HG38, split = ":") %>%
      lapply(function(x){data.frame(chr = x[1], HG38.pos = as.integer(x[2]))}) %>%
      do.call(what = rbind) %>% cbind(df[, -2])
    colnames(df)[3:4] <- c("HG19.pos", "ATAC-peak.HG38")
    df$gene <- sapply(strsplit(df$Name, split = "/"), function(x){
      if(length(x) > 1) return("multiple") else return(x)
    }) %>% gsub(pattern = "(((\\+)|(-))\\d+)+$", replacement = "")
    return(df)
  }) %>% do.call(what = rbind)
}, simplify = F)
all(SAMPLE$SAMPLE1[, -(7:8)] == SAMPLE$SAMPLE2[, -(7:8)]) # TRUE
SAMPLE <- cbind(SAMPLE$SAMPLE1[, -(7:8)], SAMPLE$SAMPLE1[, 7:8], SAMPLE$SAMPLE2[, 7:8]) %>%
  subset(gene %in% c("MIXL1", "EOMES", "SOX17", "GATA6"))
colnames(SAMPLE)[8:11] <- c("rep1.NEG", "rep1.POS", "rep2.NEG", "rep2.POS")

all <- openxlsx::read.xlsx(
  xlsxFile = "TierII-off_target-filter-gRNA_list-Renhe.xlsx",
  sheet = "Original", colNames = T, rowNames = F, check.names = F
)
all(SAMPLE$Name %in% all$Name) # TRUE

filtered <- openxlsx::read.xlsx(
  xlsxFile = "TierII-off_target-filter-gRNA_list-Renhe.xlsx",
  sheet = "off-target filtered", colNames = T, rowNames = F, check.names = F
)
SAMPLE <- subset(SAMPLE, Name %in% filtered$Name) # 3679/4289


# --------------------------------------------------
get.size.factos <- function(counts){
  counts <- as.matrix(counts)
  geom.mean <- sqrt(as.numeric(counts[, 1]) * as.numeric(counts[, 2]))
  return(data.frame(
    POS = median(counts[, 1] / geom.mean),
    NEG = median(counts[, 2] / geom.mean)
  ))
}
normalize.counts <- function(counts, size.factos){
  data.frame(
    rep1.NEG = round(counts$rep1.NEG / size.factos$rep1.NEG),
    rep1.POS = round(counts$rep1.POS / size.factos$rep1.POS),
    rep2.NEG = round(counts$rep2.NEG / size.factos$rep2.NEG),
    rep2.POS = round(counts$rep2.POS / size.factos$rep2.POS)
  )
}

pseudocount <- 1000
zs.thr <- 1
idr.thr <- 0.05
max.dist <- 100
merged.pos <- c()
SAMPLE <- lapply(c("MIXL1", "EOMES", "SOX17", "GATA6"), function(g){
  SAMPLE <- subset(SAMPLE, gene == g)
  SAMPLE[, c("rep1.POS", "rep1.NEG", "rep2.POS", "rep2.NEG")] <- pseudocount + SAMPLE[, c("rep1.POS", "rep1.NEG", "rep2.POS", "rep2.NEG")]
  size.factos <- cbind(
    rep1 = get.size.factos(counts = SAMPLE[, c("rep1.POS", "rep1.NEG")]),
    rep2 = get.size.factos(counts = SAMPLE[, c("rep2.POS", "rep2.NEG")])
  )
  SAMPLE <- cbind(SAMPLE, norm = normalize.counts(SAMPLE, size.factos))
  SAMPLE$LFC1 <- log2((SAMPLE$norm.rep1.POS + 1) / (SAMPLE$norm.rep1.NEG + 1))
  SAMPLE$LFC2 <- log2((SAMPLE$norm.rep2.POS + 1) / (SAMPLE$norm.rep2.NEG + 1))
  SAMPLE$z.score1 <- scale(SAMPLE$LFC1)[, 1]
  SAMPLE$z.score2 <- scale(SAMPLE$LFC2)[, 1]
  SAMPLE$IDR <- est.IDR(
    x = as.matrix(abs(SAMPLE[, c("z.score1", "z.score2")])),
    mu = 0.9, sigma = 1, rho = 0.7, p = 0.4
  )$IDR

  SAMPLE$sig <- SAMPLE$IDR <= idr.thr & SAMPLE$z.score1 >= zs.thr & SAMPLE$z.score2 >= zs.thr
  SAMPLE$near.sig <- sapply(1:nrow(SAMPLE), function(i){
    any(abs(SAMPLE$HG19.pos[-i] - SAMPLE$HG19.pos[i]) <= max.dist & SAMPLE$sig[-i])
  })
  SAMPLE$hit <- ifelse(SAMPLE$sig & SAMPLE$near.sig, "Yes", "No")
  write.csv(SAMPLE, file = paste0("SAMPLE/", g, ".csv"), row.names = F)

  d <- subset(SAMPLE, hit == "Yes")
  d$group <- NA; G <- 0
  merged.pos <<- sapply(1:nrow(d), function(i){
    if(is.na(d$group[i])){
      ri <- abs(d$HG19.pos - d$HG19.pos[i]) <= max.dist
      d$group[ri] <<- G; G <<- G + 1
      gr <- GRanges(
        seqinfo = hg19.info, seqnames = unique(d$chr[ri]),
        ranges = IRanges(start = min(d$HG19.pos[ri]), end = max(d$HG19.pos[ri]))
      )
      return(gr)
    } else NULL
  }) %>% do.call(what = c) %>% c(merged.pos)

  return(SAMPLE)
}) %>% do.call(what = rbind)
write.csv(SAMPLE, file = "SAMPLE/all-results.csv", row.names = F)
saveRDS(SAMPLE, file = "SAMPLE/all-results.rds")

merged.pos <- reduce(merged.pos + 75)
saveRDS(merged.pos, "SAMPLE/hits_non-overlapping_merged-pos_ext150bp.rds")
data.frame(merged.pos)[, 1:3] %>%
  cbind(ID = paste0("p", 1:length(merged.pos)), NA, strand = "*") %>%
  write.table(file = "SAMPLE/hits_non-overlapping_merged-pos_ext150bp.bed",
              col.names = F, row.names = F, quote = F, sep = "\t")

ggsave(
  filename = "SAMPLE/summary.pdf", width = 12, height = 10,
  plot = ggplot(mapping = aes(x = z.score1, y = z.score2, size = -log10(IDR))) +
    geom_point(data = subset(SAMPLE, !sig), color = "grey") +
    geom_point(data = subset(SAMPLE, sig), aes(color = gene), alpha = 0.8) +
    labs(x = "rep1 z-score", y = "rep2 z-score", color = "", shape = "",
         title = paste("IDR{|z-score1|,|z-score2|} <=", idr.thr,
                       "& z-score1 >=", zs.thr, "& z-score2 >=", zs.thr)) +
    scale_x_continuous(limits = c(-3, 8)) +
    scale_y_continuous(limits = c(-3, 8)) +
    scale_color_manual(values = c(
      EOMES = "#d93731", GATA6 = "#e16d38",
      MIXL1 = "#5d338b", SOX17 = "#40904e"
    )) + theme_bw() + labs(x = "Z-score (Rep. 1)", y = "Z-score (Rep. 2)")
)

SAMPLE.gr <-   GRanges(
  seqinfo = hg19.info, seqnames = SAMPLE$chr,
  ranges = IRanges(start = SAMPLE$HG19.pos, width = 1, names = SAMPLE$Name),
  score = (SAMPLE$LFC1 + SAMPLE$LFC2) / 2
  # score = -log10(SAMPLE$IDR)
)
export.bw(SAMPLE.gr, con = "SAMPLE/all_average_LFC.bw")

SAMPLE.hits <- subset(SAMPLE, hit == "Yes")
SAMPLE.hits.gr <-   GRanges(
  seqinfo = hg19.info, seqnames = SAMPLE.hits$chr,
  ranges = IRanges(start = SAMPLE.hits$HG19.pos, width = 1, names = SAMPLE.hits$Name),
  score = (SAMPLE.hits$LFC1 + SAMPLE.hits$LFC2) / 2
  # score = -log10(SAMPLE.hits$IDR)
)
export.bw(SAMPLE.hits.gr, con = "SAMPLE/hits_average_LFC.bw")


# --------------------------------------------------
long.peaks <- "../../XX/ATACpeaks.bed" %>%
  read.table(header = F, sep = "\t", col.names = c("chr", "start", "end")) %>%
  makeGRangesFromDataFrame(seqinfo = hg19.info, starts.in.df.are.0based = T) %>%
  subsetByOverlaps(ranges = c(SAMPLE.gr, merged.pos), type = "any") # 160

neg.regions <- lapply(1:length(long.peaks), function(pi){
  peak <- long.peaks[pi]
  ov.hits <- subsetByOverlaps(merged.pos, ranges = peak, type = "any")
  L <- length(ov.hits)
  if(L > 0){
    ov.hits <- ov.hits[order(start(ov.hits), decreasing = F)]
    if(start(ov.hits)[1] < start(peak)){
      warning("hit before peak ", peak)
      start(peak) <- start(ov.hits)[1] - 2
    }
    if(end(peak) < end(ov.hits)[L]){
      warning("hit after peak ", peak)
      end(peak) <- end(ov.hits)[L] + 2
    }
    lapply(1:(L + 1), function(i){
      if(i == 1){ S <- start(peak)
      } else S <- end(ov.hits)[i - 1] + 2
      if(i == (L + 1)){ E <- end(peak)
      } else E <- start(ov.hits)[i] - 2
      if(S <= E){
        GRanges(
          seqinfo = hg19.info, seqnames = seqnames(peak),
          ranges = IRanges(start = S, end = E)
        )
      } else NULL
    }) %>% do.call(what = c)
  } else return(peak)
}) %>% do.call(what = c)
neg.regions <- neg.regions[width(neg.regions) >= 150]
saveRDS(neg.regions, file = "SAMPLE/hits_non-overlapping_neg_min150bp.rds")
data.frame(neg.regions)[, 1:3] %>%
  cbind(ID = paste0("n", 1:length(neg.regions)), NA, strand = "*") %>%
  write.table(file = "SAMPLE/hits_non-overlapping_neg_min150bp.bed",
              col.names = F, row.names = F, quote = F, sep = "\t")


# --------------------------------------------------
comp.regions <- read.table(
  file = "competent-regions.bed", header = F, sep = "\t",
  col.names = c("chr", "start", "end", "name", "", "")
)[, 1:4] %>% makeGRangesFromDataFrame(seqinfo = hg19.info, starts.in.df.are.0based = F)

incomp.regions <- read.table(
  file = "incompetent-regions.bed", header = F, sep = "\t",
  col.names = c("chr", "start", "end", "name", "", "")
)[, 1:4] %>% makeGRangesFromDataFrame(seqinfo = hg19.info, starts.in.df.are.0based = F)

TADs <- read.csv("../TADs-to-filter.csv", header = T, row.names = 1)$region %>%
  strsplit(split = ":|-") %>%
  sapply(function(x){ c(chr = x[1], start = x[2], end = x[3]) }) %>%
  t() %>% data.frame() %>%
  makeGRangesFromDataFrame(seqinfo = hg19.info, starts.in.df.are.0based = F)
saveRDS(TADs, file = "../TADs-to-filter.rds")
incomp.regions <- subsetByOverlaps(incomp.regions, ranges = TADs, type = "within") # 62

data.frame(incomp.regions)[, 1:3] %>%
  cbind(ID = paste0("n", 1:length(incomp.regions)), NA, strand = "*") %>%
  write.table(file = "incompetent-regions-within-TADs.bed",
              col.names = F, row.names = F, quote = F, sep = "\t")

