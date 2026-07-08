# -------------------------------------------------------------------------
library(data.table)
library(metaseqR2)
library(idr)
library(plotly)
library(ggplot2)
library(gridExtra)
library(rtracklayer)

read.sam.table <- function(fn, fill){
  sam.table <- suppressMessages(fread(paste("cut -f1-11", fn), skip = 3, header = F, fill = fill)[, c(2:4, 6, 10)])
  colnames(sam.table) <- c("FLAG", "RNAME", "POS", "CIGAR", "SEQ")
  return(sam.table)
}

get.genomic.ranges <- function(sam.table, b, names){
  return(data.frame(chr = "chr13", pos = b + sam.table$POS - 1, row.names = names, stringsAsFactors = F))
}

get.sgRNA.table <- function(scr.results, index, pdx1.gr){
  scr.results <- data.frame(scr.results, stringsAsFactors = F)
  st <- cbind(
    rbind(
      pdx1.gr$intron.exon, pdx1.gr$upstream, pdx1.gr$downstream,
      data.frame(chr = gsub("^SafeHarbor", "", scr.results$sgRNA[index$safe.harbor]), pos = NA),
      data.frame(chr = rep(NA, sum(index$non.targeting)), pos = NA)
    ),
    rbind(
      scr.results[index$intron | index$exon,], scr.results[index$upstream,], scr.results[index$downstream,],
      scr.results[index$safe.harbor,], scr.results[index$non.targeting,]
    )
  )
  sh.idx <- grepl("SafeHarbor", st$sgRNA)
  rn <- st$sgRNA
  rn[sh.idx] <- paste0(rn[sh.idx], "-", 1:sum(sh.idx))
  rownames(st) <- rn
  st <- st[order(st$pos, rownames(st)),]
  return(st)
}

run.mageck <- function(sgRNAs, results.path, paired = F, gene.test.fdr.threshold = 0.25, gene.lfc.method = "median"){
  dir.create(paste0(results.path, "/mageck/"))
  cntID_file <- paste0(results.path, "/mageck/control_IDs.txt")
  count_file <- paste0(results.path, "/mageck/sgrna_count.txt")
  
  cnt <- grepl("NonTargetingControlGuide", sgRNAs$sgRNA)
  shr <- grepl("SafeHarbor", sgRNAs$sgRNA)
  write(sgRNAs$sgRNA[cnt], file = cntID_file, ncolumns = 1)
  write.table(rbind(
    data.frame(sgRNA = sgRNAs$sgRNA, gene = paste(sgRNAs$chr, sgRNAs$pos, sep = "_"), sgRNAs[, 5:8])[!cnt & !shr,],
    data.frame(sgRNA = sgRNAs$sgRNA, gene = "None", sgRNAs[, 5:8])[cnt,]
  ), file = count_file, row.names = F, col.names = T, quote = F, sep = "\t")
  
  cmd <- paste("mageck test -k", count_file, "--control-sgrna", cntID_file, "-n", paste0(results.path, "/mageck/MAGeCK"), "--pdf-report",
               "-t rep1.GFPplus,rep2.GFPplus -c rep1.GFPminus,rep2.GFPminus")
  if(paired) cmd <- paste(cmd, "--paired")
  if(gene.test.fdr.threshold != 0.25) cmd <- paste(cmd, "--gene-test-fdr-threshold", gene.test.fdr.threshold)
  if(gene.lfc.method != "median") cmd <- paste(cmd, "--gene-lfc-method", gene.lfc.method)
  cat(cmd, "\n")
}

size.factos <- function(counts){
  geom.mean <- sqrt(as.numeric(counts[, 1]) * as.numeric(counts[, 2]))
  return(data.frame(
    GFPplus = median(counts[, 1] / geom.mean),
    GFPminus = median(counts[, 2] / geom.mean)
  ))
}

normalize.counts <- function(sgRNAs, size.factos){
  sgRNAs$rep1.GFPplus <- round(sgRNAs$rep1.GFPplus / size.factos$rep1.GFPplus)
  sgRNAs$rep1.GFPminus <- round(sgRNAs$rep1.GFPminus / size.factos$rep1.GFPminus)
  sgRNAs$rep2.GFPplus <- round(sgRNAs$rep2.GFPplus / size.factos$rep2.GFPplus)
  sgRNAs$rep2.GFPminus <- round(sgRNAs$rep2.GFPminus / size.factos$rep2.GFPminus)
  return(sgRNAs)
}

Log2.Fold.Change <- function(sg){
  return(log2((sg[, 1] + 1) / (sg[, 2] + 1)))
}

average.LFC <- function(targets, N){
  if(N == 1){ w <- cbind(targets[, 1:2], window.span = 20, targets[, 9:10])
  } else{
    w <- data.frame(t(sapply(1:(nrow(targets) - N + 1), function(i){
      j <- i + N - 1
      return(c(pos = targets$pos[i], window.span = targets$pos[j] - targets$pos[i] + 20,
               rep1.LFC = mean(targets$rep1.LFC[i:j]), rep2.LFC = mean(targets$rep2.LFC[i:j])))
    })))
    w <- cbind(chr = "chr13", w)
  }
  names <- paste0("chr13:", w$pos + 1, "-", w$pos + w$window.span)
  names[duplicated(names)] <- paste0(names[duplicated(names)], "-2")
  rownames(w) <- names
  return(w)
}

window.scores <- function(targets, controls, N){
  rep1.mean <- mean(controls$rep1.LFC); rep2.mean <- mean(controls$rep2.LFC)
  rep1.stdv <- sd(controls$rep1.LFC); rep2.stdv <- sd(controls$rep2.LFC)
  if(N == 1){ w <- cbind(targets[, 1:2], window.span = 20, targets[, 9:10])
  } else{
    w <- data.frame(t(sapply(1:(nrow(targets) - N + 1), function(i){
      j <- i + N - 1
      rep1.LFC <- mean(targets$rep1.LFC[i:j]); rep2.LFC <- mean(targets$rep2.LFC[i:j])
      return(c(pos = targets$pos[i], window.span = targets$pos[j] - targets$pos[i] + 20,
               rep1.LFC = rep1.LFC, rep1.score = (rep1.LFC - rep1.mean) / rep1.stdv,
               rep2.LFC = rep2.LFC, rep2.score = (rep2.LFC - rep2.mean) / rep2.stdv))
    })))
    w <- cbind(chr = "chr13", w)
  }
  w$LFC <- (w$rep1.LFC + w$rep2.LFC) / 2
  names <- paste0("chr13:", w$pos + 1, "-", w$pos + w$window.span)
  names[duplicated(names)] <- paste0(names[duplicated(names)], "-2")
  rownames(w) <- names
  return(w)
}

merge.regions <- function(regions, max.span, method = "dist"){
  merged <- data.frame(check.names = F, `start-position` = regions$pos,
                       `end-position` = regions$pos + regions$window.span,
                       LFC = sig.regions$LFC, depth = 1)
  n <- nrow(merged)
  while(n > 1){
    span <- merged$`end-position`[2:n] - merged$`start-position`[1:(n - 1)]
    dist <- merged$`start-position`[2:n] - merged$`end-position`[1:(n - 1)]
    if(min(span) > max.span) break()
    if(method == "span"){ mi <- order(span, dist)[1]
    } else if(method == "dist"){
      do <- order(dist, span)
      mi <- do[which(span[do] < max.span)[1]]
    }
    merged$`end-position`[mi] <- merged$`end-position`[mi + 1]
    merged[mi, c("LFC", "depth")] <- colSums(merged[mi:(mi + 1), c("LFC", "depth")])
    merged <- merged[-(mi + 1),]
    n <- nrow(merged)
  }
  merged$LFC <- merged$LFC / merged$depth
  names <- paste0("chr13:", merged$`start-position` + 1, "-", merged$`end-position`)
  names[duplicated(names)] <- paste0(names[duplicated(names)], "-2")
  rownames(merged) <- names
  return(merged)
}

point.scores <- function(df, positions){
  pps <- mclapply(positions, function(pos){
    ov <- subset(df, start <= pos & end >= pos)
    if(nrow(ov) == 0){ return(NULL)
    } else return(c(pos = pos, mean.score = mean(ov$score), min.IDR.score = ov$score[which.min(ov$IDR)]))
  }, mc.cores = detectCores())
  return(data.frame(t(data.frame(pps[!sapply(pps, is.null)], check.names = F)),
                    row.names = NULL, check.names = F))
}

plot.both.LFCs <- function(sgRNAs, mt){
  tg <- subset(sgRNAs, !is.na(pos))
  if(is.null(tg$rep1.LFC)){ LFC1 <- Log2.Fold.Change(tg[, c(5, 7)])
  } else LFC1 <- tg$rep1.LFC
  if(is.null(tg$rep2.LFC)){ LFC2 <- Log2.Fold.Change(tg[, c(6, 8)])
  } else LFC2 <- tg$rep2.LFC
  m <- list(size = 5, opacity = 0.5)
  plot_ly(type = "scatter", mode = "markers", data = tg) %>%
    add_trace(text = rownames(tg), x = ~pos, y = LFC1, name = "replicate 1", marker = m) %>%
    add_trace(text = rownames(tg), x = ~pos, y = LFC2, name = "replicate 2", marker = m) %>%
    layout(colorway = c("#f4776e", "#1dc1c6"), title = mt,
           yaxis = list(title = "LogFC", zeroline = F, showline = F),
           xaxis = list(title = "start position on chr13", zeroline = F, showline = F))
}

plot.gRNAs <- function(df, title = "", color = NULL, name = NULL){
  p <- plot_ly(type = "scatter", mode = "markers", data = df)
  if(!is.null(color)) p <- p %>% add_trace(showlegend = F, x = ~pos, y = ~LFC, text = rownames(df),
                                           marker = list(color = color, colorscale = "Reds",
                                                         colorbar = list(title = "", len = 0.7)))
  if(!is.null(name)){
    p <- p %>% add_trace(showlegend = T, x = ~pos, y = ~LFC, text = rownames(df),
                         name = name, legendgroup = name)
  }
  p %>% layout(title = title, yaxis = list(title = "LogFC"),
               xaxis = list(title = "start position on chr13"))
}

plot.progress <- function(df){
  N <- as.integer(rownames(df))
  
  add.trace <- function(gp, c = T){
    if(c){ gp <- gp +
      geom_line(color = "cyan4") +
      geom_point(color = "cyan4")
    } else gp <- gp + geom_line() + geom_point()
    return(gp +
             scale_x_continuous(breaks = c(1, seq(0, 50, 5)[-1])) +
             theme_bw() +
             labs(x = "N"))
  }
  
  br1 <- c(20, seq(100, 5000, 100))
  br2 <- seq(round(min(df$cor), 2), 1, 0.01)
  suppressWarnings(grid.arrange(
    ggplot(mapping = aes(x = rep(N, 3), y = c(df$med.sp, df$mean.sp, df$q80.sp),
                         color = c(rep("median", nrow(df)), rep("mean", nrow(df)),
                                   rep("80% quantile", nrow(df))))) %>%
      add.trace(F) +
      labs(y = "window span (base pairs)", color = "") +
      scale_y_continuous(breaks = br1, labels = format(br1, width = 4)) +
      theme(legend.position = c(0.1, 0.7)),
    ggplot(mapping = aes(x = N, y = df$cor)) %>%
      add.trace() +
      labs(y = "correlation") +
      scale_y_continuous(breaks = br2, labels = format(br2, width = 4))
  ))
}

plot.peak <- function(targets, controls, windows, N, rep){
  sw <- which.min(windows[, paste0("rep", rep, ".FDR")])
  w.idx <- sw:(sw + N - 1)
  stats <- paste0("rep", rep, c(".p-value", ".FDR"))
  clt <- "non-targetting control guides"
  tlt <- paste0("The most significant ", N, "-guide window in replicate ", rep)
  ltm <- 1:2; names(ltm) <- c(tlt, clt)
  st <- paste0(tlt, " (", gsub("-2$", "", rownames(windows)[sw]), "): ",
               paste(paste(stats, signif(windows[sw, stats], 4), sep = " = "), collapse = ", "))
  mt <- paste("Distribution of LFC in", nrow(controls), "control guides and",
              N, "guides from the most significant window.")
  ggplot() +
    geom_density(aes(x = targets$rep1.LFC[w.idx], color = "replicate 1", linetype = tlt)) +
    geom_density(aes(x = targets$rep2.LFC[w.idx], color = "replicate 2", linetype = tlt)) +
    geom_density(aes(x = controls$rep1.LFC, color = "replicate 1", linetype = clt)) +
    geom_density(aes(x = controls$rep2.LFC, color = "replicate 2", linetype = clt)) +
    labs(title = mt, subtitle = st, x = "LogFC", color = "", linetype = "") +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_linetype_manual(values = ltm)
}

plot.idr1 <- function(windows, controls){
  gp <- ggplot(windows, aes(color = IDR)) +
    theme_bw() +
    scale_color_gradient(low = "green", high = "red")
  grid.arrange(
    gp +
      geom_point(aes(x = abs(rep1.LFC), y = abs(rep2.LFC)), size = 1, show.legend = F) +
      labs(x = "|LFC| for rep1", y = "|LFC| for rep2"),
    gp +
      geom_point(aes(x = rank(-abs(rep1.score)), y = rank(-abs(rep2.score))), size = 1) +
      labs(x = "rank for rep1", y = "rank for rep2"),
    ncol = 2, widths = c(1, 1.2)
  )
}

plot.idr2 <- function(windows){
  br <- c(0, 0.00001, 0.05, 0.1, 0.3, 0.5, 0.8, 1)
  return(ggplot(windows, aes(x = pos, y = LFC, color = IDR)) +
           labs(x = "start position on chr13", y = "LogFC") +
           geom_point(size = 1, alpha = 0.8) +
           theme_bw() +
           scale_colour_gradientn(values = br, colours = viridis::viridis(n = length(br)),
                                  breaks = br, guide = guide_colorbar(barheight = 30)))
}

plot.score <- function(windows){
  col <- c("#173c68", "#173c68", "#1c5785", "#1c5785", "#ffffff", "#b32339", "#b32339", "#690822", "#690822")
  gp <- ggplot(windows, aes(x = pos, y = LFC)) +
    labs(x = "start position on chr13", y = "LogFC") +
    theme_bw() +
    scale_colour_gradient2(high = "red", low = "blue", mid = "#dbdbdb", midpoint = 0)
  grid.arrange(gp +
                 geom_point(aes(color = rep1.score), size = 1, alpha = 0.8) +
                 labs(color = "replicate 1\nz-score"),
               gp +
                 geom_point(aes(color = rep2.score), size = 1, alpha = 0.8) +
                 labs(color = "replicate 2\nz-score"),
               ncol = 1, top = "z-score of Log2-fold-change (using control guides as reference)")
}

plot.sig.regions <- function(sig.regions, span.max, idr.max, azs.min, merged.span.max, merged.line = 1){
  plot_ly() %>%
    add_segments(name = "Significant Regions", data = sig.regions$regions, text = rownames(sig.regions$regions),
                 x = sig.regions$regions$pos, xend = sig.regions$regions$pos + sig.regions$regions$window.span,
                 yend = ~LFC, y = ~LFC) %>%
    add_segments(name = "Merged Regions", data = sig.regions$merged, text = rownames(sig.regions$merged),
                 x = ~`start-position`, xend = ~`end-position`, yend = merged.line, y = merged.line) %>%
    layout(title = paste0("for significant regions: maximum span = ", span.max,
                          "bp, maximum IDR = ", idr.max, ", and minimum |z-score| in both replicates = ", azs.min,
                          "\nfor merged regions: maximum span = ", merged.span.max, "bp"),
           xaxis = list(title = "start position on chr13", zeroline = F, showline = F),
           yaxis = list(title = "LogFC", zeroline = F, showline = F))
}

add.regions <- function(p, regions, name){
  return(add_trace(p, data = regions, x = ~pos, y = ~LFC, text = rownames(regions),
                   name = name, marker = list(size = 5, opacity = 0.5)))
}

add.merged <- function(p, segments, name, merged.line){
  return(add_segments(p, data = segments, x = ~start, xend = ~end, y = merged.line, yend = merged.line,
                      text = paste0(segments$chr, ":", segments$start + 1, "-", segments$end), name = name))
}

read.data <- function(data.path, y = NULL, color = NULL){
  if(grepl("\\.csv$", tolower(data.path))){
    df <- read.csv(data.path, header = T, row.names = 1, check.names = F, stringsAsFactors = F) %>%
      subset(!is.na(pos))
    colnames(df)[1:2] <- c("seqnames", "start")
    if("window.span" %in% colnames(df)){ df$end <- df$start + df$window.span
    } else df$end <- df$start + 20
    if("IDR" %in% colnames(df)) df$`-log10(IDR)` <- -log10(df$IDR + 1e-30)
    colnames(df) <- gsub("^LFC$", "average LFC", colnames(df))
    colnames(df) <- gsub("\\.LFC$", " LFC", colnames(df))
    colnames(df) <- gsub("\\.score$", " z-score", colnames(df))
  } else if(grepl("\\.((bw)|(bigwig))$", tolower(data.path))){
    df <- data.frame(suppressWarnings(import(BigWigFile(data.path))),
                     check.names = F, stringsAsFactors = F)
    df$`|score|` <- abs(df$score)
  } else if(grepl("\\.bed$", tolower(data.path))){
    df <- read.table(data.path, header = F, stringsAsFactors = F)[, 1:3]
    colnames(df) <- c("seqnames", "start", "end")
  } else return(NULL)
  if(!is.null(y) && !(y %in% colnames(df))){
    cat("\"", y, "\" not available for \"", data.path, "\"\n", sep = "")
    return(NULL)
  }
  if(!is.null(color) && !(color %in% colnames(df))){
    cat("\"", color, "\" not available for \"", data.path, "\"\n", sep = "")
    return(NULL)
  }
  return(df %>% subset(seqnames == "chr13"))
}

plot.reg.and.scores <- function(df, reg, reg.ov, y, color, xlim, leg.pos, xt = "", title = ""){
  df <- subset(df, start >= xlim[1] & end <= xlim[2])
  reg <- cbind(subset(reg, start >= xlim[1] & end <= xlim[2]), track = "regions")
  reg.ov <- cbind(subset(reg.ov, start >= xlim[1] & end <= xlim[2]), track = "regions")
  eb <- element_blank()
  gp1 <- ggplot(mapping = aes(x = start, xend = end)) +
    geom_segment(aes(y = 0, yend = 0), reg.ov, color = "red", size = 3) +
    geom_segment(aes(y = -1, yend = -1), reg, color = "blue", size = 3) +
    facet_grid(rows = "track") +
    scale_y_continuous(breaks = c(-1, 0)) + scale_x_continuous(limits = xlim) +
    labs(y = "", x = "", title = "") +
    theme_classic() + theme(axis.line.x = eb, axis.ticks.x = eb, axis.text.x = eb)
  gp2 <- ggplot(mapping = aes(x = start, xend = end)) +
    geom_segment(aes(y = df[, y], yend = df[, y], color = df[, color]), df, size = 2) +
    facet_grid(rows = "track") +
    scale_y_continuous(limits = c(min(c(-1.1, df[, y])), max(c(1.1, df[, y]))), breaks = c(-1, 0, 1)) +
    scale_x_continuous(limits = xlim) +
    scale_color_gradient(high = "#e00000", low = "#d1cdcd") +
    labs(y = y, x = xt, title = title, color = color) +
    theme_classic() + theme(plot.background = element_rect(size = 1, color = "grey"), legend.position = leg.pos)
  return(list(gp1, gp2))
}

score.per.point <- function(df, y, s, e){
  return(sapply(seq(s, e), function(p){
    w <- subset(df, p >= start & p <= end)
    if(nrow(w) == 0){ return(NA)
    } else if(nrow(w) == 1){ return(w[, y])
    } else if("IDR" %in% colnames(w)){ return(w[which.min(w$IDR), y])
    } else return(median(w[, y]))
  }))
}

plot.scores.and.cor <- function(df, reg, y, color, xlim, leg.pos,
                                pps = NULL, hm = NULL, hm.color = "", xt = "", title = "", ...){
  df <- subset(df, start >= xlim[1] & end <= xlim[2])
  reg <- cbind(reg, track = "region")
  eb <- element_blank()
  gp1 <- ggplot() + theme_classic()
  
  if(!is.null(hm)){
    pos <- reg$start:reg$end
    max.y <- sapply(pos, function(p){
      y <- subset(hm, p >= start & p <= end)[, hm.color]
      if(length(y) > 0){ return(max(y))
      } else return(NA)
    })
    gp1 <- gp1 +
      geom_raster(aes(x = pos, y = -1, fill = score),
                  cbind(na.omit(data.frame(pos = pos, score = max.y)), track = "region"), ...) +
      labs(fill = paste("max", hm.color, "per point")) + theme(legend.position = leg.pos)
  } else gp1 <- gp1 +
    geom_segment(aes(x = start, xend = end, y = -1, yend = -1), reg, color = "blue", size = 3)
  
  gp1 <- gp1 +
    facet_grid(rows = "track") +
    scale_y_continuous(breaks = c(-1, -2)) + scale_x_continuous(limits = xlim) +
    labs(y = "", x = "", title = "") +
    theme(plot.background = element_rect(size = 1, color = "grey"),
          axis.line.x = eb, axis.ticks.x = eb, axis.text.x = eb)
  
  gp2 <- ggplot() +
    geom_segment(aes(x = start, xend = end, y = df[, y], yend = df[, y], color = df[, color]), df, size = 2)
  if(!is.null(pps)) gp2 <- gp2 + geom_point(aes(x = pps$pos, y = pps[, 2]))
  gp2 <- gp2 +
    facet_grid(rows = "track") +
    scale_y_continuous(limits = c(min(c(-1.1, df[, y])), max(c(1.1, df[, y]))), breaks = c(-1, 0, 1)) +
    scale_x_continuous(limits = xlim) +
    labs(y = y, x = xt, title = title, color = color) +
    theme_classic() + theme(plot.background = element_rect(size = 1, color = "grey"), legend.position = leg.pos)
  return(list(gp1, gp2))
}


# -------------------------------------------------------------------------
sgRNA.lib.path <- "sgRNA-library/"
scr.res.path <- "XX/"
results.path <- "YY/"
min.count <- 20

# read screening results:
scr.results <- fread(scr.res.path)
index <- data.frame(upstream = grepl("PDX-1Up", scr.results$sgRNA),
                    downstream = grepl("PDX-1Down", scr.results$sgRNA),
                    intron = grepl("PDX-1Intron", scr.results$sgRNA),
                    safe.harbor = grepl("SafeHarbor", scr.results$sgRNA),
                    non.targeting = grepl("NonTargetingControlGuide", scr.results$sgRNA))
index$exon <- !index$upstream &
  !index$downstream &
  !index$intron &
  !index$safe.harbor &
  !index$non.targeting

# read sgRNA library:
pdx1.gr <- list(
  intron.exon = get.genomic.ranges(read.sam.table(paste0(sgRNA.lib.path, "/PDX1.tiling.intron.exon.sam"), T),
                                   28494168 - 20, scr.results$sgRNA[index$intron | index$exon]),
  upstream = get.genomic.ranges(read.sam.table(paste0(sgRNA.lib.path, "/PDX1.tiling.upstream.sam"), F),
                                28494276 - 50000 - 500, scr.results$sgRNA[index$upstream]),
  downstream = get.genomic.ranges(read.sam.table(paste0(sgRNA.lib.path, "/PDX1.tiling.downstream.sam"), F),
                                  28498838 - 500, scr.results$sgRNA[index$downstream])
)

# combine raw counts with gRNA library information:
sgRNAs <- get.sgRNA.table(scr.results, index, pdx1.gr)
colnames(sgRNAs)[5:8] <- c("rep1.GFPplus", "rep2.GFPplus", "rep1.GFPminus", "rep2.GFPminus")
write.csv(sgRNAs, quote = F, row.names = T, file = paste0(results.path, "/0-raw-counts.csv"))
plot.both.LFCs(sgRNAs, "LFC of raw counts.")

# filter out low-count guides:
discard <- apply(sgRNAs[, 5:8], 1, min) < min.count
sgRNAs <- sgRNAs[!discard,]

# calculate size factors and normalize the read counts:
cont.idx <- grepl("NonTargetingControlGuide", sgRNAs$sgRNA)
sf <- cbind(
  rep1 = size.factos(counts = as.matrix(sgRNAs[cont.idx, c(5, 7)])),
  rep2 = size.factos(counts = as.matrix(sgRNAs[cont.idx, c(6, 8)]))
)
sgRNAs <- normalize.counts(sgRNAs, size.factos = sf)

# calculate LogFC:
sgRNAs$rep1.LFC <- Log2.Fold.Change(sgRNAs[, c(5, 7)])
sgRNAs$rep2.LFC <- Log2.Fold.Change(sgRNAs[, c(6, 8)])
write.csv(sgRNAs, quote = F, row.names = T, file = paste0(results.path, "/1-sgRNA-summary_mc=", min.count, ".csv"))
plot.both.LFCs(sgRNAs, "LFC of normalized counts.")

# compare the average LFC of N-guide windows between the replicates:
targets <- sgRNAs[!is.na(sgRNAs$pos),]
targets <- targets[order(targets$pos),]
Ns <- c(1:30, 40, 50); n <- length(Ns)
df <- as.data.frame(t(as.data.frame(sapply(Ns, function(N){
  windows <- average.LFC(targets, N)
  return(c(median(windows$window.span), mean(windows$window.span),
           quantile(windows$window.span, 0.8),
           cor(windows$rep1.LFC, windows$rep2.LFC)))
}), row.names = c("med.sp", "mean.sp", "q80.sp", "cor"))), row.names = Ns)
plot.progress(df)


# -------------------------------------------------------------------------
sgRNA.sumr.path <- "XX/1-sgRNA-summary_mc=20.csv"
results.path <- "YY/results"
N <- 20 # for dCas9-KRAB 

sgRNAs <- read.csv(sgRNA.sumr.path, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
targets <- sgRNAs[!is.na(sgRNAs$pos),]
targets <- targets[order(targets$pos),]
controls <- subset(sgRNAs, grepl("NonTargetingControlGuide", sgRNA))

# average LFC for N-guide windows and calculate z-scores:
windows <- window.scores(targets, controls, N)

# evaluate irreproducibility:
windows$IDR <- est.IDR(abs(windows[, c("rep1.score", "rep2.score")]),
                       mu = 2, sigma = 0.6, rho = 0.6, p = 0.4)$IDR
plot.idr1(windows, controls)

write.csv(windows, file = paste0(results.path, "/2-window-summary_N=", N, ".csv"), quote = F, row.names = T)
ggsave(filename = paste0(results.path, "/", N, "-guide_windows-IDR.pdf"),
       plot = plot.idr2(windows), width = 20, height = 10)
ggsave(filename = paste0(results.path, "/", N, "-guide_windows-z-score.pdf"),
       plot = plot.score(windows), width = 20, height = 20)

plot.gRNAs(windows, title = paste0("span of ", N, "-guide window"),
           name = cut(windows$window.span, c(0, 100, 200, 500, 1000, max(windows$window.span))))

plot.gRNAs(windows, title = "-log(IDR)", color = -log(windows$IDR + 0.0001))
plot.gRNAs(windows, title = "IDR", name = cut(windows$IDR, c(0, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)))

plot.gRNAs(windows, title = "|z-score| in replicate 1", color = abs(windows$rep1.score))
plot.gRNAs(windows, title = "|z-score| in replicate 2", color = abs(windows$rep2.score))


# -------------------------------------------------------------------------
wd.sumr.path <- "AA/results/2-window-summary_N=10.csv"
results.path <- "AA/results"
span.max <- 100
idr.max <- 0.001 # (or 0.0001 for dCas9-KRAB)
azs.min <- 3 # for dCas9-KRAB
merged.span.max <- 200

windows <- read.csv(wd.sumr.path, header = T, row.names = 1, stringsAsFactors = F, check.names = F)
plot.gRNAs(windows, title = paste0("maximum span = ", span.max, "bp, maximum IDR = ", idr.max,
                                   ", and minimum |z-score| in both replicates = ", azs.min),
           name = ifelse(windows$window.span <= span.max &
                           windows$IDR <= idr.max &
                           abs(windows$rep1.score) >= azs.min &
                           abs(windows$rep2.score) >= azs.min,
                         "significant", "insignificant"))

# select significant regions and merge them:
sig.regions <- subset(windows, window.span <= span.max &
                        IDR <= idr.max &
                        abs(rep1.score) >= azs.min &
                        abs(rep2.score) >= azs.min)
merged.sig <- merge.regions(sig.regions, merged.span.max, method = "dist")

write.csv(sig.regions, quote = F, row.names = T,
          file = paste0(results.path, "/3-significant-regions_", "sm=", span.max, "_im=", idr.max, "_zm=", azs.min, ".csv"))
write.table(cbind("chr13", merged.sig), row.names = F, col.names = F, quote = F,
            file = paste0(results.path, "/4-merged-significant-regions_msm=", merged.span.max, ".bed"))

plot.sig.regions(list(regions = sig.regions, merged = merged.sig),
                 span.max, idr.max, azs.min, merged.span.max, merged.line = 0.1)