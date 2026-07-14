#!/usr/bin/env Rscript

library(sceptre)

# ============================================================
# EDIT THE CELL RANGER MATRIX LIST, INPUT TABLE, AND OUTPUT PATHS
# ============================================================
matrix_dirs_file <- "/path/to/cellranger_results/sceptre_matrix_dirs.txt"
input_dir <- "/path/to/cellranger_to_sceptre/input"
output_dir <- "/path/to/sceptre_results"
# ============================================================

matrix_dirs <- readLines(matrix_dirs_file)
grna_targets <- read.delim(file.path(input_dir, "grna_targets.tsv"))
discovery_pairs <- read.delim(file.path(input_dir, "discovery_pairs.tsv"))
positive_control_pairs <- read.delim(
    file.path(input_dir, "positive_control_pairs.tsv")
)

obj <- import_data_from_cellranger(
    directories = matrix_dirs,
    moi = "high",
    grna_target_data_frame = grna_targets
)

obj <- set_analysis_parameters(
    obj,
    discovery_pairs = discovery_pairs,
    positive_control_pairs = positive_control_pairs,
    side = "right",
    resampling_mechanism = "permutations",
    grna_integration_strategy = "union",
    resampling_approximation = "skew_normal",
    multiple_testing_method = "BH",
    multiple_testing_alpha = 0.1,
    formula_object = ~
        log(response_n_nonzero) +
        log(response_n_umis) +
        log(grna_n_nonzero + 1) +
        log(grna_n_umis + 1) +
        response_p_mito
)

obj <- assign_grnas(obj, method = "mixture")
obj <- run_qc(obj)
obj <- run_calibration_check(obj)
obj <- run_power_check(obj)
obj <- run_discovery_analysis(obj)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_outputs_to_directory(obj, output_dir)
saveRDS(obj, file.path(output_dir, "sceptre_object.rds"))

# Generate annotated discovery results
results <- get_result(obj, "run_discovery_analysis")
parts <- strsplit(results$grna_target, "_")

annotations <- data.frame(
    region_type = sapply(parts, \(x) paste(head(x, -3), collapse = "_")),
    target_chr = sapply(parts, \(x) tail(x, 3)[1]),
    target_start = as.integer(sapply(parts, \(x) tail(x, 2)[1])),
    target_end = as.integer(sapply(parts, \(x) tail(x, 1)))
)

key <- paste(results$response_id, results$grna_target)
discovery_key <- paste(discovery_pairs$response_id, discovery_pairs$grna_target)

annotated_results <- data.frame(
    response_id = results$response_id,
    grna_target = results$grna_target,
    annotations,
    design_target_gene_symbol =
        discovery_pairs$design_target_gene_symbol[
            match(key, discovery_key)
        ],
    n_nonzero_trt = results$n_nonzero_trt,
    n_nonzero_cntrl = results$n_nonzero_cntrl,
    pass_qc = results$pass_qc,
    p_value = results$p_value,
    log_2_fold_change = results$log_2_fold_change,
    significant = results$significant
)

write.csv(
    annotated_results,
    file.path(output_dir, "sceptre_discovery_results_annotated.csv"),
    row.names = FALSE
)
