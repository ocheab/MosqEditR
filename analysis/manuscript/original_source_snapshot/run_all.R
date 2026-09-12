source("R/00_setup.R")
steps <- c("01_acquire_genome_annotation.R","02_acquire_mozatlas_expression.R","03_acquire_orthology_phenotypes.R","04_acquire_ag1000g_variation.R","05_crispr_editability.R","06_build_feature_matrix.R","07_model_training.R","08_validation_and_ranking.R","09_explainability.R","10_figures_tables.R","11_update_manuscript_tables.R","99_capture_session.R")
for (s in steps) {
  message("\n===== Running ", s, " =====")
  source(file.path("R", s), local=new.env(parent=globalenv()))
}

