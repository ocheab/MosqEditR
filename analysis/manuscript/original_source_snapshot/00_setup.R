options(stringsAsFactors = FALSE)
set.seed(20260911)
required_cran <- c("renv","here","data.table","dplyr","tidyr","stringr","purrr","readr","readxl","writexl","jsonlite","digest","curl","httr2","ggplot2","patchwork","ranger","glmnet","xgboost","pROC","yardstick","rsample","recipes","parsnip","workflows","tune","dials","vip","DALEX","fastshap","igraph","WGCNA","future","furrr","reticulate")
required_bioc <- c("GEOquery","Biobase","limma","biomaRt","AnnotationDbi")
if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org")
if (!file.exists("renv.lock") || length(renv::dependencies()$Package)==0) renv::init(bare=TRUE)
missing_cran <- required_cran[!vapply(required_cran, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
missing_bioc <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask=FALSE, update=FALSE)
dir.create("data_raw", showWarnings=FALSE, recursive=TRUE)
dir.create("data_processed", showWarnings=FALSE, recursive=TRUE)
dir.create("results/tables", showWarnings=FALSE, recursive=TRUE)
dir.create("results/candidate_lists", showWarnings=FALSE, recursive=TRUE)
dir.create("logs", showWarnings=FALSE, recursive=TRUE)
writeLines(capture.output(Sys.time(), sessionInfo()), "logs/setup_environment.txt")

