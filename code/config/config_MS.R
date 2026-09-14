# 疾病网络：MS 全部病灶（MS vs Control 设计）
CFG$quick_test       <- FALSE
CFG$tag              <- "MSvC_MS"
CFG$condition_col    <- "condition"
CFG$control_label    <- "MS"
CFG$extra_filter_col   <- NULL
CFG$extra_filter_value <- NULL
CFG$qc_minLibSize    <- 500
CFG$max_genes        <- 800
CFG$nc_nNet          <- 10      # 官方默认
CFG$nc_nCells        <- 500     # 官方默认
CFG$td_K             <- 3
CFG$nc_nComp         <- 3
CFG$run_stability    <- FALSE
CFG$run_enrichment   <- FALSE
CFG$plot_annotate    <- FALSE
CFG$out_dir          <- "outputs/results/vko_ms_microglia"
