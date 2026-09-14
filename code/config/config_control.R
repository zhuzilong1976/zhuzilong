# 对照网络：Control 组小胶质（MS vs Control 设计）
CFG$quick_test       <- FALSE
CFG$tag              <- "MSvC_Control"
CFG$condition_col    <- "condition"
CFG$control_label    <- "Control"
CFG$extra_filter_col   <- NULL
CFG$extra_filter_value <- NULL
CFG$qc_minLibSize    <- 500     # 默认 1000 会把 3 个对照样本几乎清空
CFG$max_genes        <- 800
CFG$nc_nNet          <- 10      # 官方默认：网络数不足会稀释扰动信号
CFG$nc_nCells        <- 500     # 官方默认：每网络细胞数
CFG$td_K             <- 3
CFG$nc_nComp         <- 3
CFG$run_stability    <- FALSE
CFG$run_enrichment   <- FALSE
CFG$plot_annotate    <- FALSE
CFG$out_dir          <- "outputs/results/vko_ms_microglia"
