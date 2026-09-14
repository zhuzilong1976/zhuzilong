# 深度匹配敏感性分析：Control 网络（限制到 800-6000 counts 的细胞）
CFG$quick_test       <- FALSE
CFG$tag              <- "MSvC_Control_depthmatched"
CFG$condition_col    <- "condition"
CFG$control_label    <- "Control"
CFG$extra_filter_col   <- NULL
CFG$extra_filter_value <- NULL
CFG$qc_minLibSize    <- 1500    # 深度窗口下界（经窗口扫描确定：深度比 2.07 -> 1.20）
CFG$qc_maxLibSize    <- 5000    # 深度窗口上界
CFG$max_genes        <- 800
CFG$nc_nNet          <- 10
CFG$nc_nCells        <- 500
CFG$td_K             <- 3
CFG$nc_nComp         <- 3
CFG$run_stability    <- FALSE
CFG$run_enrichment   <- FALSE
CFG$plot_annotate    <- FALSE
CFG$out_dir          <- "outputs/results/vko_ms_depthmatched"
