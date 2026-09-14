# 端到端验证参数（真实数据 GSE279180 小胶质 Control 亚组）
# 注意：这是"跑通流程"的小参数，不是正式分析参数
CFG$quick_test       <- TRUE
CFG$quick_max_cells  <- 200
CFG$quick_max_genes  <- 400
CFG$nc_nNet          <- 5
CFG$nc_nCells        <- 160
CFG$run_stability    <- FALSE   # 真实数据上先只验证主流程；稳定性逻辑已在合成数据验证
CFG$run_enrichment   <- FALSE
CFG$plot_annotate    <- FALSE
CFG$out_dir          <- "work/test_output_real"
CFG$tag              <- "MSmg_Control_smoke"
