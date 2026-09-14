#!/usr/bin/env Rscript
# =============================================================================
# scTenifoldKnk 面板虚拟敲除（虚拟敲除 / virtual knockout）分析模板
# 适用：神经免疫疾病 scRNA-seq / snRNA-seq，单细胞类型 × 多候选基因面板
#
# 设计要点
#   1) 一次性建网 + 面板多次扰动（transcriptomeWide = TRUE）——省去重复建网开销
#   2) 与官方 dRegulation() 完全一致的 Z / FC / p / FDR 换算
#   3) 多种子重复的稳定性分析（Spearman rho + Top-N 一致性）
#   4) 阴性对照零分布
#   5) 结果表、图、参数、sessionInfo 全量落盘（可复现）
#
# 使用：先安装依赖 → 修改下面 CFG 与 GENE_PANEL → Rscript run_vko_panel.R
# 首次务必把 CFG$quick_test 设为 TRUE 跑通小规模，再跑全量。
# =============================================================================

## ---------------------------------------------------------------------------
## 0. 依赖
## ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  ok <- requireNamespace("scTenifoldKnk", quietly = TRUE)
  if (!ok) stop("请先安装：install.packages('scTenifoldKnk')  # CRAN v1.1+")
  library(scTenifoldKnk)
  library(Matrix)
  if (!requireNamespace("scTenifoldNet", quietly = TRUE)) stop("缺少 scTenifoldNet")
})

## ---------------------------------------------------------------------------
## 1. 配置区（按你的项目修改）
## ---------------------------------------------------------------------------
CFG <- list(
  ## 输入
  counts_file   = "work/data/counts.mtx",       # GSE279180 人 MS 病灶小胶质细胞
  genes_file    = "work/data/genes.txt",        # 若 counts 为 mtx，需提供行名
  cells_file    = "work/data/barcodes.txt",
  meta_file     = "work/data/metadata.csv",     # 含 barcode, celltype, condition, donor, lesion_type, subtype
  celltype_col  = "celltype",
  target_celltype = "Microglia",
  condition_col = "condition",
  control_label = "Control",                    # 第一套网络：正常对照；改成 "MS" 则建疾病网络
  donor_col     = "donor",
  ## 可选的额外筛选：例如做"慢性非活动病灶 CI vs 慢性活动病灶 CA"时设
  ## extra_filter_col = "lesion_type",  extra_filter_value = "CI"
  extra_filter_col   = NULL,
  extra_filter_value = NULL,

  ## 输出
  out_dir       = "outputs/results/vko_ms_microglia",
  tag           = "MS_microglia_Control",

  ## 质控（官方 scQC 参数）
  qc_minLibSize = 1000,
  qc_maxLibSize = NULL,        # 文库量上限（做深度匹配敏感性分析时使用）
  qc_minPCT     = 0.05,
  qc_maxMTratio = 0.10,
  drop_ribo_mt  = TRUE,        # 过滤 ^Rp[0-9]|^Rpl|^Rps|^Mt- （官方脚本做法）
  shared_gene_set = TRUE,      # 基因过滤基于"全部目标细胞"（两条件合并），保证两套网络基因集一致

  ## 建网参数
  nc_nNet   = 10,
  nc_nCells = 500,
  nc_nComp  = 3,
  nc_q      = 0.90,
  nc_lambda = 0,
  td_K      = 3,
  ma_nDim   = 2,
  empirical_null = TRUE,       # 用 locfdr 经验零分布（需安装 locfdr）

  ## 面板设计
  max_genes      = 3000,       # 每个网络保留多少基因（按方差取 top，且强制保留面板基因）
  max_cells_per_donor = NULL,  # 单供体细胞数上限（防止某个供体主导网络）
  top_n_per_gene = 200,        # 每个基因输出的 Top 扰动基因数
  run_stability  = TRUE,       # 稳定性分析（耗时 ≈ n_rep 次建网！）
  n_rep          = 10,
  n_cells_rep    = 500,
  run_enrichment = TRUE,
  plot_annotate  = FALSE,      # plotKO 的富集注释需要联网访问 Enrichr，离线时设为 FALSE

  ## 快速试跑（用少量细胞 + 少量基因先验证流程）
  quick_test      = TRUE,
  quick_max_cells = 300,
  quick_max_genes = 1000,

  seed = 1,
  nCores = max(1, parallel::detectCores() - 1)
)

## 稀疏矩阵上的方差：Var(x) = E[x^2] - E[x]^2（避免把稀疏矩阵转成稠密）
row_variance <- function(M) {
  m1 <- Matrix::rowMeans(M)
  m2 <- Matrix::rowMeans(M * M)
  as.numeric(m2 - m1^2)
}

## 参数覆盖：命令行第一个参数可指定覆盖文件，默认 work/vko_overrides.R
## 用法：Rscript outputs/scripts/run_vko_panel.R work/config_ms.R
## 这样可以把"对照网络"和"疾病网络"写成两份配置并行运行
args_cli <- commandArgs(trailingOnly = TRUE)
ovr_file <- if (length(args_cli) >= 1 && nzchar(args_cli[1])) args_cli[1] else "work/vko_overrides.R"
if (file.exists(ovr_file)) {
  message("应用参数覆盖：", ovr_file)
  source(ovr_file)
}

## 候选基因面板：20 个候选 + 4 个阴性对照（阴性对照应在本细胞类型中低表达/生物学无关）
GENE_PANEL <- c(
  ## --- 候选基因：MS/AD 髓系风险与神经免疫核心节点（均已确认存在于 GSE279180 矩阵）---
  "TREM2", "TYROBP", "CD33", "APOE", "BIN1",
  "PICALM", "CR1", "MS4A4A", "SPI1", "IRF8",
  "STAT3", "NFKB1", "C1QA", "C3", "CX3CR1",
  "P2RY12", "TMEM119", "BTK", "CSF1R", "IL10RA",
  ## --- 阴性对照：在本数据中表达率 5%-30%、但与 MS 小胶质生物学无关 ---
  "ACTA2",   # 平滑肌，6.4%
  "TSHR",    # 促甲状腺激素受体，8.0%
  "PECAM1",  # 内皮，30.0%（注意可能部分来自血管环境 RNA）
  "YKT6"     # 低表达且方差极低，5.0%
)
NEGATIVE_CONTROLS <- tail(GENE_PANEL, 4)

set.seed(CFG$seed)
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
log_con <- file(file.path(CFG$out_dir, paste0(CFG$tag, "_run.log")), open = "wt")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n", file = log_con, append = TRUE)
  message(msg)
}
log_msg("scTenifoldKnk 版本：", as.character(utils::packageVersion("scTenifoldKnk")))
log_msg("输出目录：", normalizePath(CFG$out_dir, mustWork = FALSE))

## ---------------------------------------------------------------------------
## 2. 读取数据
## ---------------------------------------------------------------------------
read_counts <- function(cfg) {
  ext <- tolower(tools::file_ext(cfg$counts_file))
  X <- switch(
    ext,
    "mtx" = {
      m <- Matrix::readMM(cfg$counts_file)
      rownames(m) <- readLines(cfg$genes_file)
      colnames(m) <- readLines(cfg$cells_file)
      m
    },
    "rds" = readRDS(cfg$counts_file),
    "csv" = {
      df <- utils::read.csv(cfg$counts_file, row.names = 1, check.names = FALSE)
      Matrix::Matrix(as.matrix(df))
    },
    stop("不支持的 counts 格式：", ext)
  )
  X <- as(X, "CsparseMatrix")
  if (is.null(rownames(X)) || is.null(colnames(X))) stop("counts 必须有基因名(行)与细胞名(列)")
  if (max(X, na.rm = TRUE) < 50) warning("数值看起来不像原始 counts，scTenifoldKnk 需要未归一化的原始 counts")
  X
}

X <- read_counts(CFG)
meta <- utils::read.csv(CFG$meta_file, stringsAsFactors = FALSE)
log_msg("原始矩阵：", nrow(X), " 基因 x ", ncol(X), " 细胞")

## ---------------------------------------------------------------------------
## 2b. 共享基因集：基因过滤必须基于两条件合并的细胞，否则两次运行的网络基因集不同
## ---------------------------------------------------------------------------
cells_target <- intersect(meta$barcode[meta[[CFG$celltype_col]] == CFG$target_celltype], colnames(X))
shared_genes <- NULL
if (isTRUE(CFG$shared_gene_set)) {
  Xall <- X[, cells_target, drop = FALSE]
  rate_all <- Matrix::rowMeans(Xall != 0)
  shared_genes <- rownames(X)[rate_all > CFG$qc_minPCT]
  log_msg("共享基因集（", CFG$target_celltype, " 全部 ", length(cells_target),
          " 个细胞，表达率 >", CFG$qc_minPCT, "）：", length(shared_genes), " 个基因")
  ## 关键：基因数上限也必须基于"两条件合并"的细胞来选，否则两套网络的基因集不同、无法比较
  if (!is.null(CFG$max_genes) && length(shared_genes) > CFG$max_genes) {
    v_all <- row_variance(Xall[shared_genes, , drop = FALSE])
    names(v_all) <- shared_genes
    top <- names(sort(v_all, decreasing = TRUE))[seq_len(CFG$max_genes)]
    shared_genes <- union(top, intersect(GENE_PANEL, shared_genes))
    log_msg("共享基因集按方差裁剪至 ", length(shared_genes), " 个（含全部面板基因）")
  }
  rm(Xall); gc()
}

## ---------------------------------------------------------------------------
## 3. 选择细胞（细胞类型 + 条件）
## ---------------------------------------------------------------------------
sel <- meta[[CFG$celltype_col]] == CFG$target_celltype
if (!is.null(CFG$condition_col)) {
  sel <- sel & meta[[CFG$condition_col]] == CFG$control_label
}
if (!is.null(CFG$extra_filter_col)) {
  sel <- sel & meta[[CFG$extra_filter_col]] == CFG$extra_filter_value
}
cells <- intersect(meta$barcode[sel], colnames(X))
if (length(cells) < 200) stop("可用细胞数过少（", length(cells), "），至少需要数百个细胞")
X <- X[, cells, drop = FALSE]
meta_sub <- meta[match(cells, meta$barcode), ]
log_msg("目标细胞（", CFG$target_celltype, " / ", CFG$control_label, "）：",
        ncol(X), " 细胞，来自 ", length(unique(meta_sub[[CFG$donor_col]])), " 个供体")

## 供体均衡：限制单个供体贡献的细胞数，避免某个供体主导网络
if (!is.null(CFG$max_cells_per_donor)) {
  set.seed(CFG$seed)
  keep_cells <- unlist(lapply(split(meta_sub$barcode, meta_sub[[CFG$donor_col]]), function(b) {
    if (length(b) > CFG$max_cells_per_donor) sample(b, CFG$max_cells_per_donor) else b
  }), use.names = FALSE)
  keep_cells <- intersect(keep_cells, colnames(X))
  log_msg("供体均衡：每供体上限 ", CFG$max_cells_per_donor, " 细胞，保留 ", length(keep_cells), " 个细胞")
  X <- X[, keep_cells, drop = FALSE]
  meta_sub <- meta[match(keep_cells, meta$barcode), ]
}

## ---------------------------------------------------------------------------
## 4. 质控（官方 scQC 同款参数）+ 核糖体/线粒体基因过滤
## ---------------------------------------------------------------------------
X <- scTenifoldKnk::scQC(
  X,
  minLibSize = CFG$qc_minLibSize,
  removeOutlierCells = TRUE,
  ## 启用共享基因集时，此处的基因过滤交给下面的 shared_genes 处理（保证两条件一致）
  minPCT = if (is.null(shared_genes)) CFG$qc_minPCT else 0,
  maxMTratio = CFG$qc_maxMTratio,
  label = CFG$tag
)
if (!is.null(shared_genes)) {
  n_before <- nrow(X)
  X <- X[intersect(rownames(X), shared_genes), ]
  log_msg("应用共享基因集：", n_before, " -> ", nrow(X), " 个基因")
}
if (CFG$drop_ribo_mt) {
  keep <- !grepl("^Rp[[:digit:]]|^Rpl|^Rps|^Mt-", rownames(X), ignore.case = TRUE)
  log_msg("过滤核糖体/线粒体基因：", sum(!keep), " 个")
  X <- X[keep, ]
}

## 文库量上限：用于"深度匹配"的敏感性分析（避免两组测序深度差异成为混杂）
if (!is.null(CFG$qc_maxLibSize)) {
  lib <- Matrix::colSums(X)
  keep_c <- lib <= CFG$qc_maxLibSize
  log_msg("文库量上限 ", CFG$qc_maxLibSize, "（深度匹配）：保留 ",
          sum(keep_c), "/", length(keep_c), " 个细胞")
  X <- X[, keep_c, drop = FALSE]
}
log_msg("QC 后：", nrow(X), " 基因 x ", ncol(X), " 细胞")

## 快速试跑模式：限制细胞数与基因数（基因按表达方差取 top）
if (CFG$quick_test) {
  if (ncol(X) > CFG$quick_max_cells) X <- X[, sample(ncol(X), CFG$quick_max_cells)]
  if (nrow(X) > CFG$quick_max_genes) {
    v <- row_variance(X)
    names(v) <- rownames(X)
    top <- names(sort(v, decreasing = TRUE))[seq_len(CFG$quick_max_genes)]
    ## 保证面板基因不被方差筛选剔除，否则 quick_test 会因缺失候选基因而失败
    keep <- union(top, intersect(GENE_PANEL, rownames(X)))
    X <- X[keep, ]
  }
  log_msg(">>> quick_test 模式：", nrow(X), " 基因 x ", ncol(X), " 细胞")
}

## 面板基因必须在网络中；缺失的剔除（注意：阴性对照缺失说明该细胞类型不表达，属预期）
panel_requested <- GENE_PANEL
GENE_PANEL <- GENE_PANEL[GENE_PANEL %in% rownames(X)]
missing_panel <- setdiff(panel_requested, GENE_PANEL)
log_msg("面板中可用基因：", length(GENE_PANEL), " / ", length(panel_requested), " 个")
if (length(missing_panel) > 0) {
  log_msg("!! 未进入网络的面板基因：", paste(missing_panel, collapse = ", "),
          "（原因通常是表达率低于过滤阈值）")
}
if (length(GENE_PANEL) < 3) stop("面板可用基因过少，请检查基因命名（symbol 体系）")

## ---------------------------------------------------------------------------
## 5. 距离 → 统计量（与 scTenifoldKnk::dRegulation() 内部实现一致）
## ---------------------------------------------------------------------------
dr_from_distances <- function(d, empirical_null = FALSE) {
  stopifnot(all(is.finite(d)), all(d >= 0))
  ## Box-Cox（λ 在 [-2, 2] 上网格搜索，λ = 0 处用对数替代）
  lam <- seq(-2, 2, length.out = 1000); lam <- lam[lam != 0]
  bc <- try(MASS::boxcox(d ~ 1, plot = FALSE, lambda = lam), silent = TRUE)
  nD <- if (inherits(bc, "try-error")) d else {
    l <- bc$x[which.max(bc$y)]
    if (l < 0) 1 / (d^l) else d^l
  }
  Z  <- as.numeric(scale(nD))
  FC <- d^2 / mean(d^2)
  if (isTRUE(empirical_null) && requireNamespace("locfdr", quietly = TRUE)) {
    en <- try(locfdr::locfdr(Z, plot = 0), silent = TRUE)
    p <- if (inherits(en, "try-error")) pchisq(FC, df = 1, lower.tail = FALSE) else {
      d0 <- en$fp0["mlest", "delta"]; s0 <- en$fp0["mlest", "sigma"]
      pnorm((Z - d0) / s0, lower.tail = FALSE)
    }
  } else {
    p <- pchisq(FC, df = 1, lower.tail = FALSE)
  }
  data.frame(gene = NA_character_, distance = d, Z = Z, FC = FC,
             p.value = p, p.adj = p.adjust(p, method = "fdr"))
}

panel_stats <- function(dist_mat) {
  out <- lapply(rownames(dist_mat), function(g) {
    d <- as.numeric(dist_mat[g, ])
    if (any(!is.finite(d))) return(NULL)
    df <- dr_from_distances(d, empirical_null = CFG$empirical_null)
    df$gene <- colnames(dist_mat)
    df$ko_gene <- g
    df[order(df$p.value), ]
  })
  do.call(rbind, out)
}

## ---------------------------------------------------------------------------
## 6. 主分析：一次建网 + 面板扰动
## ---------------------------------------------------------------------------
run_panel <- function(mat, panel, tag) {
  log_msg("建网 + 面板扰动开始（", tag, "）：", nrow(mat), " 基因 x ", ncol(mat), " 细胞")
  if (ncol(mat) < 500) {
    warning("细胞数 < 500：网络稳定性会明显下降（参见 2026 年 8 方法评测的最低下限建议）")
  }
  ## 每次子采样细胞数：不超过 500，且不超过 80% 的可用细胞（避免子采样退化为同一套细胞）
  n_cells_net <- max(50, min(CFG$nc_nCells, floor(0.8 * ncol(mat))))
  t0 <- Sys.time()
  res <- scTenifoldKnk::scTenifoldKnk(
    countMatrix       = mat,
    gKO               = panel,
    transcriptomeWide = TRUE,
    qc                = FALSE,          # 已在上游用 scQC 处理
    nc_nNet           = CFG$nc_nNet,
    nc_nCells         = n_cells_net,
    nc_nComp          = CFG$nc_nComp,
    nc_q              = CFG$nc_q,
    nc_lambda         = CFG$nc_lambda,
    td_K              = CFG$td_K,
    ma_nDim           = CFG$ma_nDim,
    nCores            = CFG$nCores
  )
  log_msg("完成，用时 ", round(difftime(Sys.time(), t0, units = "mins"), 1), " 分钟")

  dm <- res$perturbationDistances
  saveRDS(res, file.path(CFG$out_dir, paste0(tag, "_raw_result.rds")))
  utils::write.csv(as.data.frame(dm), file.path(CFG$out_dir, paste0(tag, "_distanceMatrix.csv")))

  ## 诊断：出度为 0 的基因，其"虚拟敲除"是空操作（网络中该行本来全为 0），结果不可解释
  wt <- res$tensorNetworks$WT
  od <- Matrix::rowSums(wt != 0)
  zero_od <- names(od)[od == 0]
  log_msg("网络诊断：", nrow(wt), " 基因中出度为 0 的共 ", length(zero_od), " 个")
  panel_zero <- intersect(panel, zero_od)
  if (length(panel_zero) > 0) {
    log_msg("!! 面板内出度为 0 的基因（敲除无效，需从结果中剔除或增大网络）：",
            paste(panel_zero, collapse = ", "))
  }
  utils::write.csv(data.frame(gene = names(od), outdegree = as.integer(od)),
                   file.path(CFG$out_dir, paste0(tag, "_outdegree.csv")), row.names = FALSE)

  st <- panel_stats(dm)
  utils::write.csv(st, file.path(CFG$out_dir, paste0(tag, "_panel_allStats.csv")), row.names = FALSE)

  ## 每个敲除基因的 Top-N 扰动基因
  top <- do.call(rbind, lapply(split(st, st$ko_gene), function(x) head(x, CFG$top_n_per_gene)))
  utils::write.csv(top, file.path(CFG$out_dir, paste0(tag, "_panel_topHits.csv")), row.names = FALSE)

  ## 面板汇总：显著基因数、最小 FDR、Z 分布
  summ <- do.call(rbind, lapply(split(st, st$ko_gene), function(x) data.frame(
    ko_gene = x$ko_gene[1],
    n_sig_FDR05 = sum(x$p.adj < 0.05),
    n_sig_FDR01 = sum(x$p.adj < 0.01),
    min_p_adj   = min(x$p.adj),
    mean_Z      = mean(x$Z),
    sd_Z        = sd(x$Z),
    is_negative_control = x$ko_gene[1] %in% NEGATIVE_CONTROLS
  )))
  summ <- summ[order(-summ$n_sig_FDR05), ]
  utils::write.csv(summ, file.path(CFG$out_dir, paste0(tag, "_panel_summary.csv")), row.names = FALSE)
  log_msg("面板汇总完成：显著基因数最多的是 ", summ$ko_gene[1], "（", summ$n_sig_FDR05[1], " 个 FDR<0.05）")
  list(dist = dm, stats = st, summary = summ)
}

main <- run_panel(X, GENE_PANEL, CFG$tag)

## ---------------------------------------------------------------------------
## 7. 阴性对照零分布
## ---------------------------------------------------------------------------
neg <- main$summary[main$summary$is_negative_control, ]
pos <- main$summary[!main$summary$is_negative_control, ]
if (nrow(neg) > 0 && nrow(pos) > 0) {
  log_msg("阴性对照显著基因数（中位数）：", stats::median(neg$n_sig_FDR05),
          "；候选基因（中位数）：", stats::median(pos$n_sig_FDR05))
  if (stats::median(neg$n_sig_FDR05) > 0.5 * stats::median(pos$n_sig_FDR05)) {
    warning("阴性对照的扰动信号与候选基因相当，提示网络/参数或细胞数存在问题，结论需谨慎")
  }
  png(file.path(CFG$out_dir, paste0(CFG$tag, "_negControl.png")), width = 1400, height = 900, res = 200)
  op <- graphics::par(mar = c(10, 4, 3, 1))
  bp <- graphics::barplot(main$summary$n_sig_FDR05, names.arg = main$summary$ko_gene, las = 2,
                          col = ifelse(main$summary$is_negative_control, "grey70", "steelblue"),
                          main = paste0("虚拟敲除面板：显著扰动基因数 (", CFG$tag, ")"),
                          ylab = "FDR < 0.05 的基因数")
  graphics::abline(h = stats::median(neg$n_sig_FDR05), lty = 2, col = "red")
  graphics::legend("topright", legend = c("候选基因", "阴性对照", "阴性对照中位数"),
                   fill = c("steelblue", "grey70", NA), border = NA, lty = c(NA, NA, 2),
                   col = c(NA, NA, "red"), bty = "n")
  graphics::par(op); dev.off()
}

## ---------------------------------------------------------------------------
## 8. 稳定性分析：多种子 × 细胞子采样，重复建网与扰动
## ---------------------------------------------------------------------------
if (CFG$run_stability) {
  n_rep <- if (CFG$quick_test) min(3, CFG$n_rep) else CFG$n_rep
  ## 子采样细胞数必须小于可用细胞数，否则"重复"会退化成同一套细胞
  n_cells_rep <- max(50, min(CFG$n_cells_rep, floor(0.8 * ncol(X))))
  log_msg("稳定性分析：", n_rep, " 次重复 × ", n_cells_rep, " 细胞（预计耗时 ≈ ", n_rep, " 次建网）")

  z_list <- lapply(seq_len(n_rep), function(r) {
    set.seed(CFG$seed + r)
    idx <- sample(ncol(X), n_cells_rep)
    rr <- scTenifoldKnk::scTenifoldKnk(
      countMatrix = X[, idx], gKO = GENE_PANEL, transcriptomeWide = TRUE, qc = FALSE,
      nc_nNet = CFG$nc_nNet, nc_nCells = max(50, min(CFG$nc_nCells, floor(0.8 * n_cells_rep))),
      nc_nComp = CFG$nc_nComp, nc_q = CFG$nc_q, nc_lambda = CFG$nc_lambda,
      td_K = CFG$td_K, ma_nDim = CFG$ma_nDim, nCores = CFG$nCores
    )
    st <- panel_stats(rr$perturbationDistances)
    ## 每个基因取 Z 向量
    z <- do.call(cbind, lapply(split(st, st$ko_gene), function(x) {
      v <- x$Z; names(v) <- toupper(x$gene); v
    }))
    rownames(z) <- toupper(rownames(z))
    z
  })

  ## 每个敲除基因：跨重复的 Z 向量 Spearman 相关 + Top-50 重叠
  stab <- do.call(rbind, lapply(GENE_PANEL, function(g) {
    M <- sapply(z_list, function(z) z[, g])
    M <- M[stats::complete.cases(M), , drop = FALSE]
    if (ncol(M) < 2 || nrow(M) < 10) {
      return(data.frame(ko_gene = g, n_genes = nrow(M), rho_mean = NA, rho_sd = NA, topK_jaccard = NA))
    }
    ## 分位数归一化（官方 STABILITY 脚本做法）；未安装 preprocessCore 时退化为原始 Z
    if (requireNamespace("preprocessCore", quietly = TRUE)) {
      Mq <- preprocessCore::normalize.quantiles(M)
      rownames(Mq) <- rownames(M)
    } else {
      Mq <- M
    }
    rho <- stats::cor(Mq, method = "sp")
    k <- min(50, nrow(Mq))
    tops <- apply(Mq, 2, function(v) names(sort(v, decreasing = TRUE))[seq_len(k)])
    jac <- mean(apply(utils::combn(seq_len(ncol(tops)), 2), 2, function(p) {
      length(intersect(tops[, p[1]], tops[, p[2]])) / length(union(tops[, p[1]], tops[, p[2]]))
    }))
    data.frame(ko_gene = g, n_genes = nrow(M),
               rho_mean = mean(rho[lower.tri(rho)], na.rm = TRUE),
               rho_sd   = stats::sd(rho[lower.tri(rho)], na.rm = TRUE),
               topK_jaccard = jac)
  }))
  utils::write.csv(stab, file.path(CFG$out_dir, paste0(CFG$tag, "_stability.csv")), row.names = FALSE)
  saveRDS(z_list, file.path(CFG$out_dir, paste0(CFG$tag, "_stability_zlist.rds")))
  log_msg("稳定性完成：rho 中位数 = ", round(stats::median(stab$rho_mean, na.rm = TRUE), 3))
}

## ---------------------------------------------------------------------------
## 9. 富集分析（对 Z 值排序做 fgsea；需要网络）
## ---------------------------------------------------------------------------
if (CFG$run_enrichment && requireNamespace("fgsea", quietly = TRUE)) {
  libs <- c(KEGG = "KEGG_2021_Human", Reactome = "Reactome_2022",
            GO_BP = "GO_Biological_Process_2021")
  top_pos <- head(main$summary$ko_gene[!main$summary$is_negative_control], 5)
  for (g in top_pos) {
    st <- main$stats[main$stats$ko_gene == g, ]
    z <- st$Z; names(z) <- toupper(st$gene)
    z <- sort(z[is.finite(z)], decreasing = TRUE)
    for (nm in names(libs)) {
      res <- try({
        p <- fgsea::gmtPathways(paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", libs[[nm]]))
        as.data.frame(fgsea::fgseaMultilevel(p, z, minSize = 5, maxSize = 500))
      }, silent = TRUE)
      if (!inherits(res, "try-error") && nrow(res) > 0) {
        res <- res[order(res$padj), ]
        utils::write.csv(res, file.path(CFG$out_dir, paste0(CFG$tag, "_enrich_", g, "_", nm, ".csv")), row.names = FALSE)
      } else {
        log_msg("富集分析跳过（", g, " / ", nm, "）：无网络或库加载失败")
      }
    }
  }
}

## ---------------------------------------------------------------------------
## 10. 子网络图（plotKO，需 enrichR 联网；失败不影响其它结果）
## ---------------------------------------------------------------------------
if (requireNamespace("igraph", quietly = TRUE)) {
  try({
    top1 <- main$summary$ko_gene[!main$summary$is_negative_control][1]
    single <- scTenifoldKnk::scTenifoldKnk(countMatrix = X, gKO = top1, qc = FALSE,
                                           nc_nNet = CFG$nc_nNet,
                                           nc_nCells = max(50, min(CFG$nc_nCells, floor(0.8 * ncol(X)))),
                                           nc_nComp = CFG$nc_nComp, nc_q = CFG$nc_q,
                                           td_K = CFG$td_K, nCores = CFG$nCores)
    saveRDS(single, file.path(CFG$out_dir, paste0(CFG$tag, "_singleKO_", top1, ".rds")))
    png(file.path(CFG$out_dir, paste0(CFG$tag, "_subnetwork_", top1, ".png")),
        width = 2400, height = 2400, res = 250)
    print(scTenifoldKnk::plotKO(single, gKO = top1,
                                annotate = isTRUE(CFG$plot_annotate), fdrThreshold = 0.05))
    dev.off()
  }, silent = TRUE)
}

## ---------------------------------------------------------------------------
## 11. 可复现性记录
## ---------------------------------------------------------------------------
writeLines(capture.output(utils::sessionInfo()),
           file.path(CFG$out_dir, paste0(CFG$tag, "_sessionInfo.txt")))
writeLines(capture.output(str(CFG)), file.path(CFG$out_dir, paste0(CFG$tag, "_params.txt")))
writeLines(GENE_PANEL, file.path(CFG$out_dir, paste0(CFG$tag, "_genePanel.txt")))
log_msg("全部完成，结果目录：", normalizePath(CFG$out_dir, mustWork = FALSE))
close(log_con)

## ---------------------------------------------------------------------------
## 附：部分"敲低（knockdown）"扩展（原生只支持出边置零的敲除）
##     论文中必须注明这是对原始方法的扩展，未被官方验证。
## ---------------------------------------------------------------------------
partial_knockdown <- function(wt_network, gene, kd = 0.5, d = 2) {
  stopifnot(kd >= 0, kd <= 1)
  KO <- wt_network
  KO[gene, ] <- KO[gene, ] * (1 - kd)     # kd = 0 对照，kd = 1 等价于完全敲除
  MA <- scTenifoldNet::manifoldAlignment(wt_network, KO, d = d)
  scTenifoldKnk::dRegulation(MA, empiricalNull = CFG$empirical_null)
}
# 用法：WT <- readRDS(".../singleKO_X.rds")$tensorNetworks$WT
#       dr <- partial_knockdown(WT, "TREM2", kd = 0.5)
