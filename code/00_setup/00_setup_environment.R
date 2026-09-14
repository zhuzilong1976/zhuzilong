#!/usr/bin/env Rscript
# =============================================================================
# 00_setup_environment.R
# scTenifoldKnk 分析环境：依赖安装 + 版本核查 + 冒烟测试
#
# 用途：在开始任何真实数据分析之前，确认"这台机器能不能跑、跑得动"
# 运行：Rscript 00_setup_environment.R
# 预期：几分钟内输出一份环境报告与一次成功的虚拟敲除结果
# =============================================================================

## ---------------------------------------------------------------------------
## 配置
## ---------------------------------------------------------------------------
OPT <- list(
  ## 国内网络建议使用镜像；如需官方源改为 "https://cloud.r-project.org"
  repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
  out_dir = "outputs/results/smoke_test",
  smoke_cells = 500,        # 冒烟测试细胞数
  smoke_genes = 100,        # 冒烟测试基因数
  run_timing_test = FALSE,  # 设为 TRUE 才跑 300 细胞 × 1000 基因的耗时基准（约 3–4 分钟）
  timing_cells = 300,
  timing_genes = 1000,
  offline = NULL            # NULL = 自动探测；TRUE/FALSE 可强制
)

options(repos = c(CRAN = OPT$repos))
dir.create(OPT$out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n================ 环境报告 ================\n")
cat("R 版本：", R.version.string, "\n")
cat("平台：  ", R.version$platform, "\n")
cat("库路径：", paste(.libPaths(), collapse = " | "), "\n")
cat("CPU 核数：", parallel::detectCores(), "\n")
cat("时间：  ", format(Sys.time()), "\n\n")

if (getRversion() < "4.1") {
  warning("建议使用 R >= 4.1（当前 ", getRversion(), "），较旧版本可能在依赖编译上遇到问题")
}

## ---------------------------------------------------------------------------
## 1. 网络探测 + 依赖检查与安装
## ---------------------------------------------------------------------------
probe_network <- function() {
  f <- tempfile()
  res <- try(
    suppressWarnings(download.file(
      "https://cran.r-project.org/web/packages/scTenifoldKnk/DESCRIPTION", f,
      quiet = TRUE, method = "libcurl"
    )),
    silent = TRUE
  )
  !inherits(res, "try-error") && file.exists(f) && file.info(f)$size > 0
}

is_offline <- isTRUE(OPT$offline)
if (is.null(OPT$offline)) {
  is_offline <- !probe_network()
  cat("网络探测：", if (is_offline) "不可达 CRAN（进入离线模式）" else "正常", "\n")
}

pkgs <- c(
  ## 核心
  "scTenifoldKnk", "scTenifoldNet", "Matrix", "MASS", "cli", "igraph", "reshape2",
  ## 统计与富集
  "locfdr", "fgsea", "preprocessCore",
  ## 可选（数据处理阶段常用）
  "Seurat", "SeuratObject", "Matrix.utils"
)

installed <- rownames(utils::installed.packages())
missing <- setdiff(pkgs, installed)

if (length(missing) > 0) {
  cat("需要安装的包：", paste(missing, collapse = ", "), "\n")
  if (is_offline) {
    cat(">> 当前离线，跳过自动安装。\n")
    cat(">> 离线安装做法：用另一台联网机器下载 CRAN 二进制包（bin/windows/contrib/<R版本>/*.zip），\n")
    cat(">> 然后执行：install.packages(c('a.zip','b.zip'), repos = NULL, type = 'win.binary')\n")
    critical <- intersect(missing, c("scTenifoldKnk", "scTenifoldNet"))
    if (length(critical) > 0) stop("缺少核心包且无法联网安装：", paste(critical, collapse = ", "))
  } else {
    ## scTenifoldKnk 与 scTenifoldNet 均来自 CRAN；preprocessCore 属 Bioconductor
    cran_ok <- setdiff(missing, "preprocessCore")
    if (length(cran_ok) > 0) install.packages(cran_ok, dependencies = TRUE)
    if ("preprocessCore" %in% missing) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("preprocessCore", update = FALSE, ask = FALSE)
    }
  }
} else {
  cat("所有依赖均已安装。\n")
}

## 版本核查
cat("\n================ 关键包版本 ================\n")
for (p in c("scTenifoldKnk", "scTenifoldNet", "Matrix", "locfdr", "fgsea")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "未安装")
  cat(sprintf("  %-16s %s\n", p, v))
}

if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  stop("scTenifoldKnk 未安装成功。请检查网络/镜像，或参考 https://cran.r-project.org/package=scTenifoldKnk")
}
if (utils::packageVersion("scTenifoldKnk") < "1.1") {
  warning("本方案基于 scTenifoldKnk >= 1.1 编写（含 transcriptomeWide 与 dr_empiricalNull），请升级")
}

suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(Matrix)
})

## 导出函数核查
cat("\n================ 导出函数核查 ================\n")
ns <- getNamespaceExports("scTenifoldKnk")
for (f in c("scTenifoldKnk", "scQC", "dRegulation", "plotKO")) {
  cat(sprintf("  %-16s %s\n", f, if (f %in% ns) "OK" else "缺失（版本不匹配？）"))
}

## ---------------------------------------------------------------------------
## 2. 冒烟测试：官方负二项模拟数据上的虚拟敲除
## ---------------------------------------------------------------------------
cat("\n================ 冒烟测试 ================\n")
set.seed(1)
nCells <- OPT$smoke_cells
nGenes <- OPT$smoke_genes
X <- round(matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.98), ncol = nCells))
rownames(X) <- c(paste0("ng", seq_len(nGenes - 10)), paste0("mt-", seq_len(10)))
cat("输入矩阵：", nrow(X), "基因 x ", ncol(X), "细胞\n")

t0 <- Sys.time()
ok <- try({
  out <- scTenifoldKnk::scTenifoldKnk(
    countMatrix   = X,
    gKO           = "ng10",
    nc_nNet       = 5,
    nc_nCells     = min(300, nCells),
    td_K          = 3,
    qc_minLibSize = 30,     # 模拟数据计数很低，需放宽
    nCores        = max(1, parallel::detectCores() - 1)
  )
}, silent = TRUE)

if (inherits(ok, "try-error")) {
  cat("冒烟测试失败：\n", as.character(ok), "\n")
  quit(status = 1)
}

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
cat("完成，用时 ", elapsed, " 分钟\n")

## 输出结构核查
expected_cols <- c("gene", "distance", "Z", "FC", "p.value", "p.adj")
dr <- out$diffRegulation
cat("\n输出结构：\n")
cat("  tensorNetworks: ", paste(names(out$tensorNetworks), collapse = ", "), "\n")
cat("  diffRegulation 列：", paste(colnames(dr), collapse = ", "), "\n")
cat("  列完整性：", if (all(expected_cols %in% colnames(dr))) "OK" else "不完整（版本差异，请检查）", "\n")
cat("  显著基因数 (FDR<0.05)：", sum(dr$p.adj < 0.05), "\n")
cat("\nTop 5 扰动基因：\n")
print(head(dr[, c("gene", "distance", "Z", "p.adj")], 5))

utils::write.csv(dr, file.path(OPT$out_dir, "smoke_diffRegulation.csv"), row.names = FALSE)
saveRDS(out, file.path(OPT$out_dir, "smoke_result.rds"))

## 可选：真实规模耗时基准
if (OPT$run_timing_test) {
  cat("\n================ 耗时基准测试 ================\n")
  set.seed(1)
  nC <- OPT$timing_cells; nG <- OPT$timing_genes
  Y <- round(matrix(rnbinom(n = nG * nC, size = 5, prob = 0.5), ncol = nC))
  rownames(Y) <- c(paste0("g", seq_len(nG - 10)), paste0("MT-", seq_len(10)))
  t1 <- Sys.time()
  invisible(scTenifoldKnk::scTenifoldKnk(
    countMatrix = Y, gKO = "g100", qc_minLibSize = 30,
    nc_nNet = 10, nc_nCells = min(300, nC), td_K = 3,
    nCores = max(1, parallel::detectCores() - 1)
  ))
  mins <- round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1)
  cat(sprintf("基准：%d 细胞 x %d 基因 单基因完整流程 = %s 分钟\n", nC, nG, mins))
  cat("据此估算面板模式（建网一次 + 24 次扰动）与 10 次稳定性复现的总耗时。\n")
}

## ---------------------------------------------------------------------------
## 3. 环境报告落盘
## ---------------------------------------------------------------------------
rep_file <- file.path(OPT$out_dir, "environment_report.txt")
writeLines(
  c(
    paste("R 版本:", R.version.string),
    paste("平台:", R.version$platform),
    paste("时间:", format(Sys.time())),
    paste("CPU 核数:", parallel::detectCores()),
    paste("scTenifoldKnk:", as.character(utils::packageVersion("scTenifoldKnk"))),
    paste("scTenifoldNet:", tryCatch(as.character(utils::packageVersion("scTenifoldNet")), error = function(e) "NA")),
    paste("冒烟测试耗时(分钟):", elapsed),
    paste("冒烟测试显著基因数:", sum(dr$p.adj < 0.05)),
    "",
    capture.output(utils::sessionInfo())
  ),
  rep_file
)

cat("\n================ 结论 ================\n")
cat("环境可用。报告已写入：", normalizePath(rep_file, mustWork = FALSE), "\n")
cat("下一步：\n")
cat("  1) 准备数据（见 outputs/next_steps_CN.md 的 W1）\n")
cat("  2) 填好 outputs/templates/gene_panel.csv（W2）\n")
cat("  3) 修改 outputs/scripts/run_vko_panel.R 的配置区，先以 quick_test = TRUE 跑通\n")
