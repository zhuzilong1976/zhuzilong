# 离线安装已下载的 CRAN 二进制包（按拓扑顺序），不联网
pkgsDir <- "work/rpkgs"
order <- readLines(file.path(pkgsDir, "install_order.txt"))
order <- order[nzchar(order)]

zips <- list.files(pkgsDir, pattern = "\\.zip$", full.names = TRUE)
names(zips) <- sub("_[0-9][^_]*\\.zip$", "", basename(zips))

missing <- setdiff(order, names(zips))
if (length(missing)) cat("缺少文件:", paste(missing, collapse = ", "), "\n")

cat("库路径:", paste(.libPaths(), collapse = " | "), "\n")
cat("待安装:", length(order), "个包\n\n")

ok <- character(0); bad <- character(0)
for (p in order) {
  if (!p %in% names(zips)) next
  if (requireNamespace(p, quietly = TRUE)) { cat(sprintf("  [已存在] %s\n", p)); ok <- c(ok, p); next }
  res <- try(
    utils::install.packages(zips[[p]], repos = NULL, type = "win.binary", quiet = TRUE),
    silent = TRUE
  )
  if (inherits(res, "try-error") || !requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [失败] %s : %s\n", p, if (inherits(res, "try-error")) conditionMessage(attr(res, "condition")) else "安装后仍无法载入"))
    bad <- c(bad, p)
  } else {
    cat(sprintf("  [成功] %-18s %s\n", p, as.character(utils::packageVersion(p))))
    ok <- c(ok, p)
  }
}

cat(sprintf("\n安装成功 %d 个，失败 %d 个\n", length(ok), length(bad)))
if (length(bad)) cat("失败列表:", paste(bad, collapse = ", "), "\n")

cat("\n--- 关键包验证 ---\n")
for (p in c("scTenifoldKnk", "scTenifoldNet", "Matrix", "MASS", "cli", "enrichR", "igraph", "reshape2", "locfdr")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "缺失")
  cat(sprintf("  %-16s %s\n", p, v))
}

cat("\n--- 载入测试 ---\n")
res <- try(suppressPackageStartupMessages(library(scTenifoldKnk)), silent = TRUE)
if (inherits(res, "try-error")) {
  cat("scTenifoldKnk 载入失败:\n", conditionMessage(attr(res, "condition")), "\n")
} else {
  cat("scTenifoldKnk 载入成功，导出函数:", paste(getNamespaceExports("scTenifoldKnk"), collapse = ", "), "\n")
}
