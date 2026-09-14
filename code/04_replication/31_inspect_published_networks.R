files <- list.files("work/author_networks", pattern = "\\.RData$", full.names = TRUE)
for (f in files) {
  cat("\n================ ", basename(f), " ================\n")
  e <- new.env()
  ok <- try(load(f, envir = e), silent = TRUE)
  if (inherits(ok, "try-error")) { cat("加载失败\n"); next }
  objs <- ls(e)
  cat("对象:", paste(objs, collapse = ", "), "\n")
  for (o in objs) {
    x <- get(o, envir = e)
    cat(sprintf("  %s: class=%s", o, paste(class(x), collapse = "/")))
    if (is.list(x)) cat("  names=", paste(head(names(x), 8), collapse = ","))
    cat("\n")
    if (is.list(x) && !is.null(x$diffRegulation)) {
      dr <- x$diffRegulation
      cat(sprintf("     diffRegulation: %d 行；列 %s\n", nrow(dr),
                  paste(colnames(dr), collapse = ",")))
      cat(sprintf("     p.adj<0.05: %d；距离最大的基因: %s\n",
                  sum(dr$p.adj < 0.05), dr$gene[which.max(dr$distance)]))
    }
    if (is.list(x) && !is.null(x$WT)) {
      cat(sprintf("     WT: %d x %d\n", nrow(x$WT), ncol(x$WT)))
    }
  }
}
