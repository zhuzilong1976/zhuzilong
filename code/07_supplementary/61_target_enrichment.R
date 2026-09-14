# B 方案：基于扰动谱的靶基因集富集分析（超几何检验，无需联网）
# 思路：对每个敲除基因，取"扰动距离排名前 N 的基因"作为靶基因集，做通路富集

suppressPackageStartupMessages({ library(Matrix) })

dir <- "outputs/results/vko_ms_microglia"
TOP_N <- 100          # 每个敲除基因取距离排名前多少的基因
MIN_SET <- 5          # 通路至少含多少背景基因才纳入检验
MAX_SET <- 300

## ---- 1. 读入扰动矩阵 ----
read_mat <- function(tag) {
  as.matrix(utils::read.csv(file.path(dir, paste0(tag, "_distanceMatrix.csv")),
                            row.names = 1, check.names = FALSE))
}
ctrl <- read_mat("MSvC_Control")
ms <- read_mat("MSvC_MS")
stopifnot(identical(dim(ctrl), dim(ms)), identical(colnames(ctrl), colnames(ms)))

## 出度为 0 的基因（空操作，统计量退化）必须剔除
read_od <- function(tag) {
  o <- utils::read.csv(file.path(dir, paste0(tag, "_outdegree.csv")), stringsAsFactors = FALSE)
  stats::setNames(o$outdegree, o$gene)
}
od_ctrl <- read_od("MSvC_Control")
od_ms <- read_od("MSvC_MS")
ko_ctrl <- rownames(ctrl)[od_ctrl[rownames(ctrl)] > 0]
ko_ms <- rownames(ms)[od_ms[rownames(ms)] > 0]
cat("Control 可用敲除基因：", length(ko_ctrl), "个；剔除出度为 0 的：",
    paste(setdiff(rownames(ctrl), ko_ctrl), collapse = ", "), "\n")
cat("MS 可用敲除基因：", length(ko_ms), "个；剔除出度为 0 的：",
    paste(setdiff(rownames(ms), ko_ms), collapse = ", "), "\n\n")

## ---- 2. 读 GMT ----
read_gmt <- function(path) {
  lines <- readLines(path)
  sets <- lapply(lines, function(l) {
    f <- strsplit(l, "\t")[[1]]
    if (length(f) < 3) return(character(0))
    unique(f[-c(1, 2)])
  })
  names(sets) <- vapply(lines, function(l) strsplit(l, "\t")[[1]][1], character(1))
  sets
}
gmt <- c(read_gmt("work/gmt/KEGG_2021_Human.gmt"),
         read_gmt("work/gmt/Reactome_2022.gmt"),
         read_gmt("work/gmt/GO_Biological_Process_2021.gmt"),
         read_gmt("work/gmt/MSigDB_Hallmark_2020.gmt"))
cat("基因集总数：", length(gmt), "\n")
background <- colnames(ctrl)
cat("背景基因数（网络基因）：", length(background), "\n\n")

## ---- 3. 超几何富集 ----
enrich_one <- function(vec, ko) {
  v <- vec[names(vec) != ko]
  v <- sort(v, decreasing = TRUE)
  target <- intersect(names(v)[seq_len(min(TOP_N, length(v)))], background)
  res <- lapply(names(gmt), function(term) {
    g <- intersect(gmt[[term]], background)
    if (length(g) < MIN_SET || length(g) > MAX_SET) return(NULL)
    hit <- length(intersect(target, g))
    if (hit < 2) return(NULL)
    exp_hit <- length(target) * length(g) / length(background)
    p <- stats::phyper(hit - 1, length(g), length(background) - length(g),
                       length(target), lower.tail = FALSE)
    data.frame(ko_gene = ko, term = term, n_set = length(g), n_hit = hit,
               expected = round(exp_hit, 1), fold = round(hit / exp_hit, 2), p = p)
  })
  out <- do.call(rbind, res)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out$p.adj <- stats::p.adjust(out$p, method = "fdr")
  out[order(out$p), ]
}
run_all <- function(mat, ko_list) do.call(rbind, lapply(ko_list, function(g) enrich_one(mat[g, ], g)))

cat("===== Control 网络富集 =====\n")
ec <- run_all(ctrl, ko_ctrl)
cat("显著结果（FDR<0.05）：", sum(ec$p.adj < 0.05), "条\n\n")
cat("===== MS 网络富集 =====\n")
em <- run_all(ms, ko_ms)
cat("显著结果（FDR<0.05）：", sum(em$p.adj < 0.05), "条\n")
utils::write.csv(ec, file.path(dir, "enrichment_Control_all.csv"), row.names = FALSE)
utils::write.csv(em, file.path(dir, "enrichment_MS_all.csv"), row.names = FALSE)

## ---- 4. 每个敲除基因的显著通路 ----
show_top <- function(e, od, label, k = 3) {
  cat("\n=====", label, "：每个敲除基因的显著通路（Top", k, "）=====\n")
  sig <- e[e$p.adj < 0.05, ]
  if (nrow(sig) == 0) { cat("（无显著通路）\n"); return(invisible(NULL)) }
  for (ko in unique(sig$ko_gene)) {
    x <- head(sig[sig$ko_gene == ko, ], k)
    cat(sprintf("\n[%s] 出度 %d\n", ko, od[ko]))
    for (i in seq_len(nrow(x))) {
      cat(sprintf("   %-58s hit=%d/%d fold=%.1f FDR=%.2g\n",
                  substr(x$term[i], 1, 58), x$n_hit[i], x$n_set[i], x$fold[i], x$p.adj[i]))
    }
  }
}
show_top(ec, od_ctrl, "Control 网络")
show_top(em, od_ms, "MS 网络")

## ---- 5. 跨敲除基因的靶基因重叠 ----
targets_of <- function(mat, ko_list) {
  out <- lapply(ko_list, function(g) {
    v <- mat[g, ]
    v <- v[names(v) != g]
    names(sort(v, decreasing = TRUE))[seq_len(TOP_N)]
  })
  names(out) <- ko_list
  out
}
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
tc_list <- targets_of(ctrl, ko_ctrl)
tm_list <- targets_of(ms, ko_ms)

J <- outer(seq_along(ko_ctrl), seq_along(ko_ctrl),
           Vectorize(function(i, j) jac(tc_list[[i]], tc_list[[j]])))
dimnames(J) <- list(ko_ctrl, ko_ctrl)
utils::write.csv(round(J, 3), file.path(dir, "target_overlap_Control.csv"))

png(file.path(dir, "target_overlap_Control.png"), width = 1400, height = 1300, res = 180)
op <- graphics::par(mar = c(8, 8, 3, 1))
graphics::image(1:nrow(J), 1:ncol(J), t(J)[, ncol(J):1],
                col = grDevices::hcl.colors(50, "YlOrRd", rev = TRUE), axes = FALSE,
                xlab = "", ylab = "", main = "靶基因集重叠（Jaccard，Control 网络）")
graphics::axis(1, at = 1:nrow(J), labels = rownames(J), las = 2, cex.axis = 0.8)
graphics::axis(2, at = 1:ncol(J), labels = rev(colnames(J)), las = 2, cex.axis = 0.8)
graphics::box()
graphics::par(op)
dev.off()

cat("\n\n===== Control 网络：靶基因重叠最高的基因对 =====\n")
ut <- which(upper.tri(J), arr.ind = TRUE)
pairs_df <- data.frame(g1 = rownames(J)[ut[, 1]], g2 = colnames(J)[ut[, 2]], jaccard = J[ut])
pairs_df <- pairs_df[order(-pairs_df$jaccard), ]
print(head(pairs_df, 12), row.names = FALSE)

cat("\n===== 同一基因在 Control 与 MS 网络中的靶基因重叠 =====\n")
common <- intersect(ko_ctrl, ko_ms)
comp <- do.call(rbind, lapply(common, function(g) {
  data.frame(ko_gene = g, jaccard = jac(tc_list[[g]], tm_list[[g]]),
             outdegree_control = od_ctrl[g], outdegree_MS = od_ms[g])
}))
comp <- comp[order(-comp$jaccard), ]
print(comp, row.names = FALSE)
utils::write.csv(comp, file.path(dir, "target_overlap_control_vs_MS.csv"), row.names = FALSE)

cat("\n完成。输出：enrichment_*.csv、target_overlap_*.csv/png\n")
