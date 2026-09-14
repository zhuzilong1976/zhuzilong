# B 方案（修正版）：特异性归一化后的靶基因富集
# 问题：直接用距离排名 → 所有敲除基因的 Top 靶标高度重叠（Jaccard 0.98），因为
#       某些低连接度基因在每次敲除中都排前，属于"共同结构"而非基因特异性效应。
# 修正：对每个靶基因，跨全部敲除基因做 z 标准化（d - mean)/sd，再按 z 排名。

suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"
TOP_N <- 60
MIN_SET <- 4
MAX_SET <- 300

read_mat <- function(tag) {
  as.matrix(utils::read.csv(file.path(dir, paste0(tag, "_distanceMatrix.csv")),
                            row.names = 1, check.names = FALSE))
}
read_od <- function(tag) {
  o <- utils::read.csv(file.path(dir, paste0(tag, "_outdegree.csv")), stringsAsFactors = FALSE)
  stats::setNames(o$outdegree, o$gene)
}

ctrl <- read_mat("MSvC_Control"); ms <- read_mat("MSvC_MS")
od_ctrl <- read_od("MSvC_Control"); od_ms <- read_od("MSvC_MS")
ko_ctrl <- rownames(ctrl)[od_ctrl[rownames(ctrl)] > 0]
ko_ms <- rownames(ms)[od_ms[rownames(ms)] > 0]

read_gmt <- function(path) {
  lines <- readLines(path)
  sets <- lapply(lines, function(l) {
    f <- strsplit(l, "\t")[[1]]
    if (length(f) < 3) return(character(0)); unique(f[-c(1, 2)])
  })
  names(sets) <- vapply(lines, function(l) strsplit(l, "\t")[[1]][1], character(1))
  sets
}
gmt <- c(read_gmt("work/gmt/KEGG_2021_Human.gmt"), read_gmt("work/gmt/Reactome_2022.gmt"),
         read_gmt("work/gmt/GO_Biological_Process_2021.gmt"), read_gmt("work/gmt/MSigDB_Hallmark_2020.gmt"))

## ---- 特异性归一化：对每个靶基因，跨全部敲除基因做 z 标准化 ----
specificity_z <- function(mat, ko_list) {
  M <- mat[ko_list, , drop = FALSE]
  mu <- colMeans(M)
  s <- apply(M, 2, stats::sd)
  s[s == 0] <- NA
  Z <- sweep(sweep(M, 2, mu, "-"), 2, s, "/")
  Z[!is.finite(Z)] <- 0
  Z
}

background <- colnames(ctrl)

enrich_z <- function(Z, ko, gmt_sets) {
  v <- Z[ko, ]; v <- v[names(v) != ko]
  target <- intersect(names(sort(v, decreasing = TRUE))[seq_len(TOP_N)], background)
  res <- lapply(names(gmt_sets), function(term) {
    g <- intersect(gmt_sets[[term]], background)
    if (length(g) < MIN_SET || length(g) > MAX_SET) return(NULL)
    hit <- length(intersect(target, g)); if (hit < 2) return(NULL)
    e <- length(target) * length(g) / length(background)
    data.frame(ko_gene = ko, term = term, n_set = length(g), n_hit = hit,
               fold = round(hit / e, 2),
               p = stats::phyper(hit - 1, length(g), length(background) - length(g),
                                 length(target), lower.tail = FALSE))
  })
  out <- do.call(rbind, res)
  if (is.null(out)) return(NULL)
  out$p.adj <- stats::p.adjust(out$p, method = "fdr")
  out[order(out$p), ]
}

Zc <- specificity_z(ctrl, ko_ctrl)
Zm <- specificity_z(ms, ko_ms)

cat("===== 特异性 z 值分布 =====\n")
cat(sprintf("Control: z>2 的靶标数（每基因中位）= %.0f ；z>2 占比 %.1f%%\n",
            median(apply(Zc, 1, function(x) sum(x > 2))),
            100 * mean(Zc > 2)))
cat(sprintf("MS:      z>2 的靶标数（每基因中位）= %.0f ；z>2 占比 %.1f%%\n\n",
            median(apply(Zm, 1, function(x) sum(x > 2))),
            100 * mean(Zm > 2)))

ec <- do.call(rbind, lapply(ko_ctrl, function(g) enrich_z(Zc, g, gmt)))
em <- do.call(rbind, lapply(ko_ms, function(g) enrich_z(Zm, g, gmt)))
cat("Control 富集显著条目（FDR<0.05）：", sum(ec$p.adj < 0.05), " / ", nrow(ec), "\n")
cat("MS      富集显著条目（FDR<0.05）：", sum(em$p.adj < 0.05), " / ", nrow(em), "\n")
utils::write.csv(ec, file.path(dir, "enrichment_specific_Control.csv"), row.names = FALSE)
utils::write.csv(em, file.path(dir, "enrichment_specific_MS.csv"), row.names = FALSE)

show <- function(e, od, label) {
  cat("\n=====", label, "：特异性靶标的显著通路 =====")
  sig <- e[e$p.adj < 0.05, ]
  if (nrow(sig) == 0) { cat("\n（无 FDR<0.05 条目）\n"); return(invisible()) }
  for (ko in unique(sig$ko_gene)) {
    x <- head(sig[sig$ko_gene == ko, ], 5)
    cat(sprintf("\n[%s] 出度 %d\n", ko, od[ko]))
    for (i in seq_len(nrow(x)))
      cat(sprintf("   %-58s hit=%d/%d fold=%.1f FDR=%.3g\n", substr(x$term[i], 1, 58),
                  x$n_hit[i], x$n_set[i], x$fold[i], x$p.adj[i]))
  }
}
show(ec, od_ctrl, "Control 网络")
show(em, od_ms, "MS 网络")

## ---- 特异性归一化后的靶标重叠 ----
targets_of <- function(Z, ko_list) {
  out <- lapply(ko_list, function(g) {
    v <- Z[g, ]; v <- v[names(v) != g]
    names(sort(v, decreasing = TRUE))[seq_len(TOP_N)]
  })
  names(out) <- ko_list; out
}
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))

for (nm in c("Control", "MS")) {
  Z <- if (nm == "Control") Zc else Zm
  kl <- if (nm == "Control") ko_ctrl else ko_ms
  tl <- targets_of(Z, kl)
  J <- outer(seq_along(kl), seq_along(kl), Vectorize(function(i, j) jac(tl[[i]], tl[[j]])))
  dimnames(J) <- list(kl, kl)
  diag(J) <- NA
  cat(sprintf("\n%s 网络：靶标 Jaccard 中位数 = %.3f（最大 %.3f）\n",
              nm, median(J, na.rm = TRUE), max(J, na.rm = TRUE)))
  utils::write.csv(round(J, 3), file.path(dir, paste0("target_overlap_specific_", nm, ".csv")))
}

cat("\n完成。输出：enrichment_specific_*.csv、target_overlap_specific_*.csv\n")
