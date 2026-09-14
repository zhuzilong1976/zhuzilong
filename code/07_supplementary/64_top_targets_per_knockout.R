# 查看特异性归一化后，各敲除基因的 top 靶基因
suppressPackageStartupMessages({ library(Matrix) })
dir <- "outputs/results/vko_ms_microglia"

read_mat <- function(tag) as.matrix(utils::read.csv(file.path(dir, paste0(tag, "_distanceMatrix.csv")),
                                                    row.names = 1, check.names = FALSE))
read_od <- function(tag) {
  o <- utils::read.csv(file.path(dir, paste0(tag, "_outdegree.csv")), stringsAsFactors = FALSE)
  stats::setNames(o$outdegree, o$gene)
}
ctrl <- read_mat("MSvC_Control"); ms <- read_mat("MSvC_MS")
od_ctrl <- read_od("MSvC_Control"); od_ms <- read_od("MSvC_MS")

spec_z <- function(mat, ko_list) {
  M <- mat[ko_list, , drop = FALSE]
  mu <- colMeans(M); s <- apply(M, 2, stats::sd); s[s == 0] <- NA
  Z <- sweep(sweep(M, 2, mu, "-"), 2, s, "/"); Z[!is.finite(Z)] <- 0; Z
}
Zc <- spec_z(ctrl, rownames(ctrl)[od_ctrl[rownames(ctrl)] > 0])
Zm <- spec_z(ms, rownames(ms)[od_ms[rownames(ms)] > 0])

cat("===== 每个候选基因的 Top 12 特异性靶基因（z 值排名）=====\n")
panel <- c("TREM2","TYROBP","CD33","APOE","C1QA","C3","CSF1R","BTK","NFKB1","SPI1","P2RY12","CX3CR1")
for (g in panel) {
  for (nm in c("Control", "MS")) {
    Z <- if (nm == "Control") Zc else Zm
    if (!g %in% rownames(Z)) { cat(sprintf("\n[%s | %s] 不在该网络\n", g, nm)); next }
    v <- Z[g, ]; v <- v[names(v) != g]
    top <- sort(v, decreasing = TRUE)[1:12]
    cat(sprintf("\n[%s | %s 网络] 出度 %d\n   %s\n", g, nm,
                if (nm == "Control") od_ctrl[g] else od_ms[g],
                paste(sprintf("%s(%.1f)", names(top), top), collapse = "  ")))
  }
}

cat("\n\n===== 反向检查：阴性对照基因的 Top 靶标 =====\n")
for (g in c("ACTA2","TSHR","PECAM1","YKT6")) {
  for (nm in c("Control","MS")) {
    Z <- if (nm == "Control") Zc else Zm
    if (!g %in% rownames(Z)) next
    v <- Z[g, ]; v <- v[names(v) != g]
    top <- sort(v, decreasing = TRUE)[1:8]
    cat(sprintf("[%s | %s] %s\n", g, nm, paste(names(top), collapse = ", ")))
  }
}
