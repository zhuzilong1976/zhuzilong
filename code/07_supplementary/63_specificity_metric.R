# 与出度无关的指标：每个敲除基因的"特异性扰动强度"
#   max_z    = 该基因敲除后，靶基因中最大的特异性 z 值
#   n_z2     = z > 2 的靶基因数
#   frac_pos = z > 0 的靶基因比例（衡量该敲除是否整体强于平均）

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
neg <- c("ACTA2", "TSHR", "PECAM1", "YKT6")

metrics <- function(mat, od) {
  ko_list <- rownames(mat)[od[rownames(mat)] > 0]
  Z <- spec_z(mat, ko_list)
  do.call(rbind, lapply(rownames(mat), function(g) {
    isneg <- g %in% neg
    if (!g %in% ko_list) {
      return(data.frame(ko_gene = g, outdegree = od[g], max_z = NA, n_z2 = NA,
                        frac_pos = NA, is_negative_control = isneg,
                        note = "出度为 0，空操作"))
    }
    v <- Z[g, ]; v <- v[names(v) != g]
    data.frame(ko_gene = g, outdegree = od[g],
               max_z = round(max(v), 2), n_z2 = sum(v > 2),
               frac_pos = round(mean(v > 0), 3), is_negative_control = isneg, note = "")
  }))
}

for (nm in c("Control", "MS")) {
  mat <- if (nm == "Control") ctrl else ms
  od <- if (nm == "Control") od_ctrl else od_ms
  m <- metrics(mat, od)
  m <- m[order(-m$max_z), ]
  cat(sprintf("\n===== %s 网络：特异性扰动强度（按 max_z 排序）=====\n", nm))
  print(m, row.names = FALSE)
  utils::write.csv(m, file.path(dir, paste0("specificity_metric_", nm, ".csv")), row.names = FALSE)

  a <- m[!m$is_negative_control & !is.na(m$max_z), ]
  b <- m[m$is_negative_control & !is.na(m$max_z), ]
  cat(sprintf("\n候选基因 max_z：中位 %.2f，范围 %.2f–%.2f\n",
              median(a$max_z), min(a$max_z), max(a$max_z)))
  cat(sprintf("阴性对照 max_z：%s\n", paste(sprintf("%s=%.2f", b$ko_gene, b$max_z), collapse = ", ")))
}
