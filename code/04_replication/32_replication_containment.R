# 跨数据集复现：在作者发表的 10 个独立网络上检验包含关系
# 命题：差异调控基因 DR ⊆ {被敲除基因} ∪ {其直接出边靶标}
suppressPackageStartupMessages({ library(Matrix) })

files <- list.files("work/author_networks", pattern = "\\.RData$", full.names = TRUE)
dataset_label <- c(
  "GSM3477499.RData" = "HNF4G / mouse liver",
  "GSM3716703.RData" = "NKX2-1 / human lung",
  "GSM4116571.RData" = "DMD / human muscle",
  "MALAT1.RData" = "MALAT1 / human",
  "Preenterocytes.RData" = "AHR / mouse intestine",
  "SRS3059998.RData" = "MECP2 / mouse brain (rep 1)",
  "SRS3059999.RData" = "MECP2 / mouse brain (rep 2)",
  "SRS3161261.RData" = "CFTR / human lung AT2 (rep 1)",
  "SRS4245406.RData" = "CFTR / human lung AT2 (rep 2)",
  "SRS3161261_Pulmonary alveolar type II cellsAkap7.RData" = "AKAP7 (negative control) / human lung AT2"
)

res <- list()
for (f in files) {
  nm <- basename(f)
  e <- new.env(); load(f, envir = e)
  obj <- get(ls(e)[1], envir = e)
  W <- as.matrix(obj$WT); K <- as.matrix(obj$KO)
  DR <- obj$diffRegulation

  ## 用数据本身确定被敲除基因：WT 该行有边、KO 该行全为 0
  od_wt <- Matrix::rowSums(W != 0)
  od_ko <- Matrix::rowSums(K != 0)
  ko_genes <- intersect(rownames(W)[od_wt > 0 & od_ko == 0], rownames(W))

  ## 直接靶标：被敲除基因在 WT 中的非零列
  targets <- unique(unlist(lapply(ko_genes, function(g) colnames(W)[W[g, ] != 0])))
  allowed <- union(ko_genes, targets)

  dr <- DR$gene[DR$p.adj < 0.05]
  dr_in <- intersect(dr, rownames(W))
  inside <- dr_in %in% allowed

  res[[nm]] <- data.frame(
    file = nm,
    dataset = dataset_label[nm],
    ko_genes = paste(ko_genes, collapse = "+"),
    network_genes = nrow(W),
    ko_outdegree = paste(od_wt[ko_genes], collapse = "+"),
    n_DR = length(dr_in),
    n_DR_inside = sum(inside),
    pct_inside = round(100 * mean(inside), 1),
    bound_ok = length(dr_in) <= sum(od_wt[ko_genes]) + length(ko_genes),
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, res)
out <- out[order(-out$pct_inside), ]
cat("\n================ 跨数据集包含关系检验 ================\n")
print(out[, c("dataset", "ko_genes", "network_genes", "ko_outdegree", "n_DR",
              "n_DR_inside", "pct_inside")], row.names = FALSE)

cat(sprintf("\n汇总：%d 个独立网络；DR ⊆ {KO} ∪ 直接靶标 完全成立的有 %d 个（%.0f%%）\n",
            nrow(out), sum(out$pct_inside == 100), 100 * mean(out$pct_inside == 100)))
cat(sprintf("所有网络合计：%d 个 DR 基因，其中 %d 个落在'KO ∪ 直接靶标'内（%.1f%%）\n",
            sum(out$n_DR), sum(out$n_DR_inside), 100 * sum(out$n_DR_inside) / sum(out$n_DR)))
cat(sprintf("|DR| ≤ 出度+1 成立：%d / %d\n", sum(out$bound_ok), nrow(out)))

dir.create("outputs/results/replication", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, "outputs/results/replication/containment_across_datasets.csv",
                 row.names = FALSE)
