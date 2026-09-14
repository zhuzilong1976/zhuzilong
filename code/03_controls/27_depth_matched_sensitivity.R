# 深度匹配敏感性分析：比较"深度匹配后"与"主分析"的两条件比较结果
suppressPackageStartupMessages({ library(Matrix); library(scTenifoldKnk) })
dm <- "outputs/results/vko_ms_depthmatched"
main <- "outputs/results/vko_ms_microglia"

run_comp <- function(A, B, label) {
  g <- intersect(rownames(A), rownames(B))
  A <- A[g, g]; B <- B[g, g]
  set.seed(1)
  MA <- scTenifoldNet::manifoldAlignment(A, B, d = 2)
  DR <- scTenifoldKnk::dRegulation(MA)
  DR <- DR[order(DR$p.value), ]
  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("基因数 %d；p<0.05：%d（%.1f%%）；FDR<0.05：%d\n",
              nrow(DR), sum(DR$p.value < 0.05), 100 * mean(DR$p.value < 0.05),
              sum(DR$p.adj < 0.05)))
  cat("Top 12：", paste(DR$gene[1:12], collapse = ", "), "\n")
  a <- as.matrix(A); b <- as.matrix(B)
  cat(sprintf("两网络相关性（Spearman）：%.3f\n",
              stats::cor(as.vector(a), as.vector(b), method = "spearman")))
  DR
}

rc <- readRDS(file.path(dm, "MSvC_Control_depthmatched_raw_result.rds"))
rm_ <- readRDS(file.path(dm, "MSvC_MS_depthmatched_raw_result.rds"))
DRdm <- run_comp(as.matrix(rc$tensorNetworks$WT), as.matrix(rm_$tensorNetworks$WT),
                 "深度匹配（800-6000 counts）")

obs <- readRDS(file.path(main, "condition_comparison.rds"))
DRmain <- run_comp(obs$A, obs$B, "主分析（全部细胞）")

meta <- utils::read.csv("work/data/metadata.csv", stringsAsFactors = FALSE)
X <- Matrix::readMM("work/data/counts.mtx")
lib <- Matrix::colSums(X)
names(lib) <- readLines("work/data/barcodes.txt")
ct <- meta$barcode[meta$condition == "Control"]
ms <- meta$barcode[meta$condition == "MS"]
inwin <- function(cells) cells[lib[cells] >= 800 & lib[cells] <= 6000]
cat("\n===== 深度窗口内的细胞（QC 前）=====\n")
cat(sprintf("Control: %d 个，中位 %s counts\n", length(inwin(ct)),
            format(round(median(lib[inwin(ct)])), big.mark = ",")))
cat(sprintf("MS:      %d 个，中位 %s counts\n", length(inwin(ms)),
            format(round(median(lib[inwin(ms)])), big.mark = ",")))
cat(sprintf("深度比（MS/Control 中位）：%.2f（主分析为 %.2f）\n",
            median(lib[inwin(ms)]) / median(lib[inwin(ct)]),
            median(lib[ms]) / median(lib[ct])))

panel <- c("TREM2","TYROBP","CD33","APOE","BIN1","PICALM","CR1","MS4A4A","SPI1","IRF8",
           "STAT3","NFKB1","C1QA","C3","CX3CR1","P2RY12","TMEM119","BTK","CSF1R","IL10RA")
cmp <- data.frame(gene = panel,
                  p_matched = DRdm$p.value[match(panel, DRdm$gene)],
                  p_main = DRmain$p.value[match(panel, DRmain$gene)])
cmp$rank_matched <- seq_len(nrow(DRdm))[match(panel, DRdm$gene)]
cmp$rank_main <- seq_len(nrow(DRmain))[match(panel, DRmain$gene)]
cmp <- cmp[order(cmp$p_matched), ]
cat("\n===== 面板基因在两次分析中的表现 =====\n")
print(cmp, row.names = FALSE, digits = 3)
cat(sprintf("\n面板基因 p 值中位数：深度匹配 %.3f ｜ 主分析 %.3f\n",
            median(cmp$p_matched, na.rm = TRUE), median(cmp$p_main, na.rm = TRUE)))

utils::write.csv(cmp, file.path(dm, "panel_depthmatched_vs_main.csv"), row.names = FALSE)
