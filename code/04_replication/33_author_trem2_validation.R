# 外部验证：用作者自己发表的 TREM2 结果检验"显著基因数 = 出度 + 1"
suppressPackageStartupMessages({ library(Matrix) })

env <- new.env()
load("work/author/GSE130626.RData", envir = env)
cat("对象：", paste(ls(env), collapse = ", "), "\n")
obj <- get(ls(env)[1], envir = env)
cat("类型：", class(obj), "；包含：", paste(names(obj), collapse = ", "), "\n\n")

DR <- obj$diffRegulation
cat("diffRegulation：", nrow(DR), "个基因\n")
cat("列：", paste(colnames(DR), collapse = ", "), "\n")
n_sig <- sum(DR$p.adj < 0.05)
cat(sprintf("p.adj<0.05 的基因数 = %d（%.2f%%）\n", n_sig, 100 * n_sig / nrow(DR)))
cat(sprintf("其中距离最大的是：%s\n\n", DR$gene[which.max(DR$distance)]))

WT <- if (!is.null(obj$tensorNetworks$WT)) obj$tensorNetworks$WT else obj$WT
WT <- as.matrix(WT)
cat("网络规模：", nrow(WT), "x", ncol(WT), "\n")
cat("非零边：", format(sum(WT != 0), big.mark = ","),
    sprintf("（密度 %.1f%%）\n", 100 * mean(WT != 0)))

od <- Matrix::rowSums(WT != 0)
gKO <- "Trem2"
if (!gKO %in% names(od)) gKO <- grep("^Trem2$|^TREM2$", names(od), value = TRUE)[1]
cat(sprintf("\n被敲除基因 %s 的出度 = %d\n", gKO, od[gKO]))
cat(sprintf("显著基因数 = %d\n", n_sig))
cat(sprintf("出度 + 1 = %d\n", od[gKO] + 1))
cat(sprintf("差值 = %d\n", n_sig - (od[gKO] + 1)))

## 显著性基因是否是 Trem2 的直接靶标
targets <- names(od)[WT[gKO, ] != 0]
sig_genes <- DR$gene[DR$p.adj < 0.05]
cat(sprintf("\nTrem2 的靶基因数（出度）= %d\n", length(targets)))
cat(sprintf("显著基因中有 %d 个是 Trem2 的直接靶标（%.1f%%）\n",
            length(intersect(sig_genes, targets)),
            100 * length(intersect(sig_genes, targets)) / length(sig_genes)))
cat(sprintf("靶标中有 %d 个被判定显著（%.1f%%）\n",
            length(intersect(sig_genes, targets)),
            100 * length(intersect(sig_genes, targets)) / length(targets)))

## 出度分布（作者的网络）
cat(sprintf("\n全网出度分布：中位 %d，最小 %d，最大 %d；出度为 0 的基因 %d 个\n",
            median(od), min(od), max(od), sum(od == 0)))

## 关键检验：显著基因数是否等于某些基因的出度+1
cat("\n===== 关键检验 =====\n")
cat(sprintf("若规律成立，应有一个基因的出度 ≈ %d（显著数 - 1）\n", n_sig - 1))
cat(sprintf("Trem2 的出度 = %d，与 %d 的差距 = %d（相对差 %.1f%%）\n",
            od[gKO], n_sig - 1, (n_sig - 1) - od[gKO],
            100 * abs((n_sig - 1) - od[gKO]) / (n_sig - 1)))

saveRDS(list(DR = DR, od = od, n_sig = n_sig, gKO = gKO),
        "work/author/author_validation.rds")
