#!/usr/bin/env Rscript
# DESeq2 차등발현 분석: rel_mutant vs WT
# 출력: 정규화 count, DE 결과표, PCA, MA, volcano, 상위 유전자 heatmap, (선택) KEGG 농축

suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap)
})

res_dir    <- Sys.getenv("RESULTS_DIR", "results")
count_file <- file.path(res_dir, "counts", "counts.txt")
sheet_file <- "samplesheet.csv"
out_dir    <- file.path(res_dir, "deseq2")
cat("결과 폴더:", res_dir, "\n")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ref_level  <- "WT"          # 기준 그룹
padj_cut   <- 0.05
lfc_cut    <- 1              # |log2FC| >= 1 (2배 이상 변화)

# ---- 1. 데이터 읽기 ----
fc <- read.delim(count_file, comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("\\.bam$", "", basename(colnames(counts)))

sheet <- read.csv(sheet_file, stringsAsFactors = FALSE)
stopifnot(!any(sheet$condition == "UNKNOWN"))
rownames(sheet) <- sheet$run
sheet <- sheet[colnames(counts), , drop = FALSE]
sheet$condition <- relevel(factor(sheet$condition), ref = ref_level)
# 시점(dawn/dusk) 열이 있으면 공변량으로 보정: ~ time + condition
has_time <- "time" %in% colnames(sheet)
if (has_time) {
  stopifnot(!any(sheet$time == "UNKNOWN"))
  sheet$time <- relevel(factor(sheet$time), ref = "dawn")
  design_formula <- ~ time + condition
} else {
  design_formula <- ~ condition
}
cat("샘플 수:", ncol(counts), "\n"); print(table(sheet$condition, if (has_time) sheet$time else "all"))

# ---- 2. DESeq2 ----
dds <- DESeqDataSetFromMatrix(counts, colData = sheet, design = design_formula)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]          # 저발현 유전자 제거
dds <- DESeq(dds)
coef_name <- paste0("condition_", levels(sheet$condition)[2], "_vs_", ref_level)
stopifnot(coef_name %in% resultsNames(dds))
cat("검정 계수:", coef_name, "\n")
res <- results(dds, name = coef_name, alpha = padj_cut)
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")

res_df <- as.data.frame(res)
res_df$log2FC_shrunk <- res_shr$log2FoldChange
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file.path(out_dir, "DE_results.csv"), row.names = FALSE)

sig <- subset(res_df, !is.na(padj) & padj < padj_cut & abs(log2FoldChange) >= lfc_cut)
cat(sprintf("유의 유전자: %d개 (up %d / down %d), 검정된 유전자 %d개 중\n",
            nrow(sig), sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0),
            sum(!is.na(res_df$padj))))
write.csv(sig, file.path(out_dir, "DE_significant.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
write.csv(counts(dds, normalized = TRUE), file.path(out_dir, "normalized_counts.csv"))

# ---- 3. 그림 ----
pca_groups <- if (has_time) c("condition", "time") else "condition"
pca_df <- plotPCA(vsd, intgroup = pca_groups, returnData = TRUE)
pv <- round(100 * attr(pca_df, "percentVar"))
p <- ggplot(pca_df, aes(PC1, PC2, colour = condition,
                        shape = if (has_time) time else condition)) +
  geom_point(size = 3) + labs(shape = if (has_time) "time" else "condition",
  x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)")) +
  theme_bw() + ggtitle("PCA (VST)")
ggsave(file.path(out_dir, "PCA.png"), p, width = 5, height = 4, dpi = 200)

png(file.path(out_dir, "MA.png"), width = 1000, height = 800, res = 150)
plotMA(res_shr, ylim = c(-5, 5), main = "MA plot (apeglm shrunk)")
dev.off()

res_df$status <- "NS"
res_df$status[res_df$gene %in% sig$gene & res_df$log2FoldChange > 0] <- "Up"
res_df$status[res_df$gene %in% sig$gene & res_df$log2FoldChange < 0] <- "Down"
v <- ggplot(subset(res_df, !is.na(padj)),
            aes(log2FoldChange, -log10(padj), colour = status)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_colour_manual(values = c(Down = "#2b6cb0", NS = "grey70", Up = "#c53030")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2) +
  geom_hline(yintercept = -log10(padj_cut), linetype = 2) +
  theme_bw() + ggtitle(coef_name)
ggsave(file.path(out_dir, "volcano.png"), v, width = 6, height = 5, dpi = 200)

top <- head(sig$gene, 40)
if (length(top) >= 2) {
  pheatmap(assay(vsd)[top, ], scale = "row",
           annotation_col = sheet[, pca_groups, drop = FALSE],
           filename = file.path(out_dir, "heatmap_top40.png"), width = 6, height = 8)
}

# ---- 4. (선택) KEGG 농축 분석 ----
# KEGG organism code 'syf' = S. elongatus PCC 7942 로 알려져 있음 — 실행 전 KEGG에서 확인.
# 인터넷 연결 필요. 실패해도 위 결과는 이미 저장됨.
tryCatch({
  suppressPackageStartupMessages(library(clusterProfiler))
  universe <- res_df$gene[!is.na(res_df$padj)]
  for (dir in c("Up", "Down")) {
    genes <- res_df$gene[res_df$status == dir]
    if (length(genes) < 5) next
    ek <- enrichKEGG(gene = genes, organism = "syf", keyType = "kegg",
                     universe = universe, pvalueCutoff = 0.05)
    write.csv(as.data.frame(ek), file.path(out_dir, paste0("KEGG_", dir, ".csv")),
              row.names = FALSE)
    cat(sprintf("KEGG %s: 유의 경로 %d개\n", dir, nrow(as.data.frame(ek))))
  }
}, error = function(e) message("[skip] KEGG 농축 분석 실패: ", conditionMessage(e)))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("[OK] 결과:", out_dir, "\n")
