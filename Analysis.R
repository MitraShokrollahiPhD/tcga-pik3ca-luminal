library(cBioPortalData)
library(MultiAssayExperiment)
library(tidyverse)
library(limma)
library(STRINGdb)
library(survival)
library(survminer)

dir.create("results", showWarnings = FALSE)

# 1. Load data
brca <- cBioDataPack("brca_tcga_pan_can_atlas_2018", ask = FALSE)
mut  <- brca[["mutations"]]
expr <- brca[["mrna_seq_v2_rsem"]]
clin <- as.data.frame(colData(brca))

# 2. Mutation table
mut_gr <- unlist(as(mut, "GRangesList"))
mut_df <- data.frame(
  sample = sub("\\..*$", "", names(mut_gr)),
  gene   = sub("^[^.]*\\.", "", names(mut_gr)),
  vclass = as.character(mcols(mut_gr)$Variant_Classification),
  change = as.character(mcols(mut_gr)$HGVSp_Short)
)

# 3. PIK3CA-mutant tumors (non-silent mutations)
nonsilent <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
               "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
pik3ca_mut_samples <- mut_df %>%
  filter(gene == "PIK3CA", vclass %in% nonsilent) %>%
  pull(sample) %>% unique()

# 4. Clean expression data
expr_mat <- if (is.matrix(expr)) expr else assay(expr)
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", ]
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
keep_samples <- intersect(colnames(expr_mat), colnames(mut))
expr_log <- log2(expr_mat[, keep_samples] + 1)

# 5. Sample table: luminal tumors only
info <- data.frame(sample = keep_samples) %>%
  mutate(pik3ca = ifelse(sample %in% pik3ca_mut_samples, "Mutant", "WT")) %>%
  left_join(dplyr::select(clin, SAMPLE_ID, SUBTYPE),
            by = c("sample" = "SAMPLE_ID")) %>%
  filter(SUBTYPE %in% c("BRCA_LumA", "BRCA_LumB"))
info$pik3ca  <- factor(info$pik3ca, levels = c("WT", "Mutant"))
info$SUBTYPE <- factor(info$SUBTYPE)

# 6. Checks
cat("Luminal tumors:", nrow(info), "\n")
print(table(info$pik3ca, info$SUBTYPE))
cat("Most common PIK3CA changes:\n")
print(head(sort(table(mut_df$change[mut_df$gene == "PIK3CA"]), decreasing = TRUE), 6))

# 7. Differential expression (limma), adjusting for luminal subtype
expr_use <- expr_log[, info$sample]
expr_use <- expr_use[rowMeans(expr_use > 1) >= 0.2, ]

design <- model.matrix(~ pik3ca + SUBTYPE, data = info)
fit <- eBayes(lmFit(expr_use, design), trend = TRUE)
res <- topTable(fit, coef = "pik3caMutant", number = Inf, sort.by = "P")
res$gene <- rownames(res)
write.csv(res, "results/pik3ca_limma_results.csv", row.names = FALSE)

# 8. Volcano plot (|log2FC| > 0.5, since PIK3CA effects are modest)
res$status <- case_when(
  res$adj.P.Val < 0.05 & res$logFC >  0.5 ~ "Up in PIK3CA mutant",
  res$adj.P.Val < 0.05 & res$logFC < -0.5 ~ "Down in PIK3CA mutant",
  TRUE ~ "Not significant"
)
volcano <- ggplot(res, aes(logFC, -log10(adj.P.Val), colour = status)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c("Up in PIK3CA mutant" = "firebrick",
                                 "Down in PIK3CA mutant" = "steelblue",
                                 "Not significant" = "grey70")) +
  labs(title = "PIK3CA mutant vs wild-type (luminal tumors, adjusted for subtype)",
       x = "log2 fold change", y = "-log10 adjusted p-value") +
  theme_minimal()
print(volcano)
ggsave("results/volcano_pik3ca.png", volcano, width = 7, height = 5, dpi = 300)

# 9. Checks
cat("Genes tested:", nrow(res), "\n")
print(table(res$status))
print(head(res[, c("gene", "logFC", "adj.P.Val")], 10))

cat("Genes tested:", nrow(res), "\n")
print(table(res$status))
print(head(res[, c("gene", "logFC", "adj.P.Val")], 15))

# 10. Where do PIK3CA and IRS2 rank?
print(res[res$gene %in% c("PIK3CA", "IRS2"), c("gene", "logFC", "adj.P.Val")])

# 11. STRING network and enrichment (top 100 genes by adjusted p-value)
dir.create("string_cache", showWarnings = FALSE)
string_db <- STRINGdb$new(version = "12.0", species = 9606,
                          score_threshold = 400, network_type = "full",
                          input_directory = "string_cache")
top_genes <- res %>% dplyr::slice_head(n = 100) %>%
  dplyr::select(gene, logFC, adj.P.Val)
mapped <- string_db$map(top_genes, "gene", removeUnmappedRows = TRUE)
cat("Genes mapped to STRING:", nrow(mapped), "of 100\n")

all_mapped <- string_db$map(res[, "gene", drop = FALSE], "gene",
                            removeUnmappedRows = TRUE)
invisible(string_db$set_background(all_mapped$STRING_id))

string_db$plot_network(mapped$STRING_id, add_summary = FALSE)
png("results/string_network_pik3ca.png", width = 2000, height = 2000, res = 300)
string_db$plot_network(mapped$STRING_id, add_summary = FALSE)
dev.off()

enrich <- string_db$get_enrichment(mapped$STRING_id)
enrich <- enrich[order(enrich$fdr), ]
write.csv(enrich, "results/string_enrichment_pik3ca.csv", row.names = FALSE)
print(head(enrich[, c("category", "description", "number_of_genes", "fdr")], 10))

print(res[res$gene %in% c("PIK3CA", "IRS2"), c("gene", "logFC", "adj.P.Val")])
print(head(enrich[, c("category", "description", "number_of_genes", "fdr")], 10))

# 12. Enrichment using all significant genes (adj. p < 0.05, |log2FC| > 0.5)
sig_genes <- res %>%
  dplyr::filter(status != "Not significant") %>%
  dplyr::select(gene, logFC, adj.P.Val)
cat("Significant genes:", nrow(sig_genes), "\n")

mapped_sig <- string_db$map(sig_genes, "gene", removeUnmappedRows = TRUE)
enrich_sig <- string_db$get_enrichment(mapped_sig$STRING_id)
cat("Enriched terms:", nrow(enrich_sig), "\n")

if (nrow(enrich_sig) > 0) {
  enrich_sig <- enrich_sig[order(enrich_sig$fdr), ]
  write.csv(enrich_sig, "results/string_enrichment_pik3ca.csv", row.names = FALSE)
  print(head(enrich_sig[, c("category", "description", "number_of_genes", "fdr")], 10))
}

print(res[res$gene %in% c("PIK3CA", "IRS2"), c("gene", "logFC", "adj.P.Val")])

# 13. Specific pathway terms only
print(table(enrich_sig$category))

pathways <- enrich_sig %>%
  dplyr::filter(category %in% c("Process", "KEGG", "WikiPathways", "RCTM"))
cat("Pathway/process terms:", nrow(pathways), "\n")
print(head(pathways[, c("category", "description", "number_of_genes", "fdr")], 15))

# 14. Survival: PIK3CA mutant vs WT in luminal tumors
surv_df <- clin %>%
  dplyr::select(SAMPLE_ID, OS_STATUS, OS_MONTHS) %>%
  dplyr::inner_join(info, by = c("SAMPLE_ID" = "sample")) %>%
  dplyr::mutate(time  = as.numeric(OS_MONTHS),
                event = ifelse(grepl("^1", OS_STATUS), 1, 0)) %>%
  dplyr::filter(!is.na(time), time > 0)

cat("Tumors in survival analysis:", nrow(surv_df), "\n")
print(table(surv_df$OS_STATUS))

fit_km <- survfit(Surv(time, event) ~ pik3ca, data = surv_df)
km_plot <- ggsurvplot(fit_km, data = surv_df, pval = TRUE, risk.table = TRUE,
                      palette = c("steelblue", "firebrick"),
                      legend.labs = c("WT", "PIK3CA mutant"),
                      xlab = "Months", ylab = "Overall survival")
print(km_plot)
png("results/km_pik3ca.png", width = 2000, height = 1800, res = 300)
print(km_plot)
dev.off()

cox <- coxph(Surv(time, event) ~ pik3ca + SUBTYPE, data = surv_df)
print(summary(cox))

# 15. Record package versions
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")