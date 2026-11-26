cd /scratch/anneduffels/proteomics/data

cat << 'EOF' > correlate_stat5a_tss_expr.R
library(tidyverse)

merged <- read.csv("prot_rna_merged.csv", stringsAsFactors = FALSE)
tss <- read.delim(
  "/scratch/Omics_in_practice/data/TSS_genes/hg38_TSS_v50.bed",
  header = FALSE, stringsAsFactors = FALSE
)
colnames(tss)[1:4] <- c("chr", "tss_start", "tss_end", "gene_name")

peaks <- read.delim(
  "/scratch/Omics_in_practice/data/TFOI_beds_9_column/STAT5A_ENCFF686RKS.bed",
  header = FALSE, skip = 1, stringsAsFactors = FALSE
)
colnames(peaks)[1:7] <- c("chr", "peak_start", "peak_end", "name_col", "score", "strand", "signalValue")
peaks <- peaks %>% mutate(peak_mid = (peak_start + peak_end) / 2)

# Koppel elke STAT5A-piek aan het dichtstbijzijnde TSS op hetzelfde chromosoom
tss_sub <- tss %>%
  filter(gene_name %in% merged$name) %>%
  distinct(gene_name, .keep_all = TRUE)

peak_nearest <- peaks %>%
  inner_join(tss_sub, by = "chr", relationship = "many-to-many") %>%
  mutate(dist_to_tss = abs(peak_mid - tss_start)) %>%
  group_by(chr, peak_start, peak_end) %>%
  slice_min(order_by = dist_to_tss, n = 1, with_ties = FALSE) %>%
  ungroup()

# Genen met een STAT5A-piek binnen 5000 bp van TSS
proximal_genes <- peak_nearest %>%
  filter(dist_to_tss <= 5000) %>%
  group_by(gene_name) %>%
  summarise(
    n_peaks = n(),
    max_signal = max(signalValue, na.rm = TRUE),
    min_dist = min(dist_to_tss, na.rm = TRUE),
    .groups = "drop"
  )

df_tfoi <- merged %>%
  mutate(
    has_stat5a_tss_5kb = name %in% proximal_genes$gene_name
  ) %>%
  left_join(proximal_genes, by = c("name" = "gene_name"))

cat("Genen met STAT5A-piek binnen 5 kb van TSS:", sum(df_tfoi$has_stat5a_tss_5kb), "\n")
print(
  df_tfoi %>%
    group_by(has_stat5a_tss_5kb) %>%
    summarise(
      n = n(),
      median_prot_log2FC = median(DMSO_drug_vs_DMSO_ratio, na.rm = TRUE),
      median_log2_tpm = median(log2_tpm1, na.rm = TRUE),
      .groups = "drop"
    )
)

wt_prot <- wilcox.test(DMSO_drug_vs_DMSO_ratio ~ has_stat5a_tss_5kb, data = df_tfoi)
wt_rna  <- wilcox.test(log2_tpm1 ~ has_stat5a_tss_5kb, data = df_tfoi)
cat(sprintf("Wilcoxon prot log2FC (STAT5A TSS +-5kb vs overig): p = %.3e\n", wt_prot$p.value))
cat(sprintf("Wilcoxon log2(TPM+1) (STAT5A TSS +-5kb vs overig): p = %.3e\n", wt_rna$p.value))

# Correlatie tussen STAT5A peak signalValue en expressie voor nabije TSS-genen
near_sig <- df_tfoi %>% filter(!is.na(max_signal))
if (nrow(near_sig) >= 5) {
  ct_rna <- cor.test(log10(near_sig$max_signal + 1), near_sig$log2_tpm1, method = "spearman")
  ct_prot <- cor.test(log10(near_sig$max_signal + 1), near_sig$DMSO_drug_vs_DMSO_ratio, method = "spearman")
  cat(sprintf("Spearman STAT5A signal vs log2(TPM+1) (n=%d): rho = %.4f (p = %.3e)\n",
              nrow(near_sig), ct_rna$estimate, ct_rna$p.value))
  cat(sprintf("Spearman STAT5A signal vs prot log2FC (n=%d): rho = %.4f (p = %.3e)\n",
              nrow(near_sig), ct_prot$estimate, ct_prot$p.value))
}

p <- ggplot(df_tfoi, aes(x = has_stat5a_tss_5kb, y = DMSO_drug_vs_DMSO_ratio, fill = has_stat5a_tss_5kb)) +
  geom_boxplot(outlier.alpha = 0.25) +
  coord_cartesian(ylim = c(-4, 4)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Proteomics fold-change by STAT5A TSS proximity (+-5 kb)",
    subtitle = sprintf("Wilcoxon p = %.2e | STAT5A protein log2FC = -0.94 (adj. p = 7.74e-07)", wt_prot$p.value),
    x = "STAT5A peak within +-5 kb of TSS",
    y = "Proteomics log2FC (DMSO_drug vs. DMSO)"
  ) +
  theme(legend.position = "none")

ggsave("stat5a_tss_expression_correlation.pdf", p, width = 7, height = 5)
cat("Opgeslagen: stat5a_tss_expression_correlation.pdf\n")
EOF

Rscript correlate_stat5a_tss_expr.R
