cd /scratch/anneduffels/proteomics/data

cat << 'EOF' > fix_tss_plot.R
library(tidyverse)

merged <- read.csv("prot_rna_merged.csv", stringsAsFactors = FALSE)
tss <- read.delim("/scratch/Omics_in_practice/data/TSS_genes/hg38_TSS_v50.bed", header = FALSE, stringsAsFactors = FALSE)
colnames(tss)[1:4] <- c("chr", "tss_start", "tss_end", "gene_name")

peaks <- read.delim("/scratch/Omics_in_practice/data/TFOI_beds_9_column/STAT5A_ENCFF686RKS.bed", header = FALSE, skip = 1, stringsAsFactors = FALSE)
colnames(peaks)[1:7] <- c("chr", "peak_start", "peak_end", "name_col", "score", "strand", "signalValue")
peaks <- peaks %>% mutate(peak_mid = (peak_start + peak_end) / 2)

tss_sub <- tss %>% filter(gene_name %in% merged$name) %>% distinct(gene_name, .keep_all = TRUE)

proximal_genes <- peaks %>%
  inner_join(tss_sub, by = "chr", relationship = "many-to-many") %>%
  mutate(dist_to_tss = abs(peak_mid - tss_start)) %>%
  group_by(chr, peak_start, peak_end) %>%
  slice_min(order_by = dist_to_tss, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(dist_to_tss <= 5000) %>%
  pull(gene_name) %>%
  unique()

df_tfoi <- merged %>%
  mutate(group = ifelse(name %in% proximal_genes, "STAT5A peak (±5 kb)", "No STAT5A peak (±5 kb)"))

counts_grp <- df_tfoi %>% count(group)
wt_prot <- wilcox.test(DMSO_drug_vs_DMSO_ratio ~ group, data = df_tfoi)

p <- ggplot(df_tfoi, aes(x = group, y = DMSO_drug_vs_DMSO_ratio, fill = group)) +
  geom_boxplot(outlier.alpha = 0.25) +
  scale_x_discrete(labels = setNames(paste0(counts_grp$group, "\n(n = ", counts_grp$n, ")"), counts_grp$group)) +
  scale_fill_manual(values = c("No STAT5A peak (±5 kb)" = "#f8766d", "STAT5A peak (±5 kb)" = "#00bfc4")) +
  coord_cartesian(ylim = c(-4, 4)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Proteomics fold-change by STAT5A TSS proximity (±5 kb)",
    subtitle = sprintf("Wilcoxon p = %.2e | STAT5A protein log2FC = -0.94 (adj. p = 7.74e-07)", wt_prot$p.value),
    x = NULL,
    y = "Proteomics log2FC (DMSO_drug vs. DMSO)",
    fill = "TSS status"
  )

ggsave("stat5a_tss_expression_correlation.pdf", p, width = 7.5, height = 5)
cat("Opgeslagen: stat5a_tss_expression_correlation.pdf\n")
EOF

Rscript fix_tss_plot.R
