cd /scratch/anneduffels/proteomics/data

cat << 'EOF' > correlate_histones_stat5a_expr.R
library(tidyverse)

merged <- read.csv("prot_rna_merged.csv", stringsAsFactors = FALSE) %>%
  distinct(name, .keep_all = TRUE)

tss <- read.delim(
  "/scratch/Omics_in_practice/data/TSS_genes/hg38_TSS_v50.bed",
  header = FALSE, stringsAsFactors = FALSE
)
colnames(tss)[1:4] <- c("chr", "tss_start", "tss_end", "gene_name")

tss_sub <- tss %>%
  filter(gene_name %in% merged$name) %>%
  distinct(gene_name, .keep_all = TRUE)

read_narrow <- function(path, prefix) {
  df <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df)[1:7] <- c("chr", "start", "end", "peak_id", "score", "strand", "signalValue")
  df %>% mutate(mid = (start + end) / 2, mark = prefix)
}

h27ac <- read_narrow("/scratch/anneduffels/chipseq_project/module_b/H3K27ac.bed", "H3K27ac")
h4me3 <- read_narrow("/scratch/anneduffels/chipseq_project/module_b/H3K4me3.bed", "H3K4me3")
stat5a <- read.delim(
  "/scratch/Omics_in_practice/data/TFOI_beds_9_column/STAT5A_ENCFF686RKS.bed",
  header = FALSE, skip = 1, stringsAsFactors = FALSE
)
colnames(stat5a)[1:7] <- c("chr", "start", "end", "peak_id", "score", "strand", "signalValue")
stat5a <- stat5a %>% mutate(mid = (start + end) / 2, mark = "STAT5A")

tss_sig <- function(peaks_df, prefix, win = 5000) {
  peaks_df %>%
    inner_join(tss_sub, by = "chr", relationship = "many-to-many") %>%
    mutate(dist = abs(mid - tss_start)) %>%
    filter(dist <= win) %>%
    group_by(gene_name) %>%
    summarise(
      !!paste0(prefix, "_peaks") := n(),
      !!paste0(prefix, "_max_sig") := max(signalValue, na.rm = TRUE),
      .groups = "drop"
    )
}

s_27ac <- tss_sig(h27ac, "H3K27ac", win = 5000)
s_4me3 <- tss_sig(h4me3, "H3K4me3", win = 5000)
s_5a   <- tss_sig(stat5a, "STAT5A", win = 5000)

df_all <- merged %>%
  left_join(s_27ac, by = c("name" = "gene_name")) %>%
  left_join(s_4me3, by = c("name" = "gene_name")) %>%
  left_join(s_5a,   by = c("name" = "gene_name")) %>%
  mutate(
    H3K27ac_log10sig = log10(replace_na(H3K27ac_max_sig, 0) + 1),
    H3K4me3_log10sig = log10(replace_na(H3K4me3_max_sig, 0) + 1),
    STAT5A_log10sig  = log10(replace_na(STAT5A_max_sig, 0) + 1),
    has_STAT5A       = !is.na(STAT5A_max_sig)
  )

write.csv(df_all, "histone_stat5a_expr_merged.csv", row.names = FALSE)

# Kwantitatieve correlaties (Spearman en Pearson)
corr_pairs <- list(
  c("H3K27ac_log10sig", "log2_tpm1", "H3K27ac signal vs RNA log2(TPM+1)"),
  c("H3K4me3_log10sig", "log2_tpm1", "H3K4me3 signal vs RNA log2(TPM+1)"),
  c("STAT5A_log10sig",  "log2_tpm1", "STAT5A signal vs RNA log2(TPM+1)"),
  c("H3K27ac_log10sig", "STAT5A_log10sig", "H3K27ac signal vs STAT5A signal"),
  c("H3K4me3_log10sig", "STAT5A_log10sig", "H3K4me3 signal vs STAT5A signal"),
  c("H3K27ac_log10sig", "DMSO_drug_vs_DMSO_ratio", "H3K27ac signal vs Prot log2FC")
)

for (pr in corr_pairs) {
  x <- df_all[[pr[1]]]
  y <- df_all[[pr[2]]]
  sp <- cor.test(x, y, method = "spearman", exact = FALSE)
  pe <- cor.test(x, y, method = "pearson")
  cat(sprintf("%-42s | Spearman rho = %+.4f (p = %.2e) | Pearson r = %+.4f (p = %.2e)\n",
              pr[3], sp$estimate, sp$p.value, pe$estimate, pe$p.value))
}

sub_27ac <- df_all %>% filter(H3K27ac_log10sig > 0)
sub_5a   <- df_all %>% filter(STAT5A_log10sig > 0)
sp_27ac_rna <- cor.test(sub_27ac$H3K27ac_log10sig, sub_27ac$log2_tpm1, method = "spearman", exact = FALSE)
sp_5a_27ac  <- cor.test(sub_5a$STAT5A_log10sig, sub_5a$H3K27ac_log10sig, method = "spearman", exact = FALSE)

p1 <- ggplot(sub_27ac, aes(x = H3K27ac_log10sig, y = log2_tpm1)) +
  geom_point(alpha = 0.25, size = 1, color = "#1b9e77") +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  theme_bw(base_size = 11) +
  labs(
    title = "H3K27ac TSS signal vs. RNA expression",
    subtitle = sprintf("n=%d | Spearman rho=%.3f (p=%.2e)", nrow(sub_27ac), sp_27ac_rna$estimate, sp_27ac_rna$p.value),
    x = "log10(H3K27ac max signal + 1)",
    y = "RNA log2(TPM + 1)"
  )

p2 <- ggplot(sub_5a, aes(x = STAT5A_log10sig, y = H3K27ac_log10sig)) +
  geom_point(alpha = 0.4, size = 1.3, color = "#7570b3") +
  geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
  theme_bw(base_size = 11) +
  labs(
    title = "STAT5A TSS signal vs. H3K27ac signal",
    subtitle = sprintf("n=%d | Spearman rho=%.3f (p=%.2e)", nrow(sub_5a), sp_5a_27ac$estimate, sp_5a_27ac$p.value),
    x = "log10(STAT5A max signal + 1)",
    y = "log10(H3K27ac max signal + 1)"
  )

p3 <- ggplot(df_all, aes(x = has_STAT5A, y = H3K27ac_log10sig, fill = has_STAT5A)) +
  geom_boxplot(outlier.alpha = 0.15) +
  scale_x_discrete(labels = c("No STAT5A\n(±5 kb)", "STAT5A peak\n(±5 kb)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "H3K27ac signal by STAT5A TSS binding",
    x = NULL,
    y = "log10(H3K27ac max signal + 1)"
  )

p4 <- ggplot(df_all, aes(x = has_STAT5A, y = log2_tpm1, fill = has_STAT5A)) +
  geom_boxplot(outlier.alpha = 0.15) +
  scale_x_discrete(labels = c("No STAT5A\n(±5 kb)", "STAT5A peak\n(±5 kb)")) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "RNA expression by STAT5A TSS binding",
    x = NULL,
    y = "RNA log2(TPM + 1)"
  )

pdf("histone_stat5a_expression_correlation.pdf", width = 10, height = 7.5)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
cat("Opgeslagen: histone_stat5a_expr_merged.csv en histone_stat5a_expression_correlation.pdf\n")
EOF

Rscript correlate_histones_stat5a_expr.R
