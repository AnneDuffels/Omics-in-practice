cat << 'EOF' > /scratch/anneduffels/chipseq_project/correlation_analysis.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 1. Load BRB-seq differential expression table
de_path = "/scratch/anneduffels/brbseq/results/de_K562_8d_vs_0h.tsv"
de_df = pd.read_csv(de_path, sep="\t")
de_df["gene_upper"] = de_df["Gene"].astype(str).str.upper().str.strip()

# 2. Load gene annotation (GTF) to extract TSS coordinates
gtf_path = "/scratch/anneduffels/chipseq_project/genomes/GRCh38/GRCh38.annotation.gtf"
gene_records = []
with open(gtf_path) as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chrom, start, end, strand, attributes = (
            fields[0],
            int(fields[3]),
            int(fields[4]),
            fields[6],
            fields[8],
        )
        gene_name = None
        for attr in attributes.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                gene_name = attr.split('"')[1] if '"' in attr else attr.split()[-1]
                break
        if gene_name:
            tss = start if strand == "+" else end
            gene_records.append((chrom, tss, gene_name.upper().strip()))

genes_df = pd.DataFrame(gene_records, columns=["chr", "tss", "gene_upper"]).drop_duplicates("gene_upper")

# 3. Load MACS2 peaks (switch sample01 to sample02 or sample03 if needed)
peaks_path = "/scratch/anneduffels/chipseq_project/results/macs2/GRCh38-sample01_peaks.xls"
peaks_df = pd.read_csv(peaks_path, sep="\t", comment="#")
peaks_df["peak_mid"] = (peaks_df["start"] + peaks_df["end"]) // 2

# Assign each peak to the nearest gene TSS on the same chromosome
annotated_chunks = []
for chrom, peak_sub in peaks_df.groupby("chr"):
    gene_sub = genes_df[genes_df["chr"] == chrom]
    if gene_sub.empty:
        continue
    gene_positions = gene_sub["tss"].values[:, None]
    peak_positions = peak_sub["peak_mid"].values[None, :]
    distances = np.abs(gene_positions - peak_positions)
    nearest_idx = np.argmin(distances, axis=0)
    min_distances = distances[nearest_idx, np.arange(len(peak_sub))]

    annotated_sub = peak_sub.copy()
    annotated_sub["gene_upper"] = gene_sub["gene_upper"].values[nearest_idx]
    annotated_sub["tss_distance"] = min_distances
    annotated_chunks.append(annotated_sub)

peaks_annotated = pd.concat(annotated_chunks, ignore_index=True)

# Keep peaks within 50 kb of a TSS and take the maximum fold enrichment per gene
nearby_peaks = peaks_annotated[peaks_annotated["tss_distance"] <= 50000]
chip_by_gene = (
    nearby_peaks.groupby("gene_upper", as_index=False)["fold_enrichment"]
    .max()
    .rename(columns={"fold_enrichment": "max_fold_enrichment"})
)

# 4. Merge ChIP enrichment with expression data
merged_df = pd.merge(de_df, chip_by_gene, on="gene_upper", how="inner").dropna(
    subset=["log2FC", "max_fold_enrichment"]
)
merged_df["log2_chip_enrichment"] = np.log2(merged_df["max_fold_enrichment"].clip(lower=1e-3))

# Calculate correlations without scipy
spearman_r = merged_df["log2FC"].corr(merged_df["log2_chip_enrichment"], method="spearman")
pearson_r = merged_df["log2FC"].corr(merged_df["log2_chip_enrichment"], method="pearson")
print(f"Matched genes: {len(merged_df)}")
print(f"Spearman r: {spearman_r:.4f}")
print(f"Pearson r:  {pearson_r:.4f}")

# Inspect STAT5A specifically
stat5a_df = merged_df[merged_df["gene_upper"] == "STAT5A"]
if not stat5a_df.empty:
    print("\nSTAT5A summary:")
    print(stat5a_df[["Gene", "log2FC", "padj", "max_fold_enrichment", "log2_chip_enrichment"]])
else:
    print("\nSTAT5A not found within 50 kb of a peak.")

# 5. Create and save scatter plot with linear trendline
plt.figure(figsize=(8, 6))
plt.scatter(
    merged_df["log2_chip_enrichment"],
    merged_df["log2FC"],
    alpha=0.35,
    s=14,
    label="All matched genes",
)

slope, intercept = np.polyfit(merged_df["log2_chip_enrichment"], merged_df["log2FC"], 1)
x_grid = np.linspace(
    merged_df["log2_chip_enrichment"].min(),
    merged_df["log2_chip_enrichment"].max(),
    100,
)
plt.plot(
    x_grid,
    slope * x_grid + intercept,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label=f"Linear fit (Spearman r = {spearman_r:.3f})",
)

if not stat5a_df.empty:
    plt.scatter(
        stat5a_df["log2_chip_enrichment"],
        stat5a_df["log2FC"],
        color="crimson",
        s=90,
        zorder=5,
        edgecolors="black",
        label="STAT5A (TFOI)",
    )
    for _, row in stat5a_df.iterrows():
        plt.annotate(
            "STAT5A",
            (row["log2_chip_enrichment"], row["log2FC"]),
            xytext=(6, 6),
            textcoords="offset points",
            color="crimson",
            fontweight="bold",
        )

plt.xlabel("log2(Max Fold Enrichment) [GRCh38-sample01]")
plt.ylabel("log2FC (K562 8d vs 0h)")
plt.title("Histone Modification vs. Gene Expression with STAT5A")
plt.legend()
plt.tight_layout()

output_png = "/scratch/anneduffels/chipseq_project/correlation_histone_expression.png"
plt.savefig(output_png, dpi=200)
print(f"\nSaved plot to: {output_png}")
EOF
python3 /scratch/anneduffels/chipseq_project/correlation_analysis.py
