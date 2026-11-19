cat << 'EOF' > /scratch/anneduffels/chipseq_project/run_all_samples.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

de_path = "/scratch/anneduffels/brbseq/results/de_K562_8d_vs_0h.tsv"
de_df = pd.read_csv(de_path, sep="\t")
de_df["gene_upper"] = de_df["Gene"].astype(str).str.upper().str.strip()

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

samples = [
    ("sample01", "K562_H3K27ac"),
    ("sample02", "K562_H3K27me3"),
    ("sample03", "K562_H3K4me3"),
]

fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

for ax, (sample_id, mark_label) in zip(axes, samples):
    peaks_path = f"/scratch/anneduffels/chipseq_project/results/macs2/GRCh38-{sample_id}_peaks.xls"
    peaks_df = pd.read_csv(peaks_path, sep="\t", comment="#")

    if peaks_df.empty:
        print(f"\n=== {sample_id} ({mark_label}): no peaks found ===")
        ax.set_title(f"{sample_id} ({mark_label})\nNo peaks")
        continue

    peaks_df["peak_mid"] = (peaks_df["start"] + peaks_df["end"]) // 2
    annotated_chunks = []
    for chrom, peak_sub in peaks_df.groupby("chr"):
        gene_sub = genes_df[genes_df["chr"] == chrom]
        if gene_sub.empty:
            continue
        distances = np.abs(gene_sub["tss"].values[:, None] - peak_sub["peak_mid"].values[None, :])
        nearest_idx = np.argmin(distances, axis=0)
        sub_copy = peak_sub.copy()
        sub_copy["gene_upper"] = gene_sub["gene_upper"].values[nearest_idx]
        sub_copy["tss_distance"] = distances[nearest_idx, np.arange(len(peak_sub))]
        annotated_chunks.append(sub_copy)

    peaks_annot = pd.concat(annotated_chunks, ignore_index=True)
    nearby = peaks_annot[peaks_annot["tss_distance"] <= 50000]
    chip_by_gene = (
        nearby.groupby("gene_upper", as_index=False)["fold_enrichment"]
        .max()
        .rename(columns={"fold_enrichment": "max_fold_enrichment"})
    )

    merged = pd.merge(de_df, chip_by_gene, on="gene_upper", how="inner").dropna(
        subset=["log2FC", "max_fold_enrichment"]
    )
    merged["log2_chip"] = np.log2(merged["max_fold_enrichment"].clip(lower=1e-3))

    spearman_r = merged["log2FC"].rank().corr(merged["log2_chip"].rank(), method="pearson")
    print(f"\n=== {sample_id} ({mark_label}) ===")
    print(f"Matched genes: {len(merged)} | Spearman r: {spearman_r:.4f}")

    stat5a = merged[merged["gene_upper"] == "STAT5A"]
    if not stat5a.empty:
        print(stat5a[["Gene", "log2FC", "padj", "max_fold_enrichment", "log2_chip"]])
    else:
        print("STAT5A not found within 50 kb.")

    ax.scatter(merged["log2_chip"], merged["log2FC"], alpha=0.35, s=12, label="Matched genes")
    if len(merged) > 1:
        slope, intercept = np.polyfit(merged["log2_chip"], merged["log2FC"], 1)
        x_grid = np.linspace(merged["log2_chip"].min(), merged["log2_chip"].max(), 100)
        ax.plot(x_grid, slope * x_grid + intercept, color="black", linestyle="--", linewidth=1.2)

    if not stat5a.empty:
        ax.scatter(stat5a["log2_chip"], stat5a["log2FC"], color="crimson", s=70, zorder=5, edgecolors="black")
        for _, row in stat5a.iterrows():
            ax.annotate("STAT5A", (row["log2_chip"], row["log2_chip"] * 0 + row["log2FC"]), xytext=(5, 5), textcoords="offset points", color="crimson", fontweight="bold")

    ax.set_title(f"{sample_id} ({mark_label})\nr = {spearman_r:.3f}, n = {len(merged)}")
    ax.set_xlabel("log2(Max Fold Enrichment)")

axes[0].set_ylabel("log2FC (K562 8d vs 0h)")
plt.tight_layout()
out_path = "/scratch/anneduffels/chipseq_project/correlation_all_histones.png"
plt.savefig(out_path, dpi=200)
print(f"\nSaved comparison plot to: {out_path}")
EOF
python3 /scratch/anneduffels/chipseq_project/run_all_samples.py
