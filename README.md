**MutSigCompare**
* It is lightweight Python tool for analyzing and comparing mutational signatures between two samples or groups.
* Requires only the input vcf files between the samples to be compared through a csv file, which then parses SBS mutations, extract trinucleotide contexts, and visualize the relative contributions of COSMIC mutational signatures through stacked barplots, heatmaps, and PCA.
* Useful for group-wise comparison support (e.g., tumor vs normal, treated vs untreated)

## Installation
git clone the repository

## Usage
```
python trinucleotide_parsing.py --sample samples.csv --ref DATA/reference.fna
```
samples.csv input format:
```
vcf_file,group
G1_sample1.vcf,Group1
G1_sample2.vcf,Group1
G2_sample1.vcf,Group2
G2_sample2.vcf,Group2
```
## Output
After running you get two types of output:

1)SBS of mutational signatures of each sample

2)Plots:
  + BAR_PLOT.jpeg → Stacked barplot of mutational signatures
  + HEAT_MAP.jpeg → Signature heatmap
  + Clustered_HEAT_MAP.jpeg → Clustered heatmap (based on cosine similarity)
  + PCA_plot.jpeg → PCA plot showing replicate and group relationships
