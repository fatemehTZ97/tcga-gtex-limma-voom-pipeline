# TCGA-LIHC vs. GTEx Healthy Liver: Differential Expression & Pathway Enrichment

## What this pipeline does

This pipeline performs differential expression analysis comparing TCGA liver
hepatocellular carcinoma (LIHC) tumor samples to GTEx healthy liver tissue,
using a harmonized RNA-seq resource to avoid batch artifacts from combining
separately-processed cohorts. It then runs GO and KEGG pathway enrichment on
the resulting gene list, with a specific focus on a pre-defined panel of
inflammation- and apoptosis-related candidate genes.

## Data source

- **Resource:** [UCSC Xena Toil RNA-seq Recompute](https://xenabrowser.net/datapages/), harmonized `TCGA TARGET GTEx` cohort (host: `toilHub`)
- **Expression matrix:** `TcgaTargetGtex_gene_expected_count` — RSEM expected counts, stored as log2(x+1), covering ~19,110 TCGA/GTEx/TARGET samples across ~60,498 genes
- **Sample metadata:** `TcgaTargetGTEX_phenotype`
- **Gene annotation:** `gencode.v23.annotation.gene.probemap` (Gencode v23)
- **Why this source:** TCGA and GTEx samples are reprocessed through an identical alignment/quantification pipeline in the Toil recompute, which minimizes (but does not eliminate — see Limitations) the technical differences that arise from combining independently-processed datasets from separate sources.

## Method summary

- **Tools:** R, `UCSCXenaTools`, `data.table`, `edgeR`, `limma`, `clusterProfiler`, `org.Hs.eg.db`, `pheatmap`, `EnhancedVolcano`, `ggplot2`
- **Sample selection:** TCGA-LIHC Primary + Recurrent Tumor (n=371) vs. GTEx liver Normal Tissue (n=110); TCGA's own adjacent "Solid Tissue Normal" samples are excluded from this comparison (GTEx is used as the healthy baseline instead)
- **Differential expression:** limma-voom. Expression values are back-transformed from log2(x+1) to linear counts, passed through `edgeR::calcNormFactors` (TMM normalization) and `voom()`, then fit with `lmFit`/`eBayes`. `eBayes()` is run on the full gene set (not a candidate-gene subset) so its variance-shrinkage prior is properly calibrated.
- **Batch/cohort correction:** none applied. Cohort (TCGA vs. GTEx) and disease status (tumor vs. normal) are perfectly confounded in this comparison — every TCGA sample here is a tumor and every GTEx sample is normal — so a cohort covariate is not statistically identifiable, and forcing a correction (e.g. ComBat) risks removing genuine biological signal along with any technical difference. This is treated as a stated limitation rather than something to correct away.
- **Gene set for enrichment:** adj.P.Val < 0.05 AND |log2FC| > 1 (a significance-only threshold is not informative at n=481, since most genes cross p<0.05 given the statistical power at this sample size)
- **Enrichment:** over-representation analysis via `clusterProfiler::enrichGO` (Biological Process) and `enrichKEGG`, using all tested genes as the background universe

## Key finding(s)

- The candidate inflammation genes (IL6, TNF, IL10, TLR2, TLR4) are all significantly **lower** in tumor vs. healthy liver — consistent with tumor-associated suppression of local immune/inflammatory signaling (immune evasion) rather than the more naively-expected direction.
- The candidate apoptosis genes show a mixed picture: BAX and BCL2 both significantly higher in tumor, CASP9 significantly lower, CASP3 not significant — and correspondingly, no apoptosis/cell-death pathway reaches significance in either GO or KEGG enrichment, despite individual apoptosis genes being significant. This suggests apoptosis dysregulation in HCC is gene-specific rather than pathway-wide.
- Cytokine-related pathways (e.g. KEGG "Cytokine-cytokine receptor interaction") are significantly enriched; no apoptosis, TNF-signaling, NF-κB, or Toll-like receptor pathway reaches significance.

## How to run it

1. **Requirements:** R (this pipeline was built with R 4.4.2), with `renv` for dependency management. All package installs in the script use `renv::install()`. Clone the repo, open it as an RStudio project, and run `renv::restore()` first to install the exact package versions this pipeline was built and verified with (recorded in `renv.lock`).
2. **Directory structure expected:**
   ```
   project/
   ├── data/                        (created automatically by the script; holds downloaded files)
   ├── results/                     (must exist before running the figure-saving steps)
   ├── renv.lock                    (exact package versions -- commit this to the repo)
   ├── .gitignore                   (excludes downloaded data and regenerated figures from git)
   └── tcga_gtex_lihc_analysis.R
   ```
   Raw data files and generated figures are excluded from version control via `.gitignore` (the expression matrix alone is ~1.2GB) — running the script regenerates both.
3. **Run order:** the script is written to run top-to-bottom in a single file, in the order: package setup → data download → sample selection → expression loading/formatting → differential expression → candidate gene extraction & sanity checks → genome-wide annotation → enrichment → figures.
4. **Expected runtime/resources:** the expression matrix download is ~1.2 GB compressed; a resumable download method (`curl::multi_download`) is used since the standard R downloader's 60-second timeout can prevent completion on slower connections. Expect several GB of transient RAM usage when the matrix is loaded.
5. **Outputs:** `results/kegg_barplot.png`, `results/go_bp_barplot.png`, `results/volcano_plot.png`, plus in-console tables (candidate gene DE results, GO/KEGG enrichment tables).

## Known limitations / caveats

- **No low-expression filter is applied before `voom()`** in this version of the pipeline — all 60,498 genes are tested, rather than first removing genes with negligible expression across samples.
- **GTEx is used as the "healthy" baseline**, not TCGA's own adjacent-normal tissue. This avoids the confound of TCGA adjacent-normal tissue often being cirrhotic/chronically-inflamed (not truly healthy), but introduces a different confound: cohort/disease status are perfectly aliased (see Method summary), so some portion of the observed differences may reflect residual technical variation between TCGA and GTEx processing rather than pure disease biology, despite the Toil harmonization.
- **The candidate gene panel in this script is 9 genes** (IL6, TNF, IL10, BAX, BCL2, CASP3, CASP9, TLR4, TLR2). CXCL8/IL-8 was explored separately during development but is not included in this published pipeline, since it was not run and independently verified end-to-end in the same way as the other 9 genes.
- **Bulk RNA-seq cannot resolve cell-type composition.** Differences in inflammatory gene expression between tumor and normal tissue may partly reflect differences in immune cell content between samples (e.g. resident Kupffer cell depletion in tumor tissue) rather than a cell-intrinsic transcriptional change, and this analysis cannot distinguish between those two explanations.
- **KEGG enrichment queries KEGG's servers live** at run time and requires an internet connection; results may shift slightly over time as KEGG updates its pathway definitions.
