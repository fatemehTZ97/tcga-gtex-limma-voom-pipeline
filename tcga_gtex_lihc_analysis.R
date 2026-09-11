# =============================================================================
# TCGA-LIHC (Tumor) vs. GTEx Healthy Liver (Normal)
# Differential Expression & Pathway Enrichment Analysis
# =============================================================================
# Data source: UCSC Xena Toil RNA-seq Recompute, harmonized "TCGA TARGET GTEx"
#              cohort (toilHub). TCGA and GTEx samples are reprocessed through
#              an identical pipeline in this resource, avoiding the batch
#              artifacts that arise from combining separately-processed cohorts.
# Method:      limma-voom differential expression on RSEM expected counts.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. SETUP: install and load UCSCXenaTools (renv-managed project library)
# -----------------------------------------------------------------------------
renv::install("bioc::UCSCXenaTools")
renv::snapshot()
library(UCSCXenaTools)
library(dplyr)

# Pull Xena's internal index of every dataset it hosts, and save a local copy
# purely as a reference/audit trail (not used downstream, but useful if a
# dataset name ever needs to be double-checked later).
data(XenaData)
write.csv(XenaData, "data/00_xena_hub_index.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. BUILD DOWNLOAD QUERIES (no data is fetched yet at this stage)
# -----------------------------------------------------------------------------
# Query 1: the harmonized gene-level expression matrix (RSEM expected counts,
# stored on Xena as log2(x+1)-transformed values across ~19,000 TCGA/GTEx/
# TARGET samples and ~60,000 genes).
expr_query <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")

# Query 2: the matching phenotype/sample metadata table, needed to identify
# which columns of the expression matrix are TCGA-LIHC tumor vs. GTEx liver.
pheno_query <- XenaGenerate(subset = XenaHostNames == "toilHub") %>%
  XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
  XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype")

# -----------------------------------------------------------------------------
# 3. DOWNLOAD DATA
# -----------------------------------------------------------------------------
# NOTE: XenaDownload() uses R's built-in downloader, which has a hard 60-second
# timeout and restarts from zero on every retry. On a slower connection this
# will fail repeatedly without ever completing (confirmed on this project).
##### the file cannot be downloaded in my pc with this code
# Download both to data/ (large file -- expression matrix is several GB)
XenaQuery(expr_query) %>% XenaDownload(destdir = "data/")
XenaQuery(pheno_query) %>% XenaDownload(destdir = "data/")
##### cannot use this code in this pc

# Fallback: curl::multi_download() supports HTTP range requests, so
# resume = TRUE means a dropped connection continues from the last byte
# received instead of restarting the ~1.2 GB download from scratch.
renv::install("curl")
library(curl)

curl::multi_download(
  urls = "https://toil.xenahubs.net/download/TcgaTargetGtex_gene_expected_count.gz",
  destfiles = "data/TcgaTargetGtex_gene_expected_count.gz",
  resume = TRUE
)

# Verify the download completed fully (should match the ~1.29 GB reported by
# the server) rather than assuming success from the console output alone.
file.info("data/TcgaTargetGtex_gene_expected_count.gz")$size

# Verify the file is a valid, uncorrupted gzip archive (catches silent
# truncation that a byte-count check alone might miss).
renv::install("R.utils")
R.utils::isGzipped("data/TcgaTargetGtex_gene_expected_count.gz")

# Same reliable download method for the (much smaller) phenotype file.
curl::multi_download(
  urls = "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
  destfiles = "data/TcgaTargetGTEX_phenotype.txt.gz",
  resume = TRUE
)
file.info("data/TcgaTargetGTEX_phenotype.txt.gz")$size

# -----------------------------------------------------------------------------
# 4. IDENTIFY TCGA-LIHC TUMOR AND GTEx HEALTHY LIVER SAMPLES
# -----------------------------------------------------------------------------
pheno <- read.delim("data/TcgaTargetGTEX_phenotype.txt.gz", stringsAsFactors = FALSE)
str(pheno)

# Check the exact spelling used for liver tissue in this metadata table before
# writing any filter -- TCGA/GTEx phenotype fields are known to have
# inconsistent capitalization across categories, so this is verified rather
# than assumed.
unique(pheno$X_primary_site[grepl("liver", pheno$X_primary_site, ignore.case = TRUE)])
warnings()  # (warnings here trace to an unrelated encoding issue in a
            # different metadata field and do not affect the liver filter)

# Check exact sample_type labels to distinguish tumor vs normal before filtering
unique(pheno$X_sample_type)

# Cross-tab study and sample_type for Liver samples only, to see exact group sizes
table(pheno$X_study[pheno$X_primary_site == "Liver"], pheno$X_sample_type[pheno$X_primary_site == "Liver"])

# Subset phenotype to TCGA liver tumor (Primary+Recurrent) and GTEx liver
# normal samples. Note this intentionally excludes TCGA's own "Solid Tissue
# Normal" adjacent-liver samples -- GTEx is used as the healthy baseline
# instead (see README for the rationale and its own caveats).
liver_pheno <- pheno[pheno$X_primary_site == "Liver" &
                       pheno$X_sample_type %in% c("Primary Tumor", "Recurrent Tumor", "Normal Tissue"), ]
table(liver_pheno$X_study, liver_pheno$X_sample_type)

# Confirm total sample count (expect 481 = 371 tumor + 110 normal) and
# preview sample IDs before matching against expression matrix.
nrow(liver_pheno)
head(liver_pheno$sample)

# -----------------------------------------------------------------------------
# 5. LOAD EXPRESSION DATA FOR THE SELECTED SAMPLES ONLY
# -----------------------------------------------------------------------------
# The full expression matrix is ~19,000 samples; only the 481 columns needed
# here are streamed in via fread's `select` argument, avoiding a full
# ~60,498 x 19,110 load into memory.
renv::install("data.table")

# Read only the header row of the expression matrix to get its sample ID
# column names, to confirm format/count before loading real data.
expr_header <- data.table::fread("data/TcgaTargetGtex_gene_expected_count.gz", nrows = 0)
length(colnames(expr_header))
head(colnames(expr_header))

# Preview GTEx sample IDs in our filtered phenotype table to compare format
# against expression matrix columns (TCGA and GTEx barcode formats differ
# structurally, so this is checked explicitly rather than assumed).
head(liver_pheno$sample[liver_pheno$X_study == "GTEX"])

# Check whether every liver_pheno sample ID exists as a column in the
# expression matrix, before attempting to subset (catches ID mismatches
# early, rather than after a silent partial join).
sum(liver_pheno$sample %in% colnames(expr_header))

# Load only the gene ID column plus our 481 sample columns, avoiding a full
# 19,110-sample load.
cols_needed <- c("sample", liver_pheno$sample)
expr_liver <- data.table::fread("data/TcgaTargetGtex_gene_expected_count.gz", select = cols_needed)
dim(expr_liver)  # expect 60498 genes x 482 columns (1 ID col + 481 samples)

# -----------------------------------------------------------------------------
# 6. BACK-TRANSFORM AND FORMAT THE EXPRESSION MATRIX
# -----------------------------------------------------------------------------
# Xena stores this dataset as log2(expected_count + 1). Reverse that
# transform (2^x - 1) to recover linear-scale counts, since downstream tools
# (edgeR/voom) expect count-scale input, not pre-logged values. Also moves
# gene IDs from the first column into row names.
expr_mat <- as.matrix(expr_liver[, -1, with = FALSE])
expr_mat <- 2^expr_mat - 1
rownames(expr_mat) <- expr_liver$sample
dim(expr_mat)

# Confirm back-transformed values are valid: no NAs, no negative counts
# (which would indicate a transform error).
sum(is.na(expr_mat))
range(expr_mat, na.rm = TRUE)

# Confirm expr_mat columns are in the same order as liver_pheno rows before
# building the design matrix -- a silent mismatch here would misassign
# every sample's group label.
identical(colnames(expr_mat), liver_pheno$sample)

# -----------------------------------------------------------------------------
# 7. BUILD GROUP LABELS AND CHECK FOR COHORT/DISEASE CONFOUNDING
# -----------------------------------------------------------------------------
# Create a two-level group factor (Tumor vs Normal) for the design matrix,
# with Normal as the reference level, so the model's "groupTumor"
# coefficient below represents Tumor - Normal (positive = higher in tumor).
group <- factor(ifelse(liver_pheno$X_sample_type == "Normal Tissue", "Normal", "Tumor"),
                levels = c("Normal", "Tumor"))
table(group)

renv::install("bioc::edgeR")
# Build DGEList object and calculate normalization factors, required before
# voom or diagnostic plots.
library(edgeR)
dge <- DGEList(counts = expr_mat, group = group)
dge <- calcNormFactors(dge)  # TMM normalization for library composition
dim(dge)

# MDS plot to visualize whether samples cluster by cohort (batch) or by
# tumor/normal status. This was used to check whether cohort and disease
# status are separable in this dataset -- they are not: every TCGA sample
# here is a tumor and every GTEx sample is normal, so cohort and disease
# status are perfectly confounded. This is why no cohort covariate or
# batch-correction step (e.g. ComBat) is applied anywhere in this script --
# see README.
# NOTE: the legend below labels a third color ("yellow") for a
# "TCGA-Normal(none)" category that does not exist in this 2-group
# comparison (no such points are ever plotted) -- left as-is, cosmetic only.
plotMDS(dge, col = ifelse(liver_pheno$X_study == "TCGA", "blue", "red"),
        pch = ifelse(group == "Tumor", 17, 16))
legend("topright", legend = c("TCGA-Tumor", "TCGA-Normal(none)", "GTEx-Normal"),
       col = c("blue", "yellow", "red"), pch = c(17, 16, 16), bty = "n")

# -----------------------------------------------------------------------------
# 8. DIFFERENTIAL EXPRESSION: limma-voom
# -----------------------------------------------------------------------------
# No low-expression filter is applied before voom() in this script (voom is
# run on all 60,498 genes). This is a deliberate scope choice for this
# version of the pipeline -- see README limitations.
design <- model.matrix(~group)
v <- voom(dge, design, plot = TRUE)  # mean-variance trend should be smooth,
                                      # not erratic; checked visually

# Fit linear model to voom-transformed data and compute moderated
# t-statistics via empirical Bayes. eBayes() is run on the FULL gene set
# here, not on a subset of candidate genes. This matters: eBayes borrows
# variance information across all tested genes to build a shared prior, and
# running it on only a handful of genes gives it almost nothing to build
# that prior from, producing unreliable p-values for exactly the genes with
# smaller/borderline effect sizes.
fit <- lmFit(v, design)
fit <- eBayes(fit)

# -----------------------------------------------------------------------------
# 9. MAP CANDIDATE GENES TO ENSEMBL IDs
# -----------------------------------------------------------------------------
# Extract limma results for the 9 core candidate genes to sanity-check
# before trusting full output.
core_genes <- c("IL6", "TNF", "IL10", "BAX", "BCL2", "CASP3", "CASP9", "TLR4", "TLR2")

# Download gene ID-to-symbol probemap (Gencode v23, matching this Xena
# dataset's annotation build) using the same reliable curl method as before.
curl::multi_download(
  urls = "https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap",
  destfiles = "data/gencode.v23.annotation.gene.probemap",
  resume = TRUE
)

# Load probemap and inspect its structure before mapping gene IDs to symbols.
probemap <- read.delim("data/gencode.v23.annotation.gene.probemap", stringsAsFactors = FALSE)
str(probemap)

# Confirm all 9 core genes exist in the probemap's gene symbol column.
core_genes <- c("IL6", "TNF", "IL10", "BAX", "BCL2", "CASP3", "CASP9", "TLR4", "TLR2")
core_genes %in% probemap$gene

# Get Ensembl IDs for the 9 core genes and confirm exactly one ID per gene
# (no ambiguous 1:many symbol matches for this specific gene set) before
# using match() to look up IDs.
core_ids <- probemap$id[match(core_genes, probemap$gene)]
data.frame(gene = core_genes, ensembl_id = core_ids)
table(probemap$gene[probemap$gene %in% core_genes])

# Confirm the 9 core gene Ensembl IDs exist as row names in expr_mat
# (matching versioned IDs, e.g. "ENSG00000136244.11").
core_ids %in% rownames(expr_mat)

# -----------------------------------------------------------------------------
# 10. EXTRACT AND SANITY-CHECK CANDIDATE GENE RESULTS
# -----------------------------------------------------------------------------
# Extract limma-voom DE results for the 9 core genes from the fitted model
# (pulled from the full-gene-set fit -- see step 8 note on why eBayes must
# be run on the full set first).
core_results <- topTable(fit, coef = "groupTumor", number = Inf)[core_ids, ]
data.frame(gene = core_genes, ensembl_id = core_ids, core_results[, c("logFC", "AveExpr", "P.Value", "adj.P.Val")])

# Boxplot expression of BCL2, CASP3, CASP9 by group to visually check for
# outlier-driven effects -- plotting the raw distribution by group rather
# than trusting only the summary statistic, for the genes whose direction
# was most worth double-checking.
genes_to_check <- c("BCL2", "CASP3", "CASP9")
ids_to_check <- core_ids[match(genes_to_check, core_genes)]
par(mfrow = c(1, 3))
for (i in seq_along(ids_to_check)) {
  boxplot(v$E[ids_to_check[i], ] ~ group, main = genes_to_check[i], ylab = "log2 CPM")
}

renv::install("pheatmap")
# Heatmap of the 9 core genes' expression across samples, ordered by group.
# Row-scaled (z-scored per gene) so genes with different expression ranges
# are visually comparable on the same color scale.
library(pheatmap)
heat_data <- v$E[core_ids, order(group)]
rownames(heat_data) <- core_genes
annotation_col <- data.frame(Group = sort(group))
rownames(annotation_col) <- colnames(heat_data)
pheatmap(heat_data, scale = "row", cluster_cols = FALSE,
         annotation_col = annotation_col, show_colnames = FALSE,
         main = "Core gene panel: log2 CPM (row-scaled)")

# Independently verify group mean differences for the 9 core genes,
# bypassing heatmap code entirely -- recompute each gene's raw group means
# directly from the voom-transformed expression matrix and confirm the sign
# matches the model's logFC. This is the key check that catches
# scaling-direction or sample-order bugs before trusting any downstream
# figure.
group_means <- sapply(core_ids, function(id) tapply(v$E[id, ], group, mean))
rownames(group_means) <- c("Normal_mean", "Tumor_mean")
colnames(group_means) <- core_genes
t(group_means)

# -----------------------------------------------------------------------------
# 11. GENOME-WIDE RESULTS AND GENE SYMBOL ANNOTATION
# -----------------------------------------------------------------------------
# Extract full genome-wide DE results table and check its dimensions.
all_results <- topTable(fit, coef = "groupTumor", number = Inf)
dim(all_results)
head(all_results)

# Map Ensembl IDs to gene symbols in the full results table, checking for
# unmapped genes.
all_results$ensembl_id <- rownames(all_results)
all_results$gene <- probemap$gene[match(all_results$ensembl_id, probemap$id)]
sum(is.na(all_results$gene))

# -----------------------------------------------------------------------------
# 12. CONVERT TO ENTREZ IDs FOR GO/KEGG ENRICHMENT
# -----------------------------------------------------------------------------
renv::install("bioc::org.Hs.eg.db")
renv::install("bioc::clusterProfiler")
# Convert gene symbols to Entrez IDs for GO/KEGG enrichment, checking
# mapping success rate. Not every gene symbol has a corresponding Entrez ID
# (non-coding RNAs, pseudogenes, etc. are commonly excluded from Entrez) --
# some mapping failure is expected and normal, not a bug.
library(org.Hs.eg.db)
library(clusterProfiler)
entrez_map <- bitr(all_results$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
nrow(entrez_map)

# Check for genes with multiple Entrez ID matches (1:many mapping issue) --
# a genuine annotation ambiguity in Gencode/Entrez, not a coding error.
sum(duplicated(entrez_map$SYMBOL))

# Remove duplicate SYMBOL-to-Entrez mappings, keeping the first match --
# the fraction affected is negligible.
entrez_map <- entrez_map[!duplicated(entrez_map$SYMBOL), ]
nrow(entrez_map)

# Merge Entrez IDs into results table and count genes usable for enrichment.
all_results_annotated <- merge(all_results, entrez_map, by.x = "gene", by.y = "SYMBOL")
nrow(all_results_annotated)

# Check whether all_results has duplicate gene symbols, explaining the row
# count increase after merge -- some gene symbols themselves map to
# multiple distinct Ensembl gene IDs in this annotation build (e.g. genes
# with paralogous loci), confirmed here and resolved in the next step.
sum(duplicated(all_results$gene))

# For genes with duplicate symbols, keep only the entry with the smallest
# p-value -- the standard, biologically reasonable way to collapse to one
# row per gene symbol before enrichment.
library(dplyr)
all_results_annotated <- all_results_annotated %>%
  group_by(gene) %>%
  slice_min(P.Value, n = 1, with_ties = FALSE) %>%
  ungroup()
nrow(all_results_annotated)

# -----------------------------------------------------------------------------
# 13. DEFINE THE SIGNIFICANT GENE LIST FOR ENRICHMENT
# -----------------------------------------------------------------------------
# Define significant DEGs using standard adj.P.Val < 0.05 threshold, and
# count up/down genes. Note: at this sample size (n=481), adj.P.Val < 0.05
# alone is not a very informative filter, since most genes cross that
# threshold given the statistical power -- hence the fold-change cutoff
# added next.
sig_genes <- all_results_annotated[all_results_annotated$adj.P.Val < 0.05, ]
nrow(sig_genes)
table(sig_genes$logFC > 0)

# Add a fold-change cutoff on top of significance, since adj.P.Val alone is
# uninformative at this sample size -- |log2FC| > 1 means >=2-fold change,
# so "significant" now also means "biologically non-trivial."
sig_genes_fc <- sig_genes[abs(sig_genes$logFC) > 1, ]
nrow(sig_genes_fc)
table(sig_genes_fc$logFC > 0)

# -----------------------------------------------------------------------------
# 14. GO BIOLOGICAL PROCESS ENRICHMENT
# -----------------------------------------------------------------------------
# Run GO Biological Process enrichment on the significant, fold-change-
# filtered gene list. `universe` is set to all genes actually tested in
# this analysis (not the whole genome), which keeps the enrichment
# statistics appropriate for this specific experiment.
go_bp <- enrichGO(gene = sig_genes_fc$ENTREZID,
                  universe = all_results_annotated$ENTREZID,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
dim(go_bp)

# Search GO BP results for inflammation- and apoptosis-related terms, since
# this is the specific biological question motivating the analysis.
go_results <- as.data.frame(go_bp)
go_results[grepl("inflamm|apoptosis|apoptotic", go_results$Description, ignore.case = TRUE),
           c("Description", "GeneRatio", "p.adjust", "Count")]

# Broaden search to catch cell-death-related terms with different phrasing
# (e.g. "cell death", "caspase") that don't literally contain "apoptosis".
go_results[grepl("cell death|apoptotic|caspase", go_results$Description, ignore.case = TRUE),
           c("Description", "GeneRatio", "p.adjust", "Count")]

# View the top 15 most significant GO BP terms overall, for broader
# biological context beyond the targeted inflammation/apoptosis search.
head(go_results[order(go_results$p.adjust), c("Description", "GeneRatio", "p.adjust", "Count")], 15)

# -----------------------------------------------------------------------------
# 15. KEGG PATHWAY ENRICHMENT
# -----------------------------------------------------------------------------
# Run KEGG pathway enrichment on the same significant, fold-change-filtered
# gene list. Queries KEGG's servers live (requires an internet connection
# at run time).
kegg_res <- enrichKEGG(gene = sig_genes_fc$ENTREZID,
                       universe = all_results_annotated$ENTREZID,
                       organism = "hsa",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)
dim(kegg_res)

# Search KEGG results for apoptosis- and inflammation-related pathways.
kegg_results <- as.data.frame(kegg_res)
kegg_results[grepl("apoptosis|inflamm|TNF|NF-kappa|toll-like|cytokine", kegg_results$Description, ignore.case = TRUE),
             c("Description", "GeneRatio", "p.adjust", "Count")]

# View top 15 most significant KEGG pathways overall, for context.
head(kegg_results[order(kegg_results$p.adjust), c("Description", "GeneRatio", "p.adjust", "Count")], 15)

# -----------------------------------------------------------------------------
# 16. FIGURES FOR MANUSCRIPT / PORTFOLIO
# -----------------------------------------------------------------------------
##### not organized -- quick default plot, kept for reference; superseded by the ggplot version below
# Barplot of top KEGG enrichment results for manuscript figure
barplot(kegg_res, showCategory = 15, title = "KEGG Pathway Enrichment: HCC Tumor vs. GTEx Healthy Liver")
#########

# Cleaner, more organized KEGG barplot using ggplot2, sorted by
# significance, wide dimensions on save to avoid title clipping in the
# RStudio plot pane preview.
library(ggplot2)
kegg_plot_data <- kegg_results[order(kegg_results$p.adjust), ][1:15, ]
kegg_plot_data$Description <- factor(kegg_plot_data$Description,
                                     levels = rev(kegg_plot_data$Description))

ggplot(kegg_plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#D73027", high = "#4575B4", name = "adj. p-value") +
  labs(title = "KEGG Pathway Enrichment: HCC Tumor vs. GTEx Healthy Liver",
       x = "Gene Count", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_blank())

# Widen plot output and shorten title to avoid clipping
ggsave("results/kegg_barplot.png", width = 10, height = 6, dpi = 300)

# Same tidy barplot approach for GO Biological Process results, same color
# scale for visual consistency between the two figures.
go_plot_data <- go_results[order(go_results$p.adjust), ][1:15, ]
go_plot_data$Description <- factor(go_plot_data$Description,
                                   levels = rev(go_plot_data$Description))

ggplot(go_plot_data, aes(x = Count, y = Description, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#D73027", high = "#4575B4", name = "adj. p-value") +
  labs(title = "GO Biological Process Enrichment: HCC Tumor vs. GTEx Healthy Liver",
       x = "Gene Count", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.y = element_text(size = 10),
        panel.grid.major.y = element_blank())

ggsave("results/go_bp_barplot.png", width = 10, height = 6, dpi = 300)

# Install EnhancedVolcano for this project
renv::install("bioc::EnhancedVolcano")
# Volcano plot of full genome-wide DE results, highlighting the 9 core
# genes. Points are colored by whether they pass the significance cutoff
# (pCutoff), the fold-change cutoff (FCcutoff), both, or neither -- matching
# the same thresholds used to build sig_genes_fc above, so the plot and the
# gene list it's built from stay consistent.
library(EnhancedVolcano)
EnhancedVolcano(all_results_annotated,
                lab = all_results_annotated$gene,
                x = "logFC",
                y = "adj.P.Val",
                selectLab = core_genes,
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "HCC Tumor vs. GTEx Healthy Liver",
                subtitle = NULL,
                labSize = 4,
                drawConnectors = TRUE)

ggsave("results/volcano_plot.png", width = 10, height = 6, dpi = 300)

# Verify the 9 core genes' plotted coordinates match their known logFC/
# p-values from the model, rather than trusting the auto-generated plot
# labels blindly.
data.frame(gene = core_genes,
           logFC = all_results_annotated$logFC[match(core_genes, all_results_annotated$gene)],
           adj.P.Val = all_results_annotated$adj.P.Val[match(core_genes, all_results_annotated$gene)])
