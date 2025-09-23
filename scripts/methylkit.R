# Get files and samples, treatments, etc
# Create methylkit db
# Histograms
# PCA
# Dendrograms
## --------------- ##
## Library loading ##
## --------------- ##

library(methylKit)
library(genomation)
library(optparse)
library(pheatmap)
library(ComplexHeatmap)

## --------- ##
## Functions ##
## --------- ##

source("scripts/methylkit_utils.R")

## -------------- ##
## Read arguments ##
## -------------- ##

option_list = list(

  make_option(c("-i", "--input_dir"), type = "character", default = NULL, 
              help = "Path to the input directory containing all the called methylation files from Bismark (Mandatory).", 
              metavar = "character"),
  
  make_option(c("-o", "--output_dir"), type = "character", default = NULL, 
              help = "Path to the output directory where all plots will be saved (Mandatory). ", 
              metavar = "character"),
  
  make_option(c("-t", "--threads"), type = "numeric", default = 1, 
              help = "Threads to use in parallel processing (Mandatory). ", 
              metavar = "numeric")
);
opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

input_dir = opt$input_dir
output_dir = opt$output_dir
methylkit_dbdir = file.path(output_dir, "methylkit_database")
threads = opt$threads

dir.create(methylkit_dbdir, recursive = T, showWarnings = F)

## ---------------- ##
## Find input files ##
## ---------------- ##

files = list.files(input_dir, pattern = "CpG_report.txt.gz", full.names = T, recursive = T)

## ------------------------------------- ##
## Prepare required inputs for methylKit ##
## ------------------------------------- ##

files_bname = basename(files)

files = lapply(files, print)

sample_ids = gsub("_S[0-9]{1,3}_R1_001_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz", "", files_bname)
sample_ids_ls = lapply(sample_ids, print)

treatments = gsub("_Donor[0-9]", "", sample_ids)
treatments = gsub("_[0-9]{1,2}ug", "", treatments)

treatments_mapping = data.frame(
    Treatment = unique(treatments),
    number = 0:(length(unique(treatments)) - 1)
)

treatments_numeric = sapply(treatments, function(x) {
    return(treatments_mapping$number[treatments_mapping$Treatment == x])
})

## ---------------- ##
## Read input files ##
## ---------------- ##

meths = methRead(
    files,
    sample.id = sample_ids_ls,
    pipeline = "bismarkCytosineReport",
    assembly = "hg38",
    dbtype = "tabix",
    dbdir = methylkit_dbdir,
    context = "CpG",
    treatment = treatments_numeric,
    mincov = 6
)

## --------------------------------- ##
## Plot CpG percentages and coverage ##
## --------------------------------- ##

pdf(file.path(output_dir, "CpG_methylation_percentage_histograms.pdf"))
for (i in 1:length(sample_ids_ls)) {
    getMethylationStats(meths[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

pdf(file.path(output_dir, "CpG_coverage.pdf"))
for (i in 1:length(sample_ids_ls)) {
    getCoverageStats(meths[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

## ----------------------------------------- ## 
## Filtering of PCR biases and normalization ##
## ----------------------------------------- ##

meths_filtered = filterByCoverage(meths, lo.count = 6, hi.count = 100)
meths_norm = normalizeCoverage(meths_filtered)

pdf(file.path(output_dir, "CpG_methylation_percentage_histograms_normalized.pdf"))
for (i in 1:length(sample_ids_ls)) {
    getMethylationStats(meths_norm[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

pdf(file.path(output_dir, "CpG_coverage_normalized.pdf"))
for (i in 1:length(sample_ids_ls)) {
    getCoverageStats(meths_norm[[i]], plot = TRUE, both.strands = FALSE)
}
dev.off()

## --------------------------------------------------- ##
## Unite all methylation datasets for further analysis ##
## --------------------------------------------------- ##

meths_united = unite(
    meths_norm,
    destrand = T,
    mc.cores = threads,
    save.db = TRUE,
    dbdir = methylkit_dbdir,
    suffix = "united_filtered_normalized"
)

# Save united dataset
# makeMethylDB(meths_united, dbdir = methylkit_dbdir)

## ---------------- ##
## Data exploration ##
## ---------------- ##

# Sample-sample correlation
meths_perc = percMethylation(meths_united)
meths_cor = cor(meths_perc)
heatmap = Heatmap(
    meths_cor,
    show_row_names = T,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8),
    top_annotation = HeatmapAnnotation(
        Treatment = treatments,
        Donor = paste("Donor", rep(1:4, 7), sep = " ")
    )
)
pdf(file.path(output_dir, "correlation_heatmap.pdf"), width = 14)
draw(heatmap)
dev.off()

pca = prcomp(t(meths_perc))
pca_df = as.data.frame(pca$x)
pca_df$donor = rep(1:4, 7)
pca_df$donor = factor(pca_df$donor)
pca_plot = ggplot(pca_df, aes(x = PC1, y = PC2, color = donor)) + geom_point()
ggsave(file.path(output_dir, "pca_ggplot.pdf"), pca_plot)

ggsave(file.path(outputDir, "correlation_heatmap.pdf"), heatmap_cor)

pdf(file.path(output_dir, "sample_correlation.pdf"), height = 20, width = 20)
getCorrelation(meths_united, plot = TRUE)
dev.off()

# PCA
pdf(file.path(output_dir, "PCA_normalized.pdf"))
PCASamples(meths_united, screeplot = TRUE)
PCASamples(meths_united)
dev.off()

# Hierarchical clustering
pdf(file.path(output_dir, "clustering_normalized.pdf"))
clusterSamples(meths_united, dist = "correlation", method = "ward", plot = TRUE)
dev.off()


## --------------------- ##
## Differential analysis ##
## --------------------- ##

groups = 0:6
pairs = combn(groups, 2, simplify = FALSE)
pairwise_results = list()

for (pair in pairs) {
  g1 = pair[1]
  g2 = pair[2]

  group1_name = treatments_mapping$Treatment[treatments_mapping$number == g1]
  group2_name = treatments_mapping$Treatment[treatments_mapping$number == g2]
  
  cat(sprintf("Testing treatment groups %s vs %s\n", group1_name, group2_name))
  
  res = test_two_treatments(meths_united, group1 = g1, group2 = g2, threads = threads, qvalue = 0.1)
  
  # Name results with pair info for easy access
  pair_name = paste0(group1_name, "_vs_", group2_name)
  pairwise_results[[pair_name]] = res
}

# Plot number of differentially methylated CpGs per comparison
sig_counts = sapply(pairwise_results, function(res) {
  df = as.data.frame(res)
  sum(df$qvalue < 0.05 & abs(df$meth.diff) >= 10)
})

# Convert to data frame for plotting
counts_df = data.frame(
  comparison = names(sig_counts),
  sig_sites = sig_counts
)

# Bar plot
p = ggplot(counts_df, aes(x = comparison, y = sig_sites, fill = comparison)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(
        title = "Number of Significant Differential Methylation Sites",
        x = "Pairwise Comparison",
        y = "Significant Sites (q < 0.05 & |meth.diff| >= 10%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "number_significant_diff_methylation_sites_per_comparison.pdf"), p)

## --------------------- ##
## Annotation of results ##
## --------------------- ##

# Genomation
refseq_annots = readTranscriptFeatures("refseq_hg38.bed")
annot_pairwise_results = lapply(
    pairwise_results, function(x) {
        annotateWithGeneParts(as(x, "GRanges"), refseq_annots)
    }
)

annot_stats = lapply(
    annot_pairwise_results,
    function(x) {
        getTargetAnnotationStats(x, percentage = TRUE, precedence = TRUE)
    }
)

pdf(file.path(output_dir, "diff_meth_annotations.pdf"))
for(i in 1:length(annot_pairwise_results)) {
    
    comparison = pairs[[i]]
    g1 = pair[1]
    g2 = pair[2]

    group1_name = treatments_mapping$Treatment[treatments_mapping$number == g1]
    group2_name = treatments_mapping$Treatment[treatments_mapping$number == g2]

    title = paste0(group1_name, "_vs_", group2_name)

    plotTargetAnnotation(
        annot_pairwise_results[[i]],
        precedence = TRUE,
        main = title
    )

}
dev.off()

# TO DOs
# Investigate other possibilities to remove donor effect?
# For each pairwise comparison, get differentially methylated positions and draw a heatmap
# Write results to table
# Analysis on CHH and CHG
# Plot with methylation distribution around TSS (i.e. https://www.researchgate.net/figure/Methylation-density-across-transcription-start-sites-TSS-regions-of-all-B-cinerea_fig1_305809128)