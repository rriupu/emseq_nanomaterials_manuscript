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
meths_perc = percMethylation(meths_norm)
meths_cor = cor(meths_perc)
heatmap_cor = pheatmap(meths_cor)
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






parts = strsplit(meths_united@sample.ids, "_")
donors = sapply(parts, function(x) {return(x[length(x)])})
treatments = sapply(parts, "[[", 1)

covariates = data.frame(
    donor = factor(donors)
)

diff_meth = calculateDiffMeth(
    meths_united,
    covariates = covariates,
    mc.cores = 15
)

# Genomation
refseq_annots = readTranscriptFeatures("refseq_hg38.bed")
diff_ann = annotateWithGeneParts(as(diff_meth, "GRanges"), refseq_annots)
plotTargetAnnotation(
    diff_ann,
    precedence = TRUE,
    main = "Differential methylation annotation"
)

promoters = regionCounts(meths_norm, refseq_annots$promoters)



as = assocComp(mBase = meths_united, sampleAnnotation = sample_annotation)
meths_united_no_donor_effect = removeComp(meths_united, comp = 1:4)

pdf(file.path(output_dir, "PCA_no_donor_effect.pdf"))
PCASamples(meths_united_no_donor_effect, screeplot = TRUE)
PCASamples(meths_united_no_donor_effect)
dev.off()

# Elana






# pdf("clustering2_destranded.pdf")
# clusterSamples(meths_united_no_donor_effect, dist = "correlation", method = "ward", plot = TRUE)
# dev.off()

# mat = percMethylation(meths_united)
# mat = mat / 100
# write.table(mat, "percMethyl.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

# library(ComBatMet)

# perc_meth = read.table("percMethyl.tsv", sep = "\t", header = T)

# samples = colnames(perc_meth)
# donors = rep(1:4, 7)
# treatments = c(
#     rep(0, 4),
#     rep(1, 4),
#     rep(2, 4),
#     rep(2, 4),
#     rep(3, 4),
#     rep(4, 4),
#     rep(5, 4)
# ) 

# adj_bv_mat = ComBat_met(perc_meth, batch = donors, group = treatments, full_mod = TRUE)

# write.table(adj_bv_mat, "adj_percMeth.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

# perc_meth = read.table("adj_percMeth.tsv", sep = "\t", header = T)
# reconstructed_meth_after_combatmet = reconstruct(perc_meth, meths_united)

# pdf("clustering_combatmet.pdf")
# clusterSamples(reconstructed_meth_after_combatmet, dist = "correlation", method = "ward", plot = TRUE)
# dev.off()