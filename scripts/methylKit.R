# Get files and samples, treatments, etc
# Create methylkit db
# Histograms
# PCA
# Dendrograms



library(methylKit)

files = list.files("/input/methylCpG_db/", pattern="txt.bgz", full.names=T)
files = files[!grepl(".tbi", files)]
files = files[!grepl("methylBase", files)]

files_bname = basename(files)

files = lapply(files, print)

sample_ids = gsub(".txt.bgz", "", files_bname)
sample_ids = lapply(sample_ids, print)
treatments = c(
    rep("CuONPs", 4),
    rep("Epoxy_Abrad", 4),
    rep("Epoxy_Incinerated", 4),
    rep("FLG", 4),
    rep("FLG-Epoxy_Abrad", 4),
    rep("FLG-Epoxy_Incinerated", 4),
    rep("Vehicle_Control", 4)
)

treatments_mapping = data.frame(
    Treatment = unique(treatments),
    number = 0:(length(unique(treatments)) - 1)
)

treatments_numeric = sapply(treatments, function(x) {
    return(treatments_mapping$number[treatments_mapping$Treatment == x])
})

meths = methRead(
    files,
    sample.id = sample_ids,
    assembly = "hg38",
    dbtype = "tabix",
    context = "CpG",
    treatment = treatments_numeric
)

meths_united = unite(meths, destrand = T, save.db = F, mc.cores = 20)

pdf("sample_correlation.pdf")
getCorrelation(meth,plot=TRUE)
dev.off()

pdf("PCA_destranded.pdf")
PCASamples(meths_united, screeplot = TRUE)
PCASamples(meths_united)
dev.off()

pdf("clustering_destranded.pdf")
clusterSamples(meths_united, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

parts = strsplit(meths_united@sample.ids, "_")
donors = sapply(parts, function(x) {return(x[length(x) - 1])})

sample_annotation = data.frame(donor = donors)

as = assocComp(mBase = meths_united, sampleAnnotation = sample_annotation)
meths_united_no_donor_effect = removeComp(meths_united, comp = 1:4)
pdf("clustering2_destranded.pdf")
clusterSamples(meths_united_no_donor_effect, dist = "correlation", method = "ward", plot = TRUE)
dev.off()

mat = percMethylation(meths_united)
mat = mat / 100
write.table(mat, "percMethyl.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

library(ComBatMet)

perc_meth = read.table("percMethyl.tsv", sep = "\t", header = T)

samples = colnames(perc_meth)
donors = rep(1:4, 7)
treatments = c(
    rep(0, 4),
    rep(1, 4),
    rep(2, 4),
    rep(2, 4),
    rep(3, 4),
    rep(4, 4),
    rep(5, 4)
) 

adj_bv_mat = ComBat_met(perc_meth, batch = donors, group = treatments, full_mod = TRUE)

write.table(adj_bv_mat, "adj_percMeth.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

perc_meth = read.table("adj_percMeth.tsv", sep = "\t", header = T)
reconstructed_meth_after_combatmet = reconstruct(perc_meth, meths_united)

pdf("clustering_combatmet.pdf")
clusterSamples(reconstructed_meth_after_combatmet, dist = "correlation", method = "ward", plot = TRUE)
dev.off()