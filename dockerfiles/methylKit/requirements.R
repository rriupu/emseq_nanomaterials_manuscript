# CRAN packages:
required_packages_cran = c("BiocManager", "devtools", "pheatmap")

install.packages(required_packages_cran)

# Install methylKit from source due to newer version fixing sorting problem 
# (https://github.com/al2na/methylKit/pull/365)
library(devtools)
install_github(
    "al2na/methylKit",
    ref = "destrand_per_chrom",
    build_vignettes = FALSE, 
    repos = BiocManager::repositories(),
    dependencies = TRUE
)

## Bioconductor packages:
required_packages_bioconductor = c(
  "genomation",
  "ComplexHeatmap"
)

BiocManager::install(required_packages_bioconductor)
