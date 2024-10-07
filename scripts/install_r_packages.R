# install_r_packages.R

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load BiocManager and install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install DESeq2
BiocManager::install("DESeq2")

install.packages(c("gplots","RUnit"), repos = "https://cloud.r-project.org")