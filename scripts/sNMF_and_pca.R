#################### packages
# LEA is on Bioconductor, install with BiocManager::install("LEA") if not already available
library(LEA)
library(vcfR)
library(adegenet)
library(ggplot2)

#################### variables
vcf_file  <- "Sample_X.filtered.vcf.gz"
k_range   <- 1:10
best_reps <- 10

#################### read vcf
vcf <- read.vcfR(vcf_file, verbose = FALSE)

#################### convert to geno format for LEA
# vcf2geno writes Sample_X.geno directly from the vcf, no intermediate PLINK step needed
vcf2geno(vcf_file, output.file = "Sample_X.geno")

#################### sNMF across K
project <- snmf("Sample_X.geno",
                 K = k_range,
                 repetitions = best_reps,
                 entropy = TRUE,
                 project = "new")

pdf("cross_entropy.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()

#################### pick best K
# lowest cross-entropy across repetitions at each K
ce <- sapply(k_range, function(k) min(cross.entropy(project, K = k)))
best_k <- k_range[which.min(ce)]
best_run <- which.min(cross.entropy(project, K = best_k))

#################### ancestry barplot at best K
qmatrix <- Q(project, K = best_k, run = best_run)

pdf("ancestry_barplot.pdf", width = 10, height = 4)
barplot(t(qmatrix), col = 1:best_k, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry coefficient")
dev.off()

#################### PCA
genlight <- vcfR2genlight(vcf)
pca <- glPca(genlight, nf = 2)

cluster <- apply(qmatrix, 1, which.max)

pca_df <- data.frame(PC1 = pca$scores[, 1],
                      PC2 = pca$scores[, 2],
                      cluster = factor(cluster))

p <- ggplot(pca_df, aes(PC1, PC2, color = cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = paste0("PCA (colored by sNMF K=", best_k, ")"))
