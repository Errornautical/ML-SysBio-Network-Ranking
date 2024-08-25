install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("impute")
install.packages("caret", dependencies = c("Depends", "Suggests"))
##Code Starts Here
library(GEOquery)
my_id <- "GSE54536"
gse <- getGEO(my_id)
gse <- gse[[1]]
gse 
pData(gse)
fData(gse) 
exprs(gse) 
summary(exprs(gse))
exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)
library(dplyr)
sampleInfoPD <- pData(gse)
sampleInfoPD
sampleInfoPD <- select(sampleInfoPD, title,characteristics_ch1.1)
sampleInfoPD
sampleInfoPD <- rename(sampleInfoPD,Patient = title, Group=characteristics_ch1.1)
sampleInfoPD


library(pheatmap)
corMatrix <- cor(exprs(gse),use="c")
 corMatrix
corMatrix[is.na(corMatrix)] <- 0

rownames(sampleInfoPD)
colnames(corMatrix)

rownames(sampleInfoPD) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfoPD)    
library(ggplot2)
library(ggrepel)
any(is.na(exprs(gse)))     
any(is.infinite(exprs(gse)))
exprs_no_na <- na.omit(exprs(gse))
finite_rows <- !apply(exprs_no_na, 1, function(x) any(is.infinite(x)))
exprs_clean <- exprs_no_na[finite_rows, ]
pca <- prcomp(t(exprs_clean))
cbind(sampleInfoPD, pca$x) %>%
  ggplot(aes(x = PC1, y = PC2, col = Group, label = paste("Patient", Patient))) +
  geom_point(aes(shape = Group), size = 3) +  # Set point size
  geom_text_repel() +
  scale_color_manual(values = c("diagnosis: Parkinson's disease" = "blue", "diagnosis: healthy" = "red")) +
  scale_shape_manual(values = c("diagnosis: Parkinson's disease" = 17, "diagnosis: healthy" = 16))


library(limma)
colnames(sampleInfoPD)
design <- model.matrix(~0+sampleInfoPD$Group)
design
colnames(design) <- c("Healthy","Parkinsons_Disease")
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)
contrasts <- makeContrasts(Parkinsons_Disease - Healthy, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
topTable(fit2)
decideTests(fit2)
table(decideTests(fit2))
anno <- fData(gse)
anno
anno <- select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)
full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

library(ggplot2)
library(ggrepel)

p_cutoff <- 0.05
fc_cutoff <- 1

sig_genes <- full_results %>%
  mutate(Significant = P.Value < p_cutoff & abs(logFC) > fc_cutoff,
         Direction = case_when(
           logFC > fc_cutoff ~ "Up",
           logFC < -fc_cutoff ~ "Down",
           TRUE ~ "Not_significant"
         ))

volcano_plot <- ggplot(sig_genes, aes(x = logFC, y = -log10(P.Value), color = Direction)) +
  geom_point(alpha = 0.8, size = 1.2, shape = 16) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not_significant" = "black"),
                     guide = guide_legend(override.aes = list(shape = 16, size = 5))) +
  xlab("log2(fold change)") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot") +
  theme_bw() +
  geom_hline(yintercept = c(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
  annotate("text", x = 10, y = -log10(p_cutoff) + 0.5, label = "P-value >= 0.05", color = "black", hjust = 0) +
  geom_text_repel(data = subset(sig_genes, Significant),
                  aes(label = rownames(sig_genes)[Significant]),
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  max.overlaps = Inf,
                  nudge_x = 0.3,
                  nudge_y = 0.2)

print(volcano_plot) 
################################################################################


DEGs_up=sig_genes[sig_genes$logFC >= +1,]
DEGs_down=sig_genes[sig_genes$logFC<= -1,]
write.csv(DEGs_up, file = "DEGs_up.csv", row.names = FALSE)
write.csv(DEGs_down, file = "DEGs_down.csv", row.names = FALSE)
