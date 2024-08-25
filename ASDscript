library(GEOquery)
id_my <- "GSE77103"
gseasd <- getGEO(id_my)
gseasd <- gseasd[[1]]
gseasd
pData(gseasd)
fData(gseasd)
exprs(gseasd)    
summary(exprs(gseasd)) 
exprs(gseasd) <- log2(exprs(gseasd))
boxplot(exprs(gseasd),outline=FALSE)
library(dplyr)
sampleInfoASD <- pData(gseasd)
sampleInfoASD
sampleInfoASD <- select(sampleInfoASD, title,characteristics_ch1)
sampleInfoASD
sampleInfoASD <- rename(sampleInfoASD, Patient_Sample = title, Group = characteristics_ch1)
sampleInfoASD

library(pheatmap)
corMatrixasd <- cor(exprs(gseasd),use="c")
corMatrixasd
corMatrixasd[is.na(corMatrixasd)] <- 0

rownames(sampleInfoASD)
colnames(corMatrixasd)

rownames(sampleInfoASD) <- colnames(corMatrixasd)
pheatmap(corMatrixasd,
         annotation_col=sampleInfoASD)



library(ggplot2)
library(ggrepel)
any(is.na(exprs(gseasd)))     
any(is.infinite(exprs(gseasd)))
exprs_no_na_asd <- na.omit(exprs(gseasd))
inite_rows <- !apply(exprs_no_na_asd, 1, function(x) any(is.infinite(x)))
exprs_clean <- exprs_no_na_asd[inite_rows, ]
pcaASD <- prcomp(t(corMatrixasd))
cbind(sampleInfoASD, pcaASD$x) %>%
  ggplot(aes(x = PC1, y = PC2, col = Group, label = paste("Patient_Sample", Patient_Sample))) +
  geom_point(aes(shape = Group), size = 3) +  # Set point size
  geom_text_repel() +
  scale_color_manual(values = c("diagnosis: autism" = "blue", "diagnosis: healthy" = "red")) +
  scale_shape_manual(values = c("diagnosis: autism" = 17, "diagnosis: healthy" = 16))


library(ggplot2)
library(ggrepel)



library(limma)
colnames(sampleInfoASD)
designasd <- model.matrix(~0+sampleInfoASD$Group)
designasd
colnames(designasd) <- c("Autism","Healthy")
fitasd <- lmFit(exprs(gseasd), designasd)
head(fitasd$coefficients)
contrastsasd <- makeContrasts(Autism - Healthy, levels=designasd)
fit2asd <- contrasts.fit(fitasd, contrastsasd)
fit2asd <- eBayes(fit2asd)
topTable(fit2asd)
decideTests(fit2asd)
table(decideTests(fit2asd))
annoasd <- fData(gseasd)
annoasd
annoasd <- select(annoasd,"GENE_SYMBOL","ENSEMBL_ID")
fit2asd$genes <- annoasd
topTable(fit2asd)
full_results_asd <- topTable(fit2asd, number=Inf)
full_results_asd <- tibble::rownames_to_column(full_results_asd,"ID")



library(ggplot2)
library(ggrepel)

p_cutoff <- 0.05
fc_cutoff <- 1

sig_genes_asd <- full_results_asd %>%
  mutate(Significantasd = P.Value < p_cutoff & abs(logFC) > fc_cutoff,
         Direction = case_when(
           logFC > fc_cutoff ~ "up",
           logFC < -fc_cutoff ~ "down",
           TRUE ~ "not_significant"
         ))

volcano_plot_asd <- ggplot(sig_genes_asd, aes(x = logFC, y = -log10(P.Value), color = Direction)) +
  geom_point(alpha = 0.8, size = 1.2, shape = 16) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "not_significant" = "black"),
                     guide = guide_legend(override.aes = list(shape = 16, size = 5))) +
  xlab("log2(fold change)") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot") +
  theme_bw() +
  geom_hline(yintercept = c(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
  annotate("text", x = 5, y = -log10(p_cutoff) + 0.5, label = "Pvalue >= 0.05", color = "black", hjust = 0) +
  geom_text_repel(data = subset(sig_genes_asd, Significantasd),
                  aes(label = rownames(sig_genes_asd)[Significantasd]),
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "grey50",
                  max.overlaps = Inf,
                  nudge_x = 0.3,
                  nudge_y = 0.2)

print(volcano_plot_asd)
DEGs_up_asd=sig_genes_asd[sig_genes_asd$logFC >= +1,]
DEGs_down_asd=sig_genes_asd[sig_genes_asd$logFC<= -1,]
write.csv(DEGs_up_asd, file = "DEGs_up_asd.csv", row.names = FALSE)
write.csv(DEGs_down_asd, file = "DEGs_down_asd.csv", row.names = FALSE)
DEGs_up_asd = sig_genes_asd[sig_genes_asd$logFC > 1 & sig_genes_asd$P.Value >= 0.05, ]
DEGs_down_asd = sig_genes_asd[sig_genes_asd$logFC <= -1 & sig_genes_asd$P.Value >= 0.05, ]
