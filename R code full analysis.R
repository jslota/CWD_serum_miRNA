###Analysis of Conserved elk serum miRNAs
###2019-06-15
###Jessy Slota

library(ggplot2)
library(DESeq2)
library(reshape2)
library(VennDiagram)
library(svglite)

###
###Elk compare mouse and human annotation
###

#load elk miRNA data

#human annotation
elk_counts_hsa <- read.csv("raw data/Elk_serum_miRNA_raw_counts.csv", row.names = 1)
#(593 miRNAs via human annotation)

#mouse annotation
elk_counts_mmu <- read.csv("raw data/Elk_serum_miRNA_raw_counts_mouse_annotation.csv", row.names = 1)
#(515 miRNAs via mouse annotation)

#find matching miRNA
rownames(elk_counts_hsa) <- gsub("hsa-", "", rownames(elk_counts_hsa))
rownames(elk_counts_mmu) <- gsub("mmu-", "", rownames(elk_counts_mmu))

#common
intersect(rownames(elk_counts_hsa), rownames(elk_counts_mmu))
#367 matching miRNAs

#get matching miRNA
elk_counts_hsa <- elk_counts_hsa[intersect(rownames(elk_counts_hsa), rownames(elk_counts_mmu)),]
elk_counts_mmu <- elk_counts_mmu[intersect(rownames(elk_counts_hsa), rownames(elk_counts_mmu)),]

colnames(elk_counts_hsa) <- paste(colnames(elk_counts_hsa), "hsa", sep = "_")
colnames(elk_counts_mmu) <- paste(colnames(elk_counts_mmu), "mmu", sep = "_")

#make full dataset
elk_counts_full <- elk_counts_hsa
elk_counts_full[,71:140] <- elk_counts_mmu

elk_samples_full <- data.frame(Sample = colnames(elk_counts_full),
                               Annotation = c(rep("Human", 70), rep("Mouse", 70)),
                               row.names = colnames(elk_counts_full))

#get normalized counts
dds_elk <- DESeqDataSetFromMatrix(elk_counts_full, elk_samples_full, ~Annotation)
dds_elk <- DESeq(dds_elk)
rlog_elk <- rlog(dds_elk)
elk_counts_norm <- assay(rlog_elk)

#make expression dataset
elk_expression <- data.frame(Human = rowMeans(elk_counts_norm[,1:69]),
                             Mouse = rowMeans(elk_counts_norm[,70:138]),
                             row.names = rownames(elk_counts_norm))
elk_expression <- elk_expression[order(rowMeans(elk_expression), decreasing = TRUE),]
elk_expression$miRNA <- 1:367
elk_expression <- melt(elk_expression, id.vars = "miRNA")
colnames(elk_expression) <- c("miRNA", "Annotation", "rlog")

#make plot
ggplot(data = elk_expression, mapping = aes(x = miRNA, y = rlog, color = Annotation)) +
  geom_line() +
  scale_color_manual(name = "Annotation:",
                     values = c("cornflowerblue", "firebrick")) +
  ylab("Normalized Count") +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.key = element_rect(fill = "white", color = "black"),
        legend.margin = margin(-0.2,0,0,0,"cm"),
        legend.position = "bottom")

ggsave("figures/Elk_miRNA_exrpression.svg", width = 6, height = 3.5)

###
###Compare expression profile for elk, hamster, human and mouse serum miRNA
###

#load elk miRNA data
elk_counts <- read.csv("raw data/Elk_serum_miRNA_raw_counts.csv", row.names = 1)
head(elk_counts)
sum(colSums(elk_counts))
elk_samples <- read.csv("raw data/Elk_Sample_Data.csv", row.names = 1)
head(elk_samples)

#load hamster miRNA data
hamster_counts <- read.csv("raw data/Hamster_serum_miRNA_raw_counts.csv", row.names = 1)
head(hamster_counts)
sum(colSums(hamster_counts))
hamster_samples <- read.csv("raw data/Hamster_samples.csv")
head(hamster_samples)

#load mouse miRNA data
mouse_counts <- read.csv("raw data/HSA_Mouse_raw_miRNA_counts.csv", row.names = 1, check.names = FALSE)
head(mouse_counts)
mouse_samples <- read.csv("raw data/Mouse_sample_info.csv")
head(mouse_samples)

#make elk dataset
elk_dds <- DESeqDataSetFromMatrix(elk_counts, elk_samples, ~Condition)
elk_dds <- DESeq(elk_dds)
elk_rlog <- rlog(elk_dds)

#make pca plot for elk
pcaData <- plotPCA(elk_rlog, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#we see that CWD_59 is an outlier

#make hamster dataset
ham_dds <- DESeqDataSetFromMatrix(hamster_counts, hamster_samples, ~Condition)
ham_dds <- DESeq(ham_dds)
ham_rlog <- rlog(ham_dds)

#make pca plot for hamster
pcaData <- plotPCA(ham_rlog, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#we see that scrapie_6 is an outlier

#make mouse dataset
mmu_dds <- DESeqDataSetFromMatrix(mouse_counts, mouse_samples, ~Condition)
mmu_dds <- DESeq(mmu_dds)
mmu_rlog <- rlog(mmu_dds)

#make pca plot for mouse
pcaData <- plotPCA(mmu_rlog, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
#we see that ctrl-1 is an outlier

###Need to make a merged dataset to normalize the data together
#remove elk outlier (CWD_59)
elk_counts <- elk_counts[,colnames(elk_counts) != "CWD_59"]
colnames(elk_counts)
#remove hamster outlier (Scrapie_6)
hamster_counts <- hamster_counts[,colnames(hamster_counts) != "Scrapie_6"]
colnames(hamster_counts)
#remove mouse outlier (Ctrl_1)
mouse_counts <- mouse_counts[,colnames(mouse_counts) != "Ctrl_1"]
colnames(mouse_counts)

#load human
human_counts <- read.csv("raw data/Human_raw_counts_hsa.csv", row.names = 1)

#find common miRNA
#1 = elk
#2 = mouse
#3 = hamster
#4 = human
common_1_2 <- intersect(rownames(elk_counts), rownames(mouse_counts))
common_1_3 <- intersect(rownames(elk_counts), rownames(hamster_counts))
common_1_4 <- intersect(rownames(elk_counts), rownames(human_counts))
common_2_3 <- intersect(rownames(mouse_counts), rownames(hamster_counts))
common_2_4 <- intersect(rownames(mouse_counts), rownames(human_counts))
common_3_4 <- intersect(rownames(hamster_counts), rownames(human_counts))
common_1_2_3 <- intersect(common_1_2, rownames(hamster_counts))
common_1_2_4 <- intersect(common_1_2, rownames(human_counts))
common_1_3_4 <- intersect(common_1_3, rownames(human_counts))
common_2_3_4 <- intersect(common_2_3, rownames(human_counts))
common_1_2_3_4 <- intersect(common_1_2_3, common_2_3_4)

#make venn diagram
venn <- draw.quad.venn(area1 = nrow(elk_counts), area2 = nrow(mouse_counts),
                       area3 = nrow(hamster_counts), area4 = nrow(human_counts),
                       n12 = length(common_1_2), n13 = length(common_1_3), n14 = length(common_1_4),
                       n23 = length(common_2_3), n24 = length(common_2_4), n34 = length(common_3_4),
                       n123 = length(common_1_2_3), n124 = length(common_1_2_4),
                       n134 = length(common_1_3_4), n234 = length(common_2_3_4),
                       n1234 = length(common_1_2_3_4),
                       fill = c("Cornflowerblue", "firebrick", "darkgreen", "goldenrod"),
                       fontfamily = rep("sans", 15))

#make merged counts data
full_counts <- elk_counts[common_1_2_3_4,]
full_counts[,70:92] <- hamster_counts[common_1_2_3_4,]
full_counts[,93:107] <- mouse_counts[common_1_2_3_4,]
full_counts[,108:109] <- human_counts[common_1_2_3_4,]
full_samples <- data.frame(Sample = colnames(full_counts),
                           Species = c(rep("Elk", 69), rep("Hamster", 23), rep("Mouse", 15), rep("Human",2)))

#merged dataset
full_dds <- DESeqDataSetFromMatrix(full_counts, full_samples, ~Species)
full_dds <- DESeq(full_dds)
full_rlog <- rlog(full_dds)
full_counts <- assay(full_rlog)

#make expression profile dataset
species_expression_data <- data.frame(Elk = as.numeric(rowMeans(full_counts[,as.character(full_samples[full_samples$Species == "Elk",]$Sample)])),
                                      Hamster = as.numeric(rowMeans(full_counts[,as.character(full_samples[full_samples$Species == "Hamster",]$Sample)])),
                                      Mouse = as.numeric(rowMeans(full_counts[,as.character(full_samples[full_samples$Species == "Mouse",]$Sample)])),
                                      Human = as.numeric(rowMeans(full_counts[,as.character(full_samples[full_samples$Species == "Human",]$Sample)])),
                                      row.names = common_1_2_3_4)

species_expression_data <- species_expression_data[order(rowMeans(species_expression_data), decreasing = TRUE),]

#format for plotting
species_expression_data$miRNA <- 1:219
species_expression_data_melt <- melt(species_expression_data, id.vars = "miRNA")
colnames(species_expression_data_melt) <- c("miRNA", "Species", "rLog")

#make plot
ggplot(data = species_expression_data_melt, mapping = aes(x = miRNA, y = rLog, color = Species)) +
  geom_line() +
  scale_color_manual(name = "Species:",
                     values = c("cornflowerblue", "darkgreen", "firebrick", "goldenrod")) +
  ylab("Normalized Count") +
  guides(color = guide_legend(nrow = 2)) +
  theme(text = element_text(size = 10),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), 
        legend.margin = margin(-0.20,0,0,0, "cm"),
        legend.key = element_rect(fill = "white", color = "black"),
        legend.position = "bottom")

#save svg
ggsave("figures/Four_species_exression_profile.svg", width = 6, height = 3)


###Hemolysis in elk serum
###2019-06-17
###Jessy Slota

library(ggplot2)
library(svglite)

#load elk miRNA data
elk_counts <- read.csv("raw data/Elk_serum_miRNA_raw_counts.csv", row.names = 1)
elk_samples <- read.csv("raw data/Elk_Sample_Data.csv", row.names = 1)

#get miR-451a counts
elk_samples$`miR-451a` <- as.numeric(elk_counts["hsa-miR-451a",])

#get miR-23a-3p counts
elk_samples$`miR-23a` <- as.numeric(elk_counts["hsa-miR-23a-3p",])

#get miR-451a/miR-23a-3p ("Hemo ratios")
elk_samples$Hemo <- elk_samples$`miR-451a`/elk_samples$`miR-23a`

#save hemo data
write.csv(elk_samples, "sample info/elk_samples_hemo.csv")

###Make plot of hemolysis ratios
sample_order <- elk_samples[order(elk_samples$Hemo),]$Sample
ggplot(elk_samples, mapping = aes(x = Sample, y = Hemo)) +
  geom_bar(stat = "identity", fill = "cornflowerblue", color = "black") +
  geom_hline(yintercept = 5, color = "firebrick3") +
  scale_y_log10() +
  scale_x_discrete(limits = sample_order) +
  coord_flip() +
  #  ggtitle("Elk") +
  ylab("miR-451a/miR-23a-3p") +
  xlab("Sample") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text.y = element_text(size = 6))

ggsave("figures/Elk_hemolysis.svg", width = 4.5, height = 5)


###Analysis of miRNA biomarkers in elk serum
###2019-06-18
###Jessy Slota

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(svglite)

#Load elk data
elk_counts <- read.csv("raw data/Elk_serum_miRNA_raw_counts.csv", row.names = 1)
elk_samples <- read.csv("sample info/elk_samples_hemo.csv", row.names = 1)
rownames(elk_samples) <- elk_samples$Sample

summary(colnames(elk_counts) == elk_samples$Sample)

elk_samples$Hemolysis <- "Hemo"
#make non-hemolyzed (cleaned) dataset
clean_samples <- as.character(elk_samples[elk_samples$Hemo < 5,]$Sample)
elk_samples[clean_samples,]$Hemolysis <- "Non-Hemo"
elk_counts_clean <- elk_counts[,clean_samples]
elk_samples_clean <- elk_samples[clean_samples,]

summary(rownames(elk_samples_clean) == colnames(elk_counts_clean))
summary(elk_samples_clean$Condition)

###DE analysis full dataset
dds_full <- DESeqDataSetFromMatrix(elk_counts, elk_samples, ~Condition)
dds_full <- DESeq(dds_full)

res_full <- results(dds_full, c("Condition", "CWD", "Healthy"))
summary(res_full$padj < 0.05)

#DE miRNA in elk full dataset
DE_miRNA_full <- as.data.frame(res_full[order(res_full$padj),])
DE_miRNA_full <- na.omit(DE_miRNA_full)
DE_miRNA_full <- DE_miRNA_full[DE_miRNA_full$padj < 0.05,]
DE_miRNA_full <- DE_miRNA_full[abs(DE_miRNA_full$log2FoldChange) > 0.585,]
DE_miRNA_full <- DE_miRNA_full[DE_miRNA_full$baseMean > 25,]

write.csv(DE_miRNA_full, "results/Elk_DE_miRNA_full.csv")

###DE analysis clean dataset
dds_clean <- DESeqDataSetFromMatrix(elk_counts_clean, elk_samples_clean, ~Condition)
dds_clean <- DESeq(dds_clean)

res_clean <- results(dds_clean, c("Condition", "CWD", "Healthy"))
summary(res_clean$padj < 0.05)

#DE miRNA in clean dataset
DE_miRNA_clean <- as.data.frame(res_clean[order(res_clean$padj),])
DE_miRNA_clean <- na.omit(DE_miRNA_clean)
DE_miRNA_clean <- DE_miRNA_clean[DE_miRNA_clean$padj < 0.05,]
DE_miRNA_clean <- DE_miRNA_clean[abs(DE_miRNA_clean$log2FoldChange) > 0.585,]
DE_miRNA_clean <- DE_miRNA_clean[DE_miRNA_clean$baseMean > 25,]

write.csv(DE_miRNA_clean, "results/Res_elk_cleaned.csv")

#best biomarkers
best_elk_biomarkers <- c(intersect(rownames(DE_miRNA_full[DE_miRNA_full$log2FoldChange > 0,]),
                                   rownames(DE_miRNA_clean[DE_miRNA_clean$log2FoldChange > 0,])),
                         intersect(rownames(DE_miRNA_full[DE_miRNA_full$log2FoldChange < 0,]),
                                   rownames(DE_miRNA_clean[DE_miRNA_clean$log2FoldChange < 0,]))
)
best_elk_biomarkers

write.csv(DE_miRNA_clean[best_elk_biomarkers,], "results/Best_elk_biomarkers_cleaned.csv")

###
###Venn diagrams
###

#up-regulated venn diagram
#area1 = elk full
#area2 = elk clean
#note get rid of category labels prior to making final figure
venn <- draw.pairwise.venn(area1 = nrow(DE_miRNA_full[DE_miRNA_full$log2FoldChange > 0,]),
                           area2 = nrow(DE_miRNA_clean[DE_miRNA_clean$log2FoldChange > 0,]),
                           cross.area = length(intersect(rownames(DE_miRNA_full[DE_miRNA_full$log2FoldChange > 0,]),
                                                         rownames(DE_miRNA_clean[DE_miRNA_clean$log2FoldChange > 0,]))),
                           scaled = FALSE,
                           fill = c("cornflowerblue", "firebrick3"), fontfamily = "sans")

#save png
#ggsave("figures/Elk_up_venn.png", venn, width = 2.5, height = 2.5, dpi = 1200)
#save svg
ggsave("figures/Elk_up_venn.svg", venn, width = 2.5, height = 2.5)

#down-regulated venn diagram
#area1 = elk full
#area2 = elk clean
#note get rid of category labels prior to making final figure
venn <- draw.pairwise.venn(area1 = nrow(DE_miRNA_full[DE_miRNA_full$log2FoldChange < 0,]),
                           area2 = nrow(DE_miRNA_clean[DE_miRNA_clean$log2FoldChange < 0,]),
                           cross.area = length(intersect(rownames(DE_miRNA_full[DE_miRNA_full$log2FoldChange < 0,]),
                                                         rownames(DE_miRNA_clean[DE_miRNA_clean$log2FoldChange < 0,]))),
                           scaled = FALSE,rotation.degree = 180,
                           fill = c("cornflowerblue", "firebrick3"), fontfamily = "sans")

#save svg
ggsave("figures/Elk_down_venn.svg", venn, width = 2.5, height = 2.5)

###
###Heatmaps
###

#expression data
rlog_full <- rlog(dds_full)

###Heatmap of DE miRNA in full dataset

#full dataset#
heatmap_data <- assay(rlog_full[rownames(DE_miRNA_full),])
heatmap_data <- t((heatmap_data - rowMeans(heatmap_data)))

#annotation info
ann <- read.csv("sample info/Elk_full_annotation.csv", row.names = 1)
rownames(ann) <- gsub("CWD", "E", rownames(ann))
ann <- ann[rownames(heatmap_data),]
ann <- ann[,-c(2,4)]
ann$Genotype_132

#pick colors
ann_colors = list(
  Condition = c(CWD = "Firebrick3", Healthy = "navy"),
  Genotype_132 = c(Leu_Leu = "goldenrod", Met_Met = "forestgreen", Met_Leu = "coral3", na = "deepskyblue")
)

#Order samples based on CWD status and genotype
heatmap_order <- ann[ann$Condition == "CWD",]
heatmap_order <- heatmap_order[rownames(heatmap_order[order(heatmap_order$Genotype_132),]),]
heatmap_order2 <- ann[ann$Condition == "Healthy",]
heatmap_order2 <- heatmap_order2[rownames(heatmap_order2[order(heatmap_order2$Genotype_132),]),]
heatmap_order[36:70,] <- heatmap_order2
rm(heatmap_order2)

heatmap_data <- heatmap_data[rownames(heatmap_order),]

heat_map <- pheatmap(heatmap_data, annotation_row = ann, show_rownames = FALSE,
                     color = colorRampPalette(c("darkorchid4", "cornflowerblue", "goldenrod1","darkorange", "firebrick"))(50),
                     annotation_colors = ann_colors,
                     cluster_rows = FALSE,
                     border_color = "NA", annotation_legend = FALSE,
                     treeheight_col = 25, fontsize = 9)

#save svg
ggsave("figures/Elk_heatmap_full.svg", heat_map, width = 6.5, height = 5)

###Heatmap of best biomarkers

#full dataset#
heatmap_data <- assay(rlog_full[best_elk_biomarkers,])
heatmap_data <- t((heatmap_data - rowMeans(heatmap_data)))

#annotation info
ann <- read.csv("sample info/Elk_full_annotation.csv", row.names = 1)
rownames(ann) <- gsub("CWD", "E", rownames(ann))
ann <- ann[rownames(heatmap_data),]
ann <- ann[,-c(2,4)]
ann$Genotype_132

#pick colors
ann_colors = list(
  Condition = c(CWD = "Firebrick3", Healthy = "navy"),
  Genotype_132 = c(Leu_Leu = "goldenrod", Met_Met = "forestgreen", Met_Leu = "coral3", na = "deepskyblue")
)

#make heatmap
heatmap_data <- heatmap_data[rownames(heatmap_order),]
heat_map <- pheatmap(heatmap_data, annotation_row = ann, show_rownames = FALSE,
                     color = colorRampPalette(c("darkorchid4", "cornflowerblue", "goldenrod1","darkorange", "firebrick"))(50),
                     annotation_colors = ann_colors,
                     cluster_rows = FALSE,
                     border_color = "NA", annotation_legend = FALSE,
                     treeheight_col = 25, fontsize = 9)

#save svg
ggsave("figures/Elk_heatmap_best.svg", heat_map, width = 3.5, height = 4.5)

###
###PCA plots
###

#format expression data
rlog_full$Group <- paste(rlog_full$Condition, rlog_full$Hemolysis, sep = "-")

#DE in full dataset
pca_data <- plotPCA(rlog_full[rownames(DE_miRNA_full),], intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

summary(elk_ref_numbers$Sample == pca_data$name)
pca_data$name <- elk_ref_numbers$Carcass.number

ggplot(pca_data, aes(x = PC1, y = PC2, label = name)) +
  geom_point(aes(color = Group), size = 1) +
  #  geom_text(size = 3) +
  scale_color_manual(name = "", values = c("firebrick3", "indianred1", "navy", "cornflowerblue")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  guides(color = guide_legend(nrow = 2)) +
  theme(text = element_text(size = 10),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(fill = "NA", color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.margin = margin(-0.3,0,0,0, "cm"))

#save svg
ggsave("figures/PCA_elk_full.svg", width = 3, height = 3)

#Best biomarkers
pca_data <- plotPCA(rlog_full[best_elk_biomarkers,], intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, label = name)) +
  geom_point(aes(color = Group), size = 1) +
  #  geom_text(size = 3) +
  scale_color_manual(name = "", values = c("firebrick3", "indianred1", "navy", "cornflowerblue")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  guides(color = guide_legend(nrow = 2)) +
  theme(text = element_text(size = 10),
        panel.background = element_rect(color = "black", fill = "white"),
        panel.border = element_rect(fill = "NA", color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.key = element_rect(fill = "white"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.margin = margin(-0.3,0,0,0, "cm"))

#save svg
ggsave("figures/PCA_elk_best.svg", width = 3, height = 3)


###Analysis of biomarkers in hamster serum
###2019-06-18
###Jessy Slota

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(svglite)

#Load hamster data
hamster_counts <- read.csv("raw data/Hamster_serum_miRNA_raw_counts.csv", row.names = 1)
hamster_samples <- read.csv("sample info/hamster_samples_hemo.csv", row.names = 1)
rownames(hamster_samples) <- hamster_samples$Sample

summary(colnames(hamster_counts) == hamster_samples$Sample)

#DE analysis
dds_hamster <- DESeqDataSetFromMatrix(hamster_counts, hamster_samples, ~Condition)
dds_hamster <- DESeq(dds_hamster)
res_hamster <- results(dds_hamster, c("Condition", "Scrapie", "PBS"))

#DE miRNA in hamster
DE_miRNA_hamster <- as.data.frame(res_hamster[order(res_hamster$padj),])
DE_miRNA_hamster <- DE_miRNA_hamster[DE_miRNA_hamster$padj < 0.05,]
DE_miRNA_hamster <- DE_miRNA_hamster[abs(DE_miRNA_hamster$log2FoldChange) > 0.585,]
DE_miRNA_hamster <- na.omit(DE_miRNA_hamster)

#save hamster results
write.csv(DE_miRNA_hamster, "results/Hamster_results.csv")

###Compare results
###compare with full elk DE miRNA
#load data
Elk_full <- read.csv("results/Elk_DE_miRNA_full.csv", row.names = 1)

#10 match
intersect(rownames(DE_miRNA_hamster), rownames(Elk_full))

###Compare with best biomarkers
#load elk results for comparison
Elk_biomarkers <- read.csv("results/Best_elk_biomarkers_cleaned.csv", row.names = 1)

#compare full list
intersect(rownames(DE_miRNA_hamster), rownames(Elk_biomarkers))
#7 match

#compare taking FC into account
c(intersect(rownames(Elk_biomarkers[Elk_biomarkers$log2FoldChange > 0,]),
            rownames(DE_miRNA_hamster[DE_miRNA_hamster$log2FoldChange > 0,])),
  intersect(rownames(Elk_biomarkers[Elk_biomarkers$log2FoldChange < 0,]),
            rownames(DE_miRNA_hamster[DE_miRNA_hamster$log2FoldChange < 0,])))
#6 match

###
###Venn diagram of matching elk and hamster DE miRNA
###

#up-regulated venn diagram
#area1 = Elk
#area2 = Hamster
#note get rid of category labels prior to making final figure
plot.new()
venn <- draw.pairwise.venn(area1 = nrow(Elk_biomarkers[Elk_biomarkers$log2FoldChange > 0,]),
                           area2 = nrow(DE_miRNA_hamster[DE_miRNA_hamster$log2FoldChange > 0,]),
                           cross.area = length(intersect(rownames(Elk_biomarkers[Elk_biomarkers$log2FoldChange > 0,]),
                                                         rownames(DE_miRNA_hamster[DE_miRNA_hamster$log2FoldChange > 0,]))),
                           scaled = FALSE,
                           fill = c("cornflowerblue", "firebrick3"), fontfamily = "sans")
#save as png
#ggsave("Figures/Venn_elk_ham_up.png", venn, width = 2.5, height = 2.5, dpi = 1200)
#save as svg
ggsave("Figures/Venn_elk_ham_up.svg", venn, width = 2.5, height = 2.5)

#down-regulated venn diagram
#area1 = elk full
#area2 = elk clean
#note get rid of category labels prior to making final figure
plot.new()
venn <- draw.pairwise.venn(area1 = nrow(Elk_biomarkers[Elk_biomarkers$log2FoldChange < 0,]),
                           area2 = nrow(DE_miRNA_hamster[DE_miRNA_hamster$log2FoldChange < 0,]),
                           cross.area = length(intersect(rownames(Elk_biomarkers[Elk_biomarkers$log2FoldChange < 0,]),
                                                         rownames(DE_miRNA_hamster[DE_miRNA_hamster$log2FoldChange < 0,]))),
                           scaled = FALSE,
                           fill = c("cornflowerblue", "firebrick3"), fontfamily = "sans")
#save as png
#ggsave("Figures/Venn_elk_ham_down.png", venn, width = 2.5, height = 2.5, dpi = 1200)
#save as svg
ggsave("Figures/Venn_elk_ham_down.svg", venn, width = 2.5, height = 2.5)

###
###heatmap
###

rlog_hamster <- rlog(dds_hamster)

heatmap_data <- assay(rlog_hamster[rownames(DE_miRNA_hamster),])
heatmap_data <- t(heatmap_data - rowMeans(heatmap_data))

#annotation
ann <- as.data.frame(hamster_samples$Condition)
rownames(ann) <- hamster_samples$Sample
colnames(ann) <- "Condition"

#pick colors
ann_colors = list(
  Condition = c(Scrapie = "Firebrick3", PBS = "navy")
)

#draw heatmap
heat_map <- pheatmap(heatmap_data, annotation_row = ann,
                     color = colorRampPalette(c("darkorchid4", "cornflowerblue", "goldenrod1","darkorange", "firebrick"))(50),
                     border_color = "NA", annotation_legend = FALSE,
                     annotation_colors = ann_colors,
                     treeheight_row = 25, treeheight_col = 25,
                     fontsize = 9, show_rownames = FALSE
)

#save as svg
ggsave("figures/Hamster_heatmap.svg", heat_map, width = 5, height = 6)

###Diagnostic potential of biomarkers in elk
###(Logistic regression modeling)
###2019-06-19
####Jessy Slota

library(ROCR)
library(ggplot2)
library(DESeq2)
library(svglite)
library(RColorBrewer)
library(glmnet)
library(caret)

#load results to get biomarker lists
elk_biomarkers <- read.csv("results/Best_elk_biomarkers_cleaned.csv", row.names = 1)
hamster_biomarkers <- read.csv("results/Hamster_results.csv", row.names = 1)

biomarkers_elk <- rownames(elk_biomarkers)
biomarkers_matching <- c(intersect(rownames(elk_biomarkers[elk_biomarkers$log2FoldChange > 0,]),
                                   rownames(hamster_biomarkers[hamster_biomarkers$log2FoldChange > 0,])),
                         intersect(rownames(elk_biomarkers[elk_biomarkers$log2FoldChange < 0,]),
                                   rownames(hamster_biomarkers[hamster_biomarkers$log2FoldChange < 0,])))

#get raw data
elk_counts <- read.csv("raw data/Elk_serum_miRNA_raw_counts.csv", row.names = 1)
elk_samples <- read.csv("raw data/Elk_Sample_Data.csv", row.names = 1)

#get normalized data
dds_elk <- DESeqDataSetFromMatrix(elk_counts, elk_samples, ~Condition)
dds_elk <- DESeq(dds_elk)
rlog_elk <- rlog(dds_elk)

#Classify based on CWD
rownames(elk_samples) <- elk_samples$Sample
CWD_elk <- as.character(elk_samples[elk_samples$Condition == "CWD",]$Sample)
elk_samples$Label <- 0
elk_samples[CWD_elk,]$Label <- 1
elk_samples$Label <- as.factor(elk_samples$Label)

###
###ROC of individual miRNAs
###

expr_data <- assay(rlog_elk)[biomarkers_elk,]
expr_data <- as.data.frame(t(expr_data))
summary(rownames(expr_data) == rownames(elk_samples))

#set up dataset
pr <- prediction(expr_data[,1], elk_samples$Label)
perf <- performance(pr, measure = "tpr", x.measure = "fpr")
values <- data.frame(FPR = unlist(slot(perf, "x.values")), TPR = unlist(slot(perf, "y.values")))
values$miRNA <- colnames(expr_data)[1]
ROC_data <- values

AUC_values <- data.frame(miRNA = colnames(expr_data),
                         AUC = 0)

AUC_values[1,2] <- unlist(slot(performance(pr, "auc"), "y.values"))


#remaining upregulated miRNA
for (i in 2:7) {
  pr <- prediction(expr_data[,i], elk_samples$Label)
  perf <- performance(pr, measure = "tpr", x.measure = "fpr")
  values <- data.frame(FPR = unlist(slot(perf, "x.values")), TPR = unlist(slot(perf, "y.values")))
  values$miRNA <- colnames(expr_data)[i]
  ROC_data[(nrow(ROC_data)+1):(nrow(ROC_data)+nrow(values)),] <- values
  AUC_values[i,2] <- unlist(slot(performance(pr, "auc"), "y.values"))
}

#downregulated miRNA
for (i in 8:21) {
  pr <- prediction(-expr_data[,i], elk_samples$Label)
  perf <- performance(pr, measure = "tpr", x.measure = "fpr")
  values <- data.frame(FPR = unlist(slot(perf, "x.values")), TPR = unlist(slot(perf, "y.values")))
  values$miRNA <- colnames(expr_data)[i]
  ROC_data[(nrow(ROC_data)+1):(nrow(ROC_data)+nrow(values)),] <- values
  AUC_values[i,2] <- unlist(slot(performance(pr, "auc"), "y.values"))
}

AUC_values <- AUC_values[order(AUC_values$AUC, decreasing = TRUE),]

write.csv(AUC_values, "results/AUC_individual_miRNAs.csv", row.names = FALSE)

#Get top 10 and 5 for alternate model
top_10_miRNA <- as.character(AUC_values[1:10,]$miRNA)
top_5_miRNA <- as.character(AUC_values[1:5,]$miRNA)

#get colors
cols <- colorRampPalette(brewer.pal(11, "Spectral"))(21)

#ROC plot
ggplot(ROC_data, aes(x = FPR, y = TPR, color = miRNA)) +
  geom_line(size = 1) +
  geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  scale_color_manual(values = cols) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 7))

#save svg
ggsave("figures/ROC_plot_miRNAs.svg", width = 3, height = 3)

###
###Lasso Penalized Logistic Regression
###

#elk biomarkers
#full 21-miRNA signature
LR_data <- assay(rlog_elk)[biomarkers_elk,]
LR_data <- as.data.frame(t(LR_data))
summary(row.names(LR_data) == row.names(elk_samples))
LR_data$Prions <- elk_samples$Label

#format data
x <- model.matrix(Prions~., LR_data)
y <- LR_data$Prions

#Find best lambda using cross-validation
set.seed(123) 
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)

LR_model_min <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)

coef(LR_model_min) #keep 6...miR-186-5p, miR-500a-3p, miR-99a-5p, miR-125b-5p, miR-181c-5p, miR-181d-5p

p <- predict(LR_model_min, newx = x, type = "response")
pr <- prediction(p, y)
perf <- performance(pr, measure = "tpr", x.measure = "fpr")

#confusion matrix
results <- factor(ifelse(p > 0.5,1,0), levels = c("0","1"))
answers <- factor(y, levels = c("0","1"))
confusionMatrix(results, answers)

#          Reference
#Prediction  0  1
#         0 29  4
#         1  6 31

auc <- performance(pr, "auc")
auc <- unlist(slot(auc, "y.values"))#0.922

temp <- data.frame(FPR = unlist(slot(perf, "x.values")), TPR = unlist(slot(perf, "y.values")))
temp$model <- "lambda_min"
values <- temp

#now use lambda.1se
LR_model_1se <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.1se)

coef(LR_model_1se) #keep 4... miR-30e-3p, miR-500a-3p, miR-125b-5p, miR-181c-5p
summary(LR_model_1se)

p <- predict(LR_model_1se, newx = x, type = "response")
pr <- prediction(p, y)
perf <- performance(pr, measure = "tpr", x.measure = "fpr")

#confusion matrix
results <- factor(ifelse(p > 0.5,1,0), levels = c("0","1"))
answers <- factor(y, levels = c("0","1"))
confusionMatrix(results, answers)

#          Reference
#Prediction  0  1
#         0 29  5
#         1  6 30

auc <- performance(pr, "auc")
auc <- unlist(slot(auc, "y.values"))#0.915

temp <- data.frame(FPR = unlist(slot(perf, "x.values")), TPR = unlist(slot(perf, "y.values")))
temp$model <- "lambda_1se"
values[72:142,] <- temp

#ROC plot
ggplot(values, aes(x = FPR, y = TPR, color = model)) +
  geom_line(size = 1) +
  geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  scale_color_manual(values = c("cornflowerblue", "forestgreen")) +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        legend.key = element_rect(fill = "white"))

ggsave("figures/ROC_plot_RL_models.svg", width = 3, height = 3)

#Matching biomarkers
#6-miRNA signature
LR_data_6 <- assay(rlog_elk)[biomarkers_matching,]
LR_data_6 <- as.data.frame(t(LR_data_6))
summary(row.names(LR_data_6) == row.names(elk_samples))
LR_data_6$Prions <- elk_samples$Label

#format data
x <- model.matrix(Prions~., LR_data_6)
y <- LR_data$Prions

#Find best lambda using cross-validation
set.seed(123) 
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")
plot(cv.lasso)

#set lamda to 0 to keep all miRNAs
LR_model_six <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)

coef(LR_model_six) #keep 5...hsa-miR-103a-3p, hsa-miR-107, hsa-miR-99a-5p, hsa-miR-125b-5p, hsa-miR-100-5p

p <- predict(LR_model_six, newx = x, type = "response")
pr <- prediction(p, y)
perf <- performance(pr, measure = "tpr", x.measure = "fpr")

#confusion matrix
results <- factor(ifelse(p > 0.5,1,0), levels = c("0","1"))
answers <- factor(y, levels = c("0","1"))
confusionMatrix(results, answers)

#          Reference
#Prediction  0  1
#         0 31  6
#         1  4 29

auc <- performance(pr, "auc")
auc <- unlist(slot(auc, "y.values"))#0.879

values <- data.frame(FPR = unlist(slot(perf, "x.values")), TPR = unlist(slot(perf, "y.values")))

#ROC plot
ggplot(values, aes(x = FPR, y = TPR)) +
  geom_line(size = 1, color = "Black") +
  geom_abline(aes(slope = 1, intercept = 0), color = "red") +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"))

ggsave("figures/ROC_plot_6-miRNA.svg", width = 3, height = 3)

###
###Custom function for repeated cross validation
###
cross_validation_fun <- function(data, fold, k, l) {
  
  # False positive rate
  fpr <- NULL
  # False negative rate
  fnr <- NULL
  # Accuracy
  acc <- NULL
  #AUC
  auc <- NULL
  
  set.seed(123)
  for (i in 1:k) {
    # Train-test splitting
    # 80% of samples -> fitting
    # 20% of samples -> testing
    #i.e. 5-fold cross validation
    smp_size <- floor((1-1/fold) * nrow(data))
    index <- createDataPartition(y = data$Prions, p = (1-1/fold), list = FALSE)
    train <- data[index, ]
    test <- data[-index, ]
    
    #format for glmnet
    x <- model.matrix(Prions~., train)
    y <- train$Prions
    
    # Fitting
    model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = l)
    
    # Predict results
    p <- predict(model, newx = model.matrix(Prions~., test), type = "response")
    
    #auc
    pr <- prediction(p, test$Prions)
    a <- performance(pr, "auc")
    auc[i] <- unlist(slot(a, "y.values"))
    
    # If prob > 0.5 then 1, else 0
    results <- factor(ifelse(p > 0.5,1,0), levels = c("0","1"))
    
    # Actual answers
    answers <- factor(test$Prions, levels = c("0","1"))
    
    # Accuracy calculation
    misClasificError <- mean(answers != results)
    
    # Collecting results
    acc[i] <- 1-misClasificError
    
    # Confusion matrix
    cm <- confusionMatrix(data=results, reference=answers)
    fpr[i] <- cm$table[2]/(nrow(LR_data)-smp_size)
    fnr[i] <- cm$table[3]/(nrow(LR_data)-smp_size)
  }
  
  res <- data.frame(AUC = auc,
                    Accuracy = acc,
                    FPR = fpr,
                    FNR = fnr)
  return(res)
}

###
###Perform 5-fold cross validation on each model
###

#Find best lambda for 21 miRNA models
x <- model.matrix(Prions~., LR_data)
y <- LR_data$Prions
set.seed(123) 
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

#21 miRNAs minimum lambda
LR_model_min_res <- cross_validation_fun(LR_data, 5, 500, cv.lasso$lambda.min)
summary(LR_model_min_res)

#Get confidence interval
error <- qt(0.95,df=length(LR_model_min_res$AUC)-1)*sd(LR_model_min_res$AUC)/sqrt(length(LR_model_min_res$AUC))
mean(LR_model_min_res$AUC) - error
mean(LR_model_min_res$AUC) + error
#0.851(0.844-0.858)

#21 miRNAs 1se lambda
LR_model_1se_res <- cross_validation_fun(LR_data, 5, 500, cv.lasso$lambda.1se)
summary(LR_model_1se_res)

#Get confidence interval
error <- qt(0.95,df=length(LR_model_1se_res$AUC)-1)*sd(LR_model_1se_res$AUC)/sqrt(length(LR_model_1se_res$AUC))
mean(LR_model_1se_res$AUC) - error
mean(LR_model_1se_res$AUC) + error
#0.793(0.785-0.801)

#Find best lambda for 6 miRNA model
x <- model.matrix(Prions~., LR_data_6)
y <- LR_data$Prions
set.seed(123) 
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial")

#6 miRNAs minimum lambda
LR_model_six_res <- cross_validation_fun(LR_data_6, 5, 500, cv.lasso$lambda.min)
summary(LR_model_six_res)

#Get confidence interval
error <- qt(0.95,df=length(LR_model_six_res$AUC)-1)*sd(LR_model_six_res$AUC)/sqrt(length(LR_model_six_res$AUC))
mean(LR_model_six_res$AUC) - error
mean(LR_model_six_res$AUC) + error
#0.867(0.859-0.872)

###miRNA target prediction
###2019-06-19
###Jessy Slota

#load miRDB targets
miRDP_miR_99a <- read.csv("targets/miRDB_miR-99a.csv")
miRDP_miR_100 <- read.csv("targets/miRDB_miR-100.csv")
miRDP_miR_103a <- read.csv("targets/miRDB_miR-103a.csv")
miRDP_miR_107 <- read.csv("targets/miRDB_miR-107.csv")
miRDP_miR_125a <- read.csv("targets/miRDB_miR-125a.csv")
miRDP_miR_125b <- read.csv("targets/miRDB_miR-125b.csv")

#load targetscan targets
targetscan_miR_99_100 <- read.delim("targets/TargetScan7.2__miR-99-5p_100-5p.predicted_targets.txt")
targetscan_miR_103_107 <- read.delim("targets/TargetScan7.2__miR-103-3p_107.predicted_targets.txt")
targetscan_miR_125a_b <- read.delim("targets/TargetScan7.2__miR-125ab-5p.predicted_targets.txt")

#load miRWalk targets
miRWalk_99a <- read.csv("targets/miRWalk_miRNA_Targets_99a.csv")
miRWalk_100 <- read.csv("targets/miRWalk_miRNA_Targets_100.csv")
miRWalk_103a <- read.csv("targets/miRWalk_miRNA_Targets_103a.csv")
miRWalk_107 <- read.csv("targets/miRWalk_miRNA_Targets_107.csv")
miRWalk_125a <- read.csv("targets/miRWalk_miRNA_Targets_125a.csv")
miRWalk_125b <- read.csv("targets/miRWalk_miRNA_Targets_125b.csv")


#tarbase v8
#Tarbase <- read.delim("targets/TarBase_v8_download.txt")
#head(Tarbase)
#tarbase_99a <- Tarbase[Tarbase$mirna == "hsa-miR-99a-5p",]
#head(tarbase_99a)
#tarbase_100 <- Tarbase[Tarbase$mirna == "hsa-miR-100-5p",]
#head(tarbase_100)
#tarbase_103a <- Tarbase[Tarbase$mirna == "hsa-miR-103a-3p",]
#head(tarbase_103a)
#tarbase_107 <- Tarbase[Tarbase$mirna == "hsa-miR-107",]
#head(tarbase_107)
#tarbase_125a <- Tarbase[Tarbase$mirna == "hsa-miR-125a-5p",]
#head(tarbase_125a)
#tarbase_125b <- Tarbase[Tarbase$mirna == "hsa-miR-125b-5p",]
#head(tarbase_125b)
#Save files
#write.csv(tarbase_99a, "targets/Tarbase_targets_99a.csv")
#write.csv(tarbase_100, "targets/Tarbase_targets_100.csv")
#write.csv(tarbase_103a, "targets/Tarbase_targets_103a.csv")
#write.csv(tarbase_107, "targets/Tarbase_targets_107.csv")
#write.csv(tarbase_125a, "targets/Tarbase_targets_125a.csv")
#write.csv(tarbase_125b, "targets/Tarbase_targets_125b.csv")

#load Tarbase targets
tarbase_99a <- read.csv("targets/Tarbase_targets_99a.csv", row.names = 1)
tarbase_100 <- read.csv("targets/Tarbase_targets_100.csv", row.names = 1)
tarbase_103a <- read.csv("targets/Tarbase_targets_103a.csv", row.names = 1)
tarbase_107 <- read.csv("targets/Tarbase_targets_107.csv", row.names = 1)
tarbase_125a <- read.csv("targets/Tarbase_targets_125a.csv", row.names = 1)
tarbase_125b <- read.csv("targets/Tarbase_targets_125b.csv", row.names = 1)

#intersect
targets_99a <- intersect(miRDP_miR_99a$Gene.Symbol, targetscan_miR_99_100$Target.gene)
targets_99a <- intersect(targets_99a, unique(miRWalk_99a$genesymbol))
targets_99a <- intersect(targets_99a, unique(tarbase_99a$geneName))

targets_100 <- intersect(miRDP_miR_100$Gene.Symbol, targetscan_miR_99_100$Target.gene)
targets_100 <- intersect(targets_100, unique(miRWalk_100$genesymbol))
targets_100 <- intersect(targets_100, unique(tarbase_100$geneName))

targets_103a <- intersect(miRDP_miR_103a$Gene.Symbol, targetscan_miR_103_107$Target.gene)
targets_103a <- intersect(targets_103a, unique(miRWalk_100$genesymbol))
targets_103a <- intersect(targets_103a, unique(tarbase_103a$geneName))

targets_107 <- intersect(miRDP_miR_107$Gene.Symbol, targetscan_miR_103_107$Target.gene)
targets_107 <- intersect(targets_107, unique(miRWalk_107$genesymbol))
targets_107 <- intersect(targets_107, unique(tarbase_107$geneName))

targets_125a <- intersect(miRDP_miR_125a$Gene.Symbol, targetscan_miR_125a_b$Target.gene)
targets_125a <- intersect(targets_125a, unique(miRWalk_125a$genesymbol))
targets_125a <- intersect(targets_125a, unique(tarbase_125a$geneName))

targets_125b <- intersect(miRDP_miR_125b$Gene.Symbol, targetscan_miR_125a_b$Target.gene)
targets_125b <- intersect(targets_125b, unique(miRWalk_125b$genesymbol))
targets_125b <- intersect(targets_125b, unique(tarbase_125b$geneName))

#very similar targets between conserved miRNA
intersect(targets_99a, targets_100)
intersect(targets_103a, targets_107)
intersect(targets_125a, targets_125b)

#target lists for miRNA clusters
targets_99_100 <- unique(c(targets_99a, targets_100))
targets_103_107 <- unique(c(targets_103a, targets_107))
targets_125_a_b <- unique(c(targets_125a, targets_125b))

#save target lists
write.csv(targets_99_100, "targets/Targets_99_100.csv", quote = FALSE)
write.csv(targets_103_107, "targets/Targets_107_103.csv", quote = FALSE)
write.csv(targets_125_a_b, "targets/Targets_125_a_b.csv", quote = FALSE)

#common targets
intersect(targets_99_100, targets_103_107)
intersect(targets_99_100, targets_125_a_b)
intersect(targets_103_107, targets_125_a_b)
common_targets <- unique(c(intersect(targets_99_100, targets_103_107),intersect(targets_99_100, targets_125_a_b),intersect(targets_103_107, targets_125_a_b)))

write.csv(common_targets, "targets/common.csv", quote = FALSE)

#table for cytoscape
targets <- data.frame(miRNA = c(rep("miR-99a/100", length(targets_99_100)),
                                rep("miR-103a/107", length(targets_103_107)),
                                rep("miR-125a/b", length(targets_125_a_b))),
                      target = c(targets_99_100, targets_103_107, targets_125_a_b))

write.csv(targets, "cytoscape/miRNA_targets_combined.csv", row.names = FALSE)

#individual miRNA targets
targets <- data.frame(miRNA = c(rep("miR-99a-5p", length(targets_99a)),
                                rep("miR-100-5p", length(targets_100)),
                                rep("miR-103a-3p", length(targets_103a)),
                                rep("miR-107", length(targets_107)),
                                rep("miR-125a-5p", length(targets_125a)),
                                rep("miR-125b-5p", length(targets_125b))),
                      target = c(targets_99a, targets_100, targets_103a,
                                 targets_107, targets_125a, targets_125b))

write.csv(targets, "cytoscape/miRNA_targets_combined_individual.csv", row.names = FALSE)

#miR-125a/b targets
targets_miR_125_a_b <- setdiff(unique(c(targets_125a, targets_125b)), common_targets)
targets_miR_103_107 <- setdiff(unique(c(targets_103a, targets_107)), common_targets)
targets_miR_99a_100 <- setdiff(unique(c(targets_99a, targets_100)), common_targets)

#gene lists for functional annotation
write.csv(targets_miR_125_a_b, "targets/Cluster_125_a_b_targets.csv", quote = FALSE)
write.csv(targets_miR_103_107, "targets/Cluster_103_107_targets.csv", quote = FALSE)
write.csv(targets_miR_99a_100, "targets/Cluster_99a_100_targets.csv", quote = FALSE)
write.csv(common_targets, "targets/Cluster_multi_targets.csv", quote = FALSE)

#all targets 
all_targets <- unique(targets$target)

write.csv(all_targets, "targets/all_targets.csv", quote = FALSE)