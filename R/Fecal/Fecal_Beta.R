library(devtools)
devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)

list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read_q2metadata("replicate-metadata-gnoto.txt")
str(metadata)
levels(metadata$`donor.type`)
colnames(metadata)[2] <- "sample.type"
colnames(metadata)[12] <- "donor.type"
str(metadata)

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)



####################### Bray-Curtis #######################
donor_type_colors <- c("gray", "lightblue", "green", "black")

bc_PCoA<-read_qza("fecal-core-metrics-results/bray_curtis_pcoa_results.qza")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

metadata$replicate <- factor(metadata$replicate)
str(metadata)
summary(metadata)

#########Change PCoA values
ggplot(bc_meta, aes(x=PC1, y=PC2, color=donor.type)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (24.90%)") +
  ylab("PC2 (17.12%)") +
  scale_color_manual(values=c("Black", "Blue", "Green", "Gray"), name = "Donor Type")

my_column <- "donor.type"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~replicate) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type")
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=3, width=4.5, device="tiff")

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "donor.type"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type")
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

bc_meta$replicate <- factor(bc_meta$replicate)

ggplot(bc_meta, aes(x=PC1, y=PC2, color=donor.type, shape=replicate)) +
  geom_point() + 
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type") +
  labs(shape = "Replicate")
ggsave(paste0("output/BC-ellipse_", my_column,"-replicate.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches



####################### Weighted UniFrac #######################
### From BC ###
#donor_type_colors <- c("Black", "Blue", "Green", "Gray")
#my_column <- "donor.type"

Wuni_PCoA<-read_qza("fecal-core-metrics-results/weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

Wuni_meta$replicate <- factor(Wuni_meta$replicate)

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= replicate), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type") +
  labs(shape = "Replicate")
ggsave(paste0("output/Wuni-ellipse_", my_column,"-replicate.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

bc_dist_mat<-read_qza("fecal-core-metrics-results/bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
### Not all values were true.
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"
#All values are true.

BC_PERMANOVA_out <- adonis2(bc_dm ~ donor.type, data = metadata_sub)

write.table(PERMANOVA_out,"output/bray_curtis_donor.type_Adonis_overall.csv",sep=",", row.names = TRUE) 



################### Unweighted UniFrac ######################
### From BC ###
#donor_type_colors <- c("Black", "Blue", "Green", "Gray")
#my_column <- "donor.type"

UWuni_PCoA<-read_qza("fecal-core-metrics-results/unweighted_unifrac_pcoa_results.qza")

UWuni_meta <- UWuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

UWuni_meta$replicate <- factor(UWuni_meta$replicate)

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),UWuni_meta,mean)

ggplot(UWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UWuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UWuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type")
ggsave(paste0("output/UWuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(UWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= replicate), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UWuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UWuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type") +
  labs(shape = "Replicate")
ggsave(paste0("output/UWuni-ellipse_", my_column,"-replicate.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches



############## Jaccard Plot #######################
### From BC ###
#donor_type_colors <- c("Black", "Blue", "Green", "Gray")
#my_column <- "donor.type"

Jaccard_PCoA<-read_qza("fecal-core-metrics-results/jaccard_pcoa_results.qza")

Jaccard_meta <- Jaccard_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

Jaccard_meta$replicate <- factor(Jaccard_meta$replicate)

# Now we are going to make an ordination plot
# Would need to change PCoA values
ggplot(Jaccard_meta, aes(x=PC1, y=PC2, color=donor.type)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=c("Blue", "Black", "Green", "Gray"), name = "Donor Type")

ggplot(Jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~replicate) +
  xlab(paste0("PC1 (", round(100*Jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type")
ggsave(paste0("output/Jaccard-basic_", my_column,".tiff"), height=3, width=4.5, device="tiff") #save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Jaccard_meta,mean)
colnames(centroids)[1] <- "donor.type"

ggplot(Jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type")
ggsave(paste0("output/Jaccard-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= replicate), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=donor_type_colors, name = "Donor Type") +
  labs(shape = "Replicate")
ggsave(paste0("output/Jaccard-ellipse_", my_column,"-replicate.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches



#####PERMANOVA#####
#BC done in previous step

##Weighted UniFrac##
wuni_dist_mat<-read_qza("fecal-core-metrics-results/weighted_unifrac_distance_matrix.qza")
wuni_dm <- as.matrix(wuni_dist_mat$data) 
rownames(wuni_dm) == metadata$SampleID ## all these values need to be "TRUE"
### Not all values were true.
metadata_sub <- metadata[match(rownames(wuni_dm),metadata$SampleID),]
rownames(wuni_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"
#All values are true.

Wuni_PERMANOVA_out <- adonis2(wuni_dm ~ donor.type, data = metadata_sub)

write.table(PERMANOVA_out,"output/weighted_unifrac_donor.type_Adonis_overall.csv",sep=",", row.names = TRUE) 
#Same results as BC. 