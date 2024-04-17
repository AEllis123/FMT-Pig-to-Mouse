library(tidyverse)
library(qiime2R)
library(ggpubr)

list.files()

meta<-read_q2metadata("replicate-metadata-gnoto.txt") 
str(meta)
colnames(meta)[2] <- "sample.type"
colnames(meta)[12] <- "donor.type"
str(meta)

evenness = read_qza("cecal-core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("cecal-core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("cecal-core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("cecal-core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
colnames(faith_pd)[1] <- "Number"
colnames(faith_pd)[2] <- "SampleID"
colnames(faith_pd)[3] <- "faith_pd"
faith_pd <- faith_pd[, -which(names(faith_pd) == "Number")]

str(meta)
str(observed_features)

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
#meta = meta[,-1]
str(meta)

#Plots
hist(meta$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta$shannon_entropy, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$pielou_e, title = "Evenness")
ggqqplot(meta$observed_features, title = "Observed Features")

library("ggpubr")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
#Data is normally distributed.
shapiro.test(meta$faith_pd)
#Data is basically normally distributed. 
shapiro.test(meta$pielou_e)
#Data is normally distributed.
shapiro.test(meta$observed_features)
#Data is normally distributed. 


###################### Pielou ###########################
#Run the ANOVA and save it as an object
aov.evenness.donor.type = aov(pielou_evenness ~ donor.type, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.evenness.donor.type)

TukeyHSD(aov.evenness.donor.type)

levels(meta$donor.type)
#Re-order the groups because the default is alphabetical order
meta$donor.type.ord = factor(meta$donor.type, c("control", "donor-piglet", "donor-sow", "Piglet", "Sow"))
levels(meta$donor.type.ord)

#Plot
boxplot(pielou_evenness ~ donor.type.ord, data=meta, ylab="Pielou Evenness")

evenness_boxplot <- ggplot(meta, aes(donor.type.ord, pielou_evenness, fill = donor.type.ord)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y="Pielou's Evenness", x = "", fill="Donor Type") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/evenness_boxplot.png", evenness_boxplot, height = 3, width = 3)

evenness_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(donor.type.ord) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

evenness_se <- ggplot(evenness_summary, aes(donor.type.ord, mean_evenness, fill = donor.type.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's Evenness  ± s.e.", x = "") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/evenness_se.png", evenness_se, height = 2.5, width = 3)



###################### Faith PD ###########################
#Run the ANOVA and save it as an object
aov.faith.pd.donor.type = aov(faith_pd ~ donor.type, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.faith.pd.donor.type)

TukeyHSD(aov.faith.pd.donor.type)

#Plot
boxplot(faith_pd ~ donor.type.ord, data=meta, ylab="Faith Phylogenetic Diversity")

faith_pd_boxplot <- ggplot(meta, aes(donor.type.ord, faith_pd, fill = donor.type.ord)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y="Faith Phylogenetic Diversity", x = "", fill = "Donor Type") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/faith_pd_boxplot.png", faith_pd_boxplot, height = 3, width = 3)

faith_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(donor.type.ord) %>%   # the grouping variable
  summarise(mean_faith = mean(faith_pd),  # calculates the mean of each group
            sd_faith = sd(faith_pd), # calculates the standard deviation of each group
            n_faith = n(),  # calculates the sample size per group
            se_faith = sd(faith_pd)/sqrt(n())) # calculates the standard error of each group

faith_se <- ggplot(faith_summary, aes(donor.type.ord, mean_faith, fill = donor.type.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_faith - se_faith, ymax = mean_faith + se_faith), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith's Phylogenetic Diversity  ± s.e.", x = "") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/faith_se.png", faith_se, height = 2.5, width = 3)



###################### Observed Features ###########################
#Run the ANOVA and save it as an object
aov.obs.feat.donor.type = aov(observed_features ~ donor.type, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.obs.feat.donor.type)

TukeyHSD(aov.obs.feat.donor.type)

#Plot
boxplot(observed_features ~ donor.type.ord, data=meta, ylab="Observed Features")

obs_feat_boxplot <- ggplot(meta, aes(donor.type.ord, observed_features, fill = donor.type.ord)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y="Observed Features", x = "", fill = "Donor Type") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/obs_feat_boxplot.png", obs_feat_boxplot, height = 3, width = 3)

obs_feature_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(donor.type.ord) %>%   # the grouping variable
  summarise(mean_obs_feature = mean(observed_features),  # calculates the mean of each group
            sd_obs_feature = sd(observed_features), # calculates the standard deviation of each group
            n_obs_feature = n(),  # calculates the sample size per group
            se_obs_feature = sd(observed_features)/sqrt(n())) # calculates the standard error of each group

obs_features_se <- ggplot(obs_feature_summary, aes(donor.type.ord, mean_obs_feature, fill = donor.type.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_obs_feature - se_obs_feature, ymax = mean_obs_feature + se_obs_feature), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features  ± s.e.", x = "") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/obs_features_se.png", obs_features_se, height = 2.5, width = 3)



###################### Shannon ###########################
#Run the ANOVA and save it as an object
aov.shannon.donor.type = aov(shannon_entropy ~ donor.type, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.shannon.donor.type)

TukeyHSD(aov.shannon.donor.type)

#Plot
boxplot(shannon_entropy ~ donor.type.ord, data=meta, ylab="Shannon Entropy")

shannon_boxplot <- ggplot(meta, aes(donor.type.ord, shannon_entropy, fill = donor.type.ord)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y="Shannon Entropy", x = "", fill = "Donor Type") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/shannon_boxplot.png", shannon_boxplot, height = 3, width = 3)

shannon_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(donor.type.ord) %>%   # the grouping variable
  summarise(mean_shannon = mean(shannon_entropy),  # calculates the mean of each group
            sd_shannon = sd(shannon_entropy), # calculates the standard deviation of each group
            n_shannon = n(),  # calculates the sample size per group
            se_shannon = sd(shannon_entropy)/sqrt(n())) # calculates the standard error of each group

shannon_se <- ggplot(shannon_summary, aes(donor.type.ord, mean_shannon, fill = donor.type.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Entropy  ± s.e.", x = "") +
  scale_fill_manual(values = c("Piglet" = "gray", "Sow" = "lightblue"))
ggsave("output/shannon_se.png", shannon_se, height = 2.5, width = 3) 






