library(ggplot2)
library(reshape2)
library(scales)
library(vegan)
library(lubridate)
library(tidyverse)

library(MASS)
library(stringr)
library(forcats)
library(broom)

getwd()

setwd("/Users/zerbe/Desktop/untargeted/")

## Load the intensity file
intens3 <- read.table(file = "peak_height_pos_2020.txt" , header = TRUE)
intens3$ID <- as.character(intens3$ID)
intens3$mz <- as.character(intens3$mz)
intens3$rt <- as.character(intens3$rt)

var_match <- read.table(file = "variable with name pos 2020.txt", header = T)
var_match$variable <- as.factor(var_match$variable)

## Load the mapping file
lc_map <- read.table(file = "lcms_map_pos_2020.txt", header = TRUE, row.names = 1)
lc_map$SampleNumber <- as.factor(lc_map$SampleNumber)
lc_map$Genotype <- factor(lc_map$Genotype, levels = c("NA", "WT", "AN2"))
lc_map2 <- lc_map[match(colnames(intens3), row.names(lc_map)),] 

summary(lc_map2)
head(lc_map2)
summary(intens3)

intens3$max <- apply(intens3[4:20],1,max)

intens4 <- intens3 %>%
  mutate (blank_calc = ifelse(max >= 3*Control,"good","blank"))

intens4 <-  intens4 %>%
  filter(blank_calc == "good")

blank_features <- intens4 %>%
  filter(blank_calc == "blank")
#katie export here

intens4 <- subset( intens4, select = -blank_calc )

#make map of just samples, no standards
map_sub <- subset(lc_map2, Genotype != "NA")
map_sub <- subset(lc_map2, Water != "D")
intens4_sub <- intens4[,match(row.names(map_sub), colnames(intens4))]

summary(intens4_sub)
head(intens4_sub)

#PCoA analysis and graphing
pc <- capscale(t(log2(intens4_sub +1 )) ~ 1, distance = "bray")
pc.axes <- cbind(map_sub, scores(pc, choices = c(1:3))$sites)
head(pc.axes)

evals  <- eigenvals(pc)
evals_percent <- (evals / sum(evals)) *100
evals_percent <- round(evals_percent,digits=2)

head(pc.axes)
summary(evals_percent)

metab_plot <- ggplot(pc.axes, aes(MDS1, MDS2, color = Genotype)) + 
  theme_bw() +
  geom_point(size = 2) +
  xlab(paste("PCo1 (", evals_percent[1], " %)", sep = "")) +
  ylab(paste("PCo2 (", evals_percent[2], " %)", sep = "")) +
  scale_color_manual(values=c("#619CFF","#F8766D")) + 
  theme(text = element_text(size = 11), legend.position = "bottom", legend.box = "horizontal") +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5), shape = guide_legend(title.position = "top", title.hjust = 0.5)) 
metab_plot
  ggsave(filename = "pcoa_metab.jpg", path = NULL,
         scale = 1, width = 3.35, height = 4, units = c("in"))

#permanova analysis
adonis(t(log2(intens4_sub+1))~ Genotype,data = map_sub)
#p value less than 0.05 is significant difference at 95% confidence

#-------------- WHICH PEAKS ARE SIG DIFF, NO NORMALIZATION

rownames(intens4) <- intens4[,"ID"]
intens5 <- intens4[5:20]
intens5 <-  as.data.frame(t(intens5))

intens5 <- tibble::rownames_to_column(intens5, "SampleID")
#intens5 <- intens5[match(row.names(intens5), row.names(map_sub)),] 

readymelt <- inner_join(map_sub, intens5, by = "SampleID")

intens3_tidy1 <- melt(readymelt, id.vars = names(map_sub))
str(intens3_tidy1)

#intens3_sub_ND <- intens4_sub[,match(row.names(map_sub), colnames(intens4_sub))]
#ntens3_tidy_ND <- melt(cbind(map_sub_ND, t(intens3_sub_ND)), id.vars = names(map_sub_ND))
#str(intens3_tidy_ND)
#nrow(intens3_tidy_ND)

comp.lm_ND <- intens3_tidy1 %>% 
  group_by(variable) %>% 
  summarise(prev = sum(value > 0) / n()) %>% 
  filter(prev >= 0) %>% 
  inner_join(intens3_tidy1, by = c("variable")) %>% 
  group_by(variable) %>% 
  na.omit() %>% 
  nest() %>% 
  mutate(nb_model_ND = map(data, ~lm(value ~ Genotype, .))) %>% 
  mutate(type = "No offset") 

comp.lm_results_ND <- comp.lm_ND %>%
  unnest_legacy(map(nb_model_ND, tidy)) %>% 
  inner_join(var_match, by = "variable")
 # mutate(p.adj = p.adjust(p.value, "BH"))
head(comp.lm_results_ND)

comp.lm_results_ND_export <- comp.lm_results_ND %>% 
  filter(p.value <= 0.05) %>%
  arrange(estimate) %>% 
  filter(term == "GenotypeAN2")

comp.lm_results_ND_notsig <- comp.lm_results_ND %>% 
  filter(p.value > 0.05) %>%
  arrange(estimate) %>% 
  filter(term == "GenotypeAN2")
nrow(comp.lm_results_ND_notsig)

comp.lm_results_ND_export <- comp.lm_results_ND_export %>%
  mutate(mode = "positive")%>%
  mutate(key = paste(mode, variable, sep = "-"))

write.table(comp.lm_results_ND_export, "/Users/zerbe/Desktop/untargeted/sig diff lists/POS/comp.lm_results_ND_export.txt", sep="\t")

head(comp.lm_results_ND_export)
head(var_match)

comp.lm_results_combo <- rbind(comp.lm_results_ND_export,comp.lm_results_ND_export_neg)
write.table(comp.lm_results_combo, "/Users/zerbe/Desktop/untargeted/sig diff lists/POS/comp.lm_results_ND_combo.txt", sep="\t")

#enriched in AN2
AN2 <-  comp.lm_results_ND_export %>% 
  filter(p.value <= 0.05) %>% 
  filter(estimate > 0)

#encriched in WT
WT <-  comp.lm_results_ND_export %>% 
  filter(p.value <= 0.05) %>% 
  filter(estimate < 0)

#boxplots of significant metabolites

new_intens <- intens4 %>%
  gather(SampleID, value, -ID) %>%
  mutate(ID = as.character(ID))%>%
  mutate(value = as.numeric(value)) %>%
  inner_join(map_sub) %>%
  group_by(SampleID)

data_for_heatC <- AN2 %>% 
  inner_join(new_intens, by = c("variable" = "ID")) %>% 
  group_by(variable) #%>%

data_for_heatD <- WT %>% 
  inner_join(new_intens, by = c("variable" = "ID")) %>% 
  group_by(variable) #%>% 

  #mutate(scaled_value = (value-min(value))/(max(value)-min(value)))

summary(data_for_heatC)

plot_C <- ggplot(subset(data_for_heatC), aes(x = Genotype, y = value, fill = Genotype)) +
  geom_boxplot() +
  facet_wrap( ~ key, scales = "free", ncol = 6) +
  theme_bw() +
  labs(y = "value") +
  theme(text = element_text(size = 10), axis.ticks.x = element_blank(), legend.key = element_blank()) + 
  scale_fill_manual(values=c( "#619CFF", "#F8766D")) + 
  theme(legend.title=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank()) +
  theme(legend.justification=c(1,0), legend.position = "bottom") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
plot_C

plot_D <- ggplot(subset(data_for_heatD), aes(x = Genotype, y = value, fill = Genotype)) +
  geom_boxplot() +
  facet_wrap( ~ key, scales = "free", ncol = 6) +
  theme_bw() +
  labs(y = "value") +
  theme(text = element_text(size = 10), axis.ticks.x = element_blank(), legend.key = element_blank()) + 
  scale_fill_manual(values=c( "#619CFF", "#F8766D")) + 
  theme(legend.title=element_blank(), axis.text.x=element_blank(), axis.title.x=element_blank()) +
  theme(legend.justification=c(1,0), legend.position = "bottom") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
plot_D

