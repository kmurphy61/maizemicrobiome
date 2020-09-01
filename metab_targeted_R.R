#targeted metabolomics microbiome study 11/5/2018

#load necessary packages
library(ggplot2)
library(lsmeans)

#set the file path (where is your text file?)
getwd()
setwd("/Users/zerbe/Desktop/untargeted")

#import text file
target <- read.table(file = "peakarea_root3.txt", header = TRUE, fill = TRUE)

head(target)
summary(target)

target$pot <- as.factor(target$pot)

#Epoxydolabradiene
full_model = lm(Epoxydolabradiene ~ genotype, target)
summary(full_model)
anova(full_model)

paired2 = lsmeans(full_model,pairwise ~ genotype)$contrasts
paired2

vscontrol = lsmeans(full_model, trt.vs.ctrl ~ genotype)$contrasts
summary(vscontrol,infer = T, side = '>')

#THD
full_model = lm(Trihydroxydolabrene ~ genotype, target)
summary(full_model)
anova(full_model)

paired2 = lsmeans(full_model,pairwise ~ genotype)$contrasts
paired2

vscontrol = lsmeans(full_model, trt.vs.ctrl ~ genotype)$contrasts
summary(vscontrol,infer = T, side = '>')

#Benzoxazinoids
full_model = lm(Benzoxazinoids ~ genotype, target)
summary(full_model)
anova(full_model)

paired1 = lsmeans(full_model,pairwise ~ genotype)$contrasts
paired1

paired2 = lsmeans(full_model,pairwise ~ genotype)$contrasts
paired2
