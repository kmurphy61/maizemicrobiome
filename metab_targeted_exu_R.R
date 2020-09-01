#targeted metabolomics microbiome study 11/5/2018

#load necessary packages
library(ggplot2)
library(lsmeans)

#set the file path (where is your text file?)
getwd()
setwd("/Users/katiemurphy/Desktop/KM for R")

#import text file
target <- read.table(file = "peakarea_exu_v2.txt", header = TRUE, fill = TRUE)

head(target)
summary(target)

target$plant <- as.factor(target$plant)

#Epoxydolabradiene
full_model = lm(Epoxydolabradiene ~ genotype, target)
summary(full_model)
anova(full_model)

vscontrol = lsmeans(full_model, trt.vs.ctrl ~ genotype)$contrasts
summary(vscontrol,infer = T, side = '>')

#THD
full_model = lm(Trihydroxydolabrene ~ genotype, target)
summary(full_model)
anova(full_model)

vscontrol = lsmeans(full_model, trt.vs.ctrl ~ genotype)$contrasts
summary(vscontrol,infer = T, side = '>')

#Benzoxazinoids
full_model = lm(Benzoxazinoids ~ genotype, target)
summary(full_model)
anova(full_model)

vscontrol = lsmeans(full_model, trt.vs.ctrl ~ genotype)$contrasts
summary(vscontrol,infer = T)
