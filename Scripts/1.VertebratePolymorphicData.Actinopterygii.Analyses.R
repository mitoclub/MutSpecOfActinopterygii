rm(list=ls(all=TRUE))

if (!require(caper)) install.packages("caper")
if (!require(geiger)) install.packages("geiger")
if (!require(ggpubr)) install.packages("ggpubr")
if (!require(dplyr)) install.packages("dplyr")
library(caper)
library(geiger)
library(ggpubr)
library(dplyr)

##########Reading MutSpec DataBase with temp and mt of Actinopterigii 
mutSpec = read.table('../Data/3results/VertebratePolymorphisms.MutSpecData.txt', header = TRUE)
table(mutSpec$Class)

mutSpecCytB = mutSpec[mutSpec$Gene == "CytB",]
mutSpecAllMean = mutSpec %>% group_by(Species, Class, temperature, matur_tm); mutSpecAllMean = mutSpecAllMean %>% summarise(A_C.N=mean(A_C.N), A_G.N=mean(A_G.N), A_T.N=mean(A_T.N), C_A.N=mean(C_A.N), C_G.N=mean(C_G.N), C_T.N=mean(C_T.N), G_A.N=mean(G_A.N), G_C.N=mean(G_C.N), G_T.N=mean(G_T.N), T_A.N=mean(T_A.N), T_C.N=mean(T_C.N), T_G.N=mean(T_G.N))
table(mutSpecAllMean[!is.na(mutSpecAllMean$temperature),]$Class)

# Average MutSpec for Actinopteri
pdf('../Figures/PolymorphicData.Actinopterygii.AverageMutSpec.pdf')
averageMutSpec = mutSpecAllMean[!is.na(mutSpecAllMean$temperature),][,5:16]; summary(averageMutSpec)
averageMutSpec = averageMutSpec[!is.na(averageMutSpec$G_A.N) & !is.na(averageMutSpec$G_C.N) & !is.na(averageMutSpec$G_T.N),]; summary(averageMutSpec)
averageMutSpec = as.data.frame(apply(averageMutSpec, 2, mean)); averageMutSpec$Subs = row.names(averageMutSpec); names(averageMutSpec) = c("Freq", "Subs")
f1 = ggbarplot(averageMutSpec, x = "Subs", y = "Freq", fill = "Subs", color = "Subs",
               palette = c("#bdbdbd", "#036a5b", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#9c3d37", "#055088", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#73514f", "#bdbdbd"), 
               xlab="Substitution types", ylab="Normalised frequencies", legend = "none")  
f1
dev.off()


#only CytB
cor.test(mutSpecCytB$A_T.N,mutSpecCytB$temperature, method = 'spearman')   
cor.test(mutSpecCytB$A_G.N,mutSpecCytB$temperature, method = 'spearman')      
cor.test(mutSpecCytB$A_C.N,mutSpecCytB$temperature, method = 'spearman')      
cor.test(mutSpecCytB$T_A.N,mutSpecCytB$temperature, method = 'spearman')   
cor.test(mutSpecCytB$T_G.N,mutSpecCytB$temperature, method = 'spearman')   
cor.test(mutSpecCytB$T_C.N,mutSpecCytB$temperature, method = 'spearman')       
cor.test(mutSpecCytB$G_A.N,mutSpecCytB$temperature, method = 'spearman')    
cor.test(mutSpecCytB$G_T.N,mutSpecCytB$temperature, method = 'spearman')    
cor.test(mutSpecCytB$G_C.N,mutSpecCytB$temperature, method = 'spearman')    
cor.test(mutSpecCytB$C_A.N,mutSpecCytB$temperature, method = 'spearman')     
cor.test(mutSpecCytB$C_T.N,mutSpecCytB$temperature, method = 'spearman')     
cor.test(mutSpecCytB$C_G.N,mutSpecCytB$temperature, method = 'spearman')    
#all by mean
cor.test(mutSpecAllMean$A_T.N,mutSpecAllMean$temperature, method = 'spearman')   
cor.test(mutSpecAllMean$A_G.N,mutSpecAllMean$temperature, method = 'spearman')      
cor.test(mutSpecAllMean$A_C.N,mutSpecAllMean$temperature, method = 'spearman')      
cor.test(mutSpecAllMean$T_A.N,mutSpecAllMean$temperature, method = 'spearman')   
cor.test(mutSpecAllMean$T_G.N,mutSpecAllMean$temperature, method = 'spearman')   
cor.test(mutSpecAllMean$T_C.N,mutSpecAllMean$temperature, method = 'spearman')       
cor.test(mutSpecAllMean$G_A.N,mutSpecAllMean$temperature, method = 'spearman')    
cor.test(mutSpecAllMean$G_T.N,mutSpecAllMean$temperature, method = 'spearman')    
cor.test(mutSpecAllMean$G_C.N,mutSpecAllMean$temperature, method = 'spearman')    
cor.test(mutSpecAllMean$C_A.N,mutSpecAllMean$temperature, method = 'spearman')     
cor.test(mutSpecAllMean$C_T.N,mutSpecAllMean$temperature, method = 'spearman')     
cor.test(mutSpecAllMean$C_G.N,mutSpecAllMean$temperature, method = 'spearman') 

##Figures

pdf('../Figures/PolymorphicData.Actinopterygii.pdf')
mutSpecAllMean = mutSpecAllMean[!is.na(mutSpecAllMean$temperature),]; N = as.character(paste("N", nrow(mutSpecAllMean), sep="")) #delete NA
f1a = ggscatter(mutSpecAllMean, x = "temperature", y = "T_C.N",
                color = "#73514f", # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                yscale = "log2", xlab="Median annual water temperature, ?C", ylab="log2 AH>GH") + stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1a

f1b1 = ggscatter(mutSpecAllMean, x = "temperature", y = "G_A.N",
                 color = "#055088", # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 yscale = "log2", xlab="Median annual water temperature, ?C", ylab="log2 CH>TH")+ stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1b1

f1b2 = ggscatter(mutSpecAllMean, x = "temperature", y = "C_T.N",
                 color = "#9c3d37", # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 yscale = "log2", xlab="Median annual water temperature, ?C", ylab="log2 GH>AH")+ stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1b2
dev.off()


######### TC / AG rank corr and log2 models
mutSpecAllMean$TCdivAG=mutSpecAllMean$T_C.N/mutSpecAllMean$A_G.N
mutSpecAllMean[mutSpecAllMean$TCdivAG == "Inf" | mutSpecAllMean$TCdivAG == "NaN",]$TCdivAG = NA
cor.test(mutSpecAllMean$TCdivAG,mutSpecAllMean$temperature, method = 'spearman')
mutSpecAllMean$GAdivCT=mutSpecAllMean$G_A.N/mutSpecAllMean$C_T.N
mutSpecAllMean[mutSpecAllMean$GAdivCT == "Inf" | mutSpecAllMean$GAdivCT == "NaN" | is.na(mutSpecAllMean$GAdivCT),]$GAdivCT = NA
cor.test(mutSpecAllMean$GAdivCT,mutSpecAllMean$temperature, method = 'spearman')



pdf('../Figures/PolymorphicData.Actinopterygii.A_GdivT_C.pdf')
f1c = ggscatter(mutSpecAllMean, x = "temperature", y = "TCdivAG",
                color = "#814194", # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                yscale = "log2", xlab="Median annual water temperature, ?C", ylab="log2 A_GdivT_C")+ stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1c 

f1d = ggscatter(mutSpecAllMean, x = "temperature", y = "GAdivCT",
                color = "#8C99A6", # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                yscale = "log2", xlab="Median annual water temperature, ?C", ylab="log2 C_TdivG_A")+ stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1d
dev.off()


#####correlation of mutspec with time of maturation in fishes
#all by mean
cor.test(mutSpecAllMean$A_T.N,mutSpecAllMean$matur_tm, method = 'spearman')   
cor.test(mutSpecAllMean$A_G.N,mutSpecAllMean$matur_tm, method = 'spearman')      
cor.test(mutSpecAllMean$A_C.N,mutSpecAllMean$matur_tm, method = 'spearman')      
cor.test(mutSpecAllMean$T_A.N,mutSpecAllMean$matur_tm, method = 'spearman')   
cor.test(mutSpecAllMean$T_G.N,mutSpecAllMean$matur_tm, method = 'spearman')   
cor.test(mutSpecAllMean$T_C.N,mutSpecAllMean$matur_tm, method = 'spearman')       
cor.test(mutSpecAllMean$G_A.N,mutSpecAllMean$matur_tm, method = 'spearman')    
cor.test(mutSpecAllMean$G_T.N,mutSpecAllMean$matur_tm, method = 'spearman')    
cor.test(mutSpecAllMean$G_C.N,mutSpecAllMean$matur_tm, method = 'spearman')    
cor.test(mutSpecAllMean$C_A.N,mutSpecAllMean$matur_tm, method = 'spearman')     
cor.test(mutSpecAllMean$C_T.N,mutSpecAllMean$matur_tm, method = 'spearman')     
cor.test(mutSpecAllMean$C_G.N,mutSpecAllMean$matur_tm, method = 'spearman')


#####################################################
############Multiple models##########################

summary(lm(formula = temperature ~ scale(T_C.N) + scale(A_G.N), data = mutSpecAllMean))
summary(lm(formula = T_C.N ~ scale(temperature) * scale(matur_tm), data = mutSpecAllMean))
summary(lm(formula = T_C.N ~ scale(temperature) + scale(matur_tm), data = mutSpecAllMean))
summary(lm(formula = T_C.N ~ scale(temperature), data = mutSpecAllMean))
summary(lm(formula = A_G.N ~ scale(temperature) + scale(matur_tm), data = mutSpecAllMean))
summary(lm(formula = A_G.N ~ scale(temperature), data = mutSpecAllMean))

mutSpecAllMeanNAzeroomit = mutSpecAllMean[!is.na(mutSpecAllMean$TCdivAG),]
mutSpecAllMeanNAzeroomit = mutSpecAllMeanNAzeroomit[mutSpecAllMeanNAzeroomit$TCdivAG != 0,]
mutSpecAllMeanNAzeroomit = mutSpecAllMeanNAzeroomit[!is.na(mutSpecAllMeanNAzeroomit$GAdivCT),]
mutSpecAllMeanNAzeroomit = mutSpecAllMeanNAzeroomit[mutSpecAllMeanNAzeroomit$GAdivCT != 0,]

summary(lm(formula = log2(TCdivAG) ~ scale(temperature), data = mutSpecAllMeanNAzeroomit))
summary(lm(formula = temperature ~ scale(TCdivAG), data = mutSpecAllMeanNAzeroomit))

########################################################################################################
### PICs ###############################################################################################
########################################################################################################
tree = read.tree('../Data/1raw/mtalign.aln.treefile.rooted')

row.names(mutSpecAllMeanNAzeroomit) = mutSpecAllMeanNAzeroomit$Species
tree_pruned = treedata(tree, mutSpecAllMeanNAzeroomit, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_pruned, mutSpecAllMeanNAzeroomit, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)
data$TCdivAG = as.numeric(as.character(data$TCdivAG))
data$temperature = as.numeric(as.character(data$temperature))
data$matur_tm = as.numeric(as.character(data$matur_tm))
data_comp <- comparative.data(tree_pruned, data, Species, vcv=TRUE)

summary(pgls(TCdivAG ~ scale(temperature), data_comp, lambda="ML"))
summary(pgls(log2(TCdivAG) ~ log2(temperature), data_comp, lambda="ML"))
