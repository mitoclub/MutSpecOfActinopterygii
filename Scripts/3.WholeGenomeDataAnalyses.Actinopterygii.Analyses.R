###################################
rm(list=ls(all=TRUE))
library(ggpubr)
library(caper)
library(geiger)
library(dplyr)


### reading whole genomes database
AllSynNuc = read.table("../data/1raw/AllGenesCodonUsageNoOverlap.txt", header = TRUE, sep = '\t')
AllSynNuc = AllSynNuc[AllSynNuc$Gene != 'ND6',]

####### obtaining neutral nucleotide fractions in whole genomes
AllSynNuc = aggregate(list(AllSynNuc$NeutralA,AllSynNuc$NeutralT,AllSynNuc$NeutralG,AllSynNuc$NeutralC), by = list(AllSynNuc$Species), FUN = sum)
names(AllSynNuc) = c('Species','NeutralA','NeutralT','NeutralG','NeutralC')
AllSynNuc$FrA = AllSynNuc$NeutralA / (AllSynNuc$NeutralA + AllSynNuc$NeutralT + AllSynNuc$NeutralG + AllSynNuc$NeutralC)
AllSynNuc$FrT = AllSynNuc$NeutralT / (AllSynNuc$NeutralA + AllSynNuc$NeutralT + AllSynNuc$NeutralG + AllSynNuc$NeutralC) 
AllSynNuc$FrG = AllSynNuc$NeutralG / (AllSynNuc$NeutralA + AllSynNuc$NeutralT + AllSynNuc$NeutralG + AllSynNuc$NeutralC) 
AllSynNuc$FrC = AllSynNuc$NeutralC / (AllSynNuc$NeutralA + AllSynNuc$NeutralT + AllSynNuc$NeutralG + AllSynNuc$NeutralC) 

### merge whole genomes with temperature
tm = read.table('../Data/1raw/FishBaseTemperature.txt', header=T)
tm = tm %>% group_by(Species) %>% summarise(temperature = median(Temperature))
matur = read.table('../Data/1raw/FishBaseMaturity_Tm.txt', header=T)
matur = matur %>% group_by(Species) %>% summarise(matur_tm = median(Tm))
AllSynNuc = merge(AllSynNuc, tm, by='Species', all.x=T)
AllSynNuc = merge(AllSynNuc, matur, by='Species', all.x=T)
nrow(AllSynNuc[!is.na(AllSynNuc$temperature),])

##############Rank corr 
#Suppl. mat. 2.b
cor.test(log2(AllSynNuc$temperature),AllSynNuc$FrA, method = "spearman")
cor.test(log2(AllSynNuc$temperature),AllSynNuc$FrT, method = "spearman")
cor.test(log2(AllSynNuc$temperature),AllSynNuc$FrG, method = "spearman")
cor.test(log2(AllSynNuc$temperature),AllSynNuc$FrC, method = "spearman")

AllSynNuc$TG = AllSynNuc$FrT+AllSynNuc$FrG
AllSynNuc$AC = AllSynNuc$FrA+AllSynNuc$FrC
AllSynNuc$TG_ACSkew = (AllSynNuc$TG-AllSynNuc$AC)/(AllSynNuc$TG+AllSynNuc$AC); summary(AllSynNuc$TG_ACSkew)
AllSynNuc$AC_TGSkew = -(AllSynNuc$TG-AllSynNuc$AC)/(AllSynNuc$TG+AllSynNuc$AC); summary(AllSynNuc$AC_TGSkew) #OUR

######Figores 
pdf("../Figures/WholeGenomeData.Actinopterygii.pdf")
AllSynNuc = AllSynNuc[!is.na(AllSynNuc$temperature),]; N = as.character(paste("N", nrow(AllSynNuc), sep="")) #delete NA
ggscatter(AllSynNuc, x = "temperature", y = "FrT",
          color = "#e61a0b", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Median annual water temperature, ?C", ylab="Whole genome neutral fraction of A")+stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))

ggscatter(AllSynNuc, x = "temperature", y = "FrG",
          color = "#009414", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Median annual water temperature, ?C", ylab="Whole genome neutral fraction of C")+stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))

ggscatter(AllSynNuc, x = "temperature", y = "FrC",
          color = "#5c5c5c", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Median annual water temperature, ?C", ylab="Whole genome neutral fraction of G")+stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))

ggscatter(AllSynNuc, x = "temperature", y = "FrA",
          color = "#0918e6", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Median annual water temperature, ?C", ylab="Whole genome neutral fraction of T")+stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))

dev.off()

pdf("../Figures/WholeGenomeData.Actinopterygii.TG-SACskew.pdf", width = 7, height = 8.5)
ggscatter(AllSynNuc, x = "temperature", y = "AC_TGSkew",
          color = "#c99bc9", # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          xlab="Mean annual water temperature, ?C", ylab="STG-SAC skew", ylim=c(0.1, 0.65), xlim=c(0,30))+stat_cor(method = "spearman", label.x = 2.5, label.y = 0.64, aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
dev.off()

### LM ANALYSES:
summary(lm(FrT ~ scale(Temperature)+scale(Tm), data = SynNuc))
summary(lm(FrT ~ log2(Temperature + 2)*log2(Tm), data = SynNuc))  # keep it for presentation!!!
summary(lm(FrT ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(FrG ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # strong
summary(lm(FrA ~ log2(Temperature + 2)+log2(Tm), data = SynNuc)) # strong

#Suppl. Mat. 2?
summary(lm(AC_TGSkew ~ log2(Temperature + 2)+log2(Tm), data = SynNuc))
summary(lm(AC_TGSkew ~ scale(Temperature + 2)+scale(Tm), data = SynNuc))# ###PICS
summary(lm(AC_TGSkew ~ Temperature, data = SynNuc))
getwd()


###################################################### phylogenetic inertia analysis
#tree = read.tree('../Data/1raw/mtalign.aln.treefile.rooted')
tree = read.tree('../Data/1raw/Fishes_Polymorphic_sp.nwk')

row.names(AllSynNuc) = AllSynNuc$Species
tree_pruned = treedata(tree, AllSynNuc, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_pruned, AllSynNuc, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)
data$AC_TGSkew = as.numeric(as.character(data$AC_TGSkew))
data$temperature = as.numeric(as.character(data$temperature))
data$matur_tm = as.numeric(as.character(data$matur_tm))
data_comp <- comparative.data(tree_pruned, data, Species, vcv=TRUE)

summary(pgls(AC_TGSkew ~ log2(temperature+2), data_comp, lambda="ML"))
summary(pgls(log2(AC_TGSkew) ~ log2(temperature+2), data_comp, lambda="ML"))
summary(pgls(log2(temperature+2) ~ AC_TGSkew, data_comp, lambda="ML"))

summary(pgls(AC_TGSkew ~ scale(temperature), data_comp, lambda="ML"))
summary(pgls(AC_TGSkew ~ log2(temperature + 2) + log2(matur_tm), data_comp, lambda="ML"))
summary(pgls(AC_TGSkew ~ temperature, data_comp, lambda="ML"))
sunnary(pgls(AC_TGSkew ~ temperature + matur_tm, data_comp, lambda="ML"))

####FIgures
pdf("../Figures/WholeGenomeData.Actinopterygii.TG-SACskewWithTempAndTm.pdf", width = 7, height = 8.5)
AllSynNuc = AllSynNuc[!is.na(AllSynNuc$matur_tm),]; N = as.character(paste("N", nrow(AllSynNuc), sep=""))
medianTm = median(AllSynNuc[!is.na(AllSynNuc$matur_tm),]$matur_tm)
AllSynNuc$Longevity = "Na"
AllSynNuc[AllSynNuc$matur_tm < medianTm,]$Longevity = "ShortLived"
AllSynNuc[AllSynNuc$matur_tm >= medianTm,]$Longevity = "LongLived"

plot(AllSynNuc[AllSynNuc$Longevity == "ShortLived",]$temperature, AllSynNuc[AllSynNuc$Longevity == "ShortLived",]$AC_TGSkew, col="#4da36c", xlab="Mean annual water temperature, C?", ylab="STG-SAC skew", ylim=c(0.1, 0.65), xlim=c(0,30), pch = 16)
abline((0.331911-0.049196), 0.006172, col="#4da36c", lwd = 2)
par(new=TRUE)
plot(AllSynNuc[AllSynNuc$Longevity == "LongLived",]$temperature, AllSynNuc[AllSynNuc$Longevity == "LongLived",]$AC_TGSkew, col="#42cbf5", xlab="", ylab="", ylim=c(0.1, 0.65), xlim=c(0,30), pch = 16)
abline(0.331911, 0.006172, col="#42cbf5", lwd = 2)
legend("bottomright", legend=c( "Short time of maturation","Long time of maturation"), col=c("#4da36c","#42cbf5"), pch = 16)
dev.off()
