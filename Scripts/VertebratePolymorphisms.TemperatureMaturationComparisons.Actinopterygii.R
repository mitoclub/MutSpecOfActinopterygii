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

pdf('../Figures/VertebratePolymorphisms.MutSpecData.Actinopterygii.pdf')
mutSpecAllMean = mutSpecAllMean[!is.na(mutSpecAllMean$temperature),]; N = as.character(paste("N", nrow(mutSpecAllMean), sep="")) #delete NA
f1a = ggscatter(mutSpecAllMean, x = "temperature", y = "T_C.N",
                color = "#73514f", # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                yscale = "log2", xlab="Median annual water temperature, ?C", ylab="AH>GH") + stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1a

f1b1 = ggscatter(mutSpecAllMean, x = "temperature", y = "G_A.N",
                 color = "#055088", # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 yscale = "log2", xlab="Median annual water temperature, ?C", ylab="CH>TH")+ stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1b1

f1b2 = ggscatter(mutSpecAllMean, x = "temperature", y = "C_T.N",
                 color = "#9c3d37", # Points color, shape and size
                 add = "reg.line",  # Add regressin line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE, # Add confidence interval
                 yscale = "log2", xlab="Median annual water temperature, ?C", ylab="GH>AH")+ stat_cor(method = "spearman", aes(label = paste(..r.label.., ..p.label.., ..N.., sep = "~`,`~")))
f1b2


dev.off()


######### TC / AG rank corr and log2 models
mutSpecAllMean$TCdivAG=mutSpecAllMean$T_C.N/mutSpecAllMean$A_G.N
mutSpecAllMean[mutSpecAllMean$TCdivAG == "Inf" | mutSpecAllMean$TCdivAG == "NaN",]$TCdivAG = NA
cor.test(mutSpecAllMean$TCdivAG,mutSpecAllMean$temperature, method = 'spearman')
mutSpecAllMean$GAdivCT=mutSpecAllMean$G_A.N/mutSpecAllMean$C_T.N
mutSpecAllMean[mutSpecAllMean$GAdivCT == "Inf" | mutSpecAllMean$GAdivCT == "NaN" | is.na(mutSpecAllMean$GAdivCT),]$GAdivCT = NA
cor.test(mutSpecAllMean$GAdivCT,mutSpecAllMean$temperature, method = 'spearman')



pdf('../Figures/VertebratePolymorphisms.MutSpecData.Actinopterygii.AGTC.pdf')
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

summary(lm(formula = log2(TCdivAG) ~ scale(temperature), data = mutSpecAllMeanNAzeroomit))
summary(lm(formula = temperature ~ scale(TCdivAG), data = mutSpecAllMeanNAzeroomit))

########################################################################################################
### PICs ###############################################################################################
########################################################################################################
tree = read.tree('../Data/1raw/mtalign.aln.treefile.rooted')

row.names(mutSpecAllMeanNAzeroomit) = mutSpecAllMeanNAzeroomit$Species

tree_pruned = treedata(tree, mutSpecAllMeanNAzeroomit, sort=T, warnings=T)$phy

#   Not found in the tree and were dropped from the dataframe:
# Boops_boops
# Helicolenus_dactylopterus
# Lepidorhombus_whiffiagonis
# Lethrinus_olivaceus
# Lithognathus_mormyrus
# Lutjanus_vitta
# Pagellus_acarne
# Plagioscion_squamosissimus
# Sarda_sarda
# Scomberomorus_commerson
# Solea_solea
# Squalius_cephalus
# Squalius_pyrenaicus
# Trisopterus_minutus

data<-as.data.frame(treedata(tree_pruned, allparameters, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$TCdivAG = as.numeric(as.character(data$TCdivAG))
data$Temperature = as.numeric(as.character(data$Temperature))

data_comp <- comparative.data(tree_pruned, data, Species, vcv=TRUE)
 ##Supplementary 1d
model = pgls(TCdivAG ~ scale(Temperature), data_comp, lambda="ML")
summary(model)

# lambda [ ML]  : 0.731
# lower bound : 0.000, p = 0.03528
# upper bound : 1.000, p = 1.2987e-08
# 95.0% CI   : (0.217, 0.900)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)         5.76508    1.93865  2.9738 0.004715 **
#   scale(Temperature)  0.55015    0.53135  1.0354 0.306023   
# scale(Tm)          -0.14488    0.57032 -0.2540 0.800630   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.819 on 45 degrees of freedom
# Multiple R-squared: 0.03013,	Adjusted R-squared: -0.01297 
# F-statistic: 0.6991 on 2 and 45 DF,  p-value: 0.5024 

model = pgls(log2(TCdivAG) ~ scale(Temperature) + scale(Tm), data_comp, lambda="ML")
summary(model)

# lambda [ ML]  : 0.535
# lower bound : 0.000, p = 0.14424
# upper bound : 1.000, p = 3.3975e-12
# 95.0% CI   : (NA, 0.830)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)         1.69761    0.63241  2.6843  0.01014 *
#   scale(Temperature)  0.47952    0.21438  2.2368  0.03030 *
#   scale(Tm)          -0.08112    0.22464 -0.3611  0.71971  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.436 on 45 degrees of freedom
# Multiple R-squared: 0.1189,	Adjusted R-squared: 0.07972 
# F-statistic: 3.036 on 2 and 45 DF,  p-value: 0.05798 

summary(lm(pic(data$TCdivAG, tree_pruned) ~ pic(data$Temperature, tree_pruned) +
             pic(data$Tm, tree_pruned)))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)                        -0.961202   2.371449  -0.405    0.687
# pic(data$Temperature, tree_pruned) -0.005479   0.091820  -0.060    0.953
# pic(data$Tm, tree_pruned)           0.110200   0.127524   0.864    0.392





####### MAMMALS
GT = read.table('../../Body/1Raw/GenerationLenghtforMammals.xlsx.txt', sep = '\t', header = TRUE)
GT$Species = gsub(' ','_',GT$Scientific_name)
ALL = merge(MUT, GT)

ALL$TCAG = ALL$T_C / ALL$A_G
ALL = ALL[ALL$TCAG != Inf,]
ALL = ALL[!is.na(ALL$TCAG),]

summary(ALL$TCAG)

cor.test(ALL$TCAG, ALL$GenerationLength_d, method = "spearman")
cor.test(ALL$T_C, ALL$GenerationLength_d, method = "spearman")
cor.test(ALL$A_G, ALL$GenerationLength_d, method = "spearman")












#################################metabolic rate approximation
allparameters = merge (TemperMut, MATUTM)
allparameters$MR=(allparameters$Tm+1)^0.75
allparameters$TemperatureK = 273.15 + allparameters$Temperature
allparameters$B=allparameters$MR * exp(-1.2/((8.617*10^-5)*allparameters$TemperatureK))
cor.test(allparameters$B, allparameters$T_C, method="spearman") #rho 0.4075402 
cor.test(allparameters$Temperature, allparameters$T_C, method = "spearman")



####################################
####################################
###Full ecology for Kuptsov
####################################
###########Taxonomy###################################################################
Taxa = read.table("../../Body/1Raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); 
for (i in (1:nrow(Taxa)))  {Taxa$Species[i] = paste(unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[1],unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[2], sep = '_')}
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)
Taxa = Taxa[,-1]
TaxaMore = read.table("../../Body/1Raw/TaxaFromKostya.2NeedTaxa.tax.txt", sep = '\t',header = FALSE) 
TaxaMore$Species = ''
for (i in (1:nrow(TaxaMore)))  
{TaxaMore$Species[i] = paste(unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[1],unlist(strsplit(as.character(TaxaMore$V1[i]),split = ' '))[2], sep = '_')}
TaxaMore$Class = gsub("; Chordata;.*",'',TaxaMore$V2); 
TaxaMore$Class = gsub(".*; ",'',TaxaMore$Class); 
TaxaMore$Class = gsub('Actinopteri','Actinopterygii',TaxaMore$Class)
TaxaMore$Class = gsub("Testudines|Squamata|Crocodylia",'Reptilia',TaxaMore$Class)
table(TaxaMore$Class)
TaxaMore = TaxaMore[,-c(1,2)]
Taxa = rbind(Taxa,TaxaMore); Taxa = unique(Taxa)

FISHsystematix= merge(TEMPE, MUT, all = T)
FISHsystematix= merge(FISHsystematix, MATULM, all = T)
FISHsystematix= merge(FISHsystematix, MATUTM, all = T)
FISHsystematix$TemperatureC=FISHsystematix$Temperature
FISHsystematix = FISHsystematix[,-2]
str(FISHsystematix)
FISHsystematix= merge(FISHsystematix, Taxa, by="Species", all = T)
table(FISHsystematix$Class)

my = c("Ceratodontiformes", "Chondrichthyes", "Actinopterygii", "Cladistia")
FISHsystematix = FISHsystematix[FISHsystematix$Class %in% my,]
FISHsystematix = FISHsystematix[!is.na(FISHsystematix$A_G),]
FISHsystematix = FISHsystematix[order(FISHsystematix$A_G),]

table(FISHsystematix$Class)
write.csv(FISHsystematix, file = '../../Body/3Results/VertebratePolymorphisms.AllActinopterygiiSystematix.csv', quote = FALSE)


install.packages("xlsx")
library("xlsx")


##########Obtaining median oxygen consumption

MetabolicData = read.csv("../../Body/2Derived/FishBaseMetabolicData.csv", sep = ";")
MetabolicData$OxyCon.mg.kg.h.= gsub(",", ".", MetabolicData$OxyCon.mg.kg.h.); MetabolicData$OnyCon.at20C.= gsub(",", ".", MetabolicData$OnyCon.at20C.)
str(MetabolicData)
MetabolicData$OxyCon.mg.kg.h.= as.numeric(MetabolicData$OxyCon.mg.kg.h.); MetabolicData$OnyCon.at20C.= as.numeric(MetabolicData$OnyCon.at20C.)
str(MetabolicData)
length(unique(MetabolicData$Spece)) #206
table(MetabolicData$Applied_stress)
table(MetabolicData$Activity)


HypoxicFishes = MetabolicData[MetabolicData$Applied_stress == "hypoxia",]
summary(HypoxicFishes)
OxConAgg = aggregate(HypoxicFishes$OxyCon.mg.kg.h., by=list(HypoxicFishes$Spece), FUN=median); names(OxConAgg) <- c("Species", "OxCon")
TempMerge = merge(OxConAgg, MUT)
nrow(TempMerge)
plot(TempMerge$OxCon, TempMerge$T_C)
abline(lm(formula = TempMerge$T_C ~ TempMerge$OxCon))
NonHypoxicFishes = MetabolicData[MetabolicData$Applied_stress != "hypoxia",]
summary(NonHypoxicFishes)


MetDataRoutine = MetabolicData[MetabolicData$Activity == "routine",]
OxConAgg = aggregate(MetDataRoutine$OxyCon.mg.kg.h., by=list(MetDataRoutine$Spece), FUN=median); names(OxConAgg) <- c("Species", "OxCon")
TempMerge = merge(OxConAgg, MUT)
nrow(TempMerge)
plot(TempMerge$OxCon, TempMerge$T_C)
abline(lm(formula = TempMerge$T_C ~ TempMerge$OxCon))
