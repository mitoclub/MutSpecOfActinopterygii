rm(list=ls(all=TRUE))

if (!require(ggpubr)) install.packages("ggpubr")
if (!require(vioplot)) install.packages("vioplot")
library("ggpubr")
library(vioplot)

###########Reading mutspec data
mutSpec = read.table('../Data/3results/VertebratePolymorphisms.MutSpecData.txt', header = TRUE)

########AnAge data reading
AA = read.table("../Data/1raw/anage_data.txt", header = TRUE, sep = '\t')
AA$Species = paste(AA$Genus,AA$Species,sep = '_')
AA$temperature = AA$Temperature..K. - 273.15

BirbTemp =  read.table("../Data/1raw/body_temp.txt", header = TRUE)
BirbTemp$Class = "Aves"

ALL = merge(mutSpec, AA, all = T)
ALL = merge(ALL, BirbTemp, all = T)

####U tests
colder <- c("Actinopterygii", "Amphibia", "Reptilia")
warmer <- c("Mammalia","Aves")
wilcox.test(MUTFROMK[MUTFROMK$Class %in% colder,]$T_C, MUTFROMK[MUTFROMK$Class %in% warmer,]$T_C, paired=F) #p-value = 2.574e-13
wilcox.test(MUTFROMK[MUTFROMK$Class %in% colder,]$A_G, MUTFROMK[MUTFROMK$Class %in% warmer,]$A_G, paired=F)  #p-value < 2.2e-16
wilcox.test(MUTFROMK[MUTFROMK$Class %in% colder,]$TCdivAG, MUTFROMK[MUTFROMK$Class %in% warmer,]$TCdivAG, paired=F) #p-value < 2.2e-16
wilcox.test(Alltemp[Alltemp$Class %in% colder,]$Temperature, Alltemp[Alltemp$Class %in% warmer,]$Temperature, paired=F) #p-value < 2.2e-16

pdf('../Figures/PolymorphicData.Actinopterygii.AverageMutSpec.Violin.pdf', width = 10, height = 5.3)
ggviolin(ALL, x = "Class", y = "temperature", select = c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "Body temperature, ?C",
         order=c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), add = "boxplot", fill="Class", palette=c("#6760db", "#7849bf", "#9145c4", "#c73a69", "#c2464c"),  ylim = c(0, 50))
ALLDIV = ALL[!is.na(ALL$A_G.N),]; ALLDIV$TCdivAG=ALLDIV$T_C.N/ALLDIV$A_G.N; ALLDIV[ALLDIV$TCdivAG == "Inf" | ALLDIV$TCdivAG == "NaN" | ALLDIV$TCdivAG == 0,]$TCdivAG  = NA; ALLDIV$TCdivAG = log2(ALLDIV$TCdivAG)
ggviolin(ALLDIV, x = "Class", y = "TCdivAG", select = c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), ylab = "AH>GH div TH>CH",
         order=c("Actinopterygii", "Amphibia", "Reptilia", "Mammalia","Aves"), add = "boxplot", fill="Class", palette=c("#6760db", "#7849bf", "#9145c4", "#c73a69", "#c2464c"))
dev.off()

table(ALL[!is.na(ALL$temperature),]$Class)
table(ALLDIV[!is.na(ALLDIV$TCdivAG),]$Class)

