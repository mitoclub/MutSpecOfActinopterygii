rm(list=ls(all=TRUE))


library(dplyr)
library(tidyr)
VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')

unzip("../Data/1raw/POLARIZEDBR_DATA.zip")
List = list.files("../Data/1raw/POLARIZEDBR_DATA/")


for (i in 1:length(List)){
  infile = paste("../Data/1raw/POLARIZEDBR_DATA/",as.character(List[i]),sep='') 
  if (length(grep('POLARISED',infile)) > 0)
  {
    Species = gsub('\\..*','',as.character(List[i]))
    Gene = gsub(Species,'',as.character(List[i])); Gene = gsub('.POLARISED.txt','',Gene); Gene = gsub('\\.POLARISED.txt','',Gene); Gene = gsub('\\.','',Gene); 
    GeneSpecies = read.table(infile, header = TRUE)
    GeneSpecies = GeneSpecies[GeneSpecies$BranchPosition == 'External',]
    NumOfSeqSpXG = nrow(GeneSpecies)
    ExternalSeqsTogether = paste(GeneSpeciess$MoreShallowNodeSeq,collapse = '')
    ExternalSeqsTogether = unlist(strsplit(ExternalSeqsTogether,'')) # 5700/3
    CodonsVec = c(); StartNuc = 1
    if (length(ExternalSeqsTogether)/3 == round(length(ExternalSeqsTogether)/3))  # if divide by 3 without the rest
    {
      for (j in 1:(length(ExternalSeqsTogether)/3))
      {
        CodonsVec = c(CodonsVec,paste(ExternalSeqsTogether[StartNuc : (StartNuc+2)],collapse = ''))
        StartNuc = StartNuc+3
      }
      AllCodons = length(CodonsVec)        # 1021
      CodonsVecNeutral = CodonsVec[CodonsVec %in% VecOfSynFourFoldDegenerateSites]
      NeutralCodons = length(CodonsVecNeutral) # 1900
      data.frame(table(CodonsVecNeutral))
      
      CodonsVecNeutral = gsub("CTA|GTA|TCA|CCA|ACA|GCA|CGA|GGA",'A',CodonsVecNeutral)
      CodonsVecNeutral = gsub("CTT|GTT|TCT|CCT|ACT|GCT|CGT|GGT",'T',CodonsVecNeutral)
      CodonsVecNeutral = gsub("CTG|GTG|TCG|CCG|ACG|GCG|CGG|GGG",'G',CodonsVecNeutral)
      CodonsVecNeutral = gsub("CTC|GTC|TCC|CCC|ACC|GCC|CGC|GGC",'C',CodonsVecNeutral)
      
      Line=c(Species,Gene,length(CodonsVecNeutral[CodonsVecNeutral == 'A']),
             length(CodonsVecNeutral[CodonsVecNeutral == 'T']),
             length(CodonsVecNeutral[CodonsVecNeutral == 'G']),
             length(CodonsVecNeutral[CodonsVecNeutral == 'C']),
             AllCodons, NeutralCodons, NumOfSeqSpXG)
      if (i == 1) {Final = Line}
      if (i >  1) {Final = rbind(Final,Line)}
    }
  }
}

Final = as.data.frame(Final); names(Final)=c('Species','Gene','A','T','G','C','NumberOfAllCodons',"NumberOfFourFoldDegenCodons")
write.table(Final, "../Data/3results/VertebratePolymorphisms.Normalization.NeutralATGC.txt", quote = FALSE, row.names = FALSE)




################### merge with classes (from Taxa & MoreTaxa), average A T G C for each species, average it for classes and draw it.
#### associate species name with Class
### Taxa 1, Cut out the third word!!!!!!!!!!!!!!!!!
Taxa = read.table("../Data/1raw/TaxaFromKostya.Names.stat", sep = '\t',header = FALSE) 
Taxa$Species = gsub(";.*",'',Taxa$V1); 
for (i in (1:nrow(Taxa)))  {Taxa$Species[i] = paste(unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[1],unlist(strsplit(as.character(Taxa$Species[i]),split = ' '))[2], sep = '_')}
Taxa$Class = gsub(";Chordata;.*",'',Taxa$V1); Taxa$Class = gsub(".*;",'',Taxa$Class); table(Taxa$Class)
Taxa$Class = gsub('Actinopteri','Actinopterygii',Taxa$Class)
Taxa$Class = gsub("Testudines|Squamata|Crocodylia|Sphenodontia",'Reptilia',Taxa$Class)
length(unique(Taxa$Species)) # 1708
table(Taxa$Class)
Taxa = Taxa[,-1]

### Taxa 2, Cut out the third word!!!!!!!!!!!!!!!!!
TaxaMore = read.table("../Data/1raw/TaxaFromKostya.2NeedTaxa.tax.txt", sep = '\t',header = FALSE) 
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

###

Final = merge(Final,Taxa, all = TRUE)

NeutralATGC = Final
#####  COUNT NUCLEOTIDE CONTENT CAREFULLY (/1raw/polarizedbr_data => external + More Shallow => codons, 4fold nucl, FrA,T,G,C) => barplot?
#####  normalized average MutSpec (pie charts or 12 boxplots for each class)?


MUT = read.table("../Data/1raw/Mutational_spectra_in_Chordata_ML.txt", header = TRUE)
length(unique(MUT$Species)) # 2404

##### FILTER 1: to take only normal substitutions and filter out species with too high fraction (> 5%) of unnormal substitutions
VecOfNormalSubstitutions <- c('A_C','C_A',
                              'A_G','G_A',
                              'C_G','G_C',
                              'C_T','T_C',
                              'G_T','T_G',
                              'T_A','A_T')

SP = data.frame(table(MUT$Species)); names(SP) = c('Species','NumberOfAllSubst')
SPN = data.frame(table(MUT[MUT$Subs %in% VecOfNormalSubstitutions,]$Species)); names(SPN) = c('Species','NumberOfNormalSubst')
SP = merge(SP,SPN); SP$FractionOfNormal = SP$NumberOfNormalSubst/SP$NumberOfAllSubst

summary(SP$FractionOfNormal) # how many to delete? ask to have more than 95% of substitutions as normal
SpeciesToDelete = SP[SP$FractionOfNormal <=0.95,]$Species; length(SpeciesToDelete) # 279 - delete
MUT = MUT[!MUT$Species %in% SpeciesToDelete,]
MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

##### FILTER 2: Synonymous Substitutions
nrow(MUT) # 461215
MUT = MUT[MUT$AncestralAA == MUT$DescendantAA,]; nrow(MUT) # 395157


MUT4fold = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]; nrow(MUT4fold) # 209120
MUT4fold$MutType = 'FourFold'
MUT = MUT4fold
MUT$Number = 1

# Count number of different mutations
AGG = aggregate(MUT$Number, by = list(MUT$Species, MUT$Gene), FUN = sum); names(AGG) = c('Species','Gene','NumOfFourFoldMut')
AGG[is.na(AGG)] <- 0

MUT = merge(MUT,NeutralATGC[,c(1:6)], by = c("Species", "Gene"))
MUT = merge(MUT,AGG, by = c('Species', "Gene"))


mut_count = aggregate(MUT$Number, by=list(MUT$Species, MUT$Gene, MUT$Subs), FUN=sum); names(mut_count) = c('Species','Gene','Subs','Freq')
mut_count =
  mut_count %>% 
  pivot_wider(names_from = Subs, values_from = Freq, values_fill = 0)

MUT = merge(MUT,mut_count, by = c('Species', "Gene"))

MUT$A = as.numeric(MUT$A)
MUT$T = as.numeric(MUT$T)
MUT$G = as.numeric(MUT$G)
MUT$C = as.numeric(MUT$C)

MUT$FrA = MUT$A / (MUT$A + MUT$T + MUT$G + MUT$C)
MUT$FrT = MUT$T / (MUT$A + MUT$T + MUT$G + MUT$C)
MUT$FrG = MUT$G / (MUT$A + MUT$T + MUT$G + MUT$C)
MUT$FrC = MUT$C / (MUT$A + MUT$T + MUT$G + MUT$C)

subset(MUT, select = -c(AncestralSeqName, DescendantSeqName,Branch,MutType,Number))

#write.table(MUT, file = '../Data/3results/VertebratePolymorphisms.MutSpecData.txt', quote = FALSE)





