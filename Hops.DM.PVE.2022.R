rm(list = ls())
#"/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis"
# Load and format both pheno and geno data

library(rrBLUP)
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(EMMREML)
library(compiler) 
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

dat=read.csv("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/DM.rqtl2.geno.csv", header = T, stringsAsFactors = F)

dat[1:6,1:6]
nan=dat[,1]
nan <- as.vector(as.matrix(nan))
dat <- dat[,2:ncol(dat)]
dat=as.matrix(dat)
dat[which(dat=="AA")]=2
dat[which(dat=="BB")]=0
dat[which(dat=="AB")]=1
rownames(dat)=nan
##
nam <- data.frame(nan)
names(nam) <- "Taxa"

df2 <- as.data.frame(dat)
df2 <- df2
df2[1:6,1:6]
df2 <- as.matrix(df2)

DM.BLUPs <- read.csv("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/Phenotypic_Analysis/DM.Parent_Offspring_Pheno_BLUP.12132021.csv", header=T, stringsAsFactors = F)
head(DM.BLUPs)

pheno <- DM.BLUPs[,c(1,4,7)]
head(pheno)
colnames(pheno) <- c("Taxa", "Percent_Sporulation_BLUP", "Percent_Lesion_BLUP")
nrow(pheno)
head(pheno)


phenames <- colnames(pheno[,-1])

nam <- data.frame(rownames(df2))
names(nam) <- "Taxa"
Pheno <- merge(pheno, nam, by="Taxa")
df2.copy <- df2 
df2 <- df2[match(Pheno$Taxa, rownames(df2)),]
df2[1:6,1:6]
df2 <- data.frame(df2)
df2 <- df2
df2[1:6,1:6]
df2 <- as.matrix(df2)
class(df2) <- "numeric"
df2[1:6,1:6]
JL_RES1 <- read.csv("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/Results_DM_12102021/Summary.QTL.Result.Lesion_BLUP.12102021.csv", header = T, stringsAsFactors = F)
JL_RES1$Trait <- rep("Percent_Lesion_BLUP", nrow(JL_RES1))
JL_RES2 <- read.csv("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/Results_DM_12102021/Summary.QTL.Result.Spor_BLUP.12102021.csv", header = T, stringsAsFactors = F)
JL_RES2$Trait <- rep("Percent_Sporulation_BLUP", nrow(JL_RES2))
JL_RES <- rbind(JL_RES1, JL_RES2)
head(JL_RES)
tail(JL_RES)
JL_RES$Trait <- as.character(JL_RES$Trait)
#JL_RES$SNP <- as.character(JL_RES$QTL)
JL_RES <- JL_RES[,c(12,1,2,3)]
head(JL_RES)
names(JL_RES)[2] <- "SNP"
head(JL_RES)
phenames <- c("Percent_Sporulation_BLUP", "Percent_Lesion_BLUP")
setwd("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/Results_DM_12102021/")

for (l in 1:length(phenames)){
  print(paste("-------------- Trait being analysed: ", phenames[l], "!!!!!!!!!!!---------------", sep = ""))
  
  ExplVar200Best <- JL_RES[which(JL_RES$Trait==phenames[l]),]
  
  bSNP<-df2[,as.character(ExplVar200Best$SNP)]
  
  phdata <- data.frame(Pheno[,1], Pheno[,phenames[l]])
  colnames(phdata)[2] <- phenames[l]
  colnames(phdata)[1] <- "Taxa"
  #sP<-as.data.frame(phdata[,phenames[l]])
  sP<-phdata
  rownames(sP) <- sP$Taxa
  da<-as.data.frame(cbind(sP, bSNP))
  
  trait_QTL_Pheno <- da
  write.table(t(data.frame(c("Trait", "QTL", "Additive Effect", "Dominance Effect", "PVE"))), paste("Hops.DM.rQTL.Effects.PhenoBLUP.12102021_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  #APV is the Among population variance in accordance to Wurschum et al. 2011 Heredity
  for(i in 3:ncol(trait_QTL_Pheno)){
    
    snp <- colnames(trait_QTL_Pheno)[i]  
    print(paste("-------------- Trait being analysed: ", phenames[l], "SNP: ", snp, "!!!!!!!!!!!---------------", sep = ""))
    
    trait_QTL_Pheno_2 <- trait_QTL_Pheno[,c(1,2,i)]
    
    AA_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==2),]
    AA <- mean(AA_class[,2], na.rm=T)
    BB_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==0),]
    BB <- mean(BB_class[,2], na.rm=T)
    
    AB_class <- trait_QTL_Pheno[which(trait_QTL_Pheno[,i]==1),]
    AB <- mean(AB_class[,2], na.rm=T)
    QTL_effect <- (AA-BB)/2
    
    Dom_effect <- AB-((AA+BB)/2)
    
    #formula.single <- as.formula(paste("Cd_comb ~ ",paste(as.character(topSNP$SNP), collapse=" + "), sep=" "))
    trait_QTL_Pheno_2$QTL <- trait_QTL_Pheno_2[,3]
    #QTL <- colnames(trait_QTL_Pheno_2[3])
    fin.anova <- lm(trait_QTL_Pheno_2[,phenames[l]] ~  QTL, data=trait_QTL_Pheno_2, na.action = na.omit)
    fin.sum <- summary(fin.anova)
    
    QVar <- round((fin.sum$adj.r.squared)*100, digits=2)#Pheno[,phenames[l]]
    
    print(paste("-------------- PVE For SNP: ", snp, "; Trait: ", phenames[l], " == ", QVar, "%  !!!!!!!!!!!---------------", sep = ""))
    write.table(t(data.frame(c(phenames[l], colnames(trait_QTL_Pheno[i]), round(abs(QTL_effect[1]), 1), round(abs(Dom_effect[1]), 1), QVar[1]))), paste("Hops.DM.rQTL.Effects.PhenoBLUP.12102021_", phenames[l],"_QTL",".txt", sep=""), sep="\t", append=T, quote=F, row.names=F, col.names=F)
  }
  
  
}


source("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/calc_snp_stats.R")

df2[1:6,1:6]
df.evachan <- t(df2) # SNP stat takes marker x individuals matrix
df.evachan[1:6,1:6]

QTL.summary <- calc_snp_stats(df.evachan)
head(QTL.summary)
QTL.summary$SNP <- rownames(QTL.summary)
QTL.summary <- QTL.summary[,c(14,5,6)]

hist(QTL.summary$maf)
summary(QTL.summary)
head(JL_RES)
head(QTL.summary)
SNP.summary <- merge(QTL.summary, JL_RES)
SNP.summary$maf <- round(SNP.summary$maf, 2)
SNP.summary
write.table(SNP.summary, "DM.CIM.SNP.MAF.12152021.csv", sep=",", quote=F, row.names = F, col.names = T)

################ Narrow sense heritability

A.Amat <- df2-1
#A.mat(geno.dip3,min.MAF=NULL,max.missing=NULL,impute.method="EM", n.core=4,shrink=FALSE,return.imputed=TRUE)
names(A.Amat)

diploid.G <- A.Amat#$imputed
diploid.G[1:6,1:6]

A.Amat.P <- A.mat(diploid.G, n.core=4,shrink=FALSE,return.imputed=FALSE)
names(A.Amat.P)
A.Amat.P[1:6,1:6]

h2<-NULL
phen <- Pheno
head(phen)

for(i in 2:dim(phen)[2]){
  
  y<-phen[,i]
  y[is.na(y)]<-mean(y,na.rm=T)
  model <- mixed.solve(y,K=A.Amat.P)
  h_sq <-model$Vu/(model$Vu+model$Ve)
  
  h2<-c(h2,h_sq)
  
}

h2<-data.frame(h2,traits=c(names(phen)[-1]))
h2
h2$h2 <- round(h2$h2,2)
h2


##
893-2072-500: 205,057.63
993-2072-500: 130,316.75
093-2072-500: 66,705.00




