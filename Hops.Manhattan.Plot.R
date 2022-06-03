rm(list=ls())
setwd("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis")

get_genome_pos <- function( chr_num, position, buffer=40000000){
  # chr_num: vector with chromosome numbers
  # position: vector with positions, corresponding to each chr_num
  # buffer: space between chromosomes in bp
  # Chromosome lengths for sorghum:
  chr_length <- c(
    chromosome_1  = 128419000, # 467,464,186
    chromosome_2  = 94254000,
    chromosome_3  = 125717000,  
    chromosome_4  = 121541000,
    chromosome_5  = 110051000,
    chromosome_6  = 119981000,
    chromosome_7  = 164823000,
    chromosome_8  = 84314000,
    chromosome_9  = 135052000,
    chromosome_10 = 185377000 
  )
  chr_length <- chr_length + buffer
  
  position + sapply( chr_num, function(x) sum(c(0,chr_length[-length(chr_length)])[ 1:x ]))
}

map <- read.delim("GAPIT.Genotype.map.txt", header = T, stringsAsFactors = F)
head(map)

trt.qtl <- read.csv("DM.Stepwise.Result.12102021.BLUPs.Only.csv", header = T, stringsAsFactors = F)
head(trt.qtl)
trt.qtl <- trt.qtl[!(trt.qtl$Name%in%c("mean", "Error")),]
head(trt.qtl)
phenames <- unique(trt.qtl$Trait) 
#Remove the 0 effect QTL from Percent Lesion S2_232509284
#trt.qtl <- trt.qtl[which(trt.qtl$Name!="S2_232509284"),]


geno_data <- read.delim("GAPIT.Genotype.map.txt", header = T, stringsAsFactors = F)
geno_data$pvalue <- rep(0, nrow(geno_data))



# There is no need for imputation since our GBS data had been previously imputed
#impute=A.mat(df2,max.missing=0.5,impute.method="mean",return.imputed=T)

#a.priori <- read.csv("/Volumes/Seagate Expansion Drive/Mac Files/PostDoc/Yasser/Maize.Panzea/mQTL_GeneEnrichment.1Mb.csv", header = T, stringsAsFactors = F)
#head(a.priori)
#a.priori$genome_pos <- get_genome_pos(a.priori$Chromosome, a.priori$Position)
#head(a.priori)

pdf("DM.StepwiseRegression.QTL.Plot.01182022.pdf", 9,18)
par(mfrow=c(5,1))

for(p in 1:length(phenames)){
  
  LPBL_JL <- trt.qtl[which(trt.qtl$Trait==phenames[p]),c(2:4,9)]
  colnames(LPBL_JL) <- c("SNP", "Chromosome", "Position", "pvalue")
  print(paste("---- Plotting Manhattan Plot For: ", phenames[p], "!!!!!!!!!---------------", sep = ""))
  geno_data2 <- rbind(geno_data, LPBL_JL)
  head(LPBL_JL)
  head(geno_data2)
  
  geno_data2$pvalue <- as.numeric(as.character(geno_data2$pvalue))
  
  geno_data3 <- geno_data2[with(geno_data2, ave(pvalue, SNP, FUN=max)==pvalue),]
  
  head(geno_data3)
  str(geno_data3)
  for(i in 2:ncol(geno_data3)){
    geno_data3[,i] <- as.numeric(geno_data3[,i])
  }
  LPBL_JL <- geno_data3
  colnames(LPBL_JL) <- c("snp", "chr", "pos", "pvalue")
  LPBL_JL <- LPBL_JL[order(LPBL_JL$chr, LPBL_JL$pos),]
  
  LPBL_JL$genome_pos <- get_genome_pos(LPBL_JL$chr, LPBL_JL$pos)
  head(LPBL_JL)
  LPBL_JL$logpvalue <- -log10(LPBL_JL$pvalue)
  tail(LPBL_JL)
  LPBL_JL$logpvalue[which(LPBL_JL$logpvalue==Inf)] <- 0
  
  #a.priori.trt <- a.priori[which(a.priori$Trait==phenames[p]),]
  marker.eff <- read.delim(paste("Hops.DM.QTL.Effects.PhenoBLUPs.12102021_", phenames[p], "_QTL.txt", sep = ""), header = T)
  # Hops.DM.QTL.Effects_Percent_Sporulation_BLUP_QTL.txt
  head(marker.eff)
  marker.eff <- marker.eff[,-1]
  colnames(marker.eff)[1] <- "snp"
  mars <- c(5,5,5,2)
  cexp=2
  cexs=2
  par(mar=mars)
  vek4<-as.numeric(LPBL_JL$chr)%%2
  
  mean_genome_pos <- tapply(LPBL_JL$genome_pos, LPBL_JL$chr, function(x) (min(x)+max(x))/2)
  
  plot( I(LPBL_JL$logpvalue) ~ LPBL_JL$genome_pos, cex=2, col=c("gray45","gray12")[as.factor(vek4)], 
        ylab=expression(paste('-log'[10],(italic(p)))), cex.lab=cexs, cex.main=cexs, cex.axis=cexs, pch=20,
        xlab='', xaxt='n', main='', frame.plot = FALSE, ylim=c(0,8))
  LPBL_JL2 <- LPBL_JL[which(LPBL_JL$logpvalue!=0),]
  LPBL_JL3 <- merge(LPBL_JL2, marker.eff)
  LPBL_JL3 <- LPBL_JL3[order(LPBL_JL3$genome_pos),]
  cexp2 <- LPBL_JL3$PVE#2
  
  points(I(LPBL_JL2$logpvalue) ~ LPBL_JL2$genome_pos, cex=cexp2, col=c("gray45","gray12")[as.factor(vek4)], pch=20)
  axis(side=1, at=mean_genome_pos, labels=1:10, tick=F, cex.axis=cexs)
  #abline(h=6, lty='dashed', col="red", lwd=cexs)
  #abline(v=c(CandGene_LPBL$genome_pos), lty='dashed', col="black", lwd=1)
  
  
  #plot(marker.eff$Pos,-log10(marker.eff$u))
  map2 <- map
  colnames(map2) <- c("snp", "chr", "pos")
  marker.eff2 <- merge(marker.eff, map2, all=T)
  head(marker.eff2)
  marker.eff2 <- marker.eff2[,-c(3,4)]
  head(marker.eff2)
  marker.eff2$Additive.Effect[is.na(marker.eff2$Additive.Effect)] <- 0
  head(marker.eff2)
  #colnames(EffectsLPBL_JL)[1] <- "snp"
  #EffectsLPBL_JL$Pos_num <- 1:nrow(EffectsLPBL_JL)
  #head(EffectsLPBL_JL)
  snp_vec <- as.vector(as.matrix(marker.eff2$snp))
  
  
  marker.eff2 <- marker.eff2[order(marker.eff2$chr, marker.eff2$pos),]
  head(marker.eff2)
  
  EffectsLPBL_JL <- marker.eff2
  
  #EffectsLPBL_JL <- EffectsLPBL_JL[order(EffectsLPBL_JL$chr, EffectsLPBL_JL$pos),]
  
  EffectsLPBL_JL$genome_pos <- get_genome_pos(EffectsLPBL_JL$chr, EffectsLPBL_JL$pos)
  head(EffectsLPBL_JL)
  
  vek_jl<-as.numeric(EffectsLPBL_JL$chr)%%2
  vek2_jl<-ifelse(as.numeric(as.character(EffectsLPBL_JL[,"Additive.Effect"]))> 0,T,F)
  
  mean_genome_pos <- tapply(EffectsLPBL_JL$genome_pos, EffectsLPBL_JL$chr, function(x) (min(x)+max(x))/2)
  #EffectsLPBL_JL.nonzero <- EffectsLPBL_JL[which(EffectsLPBL_JL$Additive.Effect!=0),]
  
  head(EffectsLPBL_JL)
  par(mar=mars)
  cexs=2
  plot(x=EffectsLPBL_JL$genome_pos, y=as.numeric(as.character(EffectsLPBL_JL$Additive.Effect)), 
       type="h", xaxt='n', frame.plot = FALSE, main="", xlab="Chromosome", ylab= "Additive Effect (%)", 
       cex.lab=cexs, cex.main=cexs, cex.axis=cexs, col=c("gray45","gray12")[as.factor(vek4)], ylim=c(0,8), lwd=2) #ifelse(as.numeric(as.character(EffectsLPBL_JL$Additive.Effect)) <= 0,'black','red')
  #legend("bottomright",legend=c('Allele with positive effect relative to RTx430','Allele with negative effect relative to RTx430'),col=c("red", "black"),lty = 1,bty="n",lwd=2,cex=0.75)
  axis(side=1, at=mean_genome_pos, labels=1:10, tick=F, cex.axis=cexs)
  
}

dev.off()


#########################################################################
