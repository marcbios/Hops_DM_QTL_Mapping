# Perform Data Formatting For RQTL 
rm(list=ls())
library(stringr)
#load("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/DM.Phenotypic.BLUPs.07232021.RData")
DM.BLUPs <- read.csv("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/Phenotypic_Analysis/DM.Parent_Offspring_Pheno_BLUP.12032021.csv", header=T, stringsAsFactors = F)
head(DM.BLUPs)
DM.BLUPs <- DM.BLUPs#[,c(1,2,3,4,5,6,7)]
#str_pad(x, 8, pad = "0")
head(DM.BLUPs)
colnames(DM.BLUPs)[1] <- "id"

geno <- read.csv("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/Epistasis/DM.rqtl2.geno.csv", header = T, stringsAsFactors = F)
geno[1:6,1:6]
geno.names <- data.frame(geno$id)
head(geno.names)
colnames(geno.names) <- "id"

# Find common genotypes between phenotypic and genotypic data
com.taxa <- merge(DM.BLUPs, geno.names, by="id")
head(com.taxa)
colnames(com.taxa) <- c("id", "Spor_BLUP", "Spor_BLUP_Int", "Spor_Mean", "Lesion_BLUP", "Lesion_BLUP_Int", "Lesion_Mean")

map <- read.csv("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/Epistasis/DM.rqtl2.gmap.csv", header = T, stringsAsFactors = F)
head(map)

geno[1:6,1:6]
geno2 <- geno[match(com.taxa$id, geno$id),]
geno2[1:6,1:6]

head(map)
# match marker names between map and geno2
dim(geno2)
geno2[1:6,1:6]
rownames(geno2) <- geno2$id
geno3 <- geno2[,match(map$marker, colnames(geno2))]
dim(geno3)
geno3[1:6,1:6]
class(geno3)
geno3 <- as.matrix(geno3)
geno3[geno3=="AA"] <- "A"
geno3[geno3=="AB"] <- "H"
geno3[geno3=="BB"] <- "B"
geno3[1:6,1:6]
geno3 <- as.data.frame(geno3)
map.t <- t(map)
map.t[1:3,1:6]
colnames(map.t) <- as.character(as.vector(as.matrix(map.t[1,])))
map.t[1:3,1:6]
geno4 <- rbind(map.t, geno3)

geno4[1:6,1:6]
geno4$id <- rownames(geno4)
geno4[1:6,1:10]

geno.pheno <- merge(com.taxa, geno4, all=T)
geno.pheno[1:6,1:10]
rownames(geno.pheno) <- geno.pheno$id
geno.pheno[1:15,1:6]
geno.pheno <- as.data.frame(geno.pheno)
geno.pheno[1:6,1:6]
map.xters <- c("marker", "chr", "pos")
taxa.names <- rownames(geno.pheno)[(!rownames(geno.pheno)%in%map.xters)]
geno.pheno2 <- geno.pheno[c(map.xters, taxa.names),]
geno.pheno3 <- geno.pheno2
geno.pheno3[is.na(geno.pheno3)] <- ""
geno.pheno4 <- geno.pheno3[,-1]
geno.pheno5 <- geno.pheno4
geno.pheno5[1,] <- colnames(geno.pheno5)
geno.pheno5[1:6,1:10]
write.table(geno.pheno5, "DM.rqtl1.data.PhenoReanalysis.12102021.csv", sep=",", quote = F, row.names = F, col.names = F) # RQTL CSV format

# Results are in "/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis"

############# Perform QTL Analysis  ######################
rm(list = ls())
date="12102021"
#install.packages("qtl")
#install.packages("qtlcharts")
#install.packages("devtools")
#devtools::install_github("rcorty/vqtl")

library(qtl)
library(qtl2)
#library(vqtl)

mycross <- read.cross("csv", ".", "DM.rqtl1.data.PhenoReanalysis.12102021.csv", crosstype="f2")
mycross <- jittermap(mycross)
mycross <- calc.genoprob(mycross, step=0, map.function="haldane", stepwidth="fixed", error.prob=0.01)#error.prob=0.002,

nind(mycross)
nchr(mycross)
totmar(mycross)
nmar(mycross)
nphen <- nphe(mycross)
nperm=100 # Specify the total number of permutations

phenames(mycross)

phenotype_names <- c(phenames(mycross))
phenotype_names
#mycross.scavar <- mycross
#mycross.mqm <- mycross
setwd("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/Results_DM_12102021")
for(j in 1:nphen){
  
  print(paste("--- performing composite interval mapping --", " Trait: ", phenotype_names[j], " ---", sep="")) 
  
  out_cim <- cim(mycross, pheno.col=j, n.marcovar=5, window=10, method="hk", imp.method="imp", error.prob=0.0001, map.function="haldane")
  
  if(max(out_cim$lod) == Inf){
  next
 }else{
   
   # Perform vairance and mean QTL mapping
   #out.scavar <- scanonevar(cross = mycross.scavar, pheno.col=j, mean_covar=NULL, var_covar=NULL, maxit=25, tol=1e-6, quiet=TRUE)
   
   #out.scavar <- scanonevar(cross = mycross.scavar, mean.formula = phenotype_names[j] ~  mean.QTL.add + mean.QTL.dom, var.formula = ~ var.QTL.add + var.QTL.dom)
   
   
   # use format="allpeaks" to get summary for each of mean and variance # also consider format="tabByCol" or format="tabByChr" 
   #summary(out.scavar, format="allpeaks")
   
   # Perform permutations for the mean part of the model
   #mean.perm <- scanonevar.meanperm(cross=mycross.scavar, pheno.col=j, mean_covar=NULL, var_covar=NULL, maxit=25, tol=1e-6, n.mean.perm = 100, seed = 27517, quiet=TRUE)
   
   # Perform permutations for the variance part of the model
  # var.perm <- scanonevar.varperm(cross=mycross.scavar, pheno.col=j, mean_covar=NULL, var_covar=NULL, maxit=25, tol=1e-6, n.var.perm = 100, seed = 27517, quiet=TRUE)
   
   # Plot the meanQTL and varQTL results
  # color <- c("slateblue", "violetred")
   
  # pdf(paste("Downy_mvQTL", phenotype_names[j],".pdf", sep=""), width=15, height=5, paper='special', )
  # plot(out.scavar, lod=1:2, col=color, bandcol="white") 
   #legend("topright", lwd=2, c("mean", "variance"), col=color)
   #add.threshold(out.scavar, chr=c(1,2,3,4,5,6,7,8,9,10), perms=mean.perm, alpha=0.1, col="forestgreen", lty=1)
   #add.threshold(out.scavar, chr=c(1,2,3,4,5,6,7,8,9,10), perms=var.perm, alpha=0.05, col="oragnge", lty=1)
   #dev.off()
   
   #vQTLs_Phen=summary(out.scavar, perm=mean.perm, lodcolumn=1, alpha=0.10)
   #write.table(QTLs_Phen, paste("Downy.QTL", phenotype_names[j],"txt", sep="."), quote=F,sep="\t")
   
   
   
   ## QTL Results in txt format
    print(paste("--- Exporting composite interval mapping result for --", phenotype_names[j], " ---", sep=""))
    write.table(out_cim, paste("CIM.Markers.Output", phenotype_names[j], date, "txt", sep="."), quote=F,sep="\t")
   
    ## 100 or 1000 Permutations
    print(paste("--- performing ", nperm, " permutations --", " Trait: ", phenotype_names[j], " ---", sep=""))
    operm <- cim(mycross, pheno.col=j, method="hk", n.perm=nperm)
    
    ## Make QTL plot
    print(paste("--- Exporting composite interval mapping QTL plot for --", phenotype_names[j], " ---", sep=""))
    pdf(paste("Downy_QTL", phenotype_names[j], date, ".pdf", sep=""), width=15, height=5, paper='special')
    
    plot(out_cim, col=c("blue", "red"), main=paste("Downy_QTL", phenotype_names[j], sep=" "))
    add.cim.covar(out_cim)
    
    add.threshold(out_cim, perms=operm, alpha=0.05, col="forestgreen", lty=1)
    add.threshold(out_cim, perms=operm, alpha=0.10, col="orange", lty=1)
    add.threshold(out_cim, perms=operm, alpha=0.25, col="red", lty=1)
    add.threshold(out_cim, perms=operm, alpha=0.50, col="red", lty=2)
    add.threshold(out_cim, perms=operm, alpha=0.75, col="red", lty=3)
    add.threshold(out_cim, perms=operm, alpha=0.80, col="red", lty=3)
    dev.off()
    
    QTLs_Phen=summary(out_cim, perm=operm, lodcolumn=1, alpha=0.80)
    write.table(QTLs_Phen, paste("Downy.QTL", phenotype_names[j], date,"txt", sep="."), quote=F,sep="\t")
    
    ########################################################################################
    # Estimate the proportion of phenotypic variation explained by the QTL
    print(paste("--- estimating the proportion of variation explained by QTL for --", phenotype_names[j], " ---", sep=""))
      postn <- as.vector(QTLs_Phen$pos)
      chrom <- as.vector(QTLs_Phen$chr)
      
      mnames_Phen <- find.marker(mycross, chr=chrom, pos=postn)
      
      fake.f2 <- subset(mycross, chr=chrom)
      fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
      
      qtl.info <- makeqtl(fake.f2, chr=chrom, pos=postn, qtl.name=mnames_Phen, what="prob")
      
      # fit an additive QTL model
      if(nrow(QTLs_Phen)==1){
        lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1, method="hk")
        } else if(nrow(QTLs_Phen)==2){
          lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2, method="hk")
          } else if(nrow(QTLs_Phen)==3){
        lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2+Q3, method="hk")
        } else if(nrow(QTLs_Phen)==4){
        lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2+Q3+Q4, method="hk")
        } else if(nrow(QTLs_Phen)==5){
        lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2+Q3+Q4+Q5, method="hk")
        }else if(nrow(QTLs_Phen)==6){
          lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2+Q3+Q4+Q5+Q6, method="hk")
        }else lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7, method="hk")
      ###########################################################################################################
      
      
      ##################################################################
      ## Make effect plot for each QTL
      
      for(q in 1:length(mnames_Phen)){
        print(paste("--- creating effect plot for --", phenotype_names[j], " QTL: ", mnames_Phen[q] ,"---",sep=""))
        pdf(paste("Efffect_Plot", phenotype_names[j], mnames_Phen[q], date, "pdf", sep = "."), 6,6)
          effectplot(mycross, pheno.col = j, mnames_Phen[q])
          dev.off()
      }
      
      ##################################################################
      
      
      #lod.add <- fitqtl(fake.f2, pheno.col=j, qtl=qtl.info, formula=y~Q1+Q2, method="hk")
      summary(lod.add)
      res.df <- data.frame(lod.add$result.drop)
      res.df$QTL <- rownames(res.df)
      QTLs_Phen$QTL <- rownames(QTLs_Phen)
      Res.QTL <- merge(QTLs_Phen, res.df)
      write.table(Res.QTL, paste("Summary.QTL.Result", phenotype_names[j], date,"csv", sep="."), sep=",", quote = F, row.names = F, col.names = T)
      
      #residuals <- attr(lod.add, "residuals")
      #plot(residuals)
      
      # Perform MQM
      #augmentedcross <- mqmaugment(mycross.mqm, minprob=1.0)
      #dm.cof <- mqmautocofactors(augmentedcross, num=10, distance=5, dominance=FALSE, plot=TRUE, verbose=TRUE)
      
      #mqm_auto <- mqmscan(augmentedcross, dm.cof)
      #setcofactors <- mqmsetcofactors(augmentedcross, 5)
     # mqm_backw <- mqmscan(augmentedcross, setcofactors)
      
     # mqm_dm1 <- mqmscan(augmentedcross)
      
     # maug_min1 <- mqmaugment(multitrait, minprob=1.0)
      #mqm_min1 <- mqmscan(maug_min1)
      
      #data(multitrait)
      #str(multitrait)
      #str(mycross.mqm)
   }                         
 }
    
    
#}

## Make new effect boxplots with backtransformed boxcox phenotypes
rm(list = ls())
traits.names <- c("Percent_Sporulation.BC", "Percent_Lesion.BC")
home.dir <- "~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/Phenotypic_Analysis/"
pheno <- read.csv(paste(home.dir, "least_squared_means_dm_bc.phenotypes_06142021.csv", sep=""), header = T, stringsAsFactors = F)
#colnames(pheno)[1] <- "id"
geno <- read.csv("hops.rqtl2.geno.csv", header = T, stringsAsFactors = F)

geno.names <- data.frame(geno$id)
head(geno.names)
colnames(geno.names) <- "id"

# Find common genotypes between phenotypic and genotypic data
com.taxa <- merge(pheno, geno.names, by="id")
head(com.taxa)
rownames(com.taxa) <- com.taxa$id

geno[1:6,1:6]
geno2 <- geno[match(com.taxa$id, geno$id),]
geno2[1:6,1:6]


for (ph in traits.names) {
  trt.pve <- read.csv(paste("Summary.QTL.Result", ph, "csv", sep="."), header = T, stringsAsFactors = F)
  trt.pve$PVE <- round(trt.pve$PVE, 1)
  trt.pheno <- com.taxa[paste(ph, "Inv", sep=".")]
  
        for (q in 1:nrow(trt.pve)){
          
          pdf(paste("Effect.Boxplot", colnames(trt.pheno), (trt.pve$QTL)[q], "pdf", sep="."), 5,7)
          par(oma=c(4,4,4,1), mgp=c(2.4, 0.8, 0), las=1)
          boxplot(trt.pheno[,1] ~ geno2[,match((trt.pve$QTL)[q], colnames(geno2))], main=paste("PVE: ", (trt.pve$PVE)[q], "%", sep=""), xlab = (trt.pve$QTL)[q], ylab=gsub(".BC.Inv", "", colnames(trt.pheno)), cex=1, cex.lab=1.2, cex.axis=1, cex.main=1.5)
          dev.off()
          
          }
  
}


###################################################################################################################################



###################### Pleiotropy
#install.packages("pleiotest")
rm(list = ls())
library(pleiotest)
setwd("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis")
# Load and format both pheno and geno data
load("~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/DM.Phenotypic.BLUPs.07232021.RData")

head(DM.BLUPs)
#str_pad(x, 8, pad = "0")

geno <- read.csv("hops.rqtl2.geno.csv", header = T, stringsAsFactors = F)
geno[1:6,1:6]
geno.names <- data.frame(geno$id)
head(geno.names)
colnames(geno.names) <- "id"

# Find common genotypes between phenotypic and genotypic data
com.taxa <- merge(DM.BLUPs, geno.names, by="id")
head(com.taxa)

map <- read.csv("hops.rqtl2.gmap.csv", header = T, stringsAsFactors = F)
head(map)

geno[1:6,1:6]
geno2 <- geno[match(com.taxa$id, geno$id),]
geno2[1:6,1:6]

head(map)
# match marker names between map and geno2
dim(geno2)
rownames(geno2) <- geno2$id
geno3 <- geno2[,match(map$marker, colnames(geno2))]
dim(geno3)
geno3[1:6,1:6]
class(geno3)

geno3[geno3=="A"] <- 0
geno3[geno3=="B"] <- 2
geno3[geno3=="H"] <- 1

class(geno3)
geno3[1:6,1:6]
geno3 <- as.matrix(geno3)

head(com.taxa)
P.Spor <- com.taxa[,c(1,3)]
head(P.Spor)
P.Spor$trait <- rep("t1", nrow(P.Spor))
colnames(P.Spor)[2] <- "y"
P.Lens <- com.taxa[,c(1,5)]
P.Lens$trait <- rep("t2", nrow(P.Lens))
colnames(P.Lens)[2] <- "y"
mt.traits <- rbind(P.Spor, P.Lens)
head(mt.traits)
mt.traits <- mt.traits[,c(1,3,2)]
head(mt.traits)
str(mt.traits)
class(geno3) <- "numeric"

pleio_model.dm <- pleioR(pheno = mt.traits, geno = geno3)

head(map)
bp.positions <- map[,-1]
rownames(bp.positions) <- map$marker
pleio.mt.dm <- mt_gwas(pleio_model.dm, save_at = NULL)#"~/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/"

pleio.test.dm <- pleio_test(pleio_model.dm, loop_breaker = 1, save_at = NULL, contrast_matrices_list = NULL)
pleio.test.dm[[2]]
pleio.test.sim[[2]]
str(pleio.test.dm)
str(pleio.test.sim)
str(pleio.mt.dm)

pleio_plot(pleio.mt.dm, alpha = "bonferroni05", n_traits = 2, bp_positions = bp.positions,
            set_colors = NULL,set_text = NULL,set_plot = TRUE, chr_spacing = 1e+05)

source("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/pleio.plot.ammended.R")
pleio_plot.adj(pleio.mt.dm, alpha = "bonferroni05", n_traits = 2, bp_positions = bp.positions,
           set_colors = NULL,set_text = NULL,set_plot = TRUE, chr_spacing = 1e+05)

source("/Users/marcus.olatoye/OneDrive - USDA/Desktop/Files/Downy_Mildew/Data/Downy_Mildew/Manuscript_2021/DM_QTL_Analysis/pleio.plot.ammended.R")
pleio_plot.adj( pleio.test.dm, alpha = "bonferroni05", n_traits = 2, bp_positions = bp.positions,
  set_colors = NULL,set_text = NULL,set_plot = TRUE, chr_spacing = 1e+05)

head(pleio.test.dm)

p_values <- apply(pleio.test.dm[[1]][, 1:n_traits, drop = F], 1, max)

pleio.test.dm[[2]]

####################################
rm(list = ls())

sim1 <- pleio_simulate(n_traits = 3, n_individuals = 1e4, n_snp = 1e3, percentage_mv = 0.1)
str(sim1$pheno)

dim(sim1$pheno)
head(sim1$pheno)
#pheno.sim <- sim1$pheno[which(sim1$pheno$trait!="t3"),]
dim(sim1$geno)
sim1$geno[1:6,1:6]

pleio_model <- pleioR(pheno = sim1$pheno, geno = sim1$geno)
pleio.sim.mt <- mt_gwas(pleio_model, save_at = NULL)
pleio.test.sim<- pleio_test(pleio_model, loop_breaker = 1, save_at = NULL, contrast_matrices_list = NULL)
pleio.test.sim[[2]]
pleio_plot( pleio.test.sim, alpha = "bonferroni05", n_traits = 3, bp_positions = NULL, set_colors = NULL,set_text = NULL,set_plot = TRUE, chr_spacing = 1e+05)



pleio_model_test <- pleio_test(pleio_model)



