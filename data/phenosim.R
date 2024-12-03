# Script used to generate the quantitative and binary traits in the example data provided with the package

wk_dir <- c("C:/Users/bross/OneDrive/Haplotypes_CLQD_MLC/RegionScan_paper/RegionScan_dev/GitHub2024_RegionScan/example/")
setwd(wk_dir)

# phenotype simulation for the illustrative example
library(devtools)
load(paste0(wk_dir, "geno.rda"))
load(paste0(wk_dir, "SNPinfo.rda"))
load(paste0(wk_dir, "REGIONinfo.rda"))
geno.rda<-as.data.frame(cbind(ID=1:nrow(geno),geno))
SNPinfo<-as.data.frame(SNPinfo)
# causal SNPs (from region 8, in two different LD bins)
causal1<-"chr16.46880510.A.G" 
causal2<-"chr16.46889594.G.A"
cor(geno.rda[,c(causal1,causal2)])
subset(SNPinfo,variant==causal1)
subset(SNPinfo,variant==causal2)
y1<-0.9+0.065*geno[,causal1]+0.15*geno[,causal2]+rnorm(nrow(geno),mean=1.5, sd=1.5)
summary(glm(as.formula(paste0("y1~",causal1,"+",causal2)),data=cbind(y1,geno)))
phenocov <- as.data.frame(cbind(ID=1:nrow(geno),sim_QT=y1,sim_bin=ifelse(y1>=quantile(y1,0.25),1,0)))

# saving data
save(geno,file="geno.rda")
save(SNPinfo,file="SNPinfo.rda")
save(REGIONinfo,file="REGIONinfo.rda")
save(phenocov,file="phenocov.rda")