setwd("/Users/TomDugal/Desktop/Project_results/SISG2019Data/")
library(gdsfmt)
library(SNPRelate)
map<-read.table(file="YRI_CEU_ASW_MEX_NAM.bim",sep="\t", header=FALSE,na="NA")
head(map)
dim(map)

pop_samp=read.table(file="Population_Sample_Info.txt",header = TRUE)

bed<-"YRI_CEU_ASW_MEX_NAM.bed"
fam<-"YRI_CEU_ASW_MEX_NAM.fam"
bim<-"YRI_CEU_ASW_MEX_NAM.bim"

snpgdsBED2GDS(bed, fam, bim, "hapmapdata_2.gds") #Conversion from PLINK BED to GDS

data <- snpgdsOpen("hapmapdata_2.gds")
pruned<- snpgdsLDpruning(data, ld.threshold=0.8) #LD-based SNP pruning at threshold 0.2
pruned.id<-unlist(pruned)

sample.id<-read.gdsn(index.gdsn(data, "sample.id"))
pop=pop_samp$Population

pca <- snpgdsPCA(data, snp.id=pruned.id)
df<- data.frame(sample.id =pca$sample.id, 
                  pop = factor(pop)[match(pca$sample.id, sample.id)],
                  PC1 = pca$eigenvect[,1],
                  PC2 = pca$eigenvect[,2],
                  PC3 = pca$eigenvect[,3],
                  PC4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

write.csv(df, file="hapmap_PCA_3.csv")

pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
getwd()