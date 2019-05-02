library(dplyr)

#names.df <- read.table("pops_sample_test.txt", header = FALSE) #File I was testing the code with
names.df.2 <- read.csv("./rawdata/JordanEmicroMetadataBySample.csv")
n.df <- names.df.2[,c(1,2)] #don't need to do this step, as the function calls which columns specifically
#write.table(names.df, "pops_sample_test.txt", quote =FALSE, col.names=FALSE, row.names = FALSE)
full.data.frame <-read.table("./rawdata/JordanEmicroSNP_012.tsv", header =FALSE)

#1 pop dataframe (Labeled "Population" and "Sample")
#2 snp dataframe (sample labeled "V1")
#3 pop.size <- 10
#4 numb <- "testing"  # for labeling of the final file (probably species)
#5 maf.x <- 0.01
#6 prop <- 0.9 (proportion of intact data; i.e. 10% missing data)

#identify the parameters for all pops
minor.allele <- c(0.1) #, 0.05, 0.01)
prop.md <- c(0.5, 0.6) #, 0.7, 0.8, 0.9)
ind.pop <- c(25) #, 20, 15, 10, 6)

for (i in minor.allele) {
  for(j in prop.md) {
    for(k in ind.pop) {
choose.pop(n.df, full.data.frame, k, "species", i, j) 
#third is the label of the written file, and probably the species
    }
  }
}
#####################################################
#   function
#####################################################

choose.pop <- function(pops, snps, pop.size, numb, maf.x, prop){
  size.x <- pop.size
  count.pops<-pops %>% group_by(Population) %>% count(Population)
  keep<-count.pops %>% filter(n > size.x)
  indata <- subset(pops, Population %in% keep$Population)
  sub.10 <- indata %>%
    group_by(Population) %>%
    sample_n(size = size.x)
  outdata <- subset(pops, !Population %in% keep$Population)
  final <- rbind(as.data.frame(sub.10), as.data.frame(outdata))
  final.merge.snp <- merge(final, snps, by.x= "Sample", by.y = "V1")
  final.merge.snp <- final.merge.snp[,-2]
  #write.csv(final.merge.snp, paste0(numb, pop.size, ".csv"))
  out <- final.merge.snp[,-1]
  rownames(out) <- final.merge.snp[,1]

  #clean snps
  geno <- t(out)
  # look at counts across samples
  samplecalls<-apply(geno,1, function(x) sum(!is.na(x)));
  
  samp.thres <- prop*ncol(geno)
  
  g.snp <- geno[samplecalls > samp.thres,]
  ###### minor allele frequency ###########
  maf<-rowSums(g.snp,na.rm=T)/((ncol(g.snp)-rowSums(is.na(g.snp)))*2)
  #table(maf>mafThresh)
  maf<-g.snp[maf>maf.x,]
  #major allele freq
  jaf<-rowSums(maf,na.rm=T)/((ncol(g.snp)-rowSums(is.na(maf)))*2)
  jafThresh<-1-maf.x
  #table(jaf<jafThresh)
  snps<-maf[jaf<jafThresh,]
  ind <- ncol(snps)
  ####### write table ###########
  write.table(snps, paste0(numb,"_",ind,"_",k,"_",j,"_",i,".txt"), quote = FALSE)
}






library("reshape2")
m <-read.table("species_559_25_0.6_0.1_.txt", header =TRUE)
m.short <- t(m[1:500, 1:559])
test <- cor(m.short, use = "na.or.complete" )
test.melt <- melt(test)
head(test.melt)
test.sub <- subset(test.melt, value > 0)

test.cfk <- cor.fk(t(m.short))
