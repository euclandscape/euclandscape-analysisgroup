library(dplyr)

#names.df <- read.table("pops_sample_test.txt", header = FALSE) #File I was testing the code with
names.df.2 <- read.csv("./rawdata/JordanEmicroMetadataBySample.csv")
n.df <- names.df.2[,c(1,2)] #don't need to do this step, as the function calls which columns specifically
#write.table(names.df, "pops_sample_test.txt", quote =FALSE, col.names=FALSE, row.names = FALSE)
full.data.frame <-read.table("./rawdata/JordanEmicroSNP_012.tsv", header =FALSE)

#1 pop dataframe (Labeled "Population" and "Sample")
#2 snp dataframe (sample labeled "V1")
#3 pop.size <- 10
#4 numb <- "testing"  # for labeling of the final file
#5 maf.x <- 0.01
#6 prop <- 0.9 (proportion of intact data; i.e. 10% missing data)


choose.pop(n.df, full.data.frame, 15, "testing", 0.1, 0.6) #second is the label of the written file

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
  
  #final[1:10,1:2]
  #out[1:10,1:10]
  #full.data.frame[1:10,1:10]
  #final.merge.snp[1:10,1:10]
  
  #complete data (10% missing data)
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
  ####### write table ###########
  write.table(snps, paste(numb, size.x, prop, maf.x, ".txt", sep = "_"), quote = FALSE)
}

