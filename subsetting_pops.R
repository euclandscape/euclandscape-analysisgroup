library(dplyr)
names.df <- read.table("pops_sample_test.txt", header = FALSE) # a few of the pops need to be fixed
full.data.frame <-read.table("./rawdata/JordanEmicroSNP_012.tsv", header =FALSE) 

choose.pop(20, "testing", 0.95, 0.05, 0.9) #second is the label of the writtin file
# (number of individuals/pop, first portion for label of file, alternate allele freq, minor allele freq, proportion of present data)
#1 pop.size <- 20
#2 numb <- "testing"
#3 jaf.x <- 0.99
#4 maf.x <- 0.01
#5 prop <- 0.9


###################################################################
#         function for choosing pops
###################################################################
choose.pop <- function(pop.size, numb,jaf.x, maf.x, prop){
size.x <- pop.size
count.pops<-names.df %>% group_by(V2) %>% count(V2)
keep<-count.pops %>% filter(n > size.x)
indata <- subset(names.df, V2 %in% keep$V2)
sub.10 <- indata %>%
          group_by(V2) %>%
          sample_n(size = size.x)
outdata <- subset(names.df, !V2 %in% keep$V2)
final <- rbind(as.data.frame(sub.10), as.data.frame(outdata))
final.merge.snp <- merge(final, full.data.frame, by.x= "V1", by.y = "V1")
final.merge.snp <- final.merge.snp[,-2]
#write.csv(final.merge.snp, paste0(numb, pop.size, ".csv"))
out <- final.merge.snp[,-1]
rownames(out) <- final.merge.snp[,1]

#final[1:10,1:2]
#out[1:10,1:10]
#full.data.frame[1:10,1:10]
#final.merge.snp[1:10,1:10]

#clean snps
geno <- t(out)
# look at counts across samples
samplecalls<-apply(geno,1, function(x) sum(!is.na(x)));

samp.thres <- prop*ncol(geno)

g.snp <- geno[samplecalls > samp.thres,]
###### minor allele frequency ###########
maf<-rowSums(g.snp,na.rm=T)/((ncol(g.snp)-rowSums(is.na(g.snp)))*2)
mafThresh<-maf.x
#table(maf>mafThresh)
maf<-g.snp[maf>mafThresh,]
#major allele freq
jaf<-rowSums(maf,na.rm=T)/((ncol(g.snp)-rowSums(is.na(maf)))*2)
jafThresh<-jaf.x
#table(jaf<jafThresh)
snps<-maf[jaf<jafThresh,]
####### write table ###########
write.table(snps, paste(numb, size.x, prop, maf.x, ".txt", sep = "_"), quote = FALSE)
}
