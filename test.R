t <- '/data/CCD201-T1-DNA-rep1_S1.bam'
b <- '/data/01-201-B1-DNA-rep1_S1.bam'

source('cont.R')

panel <- read.csv2("contPanel.csv", sep=",", stringsAsFactors = F)
panel <- panel[with(panel,order(contig,position)),]
panel$contig <- paste0('chr', panel$contig)
#panel <- subset(panel, Source=='BMCMED')

res <- estCont(b, t, panel, mode='pair', percentHomo=10, minReads=10, maxContLevelGerm=5)
rbind(res$QCGerm, res$QCTumor)
