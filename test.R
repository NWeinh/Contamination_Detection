library(plyr)

source('cont.R')

# read in the samples sheet
samples <- read.csv2('sampleSheet.csv', sep=',', stringsAsFactors=F)

# read in the SNP panel
panel <- read.csv2("contPanel.csv", sep=",", stringsAsFactors = F)
panel <- panel[with(panel,order(contig,position)),]
panel$contig <- paste0('chr', panel$contig) # comment this out if your bam files don't use the chr notations for chromosomes
#panel <- subset(panel, Source=='BMCMED') # remove comment if you only want to use SNPs from this publication: 

# estimate contamination for list of samples in the sample sheet
res <- apply(samples, 1, function(x) {
  estCont(x['CONTROL'], x['TUMOR'], panel)
})
names(res) <- samples$SAMPLE_ID

# format the output for tabular display
df <- rbind(do.call(rbind, lapply(res, '[[', 1)),
            do.call(rbind, lapply(res, '[[', 2))
            )
df1 <- data.frame(Sample=rownames(df), df)
df1 <- arrange(df1, Sample)
df1