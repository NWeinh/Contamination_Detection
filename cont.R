library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)

.pileupFreq <- function(pileupres) {
  nucleotides <- levels(pileupres$nucleotide)
  res <- split(pileupres, pileupres$seqnames)
  res <- lapply(res, function (x) {split(x, x$pos)})
  res <- lapply(res, function (positionsplit) {
    nuctab <- lapply(positionsplit, function(each) {
      chr = as.character(unique(each$seqnames))
      pos = as.character(unique(each$pos))
      tabcounts <- sapply(nucleotides, function (n) {sum(each$count[each$nucleotide == n])})
      c(chr,pos, tabcounts)
    })
    nuctab <- data.frame(do.call("rbind", nuctab),stringsAsFactors=F)
    rownames(nuctab) <- NULL
    nuctab
  })
  res <- data.frame(do.call("rbind", res),stringsAsFactors=F)
  rownames(res) <- NULL
  colnames(res) <- c("seqnames","start",levels(pileupres$nucleotide))
  res[3:ncol(res)] <- apply(res[3:ncol(res)], 2, as.numeric)
  res
}

.readBam <- function(bam, vcf.ranges, min_base_quality) {
  bamfile <- as.character(bam)
  bf <- BamFile(bamfile)
  param <- ScanBamParam(which = vcf.ranges)
  flag <- scanBamFlag(isDuplicate = FALSE)
  p_param <- PileupParam(max_depth = 10000,
                         min_base_quality = min_base_quality,
                         ignore_query_Ns = FALSE)
  res <- pileup(bf, 
                scanBamParam=param, 
                scanBamFlag=flag,
                pileupParam=p_param)
  res
}

.qcTab <- function(type, ranges, tab, contPerSNP, maxContLevelGerm, minReads, percentHomo, aberrantSNP) {
  tab$REF_A <- as.vector(ranges$REF)
  tab$ALT_A <- as.vector(ranges$ALT)
  tab$REF <- NA
  tab$ALT <- NA
  tab$REF <- ifelse(ranges$REF == "T", tab$T, tab$REF)
  tab$REF <- ifelse(ranges$REF == "A", tab$A, tab$REF)
  tab$REF <- ifelse(ranges$REF == "C", tab$C, tab$REF)
  tab$REF <- ifelse(ranges$REF == "G", tab$G, tab$REF)
  tab$ALT <- ifelse(ranges$ALT == "T", tab$T, tab$ALT)
  tab$ALT <- ifelse(ranges$ALT == "A", tab$A, tab$ALT)
  tab$ALT <- ifelse(ranges$ALT == "C", tab$C, tab$ALT)
  tab$ALT <- ifelse(ranges$ALT == "G", tab$G, tab$ALT)
  
  tab$A1 <- ifelse(tab$REF_A =="A" & tab$ALT_A == "C", tab[,5], 0)
  tab$A1 <- ifelse(tab$REF_A =="A" & tab$ALT_A == "G", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="A" & tab$ALT_A == "T", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="C" & tab$ALT_A == "A", tab[,5], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="C" & tab$ALT_A == "G", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="C" & tab$ALT_A == "T", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="G" & tab$ALT_A == "A", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="G" & tab$ALT_A == "C", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="G" & tab$ALT_A == "T", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="T" & tab$ALT_A == "A", tab[,4], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="T" & tab$ALT_A == "C", tab[,3], tab$A1)
  tab$A1 <- ifelse(tab$REF_A =="T" & tab$ALT_A == "G", tab[,3], tab$A1)
  
  tab$A2 <- ifelse(tab$REF_A == "A" & tab$ALT_A == "C", tab[,6], 0)
  tab$A2 <- ifelse(tab$REF_A == "A" & tab$ALT_A == "G", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "A" & tab$ALT_A == "T", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "C" & tab$ALT_A == "A", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "C" & tab$ALT_A == "G", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "C" & tab$ALT_A == "T", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "G" & tab$ALT_A == "A", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "G" & tab$ALT_A == "C", tab[,6], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "G" & tab$ALT_A == "T", tab[,4], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "T" & tab$ALT_A == "A", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "T" & tab$ALT_A == "C", tab[,5], tab$A2)
  tab$A2 <- ifelse(tab$REF_A == "T" & tab$ALT_A == "G", tab[,4], tab$A2)
  tab <- tab[,c(1:2,11:16)]
  
  tab$P_REF <- round(100/(as.numeric(tab$REF) + as.numeric(tab$ALT)) * as.numeric(tab$REF), 0)
  
  tab$minReads <- ifelse(as.numeric(tab$REF) + as.numeric(tab$ALT) < minReads, "no", "yes")
 
  if (type=='germline') {
    tab$Homo <- ifelse(tab$P_REF == 100 | tab$P_REF <= percentHomo, "yes", "no")
    tab$Conta <- as.numeric(ifelse(tab$P_REF >= (100 - maxContLevelGerm), tab$ALT, tab$REF))
    tab$Geno_A <- ifelse(tab$P_REF >= (100 - maxContLevelGerm), tab$REF_A, tab$ALT_A)
    tab$Conta_A <- ifelse(tab$P_REF >= (100 - maxContLevelGerm), tab$ALT_A, tab$REF_A)  
    
    tab <- subset(tab, Homo == 'yes' & minReads == 'yes')
    
    tab_QC <- tab
    tab_QC$Norm <- as.numeric(ifelse(tab_QC$P_REF >= 90, tab_QC$REF, tab_QC$ALT))
    tab_QC$ProContaAlt <- 100 / (tab_QC$Norm+tab_QC$Conta + tab_QC$A1 + tab_QC$A2) * (tab_QC$Conta)
    tab_QC$ProContaA1 <- 100 / (tab_QC$Norm+tab_QC$Conta + tab_QC$A1 + tab_QC$A2) * (tab_QC$A1)
    tab_QC$ProContaA2 <- 100 / (tab_QC$Norm+tab_QC$Conta + tab_QC$A1 + tab_QC$A2) * (tab_QC$A2)
    
    QC <- list(N_HOMO = dim(tab_QC)[1],
               NumSNP_Alt = length(which(as.numeric(tab_QC[ ,12]) > contPerSNP)),
               NumSNP_A1 = length(which(as.numeric(tab_QC[ ,7]) > contPerSNP)),
               NumSNP_A2 = length(which(as.numeric(tab_QC[ ,8]) > contPerSNP)),
               Sum_ALT = sum(as.numeric(tab_QC[ ,12])),
               Sum_A1 = sum(as.numeric(tab_QC[ ,7])),
               Sum_A2 = sum(as.numeric(tab_QC[ ,8])),
               Percent_Alt = 100 / dim(tab_QC)[1] * length(which(as.numeric(tab_QC[ ,12]) > contPerSNP)),
               Percent_A1 = 100 / dim(tab_QC)[1] * length(which(as.numeric(tab_QC[ ,7]) > contPerSNP)),
               Percent_A2 = 100 / dim(tab_QC)[1] * length(which(as.numeric(tab_QC[ ,8]) > contPerSNP)),
               PercentCont_Alt = sum(tab_QC$ProContaAlt[which(as.numeric(tab_QC[ ,12]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A1 = sum(tab_QC$ProContaA1[which(as.numeric(tab_QC[ ,7]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A2 = sum(tab_QC$ProContaA2[which(as.numeric(tab_QC[ ,8]) > contPerSNP)]) / dim(tab)[1],
               Type = 'Germline'
    )
    QC$Background <- ifelse(QC$PercentCont_A1 > QC$PercentCont_A2, QC$PercentCont_A1, QC$PercentCont_A2)
    QC$ContAdj <- ifelse(QC$N_HOMO / 100 * QC$Percent_Alt > aberrantSNP, 2 * (QC$PercentCont_Alt - QC$Background), 0)
    QC$ContAdj <- round(ifelse(QC$ContAdj < 0, 0, QC$ContAdj), 2)
  }
  
  if (type=='tumor') {
    tab$Homo <- ifelse(tab$P_REF >= 90 | tab$P_REF <= 10 , "yes", "no")
    tab$Conta <- tab$ALT
    tab$Norm=tab$REF
    tab$ProContaAlt <- 100 / (tab$Norm + tab$Conta + tab$A1 + tab$A2) * (tab$Conta)
    tab$ProContaA1 <-  100 / (tab$Norm + tab$Conta + tab$A1 + tab$A2) * (tab$A1)
    tab$ProContaA2 <- 100 / (tab$Norm + tab$Conta + tab$A1 + tab$A2) * (tab$A2)
    
    QC <- list(N_HOMO = dim(tab)[1],
               NumSNP_Alt = length(which(as.numeric(tab[ ,12]) > contPerSNP)),
               NumSNP_A1 = length(which(as.numeric(tab[ ,7]) > contPerSNP)),
               NumSNP_A2 = length(which(as.numeric(tab[ ,8]) > contPerSNP)),
               Sum_Alt = sum(as.numeric(tab[ ,12])),
               Sum_A1 = sum(as.numeric(tab[ ,7])),
               Sum_A2 = sum(as.numeric(tab[ ,8])),
               Percent_Alt = 100 / dim(tab)[1] * length(which(as.numeric(tab[ ,12]) > contPerSNP)),
               Percent_A1 = 100 / dim(tab)[1] * length(which(as.numeric(tab[ ,7]) > contPerSNP)),
               Percent_A2 = 100 / dim(tab)[1] * length(which(as.numeric(tab[ ,8]) > contPerSNP)),
               PercentCont_Alt = sum(tab$ProContaAlt[which(as.numeric(tab[ ,12]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A1 = sum(tab$ProContaA1[which(as.numeric(tab[ ,7]) > contPerSNP)]) / dim(tab)[1],
               PercentCont_A2 = sum(tab$ProContaA2[which(as.numeric(tab[ ,8]) > contPerSNP)]) / dim(tab)[1],
               Type = 'Tumor'
    )
    QC$Background <- ifelse(QC$PercentCont_A1 > QC$PercentCont_A2, QC$PercentCont_A1, QC$PercentCont_A2)
    QC$ContAdj <- ifelse(QC$N_HOMO / 100 * QC$Percent_Alt > aberrantSNP, 2 * (QC$PercentCont_Alt - QC$Background), 0)
    QC$ContAdj <- round(ifelse(QC$ContAdj < 0, 0, QC$ContAdj), 2)
  }
  return(list(QC=QC, tab=tab))
}

estCont <- function(bamGermline, bamTumor, panel, mode='pair', percentHomo=10, minReads=50, maxContLevelGerm=10,
                    min_base_quality=20, contPerSNP=0, aberrantSNP=3) {
  if(mode == 'pair') {
    ## GERMLINE
    vcf.ranges=GRanges(seqnames = panel$contig, 
                       ranges = IRanges(start=panel$position, end=panel$position), 
                       strand = "*", 
                       paramRangeID = NA, 
                       REF = panel$ref_allele, 
                       ALT = panel$alt_allele, 
                       QUAL = NA, 
                       FILTER = "REJECT", 
                       SNP = panel$SNP)
    res <- .readBam(bamGermline, vcf.ranges, min_base_quality)
    tab <- .pileupFreq(res)
    
    QCGerm <- .qcTab('germline', vcf.ranges, tab, contPerSNP, maxContLevelGerm, minReads, percentHomo, aberrantSNP)
    
    ## TUMOR
    vcf.ranges2 <- GRanges(seqnames = QCGerm$tab$seqnames, 
                           ranges = IRanges(start=as.numeric(QCGerm$tab$start), end=as.numeric(QCGerm$tab$start)), 
                           strand = "*", 
                           paramRangeID = NA, 
                           REF = QCGerm$tab$Geno_A, 
                           ALT = QCGerm$tab$Conta_A, 
                           QUAL = NA, 
                           FILTER = "REJECT"
                           )
    resTumor <- .readBam(bamTumor, vcf.ranges2, min_base_quality)
    tabTumor <- .pileupFreq(resTumor)
    
    QCTumor <- .qcTab('tumor', vcf.ranges2, tabTumor, contPerSNP, maxContLevelGerm, minReads, percentHomo, aberrantSNP)
  } 
  
  if(mode == 'single') {
    print('No support for this, yet!')
  }
  return(list(QCGerm=QCGerm$QC,
              QCTumor=QCTumor$QC,
              TabGerm=QCGerm$tab,
              TabTumor=QCTumor$tab
              )
         )
}

