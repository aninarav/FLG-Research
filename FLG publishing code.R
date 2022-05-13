#clearing global environment for a clean storage
rm(list=ls())

#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
#libraries called, which have to be installed first
#Everyone except ggplot2 requires biocmanager::install() to install
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(IRanges)
library(Biostrings)
library(ggplot2)
#store file name
fl <- "FLG_SNPs.VCF"

hdr <- info(scanVcfHeader(fl))
vcf <- readVcf(fl,"hg38")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#renames the vcf files' Seqlevels so that it fits the txdb file
vcf <- renameSeqlevels(vcf, c('NC_000001.11'='chr1'))
seqlevels(txdb) <- "chr1"
#coding locations
locCod <- locateVariants(vcf,txdb,CodingVariants())

#All the locations 
locAll <- locateVariants(vcf,txdb,AllVariants())


#predict coding based on amino acid changes
AAs <- predictCoding(vcf,txdb, Hsapiens)

#finds the intersection between AAs and locCod
#ranges is a Granges object with 26 metadata columns (from both AAs and locCod)
library(IRanges)
hits <- (IRanges::findOverlaps(AAs,locCod))
ranges <- subsetByOverlaps(locCod,AAs)
idx <- unique(subjectHits(hits))
values <- DataFrame(
  paramRangeID = AAs$paramRangeID[idx],
  REF = AAs$REF[idx],
  ALT = AAs$ALT[idx],
  QUAL = AAs$QUAL[idx],
  FILTER = AAs$FILTER[idx],
  varAllele = AAs$varAllele[idx],
  CDALOC = AAs$CDSLOC[idx],
  PROTEINLOC = AAs$PROTEINLOC[idx],
  QUERYID = AAs$QUERYID[idx],
  TXID = AAs$TXID[idx],
  CDSID = AAs$CDSID[idx],
  GENEID = AAs$GENEID[idx],
  CONSEQUENCE = AAs$CONSEQUENCE[idx],
  REFCODON = AAs$REFCODON[idx],
  VARCODON = AAs$VARCODON[idx],
  REFAA = AAs$REFAA[idx],
  VARAA = AAs$VARAA[idx])

mcols(ranges) <- c(mcols(ranges),values)


#subsetting the ranges object to the proper FLG range of Chr1: 152,302,165 - 152,325,239
e <- which(as.data.frame(ranges)$end > 152325239)[1]-1
ranges <- ranges[1:e]

#this step removes the synonymous variants and keeps the AA variants that makes the impact
#this step stores the new data in r2 with a granges class and 461 entries
#levels(ranges$CONSEQUENCE)[4] is synonymous

ids = which(ranges$CONSEQUENCE == levels(ranges$CONSEQUENCE)[4])
idn <- 1:length(ranges)
idn <- idn[-ids]
r2 <- ranges[idn]


#this step creates 3 new granges object :

#1 for frameshift (rnfs)
idfs = which(r2$CONSEQUENCE == levels(r2$CONSEQUENCE)[1])
rfs <- r2[idfs]
#1 for nonsense (rnse)
idnse = which(r2$CONSEQUENCE == levels(r2$CONSEQUENCE)[2])
rnse <- r2[idnse]

#1 for non-synonymous (rnsy)
idnsy = which(r2$CONSEQUENCE == levels(r2$CONSEQUENCE)[3])
rnsy <- r2[idnsy]

#-----------------------------------------------------------------------------------#
#non-synonymous section >>> 
#ordering the rnsy in terms of comparison between reference and variable Amino Acids


data("BLOSUM62")
refAA2 <- as.character(rnsy$REFAA)
varAA2 <- as.character(rnsy$VARAA)

#creating the comparison list

SubVec <- c()
for (n in 1:length(varAA2)){
  SubVec <- append(SubVec, BLOSUM62[refAA2[n],varAA2[n]])
}

SubDat <- DataFrame(COMP = SubVec)

# adding the comparison list to the granges object rnsy2
rnsy2 <- rnsy
mcols(rnsy2) <- c(mcols(rnsy), SubDat)


#limits the non-synonymous mutations to exon 3 because the other 2 domains do not have protein domains
nsye1 <- which(as.data.frame(rnsy2)$end > 152314747)[1]-1
rnsy2 <- rnsy2[1:nsye1]

#ordering rnsy2 based on most impactful
o <- order(rnsy2$COMP,start(rnsy2))
rnsy2 <- rnsy2[o]

#picking only non-positive comparison (BLOSUM62) scores because thy make the greatest impact
insyf <- which(rnsy2$COMP < 1)
rnsyf <- rnsy2[insyf]

#indicating calcium binding site on non-synonymous
o <- order(start(rnsyf))
rnsyf <- rnsyf[o]
cbs <- append(seq(0,0,length.out = (which(start(rnsyf) > 152314622)[1]-1)),
              seq(1,1,length.out = length(which(start(rnsyf) > 152314622))))
mcols(rnsyf) <- c(mcols(rnsyf), DataFrame(Ca_bindsite = cbs))
o <- order(rnsyf$COMP)
rnsyf <- rnsyf[o]

#----------------------------------------------------------------------#
# Nonsense mutation section >>

#extract reference nucleotide sequence of exon deleted by the rnse's first and last mutations.
#the rest don't matter because the first and last RNSE mutation 
#come before the last exon group and after the second to last exon group.
#the last exon group is Chr1 : 152325189 - 152325239
# All information comes from NCBI Variation Viewer
#currently the exon 3 is chosen to display this information
refNT <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38,"chr1",152302700,152314747)
refNTAA <- translate(refNT)
#input the output of 2 lines from before (refNT) into the NCBI conserved domain database search  

#all three regions are taken because exon 1 and 2 happen before the important exon 3
#exon 3
e <- 1:(which(as.data.frame(rnse)$end > 152314747)[1]-1)
#exon 1
append(e,intersect(which(as.data.frame(rnse)$end > 152325189), which(as.data.frame(rnse)$end < 152325239) ))
#exon 2
append(e,intersect(which(as.data.frame(rnse)$end > 152315319), which(as.data.frame(rnse)$end < 152315456) ))

rnse <- rnse[e]

#--------------------------------------------------------------------------#
#frameshift section >>

#all three regions are taken because exon 1 and 2 happen before the important exon 3
#exon 3
fsid <- 1:(which(as.data.frame(rfs)$end > 152314747)[1]-1)
#exon 1
append(fsid,intersect(which(as.data.frame(rfs)$end > 152325189), which(as.data.frame(rfs)$end < 152325239) ))
#exon 2
append(fsid,intersect(which(as.data.frame(rfs)$end > 152315319), which(as.data.frame(rfs)$end < 152315456) ))

rfs <- rfs[fsid]

#----------------------------------------------------------------#
#Graphs

#figure 2
dat1 <- data.frame(Location = c(start(locAll),
                                start(locCod)),
                   Legend = append(rep("All SNPs",
                                       each = length(locAll)),
                                   rep("Coding SNPs",
                                       each = length(locCod))))

library(ggplot2)
figure2 <- ggplot(dat1, aes(x = Location, fill = Legend)) +
  geom_rect(aes(xmin = 152302700, xmax = 152314747, ymin = 0, ymax = 3000, fill = "Region : Exon 3"),alpha = 0.01) +
  geom_rect(aes(xmin = 152315319, xmax = 152315456, ymin = 0, ymax = 3000, fill = "Region : Exon 2"),alpha = 0.01) +
  geom_rect(aes(xmin = 152314622, xmax = 152314747, ymin = 0, ymax = 3000, fill = "Region : Ca (2+) Binding"),alpha = 0.75) +
  geom_rect(aes(xmin = 152325189, xmax = 152325239, ymin = 0, ymax = 3000, fill = "Region : Exon 1"),alpha = 0.01) +
  geom_histogram(bins = 100, alpha = 0.9) +
  scale_fill_manual(values = c("orange","green","red","blue4", "dodgerblue1","cadetblue2")) +
  geom_vline(xintercept = c(152302700,152314747), color = "cadetblue2",alpha = 0.25) +
  geom_vline(xintercept = c(152315319,152315456), color = "dodgerblue1", alpha = 0.25) + 
  xlab("Genomic Location (CHR1 : Nucleotides)") +
  ylab('Counts of SNPs')

#figure 3
barplot(table(rnsy2$COMP),
        xlab = "BLOSUM62 scores",
        ylab = "Count of SNPs",
        col = "pink")
abline(v = 6, col = "red",lwd= 3)

#figure 4
dat4 <- data.frame(Location = c(start(locCod),start(rnse),
                                start(rfs),start(rnsyf)),
                   Legend = c(rep("SNPs : Coding", each = length(locCod)),
                              rep("SNPs : Nonsense", each = length(rnse)),
                              rep("SNPs : Frameshift", each = length(rfs)),
                              rep("SNPs : Nonsynonymous", each = length(rnsyf))))

figure4 <- ggplot(dat4, aes(x = Location, fill = Legend)) +
  geom_rect(aes(xmin = 152302700, xmax = 152314747, ymin = 0, ymax = 200, fill = "Region spanning Exon 3"),alpha = 0.01) +
  geom_rect(aes(xmin = 152315319, xmax = 152315456, ymin = 0, ymax = 200, fill = "Region spanning Exon 2"),alpha = 0.01) +
  geom_rect(aes(xmin = 152314622, xmax = 152314747, ymin = 0, ymax = 200, fill = "Region : Ca (2+) Binding"),alpha = 0.75) +
  geom_rect(aes(xmin = 152325189, xmax = 152325239, ymin = 0, ymax = 200, fill = "Region spanning Exon 1"),alpha = 0.01) +
  geom_histogram(data = dat4[1:length(locCod),1:2], bins = 100 , alpha = 0.5) +
  geom_histogram(data = dat4[-(1:length(locCod)),1:2], position = "stack", bins = 100, alpha = 0.9) +
  scale_fill_manual(values = c("red","blue4", "dodgerblue1","cadetblue1",
                               "green", "purple", "yellow", "pink","black")) +
  geom_vline(xintercept = c(152302700,152314747), color = "cadetblue2",alpha = 0.1) +
  geom_vline(xintercept = c(152315319,152315456), color = "dodgerblue1", alpha = 0.25) +
  xlab("Genomic Location (CHR1 : Nucleotides)") +
  ylab("Counts of SNPs")



