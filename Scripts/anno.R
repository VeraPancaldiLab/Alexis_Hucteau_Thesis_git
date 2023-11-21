library(Homo.sapiens)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
Promoters_genes <- promoters(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))
AnnoBMIQ <- read.table("GitHub/Thesis_paper/Datasets/Annotation_EPIC_BMIQ.tsv", sep = "\t", header = T)

Anno_Granges <- GRanges(AnnoBMIQ[1:3])
head(Anno_Granges)
Over <- overlapsRanges(Anno_Granges, Promoters_genes@ranges)
New_anno <- data.frame(mcols(Anno_Granges@ranges[queryHits(Over),]),
                       data.frame(mcols(DMP_GRanges[subjectHits(overlaps),])))

head(Promoters_genes@ranges)
