---
title: "Change Pchic Annotation"
author: Alexis Hucteau
date: "`r Sys.Date()`"
output:
    html_document:
      toc: yes
      toc_float: yes
      theme: united
editor_options:
  markdown:
    wrap: sentence
---

```{r}
rm(list = ls())


"%ni%" <- Negate("%in%")
```

```{r}
prepare_pchic <- function(cell_lines = "all", minimum_interaction = 5){
  load("/media/alexis/DATA/pchic.RData")
  if (length(cell_lines) >= 1){
    cell_lines = c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
  }
  pchic <- data.frame(pchic[rowSums(pchic[,cell_lines] >= minimum_interaction) >= 1, 1:10]) %>% na.omit(.)
  colnames(pchic)[c(1:5, 6:10)] <- rep(c("chr", "start", "end", "ID", "Name"), 2)
  return(pchic)
}

Pchic_all_network <- prepare_pchic(c("Mon", "Mac0", "Mac1", "Mac2", "MK", "EP", "Ery"))
pchic_bed <- unique(rbind(Pchic_all_network[, c(1:4, 5)], Pchic_all_network[, c(6:9, 10)]))
```

```{r}
library("GenomicRanges")
library("AnnotationHub")
library("rtracklayer")
library(data.table)
library("BiocGenerics")
library(parallel)

# BiocManager::install("Repitools")
BiocManager::install("Repitools", lib= "/media/alexis/DATA/R/Rpackage/", force = TRUE)
```

```{r}
Ahub <- AnnotationHub()
Ahub <- subset(Ahub, species == "Homo sapiens")
unique(Ahub$dataprovider)

promoter_data <- query(Ahub, c("GRanges", "hg38"))
promoter_data$title
promoter_data

data.frame("id" = promoter_data@.db_uid, "desc" = promoter_data$title)

test <- promoter_data[["AH75192"]]
test2 <- promoters(test)
dt2 <- copy(as.data.table(test2))
# test3 <- Repitools::annoGR2DF(test2)
dt2_GR
```

```{r}
prepare_pchic <- function(cell_lines = NULL, minimum_interaction = 5, pchic = NULL){
  if (is.null(pchic)){
    load("/media/alexis/DATA/pchic.RData")
  }
  if (is.null(cell_lines)){
    cell_lines = c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
  }
  pchic <- data.frame(pchic[rowSums(pchic[,cell_lines] >= minimum_interaction) >= 1,c(1:3,5:8,10)]) %>% na.omit(.)
  pchic$baitChr <- paste0("chr", pchic$baitChr)
  pchic$oeChr <- paste0("chr", pchic$oeChr)
  return(pchic)
}
```

### Myeloid pchic

```{r}
Myelo_cell_lines <- c("Mon", "Mac1", "Mac0", "Mac2", "MK", "Ery", "EP")
HSC_pchic <- read.table("../../Datasets/HSC_Pchic_cleaned.tsv", sep = "\t", header = T)

pchic <- prepare_pchic(cell_lines = Myelo_cell_lines)
pchic_bed <- pchic
colnames(pchic_bed) <- rep(c("chr", "start", "end", "Name"), 2)
pchic_bed <- rbind(pchic_bed[,1:4], pchic_bed[,5:8]) %>%
  unique
rownames(pchic_bed) <- paste0(pchic_bed$chr, ":", pchic_bed$start, "-", pchic_bed$end)
```

```{r}
pchic_GRanges <- GRanges(seqnames = pchic_bed$chr, 
                         ranges = IRanges(start = pchic_bed$start, end = pchic_bed$end),
                         Pchic_Chrom = pchic_bed$chr, Pchic_Start = pchic_bed$start, Pchic_End = pchic_bed$end,
                         Pchic_GeneName = pchic_bed$Name)

dt2_GRanges <- GRanges(seqnames = dt2$seqnames, 
                       ranges = IRanges(start = dt2$start, end = dt2$end),
                       Chrom_dt2 = dt2$seqnames, Start_dt2 = dt2$start, End_dt2 = dt2$end, Gene_name = dt2$Gene)

Over <- findOverlaps(pchic_GRanges, dt2_GRanges)
New_anno <- data.frame(mcols(pchic_GRanges[queryHits(Over),]),
                       data.frame(mcols(dt2_GRanges[subjectHits(Over),]))) %>%
  .[c(1:3, 8)] %>% unique

Pchic_not_in_new_anno <- pchic_bed %>% dplyr::filter

New_anno$Rownames <- paste0(New_anno$Pchic_Chrom, ":", New_anno$Pchic_Start, "-", New_anno$Pchic_End)

New_anno_cleaned <- New_anno$Rownames %>% unique %>% mclapply(function(frag){
  tmp <- dplyr::filter(New_anno, Rownames == frag)
  Gene_names <- tmp$Gene_name %>% na.omit %>% unique %>% paste(collapse = ";")
  data.table("Rownames" = frag, "chr" = tmp$Pchic_Chrom, "start" = tmp$Pchic_Start[1], 
                    "end" = tmp$Pchic_End[1], "Name" = Gene_names) %>% unique

}, mc.cores = 14) %>% data.table::rbindlist()
rownames(New_anno_cleaned) <- New_anno_cleaned$Rownames

Pchic_not_in_new_anno <- pchic_bed[rownames(pchic_bed) %ni% New_anno_cleaned$Rownames,]

Gene_not_new_anno <- Pchic_not_in_new_anno$Name %>% sapply(stringr::str_split, ";") %>% unlist %>% unique
Gene_new_anno <- New_anno_cleaned$Name %>% sapply(stringr::str_split, ";") %>% unlist %>% unique

intersect(Gene_not_new_anno, Gene_new_anno)

New_anno_cleaned %>% write.table("~/GitHub/Thesis_paper/Datasets/Chromatine_part/New_annotation_fragment_to_gene_promoter.tsv", sep = "\t")
```
