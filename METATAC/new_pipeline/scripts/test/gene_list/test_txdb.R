library(GenomicRanges)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ArchR)
library(GenomicFeatures)

setwd("/data/Lab/otherwork/GeACT/ATAC/database")

# test official TxDb -----------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# gencode gene with entrez/ncbi ids
have_id <- unlist(transcriptsBy(txdb, by = "gene"))

transcripts(txdb) %>% length
have_id %>% length

# ensembl/gencode gene names with no entrez/ncbi ids
setdiff(transcripts(txdb)$tx_name, sort(have_id)$tx_name) %>% length
no_id <- transcripts(txdb)[transcripts(txdb)$tx_name %in% setdiff(transcripts(txdb)$tx_name, sort(have_id)$tx_name)]

# test direct generation of txdb ----------------------------
# not feasible due to difficulty in preprocessing of metadata
# my_txdb <- makeTxDbFromGFF("GENCODE/v26/gencode.v26.primary_assembly.annotation.gff3.gz", format = "gff3", taxonomyId = 9606, dataSource = "GENCODEv26",
#                            organism = "Homo sapiens", chrominfo = seqinfo(BSgenome.Hsapiens.UCSC.hg38))
# saveDb(my_txdb, 'txdb/txdb.gencodev26.sqlite')
# my_txdb <- loadDb(file = 'txdb/txdb.gencodev26.sqlite')



# test indirect generation of txdb from grange(by gff3) -------------
# adapted from LieberInstitute/GenomicState/R/gencode_txdb.R
chrs <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
gff3_file <- "GENCODE/v26/gencode.v26.primary_assembly.annotation.gff3.gz"
gencode_gff3 <- rtracklayer::import(con = gff3_file)

## Keep only the main chrs
gencode_gff3 <- GenomeInfoDb::keepSeqlevels(gencode_gff3, chrs,
                                           pruning.mode = "coarse")

metadata <- GenomicFeatures:::.prepareGFFMetadata(
  file = gff3_file,
  dataSource = NA, organism = "Homo sapiens",
  taxonomyId = NA, miRBaseBuild = NA, metadata = NULL
)

my_seqinfo <- GenomeInfoDb::Seqinfo(genome = "hg38") %>% GenomeInfoDb::keepSeqlevels(chrs)
gr <- GenomicFeatures:::.tidy_seqinfo(
  gr = gencode_gff3,
  circ_seqs = DEFAULT_CIRC_SEQS,
  chrominfo = my_seqinfo
)

## Prune again since GenomeInfoDb::Seqinfo() will return many seqlevels # changed
gr <- GenomeInfoDb::keepSeqlevels(gr, chrs, pruning.mode = "coarse")

# trim version number
# gr$ID <- gsub(pattern = "\\.[0-9]*", replacement = "", gr$ID)
# gr$Parent <- gsub(pattern = "\\.[0-9]*", replacement = "", gr$Parent)


my_txdb <- GenomicFeatures::makeTxDbFromGRanges(gr, metadata = metadata)

saveDb(my_txdb, 'txdb/txdb.gencodev26.sqlite')



# test ArchR with self generated txdb -----------------
# my_txdb <- loadDb(file = 'txdb/txdb.gencodev26.sqlite')
# geneAnnotation  <- createGeneAnnotation(TxDb = my_txdb, OrgDb = org.Hs.eg.db, annoStyle = "ENSEMBL")

# import blacklists
blacklist1 <- rtracklayer::import("blacklists/hg38-blacklist.v2.bed.gz")
blacklist2 <- rtracklayer::import("blacklists/mito_hg38_peaks.narrowPeak")

blacklist <- c(blacklist1, blacklist2)
blacklist <- GenomicFeatures:::.tidy_seqinfo(
  gr = blacklist,
  circ_seqs = GenomicFeatures::DEFAULT_CIRC_SEQS,
  chrominfo = GenomeInfoDb::Seqinfo(genome = "hg38")
)
blacklist <- sort(blacklist)

my_genome <- BSgenome.Hsapiens.UCSC.hg38
seqlevels(my_genome) <- chrs


chromsize <- GRanges(names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)), IRanges(1, seqlengths(BSgenome.Hsapiens.UCSC.hg38)))
chromsize <- GenomeInfoDb::keepSeqlevels(chromsize, chrs, pruning.mode = "coarse")
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Hsapiens.UCSC.hg38, blacklist = blacklist, 
                                           chromSizes = chromsize)


# test generate geneAnnotation(gene, exon, TSS) manually ----------------------
my_txdb <- loadDb("txdb/txdb.gencodev26.sqlite")

genes <- genes(my_txdb)
gene_symbol <- read_tsv("GENCODE/v26/gene_ID2Name_fixed.txt", col_names = F) %>% column_to_rownames(var = "X1")

mcols(genes)$symbol <- gene_symbol[mcols(genes)$gene_id, ]
names(genes) <- NULL
genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)

###########################
message("Getting Exons..")
exons <- unlist(GenomicFeatures::exonsBy(my_txdb, by = "tx"))
exons$tx_id <- names(exons)
mcols(exons)$gene_id <- suppressMessages(AnnotationDbi::select(my_txdb, keys = paste0(mcols(exons)$tx_id), column = "GENEID", keytype = "TXID")[, "GENEID"])
exons <- exons[!is.na(mcols(exons)$gene_id), ]
mcols(exons)$symbol <- gene_symbol[mcols(exons)$gene_id, ]
mcols(exons)$symbol[is.na(mcols(exons)$symbol)] <- paste0("NA_", mcols(exons)$gene_id)[is.na(mcols(exons)$symbol)]
names(exons) <- NULL
mcols(exons)$exon_id <- NULL
mcols(exons)$exon_name <- NULL
mcols(exons)$exon_rank <- NULL
mcols(exons)$tx_id <- NULL
exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)

###########################
message("Getting TSS..")
TSS <- unique(resize(GenomicFeatures::transcripts(my_txdb), width = 1, fix = "start"))

geneAnnotation <- SimpleList(genes = genes, exons = exons, TSS = TSS)

# save gene and genome Annotation
saveRDS(object = genomeAnnotation, file = "annotation/genomeAnnotation.rds")
saveRDS(object = geneAnnotation, file = "annotation/geneAnnotation.rds")
