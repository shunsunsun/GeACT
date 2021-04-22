library(cicero)

Args <- commandArgs()
setwd(Args[6])
dir.create("cicero",recursive = T)
cicero_data <- read.table(gzfile(Args[7]),sep='\t',header=F,stringsAsFactors = F)

# Loading data from a simple sparse matrix format
input_cds <- make_atac_cds(cicero_data, binarize = TRUE)

# Constructing cis-regulatory networks
## Running Cicero
### Create a Cicero CDS

#### dimensionality reduction by UMAP
set.seed(2020)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', preprocess_method = "LSI")
# plot_cells(input_cds)

#### access the UMAP coordinates
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

### Run Cicero

#### run_cicero: get Cicero outputs with all defaults
genome_size = read.table("~/references/human/hg38.chrom22X.sizes",
                         sep='\t',header = F)
conns <- run_cicero(cicero_cds, genome_size, sample_num = 2)
write.table(conns,gzfile("cicero/conns.txt.gz"),sep = '\t',row.names = F,col.names = F,quote = F)

## Finding cis-Co-accessibility Networks (CCANS)
CCAN_assigns <- generate_ccans(conns)
write.table(CCAN_assigns,gzfile("cicero/CCAN_assigns.txt.gz"),sep='\t',row.names = F,col.names = F,quote = F)

## Cicero gene activity scores
#### Add a column for the pData table indicating the gene if a peak is a promoter ####
# Create a gene annotation set that only marks the transcription start sites of 
# the genes. We use this as a proxy for promoters.

gene_annotation_sub <- read.table("~/references/human/gene_position.hg38.v26.bed",
                                  sep='\t',header = F,stringsAsFactors = F)
chr_ids = paste0("chr",c(1:22,"X"))
gene_annotation_sub <- subset(gene_annotation_sub,is.element(gene_annotation_sub$V1,chr_ids))
gene_annotation_sub[gene_annotation_sub$V6=="+",2] = gene_annotation_sub[gene_annotation_sub$V6=="+",2]-2000
gene_annotation_sub[gene_annotation_sub$V6=="-",3] = gene_annotation_sub[gene_annotation_sub$V6=="-",3]+2000
gene_annotation_sub <- gene_annotation_sub[,1:4]
names(gene_annotation_sub) <- c("chromosome", "start", "end", "gene")
dup_gene = read.table("~/references/human/gene_dupName.txt",header = F,stringsAsFactors = F)[,1]
gene_annotation_sub <- subset(gene_annotation_sub,!is.element(gene_annotation_sub$gene,dup_gene))

input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)

write.table(fData(input_cds),"cicero/annotate.txt",sep='\t',row.names = F,col.names = T,quote = F)

# tail(fData(input_cds))

#### Generate gene activity scores ####
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]
write.table(as.matrix(unnorm_ga),gzfile("cicero/cicero_gene_activities_unnorm.tsv.gz"),
            sep = '\t',row.names = T,col.names = T, quote = F)

# make a list of num_genes_expressed
num_genes <- pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
write.table(as.matrix(cicero_gene_activities),gzfile("cicero/cicero_gene_activities.tsv.gz"),
            sep = '\t',row.names = T,col.names = T, quote = F)