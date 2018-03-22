# filter cells

setwd("/lustre/user/tianf/06-Human_cell_atlas/test_samples")

exprMatr <- read.table("03-expression/merged/UMIcount_allGenes.txt",header = T,row.names = 1,sep = "\t", stringsAsFactors = F)
print("Before filter:")
dim(exprMatr)

## cell level
# 1
mapStat_per_cell <- read.table("02-alignment/merged/mapStat.txt",header = F, sep = "\t", row.names = 2, stringsAsFactors = F)
#dim(mapStat_per_cell)
hist(mapStat_per_cell$V7,breaks = 50, xlab = "Mapping ratio", main = NA)
# 2
expressed_genes_per_cell <- apply(exprMatr>0, 2, sum)
hist(expressed_genes_per_cell,breaks = 50)
mapStat_per_cell$Genes <- expressed_genes_per_cell[rownames(mapStat_per_cell)]
# 3
nUMI_per_cell <- colSums(exprMatr)
hist(nUMI_per_cell/1e6,breaks = 50)
mapStat_per_cell$nUMI <- nUMI_per_cell
# 4
#mito.genes <- grep(pattern = "^MT-", x = rownames(exprMatr), value = T)
percent.mito <- colSums(exprMatr[grep(pattern = "^MT-", x = rownames(exprMatr)),])/colSums(exprMatr)
hist(percent.mito, breaks = 50)
mapStat_per_cell$mito <- percent.mito
# 5
cor_cells <- cor(x = exprMatr, method = "spearman")
diag(cor_cells) <- NA
avg.cor <- rowMeans(cor_cells, na.rm = T)
hist(avg.cor, breaks = 50)
mapStat_per_cell$avg.cor <- avg.cor

mapStat_per_cell$filter <- mapStat_per_cell$V3>=5e5 & mapStat_per_cell$V7>=0.75 & mapStat_per_cell$Genes>7000 & mapStat_per_cell$Genes<12000 & percent.mito<0.1 & avg.cor>=0.775
sum(mapStat_per_cell$filter)
#mapStat_filtered <- mapStat_per_cell[mapStat_per_cell$V3>=5e5 & mapStat_per_cell$V7>=0.7 & mapStat_per_cell$Genes>=7000, ]

pdf("03-expression/merged/cell_filtering.pdf", width = 4, height = 4, useDingbats = F)
p <- ggplot(mapStat_per_cell, aes(V3/1e6, fill=filter)) + geom_histogram(binwidth = 0.1, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("Reads number / (Million)") + ylab("Sample number")

p <- ggplot(mapStat_per_cell, aes(V7, fill=filter)) + geom_histogram(binwidth = 0.01, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("Mapping ratio") + ylab("Sample number")

p <- ggplot(mapStat_per_cell, aes(nUMI/1e6, fill=filter)) + geom_histogram(binwidth = 0.02, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("UMI count / (Million)") + ylab("Sample number")

p <- ggplot(mapStat_per_cell, aes(Genes, fill=filter)) + geom_histogram(binwidth = 200, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("Dected gene count") + ylab("Sample number")

p <- ggplot(mapStat_per_cell, aes(mito, fill=filter)) + geom_histogram(binwidth = 0.005, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("Expression raio of mitochondrial genes") + ylab("Sample number")

p <- ggplot(mapStat_per_cell, aes(avg.cor, fill=filter)) + geom_histogram(binwidth = 0.005, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("Average Spearman correlation with other cells") + ylab("Sample number")

p <- ggplot(mapStat_per_cell, aes(x = V3/1e6, y = avg.cor, color=filter)) + geom_point()
p + theme_grey() + xlab("Reads number / (Million)") + ylab("Average Spearman correlation with other cells")

dev.off()

### The detected gene number in different gene annotation files
detected_genes_gencode_v26 <- apply(exprMatr>0, 2, sum)

gene_list_refGene_hg19 <- read.table("Genomes/human/hg19_refGene.list", header = F, sep = "\t", stringsAsFactors = F)[,1]
length(gene_list_refGene_hg19)
detected_genes_refGene_hg19 <- apply(exprMatr[rownames(exprMatr)%in%gene_list_refGene_hg19,]>0, 2, sum)

gene_list_refGene_hg38 <- read.table("Genomes/human/hg38_refGene.list", header = F, sep = "\t", stringsAsFactors = F)[,1]
length(gene_list_refGene_hg38)
detected_genes_refGene_hg38 <- apply(exprMatr[rownames(exprMatr)%in%gene_list_refGene_hg38,]>0, 2, sum)

detected_genes_DF <- data.frame(refGene_hg19=detected_genes_refGene_hg19, refGene_hg38=detected_genes_refGene_hg38, gencode_v26=detected_genes_gencode_v26)
detected_genes_DF <- detected_genes_DF[order(detected_genes_DF$gencode_v26),]
detected_genes_DF$cell <- rownames(detected_genes_DF)
library("reshape2")
detected_genes_DF_melt <- melt(data = detected_genes_DF, id.vars = "cell", variable.name = "annotation", value.name = "Gene")
detected_genes_DF_melt$cell <- factor(detected_genes_DF_melt$cell, levels = detected_genes_DF$cell)

pdf("03-expression/merged/numGene_across_annotation.pdf", width = 6, height = 4, useDingbats = F)
p <- ggplot(detected_genes_DF_melt, aes(x=cell, y=Gene, group=annotation, color=annotation)) + geom_line()
p + theme_bw() + theme(axis.text.x=element_blank(), axis.ticks = element_blank())

p <- ggplot(detected_genes_DF_melt, aes(Gene, fill=annotation)) + geom_histogram(binwidth = 100, alpha=.5, position="identity", color="black")
p + theme_grey() + xlab("Dected gene count") + ylab("Sample number")
dev.off()
###

# gene level
expressed_cells_per_gene <- apply(exprMatr[,mapStat_per_cell$filter]>0, 1, sum)
#hist(expressed_cells_per_gene/ncol(exprMatr),breaks = 50)
#abline(v=0.005,lty=2,col="blue")
#ncol(exprMatr)*0.005

nCountsPerGene <- apply(exprMatr, 1, sum)
hist(nCountsPerGene, breaks = 1000000, xlim = c(0,1000))

# filter
# expressed_cells_per_gene/ncol(exprMatr)>=0.005
exprMatr_filtered <- exprMatr[expressed_cells_per_gene>=3, mapStat_per_cell$filter]
print("After filter:")
dim(exprMatr_filtered)

write.table(x = mapStat_per_cell,file = "03-expression/merged/cell_filtering.txt",row.names = T,col.names = F,quote = F,sep = "\t")
write.table(x = exprMatr_filtered,file = "03-expression/merged/UMIcount_filtered.txt",row.names = T,col.names = NA,quote = F,sep = "\t")
