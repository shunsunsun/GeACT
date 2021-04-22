setwd("~/xiegroup/ATAC/datasets/Heart/20191018_1122/merge")
library(ggplot2)
library(pheatmap)
library(gridExtra)
library(GGally)

cleanFqStat = read.table('cleanFqStat.txt',header=T,row.names = 3)

mapStat = read.table('mapStat_human.txt',header=T,row.names = 1,sep='\t')
mapStat = mapStat[rownames(cleanFqStat),]
all.equal(rownames(cleanFqStat),rownames(mapStat))
all.equal(cleanFqStat$Barcode,mapStat$Reads)

# con_freq = read.csv(gzfile('con_freq_human.csv.gz'),sep=',',header=T,row.names = 1,check.names = F)
# all.equal(rownames(cleanFqStat),rownames(con_freq))

hist(mapStat$Reads,breaks = 50)
hist(log10(mapStat$Reads),breaks = 50)

hist(mapStat$Aligned_ratio,breaks = 20)

# hist(mapStat$Num_frags,breaks = 20)
# hist(log10(mapStat$Num_frags),breaks = 20)

hist(mapStat$Num_frags_decon,breaks = 20)
hist(log10(mapStat$Num_frags_decon),breaks = 50)

mapStat$con_rate = 1 - mapStat$Num_frags_decon/mapStat$Num_frags
hist(mapStat$con_rate,breaks = 50)

mapStat$confused_rate = mapStat$Num_frags_con/mapStat$Num_frags
hist(mapStat$confused_rate,breaks = 20)

mapStat$Mito_ratio = mapStat$Num_mito_frags/(mapStat$Num_mito_frags+mapStat$Num_frags)
hist(mapStat$Mito_ratio,breaks = 50)
hist(log10(mapStat$Mito_ratio),breaks=50)

# Batch = data.frame(Batch = cleanFqStat$Batch,row.names = rownames(cleanFqStat))
# Batch$Batch = gsub('PD10_XXL_15_A','',Batch$Batch)

# hist(mapStat$Frags_in_peaks,breaks = 20)
# mapStat$FRiP_C1 = mapStat$Frags_in_peaks/mapStat$Num_frags
# hist(mapStat$FRiP_C1,breaks = 20)
# hist(mapStat$Peaks_detected,breaks = 20)
# plot(mapStat$FRiP_C1~mapStat$log10_Num_frags)

nReads_l = 1e5
nReads_h = 6e6
mapping_l = 0.9
nFrags_l = 4.2
nFrags_h = 5.7
con_rate_h = 0.54
mito_h = 0.1

keep_DF <- data.frame(reads = mapStat$Reads>=nReads_l & mapStat$Reads<=nReads_h,
                      mapping_ratio = mapStat$Aligned_ratio>=mapping_l,
                      num_frags = log10(mapStat$Num_frags_decon)>=nFrags_l & log10(mapStat$Num_frags_decon)<=nFrags_h,
                      con_rate = mapStat$con_rate<=con_rate_h,
                      # FRiP = mapStat$FRiP>=0.5,
                      # peaks_detected = mapStat$Peaks_detected>=2e4,
                      Mito_ratio = mapStat$Mito_ratio <= mito_h
                      )
keep_DF=keep_DF*1
rownames(keep_DF)=rownames(mapStat)
sum(rowSums(keep_DF)==ncol(keep_DF))

mapStat$keep <- (rowSums(keep_DF)==ncol(keep_DF))
# cleanFqStat$keep <- mapStat$keep

pheatmap(keep_DF,
         show_colnames = T, show_rownames = F,
         cluster_rows = F, cluster_cols = F,
         color = c("blue","yellow"), legend = F,
         fontsize = 5)

# mean(mapStat$Reads[mapStat$keep])
# mean(mapStat$Num_frags_decon[mapStat$keep])
# mean(mapStat$Mito_ratio[mapStat$keep])
# mean(mapStat$Peaks_detected[mapStat$keep])
# mapStat$FRiP <- mapStat$Frags_in_peaks/mapStat$Num_frags
# mean(mapStat$FRiP[mapStat$keep])

pdf("./QC_plot.pdf", width = 6, height = 6, useDingbats = F, onefile = T)

print(pheatmap(t(keep_DF),
               show_colnames = F, show_rownames = T,
               cluster_rows = F, cluster_cols = F,
               color = c("blue","yellow"), legend = F,
               fontsize = 5))

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())

hist_read <- ggplot(mapStat, aes(Reads/1e6, fill=keep)) +
  geom_histogram(binwidth = 0.5, alpha=.5, position="identity", color="black") +
  geom_vline(xintercept = nReads_h/1e6, linetype = "dashed") +
  geom_vline(xintercept = nReads_l/1e6, linetype = "dashed")
print(hist_read + theme_classic() +
        xlab("Reads number / (Million)") +
        ylab("Cell number")
        # theme(panel.grid = element_blank(), legend.position = c(0.99, 0.99), legend.justification = c("right", "top"))
        )

p <- ggplot(mapStat, aes(Aligned_ratio, fill=keep)) +
  geom_histogram(binwidth = 0.01, alpha=.5, position="identity", color="black") +
  geom_vline(xintercept = mapping_l, linetype = "dashed")
print(p + theme_classic() +
        xlab("Aligned ratio") +
        ylab("Cell number")
        # theme(panel.grid = element_blank(), legend.position = c(0.99, 0.99), legend.justification = c("right", "top"))
        )

hist_Num_frags <- ggplot(mapStat, aes(log10(Num_frags_decon), fill=keep)) +
  geom_histogram(binwidth = 0.1, alpha=.5, position="identity", color="black") +
  geom_vline(xintercept = nFrags_h, linetype = "dashed") +
  geom_vline(xintercept = nFrags_l, linetype = "dashed")
print(hist_Num_frags + theme_classic() +
        xlab("log10(Number of fragments)") +
        ylab("Cell number"))

scatter <- ggplot(mapStat, aes(x = Reads/1e6, y = log10(Num_frags_decon), color=keep)) +
  geom_point() + theme_classic() +
  xlab("Reads number / (Million)") +
  ylab("log10(Number of fragments)") +
  geom_vline(xintercept = nReads_h/1e6, linetype = "dashed") +
  geom_vline(xintercept = nReads_l/1e6, linetype = "dashed") +
  geom_hline(yintercept = nFrags_h, linetype = "dashed") +
  geom_hline(yintercept = nFrags_l, linetype = "dashed")

print(grid.arrange(hist_read+ theme_classic() + xlab("") + ylab("Cell number") + guides(fill=F),
                   empty, scatter + guides(color=F),
                   hist_Num_frags+ theme_classic() + xlab("") + ylab("Cell number")+coord_flip() + guides(fill=F),
                   ncol=2, nrow=2, widths=c(3, 1.5),heights=c(1.5, 3)))

hist_con_rate <- ggplot(mapStat, aes(con_rate, fill=keep)) +
  geom_histogram(binwidth = 0.02, alpha=.5, position="identity", color="black") +
  geom_vline(xintercept = con_rate_h, linetype = "dashed")
print(hist_con_rate + theme_classic() +
        xlab("Contamination rate") +
        ylab("Cell number"))

scatter <- ggplot(mapStat, aes(x = log10(Num_frags_decon), y = con_rate, color=keep)) +
  geom_point() + theme_classic() +
  xlab("log10(Number of fragments)") +
  ylab("Contamination rate") +
  geom_vline(xintercept = nFrags_h, linetype = "dashed") +
  geom_vline(xintercept = nFrags_l, linetype = "dashed") +
  geom_hline(yintercept = con_rate_h, linetype = "dashed")

print(grid.arrange(hist_Num_frags + theme_classic() + xlab("") + ylab("Cell number") + guides(fill=F),
                   empty, scatter + guides(color=F),
                   hist_con_rate + theme_classic() + xlab("") + ylab("Cell number")+coord_flip() + guides(fill=F),
                   ncol=2, nrow=2, widths=c(3, 1.5),heights=c(1.5, 3)))

hist_mito <- ggplot(mapStat, aes(log10(Mito_ratio), fill=keep)) +
  geom_histogram(binwidth = 0.1, alpha=.5, position="identity", color="black") +
  geom_vline(xintercept = log10(mito_h), linetype = "dashed")
print(hist_mito + theme_classic() +
        xlab("log10(Mitochondrial content)") +
        ylab("Cell number"))

scatter <- ggplot(mapStat, aes(x = log10(Num_frags_decon), y = log10(Mito_ratio), color=keep)) +
  geom_point() + theme_classic() +
  xlab("log10(Number of fragments)") +
  ylab("Mitochondrial content") +
  geom_vline(xintercept = nFrags_h, linetype = "dashed") +
  geom_vline(xintercept = nFrags_l, linetype = "dashed") +
  geom_hline(yintercept = log10(mito_h), linetype = "dashed")

print(grid.arrange(hist_Num_frags+ theme_classic() + xlab("") + ylab("Cell number") + guides(fill=F),
                   empty, scatter + guides(color=F),
                   hist_mito+ theme_classic() + xlab("") + ylab("Cell number")+coord_flip() + guides(fill=F),
                   ncol=2, nrow=2, widths=c(3, 1.5),heights=c(1.5, 3)))

# hist_FRiP <- ggplot(mapStat, aes(FRiP, fill=keep)) +
#   geom_histogram(binwidth = 0.02, alpha=.5, position="identity", color="black") +
#   geom_vline(xintercept = 0.5, linetype = "dashed")
# print(hist_FRiP + theme_classic() +
#         xlab("FRiP") +
#         ylab("Cell number"))
# 
# hist_peak <- ggplot(mapStat, aes(Peaks_detected, fill=keep)) +
#   geom_histogram(binwidth = 5e3, alpha=.5, position="identity", color="black") +
#   geom_vline(xintercept = 2e4, linetype = "dashed")
# print(hist_peak + theme_classic() +
#         xlab("Number of detected peaks") +
#         ylab("Cell number"))
# 
# scatter <- ggplot(mapStat, aes(x = log10(Num_frags), y = Peaks_detected, color=keep)) +
#   geom_point() + theme_classic() +
#   xlab("log10(Number of fragments)") +
#   ylab("Number of detected peaks") +
#   geom_vline(xintercept = 4.8, linetype = "dashed") +
#   geom_hline(yintercept = 2e4, linetype = "dashed")
# 
# print(grid.arrange(hist_Num_frags+ theme_classic() + xlab("") + ylab("Cell number") + guides(fill=F),
#                    empty, scatter + guides(color=F),
#                    hist_peak+ theme_classic() + xlab("") + ylab("Cell number")+coord_flip() + guides(fill=F),
#                    ncol=2, nrow=2, widths=c(3, 1.5),heights=c(1.5, 3)))
# 
# scatter <- ggplot(mapStat, aes(x = log10(Num_frags), y = FRiP, color=keep)) +
#   geom_point() + theme_classic() +
#   xlab("log10(Number of fragments)") +
#   ylab("FRiP") +
#   geom_vline(xintercept = 4.8, linetype = "dashed") +
#   geom_hline(yintercept = 0.5, linetype = "dashed")
# 
# print(grid.arrange(hist_Num_frags+ theme_classic() + xlab("") + ylab("Cell number") + guides(fill=F),
#                    empty, scatter + guides(color=F),
#                    hist_FRiP+ theme_classic() + xlab("") + ylab("Cell number")+coord_flip() + guides(fill=F),
#                    ncol=2, nrow=2, widths=c(3, 1.5),heights=c(1.5, 3)))

# print(pheatmap(log10(con_freq+1),
#                color=colorRampPalette(c("white","red"))(200),
#                show_colnames = F, show_rownames = F,
#                cluster_rows = F, cluster_cols = F,
#                main = '20190407',cellwidth = 0.4,cellheight = 0.4,
#                annotation_row = Batch,annotation_col = Batch,
#                annotation_legend = F))

dev.off()

pdf("./QC_summary.pdf", width = 12, height = 12, useDingbats = F, onefile = T)
p <- ggpairs(data=cbind(mapStat[,c(1,6,11,14)],log10(mapStat[,c(8,13)])),columns = c(1:3,5,6), aes(color=keep))
print(p + theme(axis.text=element_text(size=8),text=element_text(size=8)))
dev.off()

mapStat$Batch = cleanFqStat$Batch
mapStat$R = NA
mapStat$C = NA
for(i in 1:nrow(mapStat))
{
  id = strsplit(rownames(mapStat)[i],'_')[[1]][length(strsplit(rownames(mapStat)[i],'_')[[1]])]
  mapStat$R[i] = as.numeric(strsplit(id,'-')[[1]][2])
  mapStat$C[i] = as.numeric(strsplit(id,'-')[[1]][1])
}
mapStat_plate = split(mapStat,mapStat$Batch)

Batches = levels(mapStat$Batch)

pdf("./QC_matrix.pdf", width = 8, height = 6, useDingbats = F, onefile = T)

for(batch in Batches)
{
  p <- ggplot(mapStat_plate[[batch]], aes(x=factor(C), y=factor(R)))
  print(p + geom_point(aes(size=Num_frags_decon,color=con_rate)) +
          scale_size_area(max_size = 10, guide = F) +
          geom_text(aes(label = round(Num_frags_decon/1e4,2)),size=3,nudge_y = -0.3) +
          scale_colour_gradient2(low = "aquamarine", mid = "dodgerblue", high = "blue",
                                 midpoint = 0.2, space = "Lab", limits = c(0,0.5),
                                 na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
          scale_x_discrete(limits=levels(factor(mapStat_plate[[batch]]$C))) +
          scale_y_discrete(limits=rev(levels(factor(mapStat_plate[[batch]]$R)))) +
          labs(title=batch, x="Column", y="Row",colour = "Contamination rate") +
          theme(plot.title = element_text(hjust = 0.5, size = 20),
                legend.position = "bottom",legend.key.width = unit(2,"cm")))

  p <- ggplot(mapStat_plate[[batch]], aes(x=factor(C), y=factor(R)))
  print(p + geom_point(aes(size=Num_frags_decon,color=Mito_ratio)) +
          scale_size_area(max_size = 10, guide = F) +
          geom_text(aes(label = round(Num_frags_decon/1e4,2)),size=3,nudge_y = -0.3) +
          scale_colour_gradient2(low = "aquamarine", mid = "dodgerblue", high = "blue",
                                 midpoint = 0.05, space = "Lab", limits = c(0,0.1),
                                 na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
          scale_x_discrete(limits=levels(factor(mapStat_plate[[batch]]$C))) +
          scale_y_discrete(limits=rev(levels(factor(mapStat_plate[[batch]]$R)))) +
          labs(title=batch, x="Column", y="Row",colour = "Mitochondrial content") +
          theme(plot.title = element_text(hjust = 0.5, size = 20),
                legend.position = "bottom",legend.key.width = unit(2,"cm")))

  p <- ggplot(mapStat_plate[[batch]], aes(x=factor(C), y=factor(R))) 
  print(p + geom_point(aes(size=Num_frags_decon,color=keep)) +
          scale_size_area(max_size = 10, guide = F) +
          geom_text(aes(label = round(Num_frags_decon/1e4,2)),size=3,nudge_y = -0.3) +
          scale_x_discrete(limits=levels(factor(mapStat_plate[[batch]]$C))) +
          scale_y_discrete(limits=rev(levels(factor(mapStat_plate[[batch]]$R)))) +
          labs(title=batch, x="Column", y="Row",colour = "keep") +
          theme(plot.title = element_text(hjust = 0.5, size = 20),
                legend.position = "bottom",legend.key.width = unit(2,"cm")))
}

dev.off()

# write.table(x = cleanFqStat,file = "./cleanFqStat_filter.txt",row.names = T,col.names = T,quote = F,sep = "\t")
write.table(x = mapStat,file = "./mapStat_human_filter.txt",row.names = T,col.names = T,quote = F,sep = "\t")
write.table(x = keep_DF,file = "./keep_DF.txt",row.names = T,col.names = T,quote = F,sep = "\t")
write.table(x = rownames(mapStat)[mapStat$keep],file = "kept_cell.txt",
            row.names = F,col.names = F,quote = F)
