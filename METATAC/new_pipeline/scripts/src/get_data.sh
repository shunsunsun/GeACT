cd /home/chenzy/Lab/otherwork/GeACT/ATAC/data/11-14w
for i in `ls`
do
	if  [[ -d $i ]]; then
		if [[ ! -e $i/RNA/19-22w ]]; then
		mkdir -p $i/RNA/19-22w
		fi
		rsync -av -e ssh chenzy@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/pooled_data/$i/03-expression/merged/cellCluster/Seurat_metaData.txt $i/RNA/19-22w/Seurat_metaData.txt
		rsync -av -e ssh chenzy@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/pooled_data/$i/03-expression/merged/filtering/UMIcount_cellFiltered.txt $i/RNA/19-22w/UMIcount_cellFiltered.txt
	fi
done
