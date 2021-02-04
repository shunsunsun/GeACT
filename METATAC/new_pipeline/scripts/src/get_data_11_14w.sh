cd /home/chenzy/Lab/otherwork/GeACT/ATAC/data/11-14w
lt=(01_stomach 02_small_intestine 03_kidney  04_lung 17_thymus)

for i in ${lt[@]}
do
	if  [[ -d $i ]]; then
		if [[ ! -e $i/RNA/11-14w ]]; then
		mkdir -p $i/RNA/11-14w
		fi
		rsync -av -e ssh chenzy@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/pooled_data_14w/$i/03-expression/merged/cellCluster/Seurat_metaData.txt $i/RNA/11-14w/Seurat_metaData.txt
		rsync -av -e ssh chenzy@162.105.250.222:/rd/user/tianf/06-Human_cell_atlas/pooled_data_14w/$i/03-expression/merged/filtering/UMIcount_cellFiltered.txt $i/RNA/11-14w/UMIcount_cellFiltered.txt
	fi
done
