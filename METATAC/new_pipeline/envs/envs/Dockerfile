From rocker/tidyverse:4.0.3

COPY installBioc.r installGithub.r /rocker_scripts/
RUN chmod 755 /rocker_scripts/installBioc.r /rocker_scripts/installGithub.r \
  && ln -nsf /rocker_scripts/installBioc.r /usr/local/bin/installBioc.r \
  && ln -nsf /rocker_scripts/installGithub.r /usr/local/bin/installGithub.r

RUN apt-get update\
  && apt-get install -y --no-install-recommends libgsl-dev libhdf5-dev libglpk-dev libxt-dev libbz2-dev liblzma-dev libboost-dev libboost-filesystem-dev \
  && apt-get clean\
  && rm -rf /var/lib/apt/lists/ \
  && installBioc.r -t 600 -e -s -m https://mirrors.tuna.tsinghua.edu.cn/bioconductor \
	multtest \
	rhdf5 \
	GenomicFeatures \
	BSgenome.Hsapiens.UCSC.hg38 chromVAR \
	ComplexHeatmap \
	motifmatchr \
	SingleCellExperiment \
	batchelor \
	DelayedMatrixStats \
	limma \
	monocle \
	TxDb.Hsapiens.UCSC.hg38.knownGene\
  && installGithub.r -t 600 -b https://mirrors.tuna.tsinghua.edu.cn/bioconductor \ 
	GreenleafLab/ArchR@release_1.0.1 \
	cole-trapnell-lab/monocle3 \
	cole-trapnell-lab/cicero-release@monocle3 \
	immunogenomics/harmony \
	satijalab/seurat-data \
	mojaveazure/seurat-disk \
	satijalab/seurat-wrappers \
  && wget https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.2.3.tar.gz \
  && install2.r --error \
    hdf5r\
    Seurat_3.2.3.tar.gz \
	cowplot \
	pheatmap \
	networkD3 \
	plyr \
	gridExtra \
    Cairo \
  && rm Seurat_3.2.3.tar.gz \
  && rm -rf /tmp/R* /tmp/downloaded_packages

VOLUME /data/Lab
