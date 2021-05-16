
ident2clgrp <- function(ident.input) {
  ident.input[grepl("^Epi", ident.input)] <- "Epithelial"
  ident.input[grepl("^Endo", ident.input)] <- "Endothelial"
  ident.input[grepl("^SM-", ident.input)] <- "Smooth muscle"
  ident.input[grepl("^SM$", ident.input)] <- "Smooth muscle"
  ident.input[grepl("^SKM$", ident.input)] <- "Skeletal muscle"
  ident.input[grepl("^Fibro", ident.input)] <- "Fibroblast"
  # immune
  ident.input[grepl("^B-", ident.input)] <- "B"
  ident.input[grepl("^Pro-B", ident.input)] <- "B"
  ident.input[grepl("^Pre-B", ident.input)] <- "B"
  ident.input[grepl("^DC/Macro", ident.input)] <- "DC/Macrophage"
  ident.input[grepl("^Mast-", ident.input)] <- "Mast"
  ident.input[grepl("^Macro", ident.input)] <- "Macrophage"
  ident.input[grepl("^Neutrophil-", ident.input)] <- "Neutrophil"
  ident.input[grepl("^NKT-", ident.input)] <- "NKT"
  ident.input[grepl("^T-", ident.input)] <- "T"
  ident.input[grepl("^Pre-T", ident.input)] <- "T"
  ident.input[grepl("^Immune-", ident.input)] <- "Immune"
  #
  ident.input[grepl("^Erythrocyte-", ident.input)] <- "Erythrocyte"
  ident.input[grepl("^CACNA1A-", ident.input)] <- "CACNA1A"
  ident.input[ident.input %in% c("PT", "LoH", "LoH-Prog", "DT", "PC-CLU", "PC-BCAT1", "Podocyte-GPC3", "Podocyte-PLA2R1")] <- "Epithelial"
  ident.input[ident.input %in% c("LSEC")] <- "Endothelial"
  ident.input[grepl("^Sertoli-", ident.input)] <- "Sertoli"
  ident.input[grepl("^Granulosa-", ident.input)] <- "Granulosa"
  # FGC
  ident.input[grepl("^SSC$", ident.input)] <- "FGC"
  return(ident.input)
}