
# ROOT DIRECTORY (to modify on your computer)
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.script  <- file.path(path.root, "04_scripts/analysis-individual/kan_2024")
path.data   <- file.path(path.root, "03_data/analysis-individual/kan_2024")

## Create phyloseq file
library(qiime2R)
library(phyloseq)
taxonomy_qza<-read_qza(file.path(path.data,"03_qiime2/taxonomy.qza"))

#taxonomic table conversion

taxonomy <- as.matrix(taxonomy_qza$data)
row.names(taxonomy) = taxonomy[,1]
taxonomy = taxonomy[,-1];taxonomy
library(dplyr)
library(tidyr)
taxonomy <- data.frame(taxonomy) %>% select(-Confidence)#tanonomy由phyloseq模式转换为普通表
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right")
taxonomy[is.na(taxonomy)] <- "Unassigned"
taxonomy <- as.matrix(taxonomy)
taxonomy<-phyloseq::tax_table(taxonomy)

otutab_qza =read_qza(file.path(path.data,"03_qiime2/otutab.qza"))
otutab <- as.matrix(otutab_qza$data, taxa_are_rows = TRUE)
otutab<-otu_table(otutab, taxa_are_rows = TRUE)
  
tree_qza=read_qza(file.path(path.data,"03_qiime2/rooted-tree.qza"))
rooted_tree <- phy_tree(tree_qza$data)

metadata_qza <- read.table(file.path(path.data,"01_Metadata/metadata.txt"), header = TRUE, sep = "\t", row.names = 1)
sample_data <- sample_data(metadata_qza)

refseq_qza=read_qza(file.path(path.data,"03_qiime2/repseq.qza"))
refseq <- refseq(refseq_qza$data)

#import data
physeq <- phyloseq(
  otutab,
  taxonomy,
  rooted_tree,
  sample_data,
  refseq
);physeq

physeq_net <- subset_taxa(physeq, Kingdom != "k_mitochondria" & Kingdom != "k_chloroplast");physeq_net

saveRDS(physeq_net,file.path(path.root,"03_data/phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
#____________________________________________________________________
