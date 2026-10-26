# # ROOT DIRECTORY (to modify on your computer)
# getwd()
# path.root <- "E:/b_myproject/kan_2024_full16s"
# path.script  <- file.path(path.root, "04_scripts/analysis-individual/kan_2024")
# path.data   <- file.path(path.root, "03_data/analysis-individual/kan_2024")
# #____________________________________________________________________
# 
# ## 创建phyloseq文件
# library(qiime2R)
# library(phyloseq)
# taxonomy_qza<-read_qza(file.path(path.data,"03_qiime2/taxnomy.qza"))
# 
# #分类表转换
# taxonomy <- tax_table(as.matrix(taxonomy_qza$data))
# row.names(taxonomy) = taxonomy[,1]
# taxonomy = taxonomy[,-1];taxonomy
# library(dplyr)
# library(tidyr)
# taxonomy <- data.frame(taxonomy) %>% select(-Confidence)#tanonomy由phyloseq模式转换为普通表
# taxonomy <- taxonomy %>%
#   separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right")
# taxonomy[is.na(taxonomy)] <- "Unassigned"
# taxonomy <- tax_table(as.matrix(taxonomy))
# 
# 
# otutab_qza =read_qza(file.path(path.data,"03_qiime2/otutab.qza"))
# otutab <- otu_table(as.matrix(otutab_qza$data), taxa_are_rows = TRUE)
# 
# tree_qza=read_qza(file.path(path.data,"03_qiime2/rooted_tree.qza"))
# rooted_tree <- phy_tree(tree_qza$data)
# 
# metadata_qza <- read.table(file.path(path.data,"01_Metadata/metadata.txt"), header = TRUE, sep = "\t", row.names = 1)
# # library(readxl)
# # metadata_qza <- read_excel(file.path(path.data,"01_Metadata/metadata.xlsx"))
# # library(tibble)
# # metadata_qza <- metadata_qza %>% 
# #   tibble::column_to_rownames("sample-id") 
# sample_data <- sample_data(metadata_qza)
# 
# refseq_qza=read_qza(file.path(path.data,"03_qiime2/repseq.qza"))
# refseq <- refseq(refseq_qza$data)
# 
# #导入数据
# physeq <- phyloseq(
#   otutab,
#   taxonomy,
#   rooted_tree,
#   sample_data,
#   refseq
# );physeq
# 
# tax_table(physeq)
# 
# physeq <- subset_taxa(physeq, tax_table(physeq)[, "Phylum"] != "Unassigned"
#                               & tax_table(physeq)[, "Phylum"] != " p__"
#                                & tax_table(physeq)[, "Kingdom"] != "k_mitochondria"
#                                & tax_table(physeq)[, "Kingdom"] != "k_chloroplast")
# 
# physeq
# 
# 
# saveRDS(physeq,file.path(path.root,"03_data/phyloseq-objects/physeq-with-phylotree/phyloseq_kan2024.rds"))
# #____________________________________________________________________
# 
# physeq_notree <- phyloseq(
#   otutab,
#   taxonomy,
#   #rooted_tree,
#   sample_data,
#   refseq
# );physeq_notree
# 
# physeq_notree <- subset_taxa(physeq_notree, tax_table(physeq_notree)[, "Phylum"] != "Unassigned"
#                       & tax_table(physeq_notree)[, "Phylum"] != " p__"
#                       & tax_table(physeq_notree)[, "Kingdom"] != "k_mitochondria"
#                       & tax_table(physeq_notree)[, "Kingdom"] != "k_chloroplast")
# 
# physeq_notree
# 
# saveRDS(physeq_notree,file.path(path.root,"03_data/phyloseq-objects/physeq-without-phylotree/phyloseq_kan2024.rds"))
# 
# ps = base::readRDS(file.path(path.root,"03_data/phyloseq-objects/physeq-without-phylotree/phyloseq_kan2024.rds"))
#____________________________________________________________________


# for network
# ROOT DIRECTORY (to modify on your computer)
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.script  <- file.path(path.root, "04_scripts/analysis-individual/kan_2024")
path.data   <- file.path(path.root, "03_data/analysis-individual/kan_2024")

## 创建phyloseq文件
library(qiime2R)
library(phyloseq)

taxonomy_qza<-read_qza(file.path(path.data,"03_qiime2/taxonomy.qza"))

#分类表转换
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
# library(readxl)
# metadata_qza <- read_excel(file.path(path.data,"01_Metadata/metadata.xlsx"))
# library(tibble)
# metadata_qza <- metadata_qza %>% 
#   tibble::column_to_rownames("sample-id") 
sample_data <- sample_data(metadata_qza)

refseq_qza=read_qza(file.path(path.data,"03_qiime2/repseq.qza"))
refseq <- refseq(refseq_qza$data)

#导入数据
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
