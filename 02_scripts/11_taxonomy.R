
# Fig2A
rm(list = ls())
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(ggplot2)
library(MicrobiotaProcess)
library(phyloseq)
library(file2meco)
library(microeco)

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
mt <- phyloseq2meco(ps)
# show the pie chart with the group mean values,phylum
t4 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 30, groupmean = "Group")
View(t4$data_abund)

View(phyloseq::tax_table(ps))


# ROOT DIRECTORY (to modify on your computer)
rm(list = ls())
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
ps = rarefy_even_depth(ps, sample.size = min(sample_sums(ps)))
sample_sums(ps)
##observed_otu转化为百分比
# ps_perc <- transform_sample_counts(ps, function(x) x / sum(x));ps_perc
ps_perc<- microbiome::transform(ps, "compositional")#“Z,” “clr,” “hellinger,"compositional",” and “shift”
sample_sums(ps_perc)

#ref:https://github.com/gmteunisse/fantaxtic
# devtools::install_github("gmteunisse/fantaxtic")
require("fantaxtic")
require("phyloseq")
library(dplyr)
library(ggplot2)

top_nested <- nested_top_taxa(ps_perc,
                              top_tax_level = "Phylum",
                              nested_tax_level = "Genus",
                              n_top_taxa = 5, 
                              n_nested_taxa = 3,
                              nested_merged_label = "other <tax>",
                              by_proportion = T,
                              #grouping = "Group",
                              FUN = mean)

sample_labels <- c("A1" = "SG1", "A2" = "SG2","A3" = "SG3","A4" = "SG4","A5" = "SG5","A6" = "SG6",
                   "B1" = "PA1", "B2" = "PA2","B3" = "PA3","B4" = "PA4","B5" = "PA5","B6" = "PA6",
                   "C1" = "SV1", "C2" = "SV2","C3" = "SV3","C4" = "SV4","C5" = "SV5","C6" = "SV6")

p1=plot_nested_bar(top_nested$ps_obj,
                   top_level = "Phylum",
                   nested_level = "Genus",
                   nested_merged_label = "NA and other <tax>",
                   legend_title = "Phylum and Genus") +
  facet_wrap(~Group,
             scales = "free_x") +
  scale_x_discrete(labels = sample_labels) +
  #labs(title = "Relative abundances of the top 5 genus for each of the top 5 phyla") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 8, 
                                  face = "bold"),
        legend.key.size = unit(10, 
                               "points"),
        axis.text.x = element_text(angle = 0, hjust = 0.5));p1


ggsave(file.path(path.plot,"06_taxnomy/06_taxnomy-4.jpeg"), plot = p1, width = 8, height = 6)
#______________________________________________________________________________________________________________________________

# 查看各个门百分比
rm(list = ls())
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(ggplot2)
library(MicrobiotaProcess)
library(phyloseq)

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
mt <- phyloseq2meco(ps)
# show the pie chart with the group mean values,phylum
t4 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 5, groupmean = "Group")
View(t4$data_abund)
#--------------------------------------------------------------------------------------------------------------------
# Fig3A
# visualize the relative abundance of top 30 genus.
rm(list = ls())
getwd()
path.root <- "E:/b_myproject/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(ggplot2)
library(MicrobiotaProcess)
library(phyloseq)

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))

#稀疏化（样本深度均一化）
ps <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123, replace = FALSE)
sample_sums(ps)

mpse <- as.MPSE(ps)

mpse %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups
    .abundance=RareAbundance,
    .group=Group
  )

p13 <- mpse %>%
  mp_plot_abundance(
    .abundance=RareAbundance,
    .group=Group,
    plot.group = T,
    taxa.class = Genus,
    topn = 30,
    relative = TRUE, # or FALSE
    geom = "flowbar"
  )+
  theme(strip.background = element_rect(colour=NA, fill="grey"),axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "right",  # 确保图例在右侧
        legend.text = element_text(size = 8),  # 调整图例字体大小
        legend.title = element_text(size = 8))
p13
ggsave(file.path(path.plot,"06_taxnomy/06_taxnomy_genus-1.jpeg"), plot = p13, width = 8, height = 6)
#--------------------------------------------------------------------------------------------------------------------

