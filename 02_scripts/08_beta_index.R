# Fig1C
rm(list = ls())
getwd()
path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
library(microbiome)
library(phyloseq)
# install file2meco package (>= 0.7.1)
library(file2meco)

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
sample_sums(ps)
# summarize_phyloseq(ps)

#标准化
# ps <- rarefy_even_depth(ps, sample.size = 1000, rngseed = 123, replace = FALSE)
ps<- microbiome::transform(ps, "compositional")#“Z,” “clr,” “hellinger,"compositional",” and “shift”
sample_sums(ps)

# 进行Adonis分析
library(vegan)
set.seed(123)  # 为了结果可重复
# pcoa <- ordinate(ps, method = "PCoA", distance = "bray") # method:"DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
# pcoa#查看

otutab<-as.data.frame(t(otu_table(ps)))
class(otutab)
metadata <- data.frame(sample_data(ps))
class(metadata)

BC.dist = vegdist(otutab, distance="bray", k=2, trymax=1000)

adonis_result<-adonis2(BC.dist ~ Group, data = metadata, by='margin', permutations = 9999);adonis_result
F_value <- adonis_result[1,4];F_value
P_value <- adonis_result[1,5];P_value

p_adjusted <- p.adjust(
  P_value,
  method = "BH"
)

p_adjusted

R2 = adonis_result[1,3];R2
#______________________________________________________________________________________________________________________________

# PCoA
mt <- phyloseq2meco(ps)
mt$cal_betadiv()
theme_set(theme_bw())
t1 <- trans_beta$new(dataset = mt, group = "Group", measure = "bray")
t1$cal_ordination(method = "PCoA")# method:"DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
# extract the axis scores
tmp <- t1$res_ordination$scores
# differential test with trans_env class
t2 <- trans_env$new(dataset = mt, add_data = tmp[, 1:2])
# 'KW_dunn' for non-parametric test
t2$cal_diff(group = "Group", method = "KW_dunn")

p1 <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))+
  annotate("label", x = Inf, y = Inf, label = paste("F =", F_value, "\np =", p_adjusted,"\nR2 =",R2), 
           hjust = 1, vjust = 1, size = 3.0, color = "black", fill = "white", linewidth = 0.5, label.r = unit(0.2, "lines")) +
  theme(plot.margin = margin(10, 10, 10, 10));p1
# groups order in p2 is same with p1; use legend.position = "none" to remove redundant legend
p2 <- t2$plot_diff(measure = "PCo1", add_sig = T) + theme_bw() + coord_flip() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
p3 <- t2$plot_diff(measure = "PCo2", add_sig = T) + theme_bw() + 
  theme(legend.position = "none", axis.title.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
# height of the upper figure and width of the right-hand figure are both 0.2-fold of the main figure
g <- p1 %>% aplot::insert_top(p2, height = 0.2) %>% aplot::insert_right(p3, width = 0.2);g
ggsave(file.path(path.plot,"04_beta_index/04_beta_index-pcoa-1.1.jpeg"), plot = g, width = 8, height = 6)
#______________________________________________________________________________________________________________________________


# Test homogeneity of multivariate dispersion
set.seed(123)

bd <- betadisper(BC.dist, metadata$Group)

# Permutation test
bd_test <- permutest(
  bd,
  permutations = 9999
)

print(bd_test)

# Descriptive statistics for each group
print(bd)



# #plot and compare the group distances
# rm(list = ls())
# getwd()
# path.root <- "D:/test/kan_2024_full16s"
# path.plot <- file.path(path.root, "05_plot")
# path.data <- file.path(path.root, "03_data")
# 
# library(microeco)
# library(magrittr)
# library(ggplot2)
# library(aplot)
# library(microbiome)
# library(phyloseq)
# # install file2meco package (>= 0.7.1)
# library(file2meco)
# 
# ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
# sample_sums(ps)
# # summarize_phyloseq(ps)
# 
# # ps <- rarefy_even_depth(ps, sample.size = 1000, rngseed = 123, replace = FALSE)
# ps<- microbiome::transform(ps, "compositional")#“Z,” “clr,” “hellinger,"compositional",” and “shift”
# sample_sums(ps)
# 
# mt <- phyloseq2meco(ps)
# mt$cal_betadiv()
# 
# t1 <- trans_beta$new(dataset = mt, group = "Group", measure = "bray")
# # calculate and plot sample distances within groups
# t1$cal_group_distance(within_group = TRUE)
# # return t1$res_group_distance
# # perform Wilcoxon Rank Sum and Signed Rank Tests
# t1$cal_group_distance_diff(method = "wilcox")
# # plot_group_order parameter can be used to adjust orders in x axis
# p4 <- t1$plot_group_distance(add = "mean")+
#   theme(legend.position = "right", 
#         axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
#         axis.title.x = element_blank());p4
# ggsave(file.path(path.plot,"04_beta_index/04_beta_index-9.jpeg"), plot = p4, width = 6, height = 6)
# #______________________________________________________________________________________________________________________________
# 
