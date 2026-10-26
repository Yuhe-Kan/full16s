# Fig1B
# ## # rarefaction visualization
# ROOT DIRECTORY (to modify on your computer)
rm(list = ls())
getwd()
path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))

library(MicrobiotaProcess)
library(ggplot2)
mpse <- as.MPSE(ps) 
mpse %<>% mp_rrarefy()

mpse %<>% 
  mp_cal_rarecurve(
    .abundance = RareAbundance,
    chunks = 100
  )

mpse %>% print(width=180)

theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  # panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 16, family="Arial"),
  axis.text = element_text(size = 13, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank(),
  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

p5 <- mpse %>% 
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve, 
    .alpha = c("Shannon","Chao1","Observe"), 
    .group = Group, 
    plot.group = TRUE
  ) +
  #scale_color_manual(values=c("#00A087FF", "#3C5488FF","grey","blue")) +
  scale_fill_manual(values=c("#FDDED7", "#A28CC2","#CCD376","#5CB0C3"),guide="none")  +
  theme+
  labs(y="Value", x="Reads Number")
p5
ggsave(file.path(path.plot,"02_rarefraction_curve/02_rarefraction_curve.jpeg"),plot=p5, width = 10, height = 3)


#------------------------------------------------------------------------
# calculate good's coverage

library(phyloseq)
library(microeco)
library(file2meco)
library(dplyr)

mt <- phyloseq2meco(ps)

t1 <- trans_alpha$new(dataset = mt, group = "Group")

otu <- as(otu_table(ps), "matrix")

if (taxa_are_rows(ps)) {
  otu <- t(otu)
}

coverage_df <- data.frame(
  Sample = rownames(otu),
  N = rowSums(otu),
  F1 = rowSums(otu == 1),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Coverage = ifelse(N > 0, 1 - F1 / N, NA_real_)
  )

metadata <- data.frame(sample_data(ps))
metadata$Sample <- rownames(metadata)

coverage_df <- coverage_df %>%
  left_join(
    metadata %>% select(Sample, Group),
    by = "Sample"
  )

print(coverage_df)

coverage_summary <- coverage_df %>%
  group_by(Group) %>%
  summarise(
    n = sum(!is.na(Coverage)),
    Mean = mean(Coverage, na.rm = TRUE),
    SD = sd(Coverage, na.rm = TRUE),
    Min = min(Coverage, na.rm = TRUE),
    Max = max(Coverage, na.rm = TRUE),
    `Mean ± SD` = sprintf("%.6f ± %.6f", Mean, SD),
    .groups = "drop"
  )

print(coverage_summary)
#------------------------------------------------------------------------
#或
#------------------------------------------------------------------------
#calculate goods' coverage
wd=/mnt/e/b_myproject/kan_2024_full16s/03_data/analysis-individual/kan_2024/03_qiime2

conda activate qiime2-amplicon-2023.5

qiime diversity alpha \
--i-table $wd/otutab.qza \
--p-metric goods_coverage \
--o-alpha-diversity $wd/coverage

qiime metadata tabulate \
--m-input-file $wd/coverage.qza \
--o-visualization $wd/coverage.qzv
#______________________________________________________________________________________________________________________________


