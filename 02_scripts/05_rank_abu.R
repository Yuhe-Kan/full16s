# Fig1A

# ROOT DIRECTORY (to modify on your computer)
getwd()
rm(list = ls())
path.root <- "E:/b_myproject/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
# ps <- transform_sample_counts(ps, function(x) x / sum(x));ps
psmelt<-psmelt(ps)

psmelt$Rank<-1:nrow(psmelt)

library(ggplot2)
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
 # panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 13, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

p <- ggplot(psmelt, aes(x = Rank, y = Abundance, color = Group)) +
  geom_line() +
  geom_point() +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(x = "Abundance Rank", y = "Reads number") +
  theme_minimal()+
  BioR.theme
p
ggsave(file.path(path.plot,"01_rank_abun/01_rank_abun.jpeg"), plot = p, width = 6, height = 4)
