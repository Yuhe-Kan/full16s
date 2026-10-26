
# Fig2B
rm(list = ls())
getwd()
path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

# library(microbiomeSeq)
library(phyloseq)
library(ggpubr)
library(dplyr)

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
ps = rarefy_even_depth(ps, sample.size = min(sample_sums(ps)))
ps<- microbiome::transform(ps, "compositional")#“Z,” “clr,” “hellinger,"compositional",” and “shift”
sample_sums(ps)

ps_phylum<-tax_glom(ps,taxrank = rank_names(ps)[2])
taxtab(ps_phylum)
ps_phylum <- subset_taxa(ps_phylum, Phylum %in% c(" p__Proteobacteria"," p__Cyanobacteria",
                                                  " p__Actinobacteria",
                                                  " p__Firmicutes"," p__Bacteroidetes"))
ps_phylum_melt<-psmelt(ps_phylum)


library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)

ps_phylum_melt <- ps_phylum_melt %>%
  mutate(
    Group = factor(Group, levels = c("CK", "SG", "PA", "SV")),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(!is.na(Group), !is.na(Abundance), !is.na(Phylum))

# Pairwise comparisons
comparisons <- list(
  c("CK", "SG"),
  c("CK", "PA"),
  c("CK", "SV"),
  c("SG", "PA"),
  c("SG", "SV"),
  c("PA", "SV")
)


kw_results <- ps_phylum_melt %>%
  group_by(Phylum) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "BH")
  )

print(kw_results)


dunn_results <- ps_phylum_melt %>%
  group_by(Phylum) %>%
  dunn_test(
    Abundance ~ Group,
    p.adjust.method = "BH"
  ) %>%
  ungroup()

print(dunn_results)

# Export statistical results
write.csv(
  kw_results,
  "phylum_Kruskal_Wallis_BH.csv",
  row.names = FALSE
)

write.csv(
  dunn_results,
  "phylum_Dunn_pairwise_BH.csv",
  row.names = FALSE
)


plot_dunn <- dunn_results %>%
  mutate(
    p.adj.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj < 0.05    ~ "*",
      TRUE            ~ "ns"
    )
  )

# Calculate bracket positions separately for each phylum
y_positions <- ps_phylum_melt %>%
  group_by(Phylum) %>%
  summarise(
    y_max = max(Abundance * 100, na.rm = TRUE),
    y_min = min(Abundance * 100, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(y_max > y_min, y_max - y_min, abs(y_max) + 1)
  )

plot_dunn <- plot_dunn %>%
  left_join(y_positions, by = "Phylum") %>%
  group_by(Phylum) %>%
  arrange(match(group1, c("CK", "SG", "PA", "SV")),
          match(group2, c("CK", "SG", "PA", "SV"))) %>%
  mutate(
    y.position = y_max + row_number() * 0.12 * y_range
  ) %>%
  ungroup()


p1 <- ggplot(
  data = ps_phylum_melt,
  aes(x = Group, y = Abundance * 100, color = Group)
) +
  geom_boxplot(
    aes(fill = Group),
    position = position_dodge(1),
    color = "black",
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "#a6cee3",
    width = 0.15,
    size = 2.5,
    alpha = 0.8
  ) +
  stat_pvalue_manual(
    plot_dunn,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 3,
    hide.ns = TRUE
  ) +
  facet_wrap(~Phylum, scales = "free_y", ncol = 3) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey"),
    strip.text = element_text(
      color = "blue",
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0,
      color = "black"
    )
  ) +
  labs(
    x = "",
    y = "Relative Abundance (%)",
    title = NULL
  ) +
  scale_fill_manual(
    values = c(
      "#a6cee3", "#bebada", "#e5c494",
      "#386cb0", "#1b9e77", "#fdc086"
    )
  ) +
  scale_color_manual(
    values = c(
      "#a6cee3", "#bebada", "#e5c494",
      "#386cb0", "#1b9e77", "#fdc086"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  coord_cartesian(clip = "off")

print(p1)

ggsave(
  "phylum_relative_abundance_Dunn_BH.pdf",
  p1,
  width = 12,
  height = 10
)

ggsave(
  "phylum_relative_abundance_Dunn_BH.tiff",
  p1,
  width = 8,
  height = 6,
  dpi = 600
)

# ggsave(file.path(path.plot,"07_diff/07_diff-phylum-1.jpeg"),plot=p1, width = 8, height = 6)
#______________________________________________________________________________________________________________________________

# Fig3B
rm(list = ls())
getwd()
path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))

# ps <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123, replace = FALSE)
# sample_sums(ps)
ps<- microbiome::transform(ps, "compositional")#“Z,” “clr,” “hellinger,"compositional",” and “shift”
sample_sums(ps)

library(phyloseq)
library(ggpubr)
library(dplyr)
ps_genus<-tax_glom(ps,taxrank = rank_names(ps)[6])
taxtab(ps_genus)
ps_genus <- subset_taxa(ps_genus, Genus %in% c(" g__Salinicoccus"," g__Kushneria"," g__Sphingomonas"," g__Delftia",
                                               " g__Rhizobium"," g__Exiguobacterium", " g__Pediococcus",
                                               " g__Flavisolibacter"," g__Microcoleus"," g__Tumebacillus"," g__Jeotgalibacillus",
                                               " g__Halomonas"," g__Archangium"," g__Limosilactobacillus"," g__Rubellimicrobium",
                                               " g__Kocuria"," g__Ramlibacter"," g__Devosia"," g__Acinetobacter"
                                               ," g__Saccharibacillus"," g__Rhodocytophaga"," g__Pseudoxanthomonas"," g__Flexibacter"," g__Paracoccus"))
ps_genus_melt<-psmelt(ps_genus)


library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)


ps_genus_melt <- ps_genus_melt %>%
  mutate(
    Group = factor(Group, levels = c("CK", "SG", "PA", "SV")),
    Abundance = as.numeric(Abundance)
  ) %>%
  filter(!is.na(Group), !is.na(Abundance), !is.na(Genus))

# Pairwise comparisons
comparisons <- list(
  c("CK", "SG"),
  c("CK", "PA"),
  c("CK", "SV"),
  c("SG", "PA"),
  c("SG", "SV"),
  c("PA", "SV")
)


kw_results <- ps_genus_melt %>%
  group_by(Genus) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "BH")
  )

print(kw_results)


dunn_results <- ps_genus_melt %>%
  group_by(Genus) %>%
  dunn_test(
    Abundance ~ Group,
    p.adjust.method = "BH"
  ) %>%
  ungroup()

print(dunn_results)

# Export statistical results
write.csv(
  kw_results,
  "genus_Kruskal_Wallis_BH.csv",
  row.names = FALSE
)

write.csv(
  dunn_results,
  "genus_Dunn_pairwise_BH.csv",
  row.names = FALSE
)


plot_dunn <- dunn_results %>%
  mutate(
    p.adj.signif = case_when(
      p.adj <= 0.0001 ~ "****",
      p.adj <= 0.001  ~ "***",
      p.adj <= 0.01   ~ "**",
      p.adj < 0.05    ~ "*",
      TRUE            ~ "ns"
    )
  )

# Calculate bracket positions separately for each genus
y_positions <- ps_genus_melt %>%
  group_by(Genus) %>%
  summarise(
    y_max = max(Abundance * 100, na.rm = TRUE),
    y_min = min(Abundance * 100, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = ifelse(y_max > y_min, y_max - y_min, abs(y_max) + 1)
  )

plot_dunn <- plot_dunn %>%
  left_join(y_positions, by = "Genus") %>%
  group_by(Genus) %>%
  arrange(match(group1, c("CK", "SG", "PA", "SV")),
          match(group2, c("CK", "SG", "PA", "SV"))) %>%
  mutate(
    y.position = y_max + row_number() * 0.12 * y_range
  ) %>%
  ungroup()


p1 <- ggplot(
  data = ps_genus_melt,
  aes(x = Group, y = Abundance * 100, color = Group)
) +
  geom_boxplot(
    aes(fill = Group),
    position = position_dodge(1),
    color = "black",
    linewidth = 0.5,
    outlier.shape = NA
  ) +
  geom_jitter(
    color = "#a6cee3",
    width = 0.15,
    size = 2.5,
    alpha = 0.8
  ) +
  stat_pvalue_manual(
    plot_dunn,
    label = "p.adj.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 3,
    hide.ns = TRUE
  ) +
  facet_wrap(~Genus, scales = "free_y", ncol = 3) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey"),
    strip.text = element_text(
      color = "blue",
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0,
      color = "black"
    )
  ) +
  labs(
    x = "",
    y = "Relative Abundance (%)",
    title = NULL
  ) +
  scale_fill_manual(
    values = c(
      "#a6cee3", "#bebada", "#e5c494",
      "#386cb0", "#1b9e77", "#fdc086"
    )
  ) +
  scale_color_manual(
    values = c(
      "#a6cee3", "#bebada", "#e5c494",
      "#386cb0", "#1b9e77", "#fdc086"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  coord_cartesian(clip = "off")

print(p1)

ggsave(
  "genus_relative_abundance_Dunn_BH.pdf",
  p1,
  width = 14,
  height = 16
)

ggsave(
  "genus_relative_abundance_Dunn_BH.tiff",
  p1,
  width = 14,
  height = 16,
  dpi = 600
)

# ggsave(file.path(path.plot,"07_diff/07_diff-genus-1.jpeg"),plot=p2, width = 14, height = 16)
#______________________________________________________________________________________________________________________________

