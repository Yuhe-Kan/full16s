# Fig4
# ROOT DIRECTORY
rm(list = ls())
getwd()

path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(coin)
library(MicrobiotaProcess)
library(ggplot2)
library(phyloseq)
library(microbiome)

# Load phyloseq object
ps <- readRDS(
  file.path(
    path.data,
    "phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"
  )
)

# Rarefy samples to an even sequencing depth
ps <- rarefy_even_depth(
  ps,
  sample.size = min(sample_sums(ps)),
  rngseed = 123,
  replace = FALSE,
  verbose = FALSE
)

# Hellinger transformation
ps <- microbiome::transform(ps, "hellinger")

# Set random seed
set.seed(1024)

# Differential analysis
deres <- MicrobiotaProcess::diff_analysis(
  obj = ps,
  classgroup = "Group",
  
  # LDA-based effect-size estimation
  mlfun = "lda",
  
  # Overall multi-group test
  firstcomfun = "kruskal.test",
  firstalpha = 0.05,
  
  # Multiple-testing adjustment
  padjust = "fdr",
  
  # Feature filtering
  filtermod = "pvalue",
  
  # Strict one-against-one comparisons
  strictmod = TRUE,
  secondcomfun = "wilcox.test",
  secondalpha = 0.05,
  
  # Subgroup settings
  subclmin = 3,
  subclwilc = TRUE,
  
  # LDA effect-size threshold
  ldascore = 3
)

# Inspect results
deres

taxa_order <- c(
  # p__;CK
  "o__Candidatus_Babeliales",
  "c__Candidatus_Babeliae",
  
  # Proteobacteria
  "o__Chromatiales",
  "g__Delftia",
  "f__Ectothiorhodospiraceae",
  "s__Herbaspirillum_huttiense",
  "g__Herbaspirillum",
  
  # Firmicutes
  "o__Lactobacillales",
  "f__Lactobacillaceae",
  "s__Pediococcus_pentosaceus",
  "g__Pediococcus",
  "g__Staphylococcus",
  "g__Streptococcus",
  "f__Streptococcaceae",
  "s__Streptococcus_parauberis",
  
  # Actinobacteria
  "c__Acidimicrobiia",
  "o__Acidimicrobiales",
  "f__Propionibacteriaceae",
  "s__Cutibacterium_acnes",
  "g__Cutibacterium",
  
  # Acidobacteria
  "o__Thermoanaerobaculales",
  "f__Thermoanaerobaculaceae",
  "c__Thermoanaerobaculia",
  "s__bacterium_3−5_YC−ZSS−LKJ−3024",
  "g__Subgroup_10",
  
  # PA; Proteobacteria
  "s__Methyloceanibacter_marginalis",
  "s__Roseomonas_vinacea",
  
  # SG
  # Proteobacteria
  "o__Oceanospirillales",
  "f__Halomonadaceae",
  "g__Kushneria",
  "g__Halomonas",
  "g__Paracoccus",
  "s__Lysobacter_bugurensis",
  # Bacteroidetes
  " f__Catalimonadaceae",
  
  #SA
  # Proteobacteria
  "g__Rubellimicrobium",
  "g__Azohydromonas",
  "s__Azohydromonas_australica",
  "g__Sorangium",
  " s__Geminicoccus_roseus",
  "g__Geminicoccus",
  "s__Sorangium_cellulosum",
  
  # Bacteroidetes
  "g__Flavisolibacter",
  "s__Flavisolibacter_ginsengisoli",
  
  # Actinobacteria
  "g__Arthrobacter",
  "s__Microbacterium_oleivorans"
)


taxa_order <- trimws(taxa_order)


es_p <- ggeffectsize(
  obj = deres,
  removeUnknown = TRUE,
  lineheight = 0.3,
  linewidth = 0.6,
  pointsize = 2.0
)


print(es_p$mapping)
print(names(es_p$data))


matched_cols <- names(es_p$data)[
  vapply(
    es_p$data,
    function(x) {
      any(trimws(as.character(x)) %in% taxa_order, na.rm = TRUE)
    },
    logical(1)
  )
]

if (length(matched_cols) != 1) {
  stop(
    paste0(
      "无法唯一识别分类名称列。候选列为：",
      paste(matched_cols, collapse = ", "),
      "。请检查 es_p$mapping 和 names(es_p$data)。"
    )
  )
}

taxa_col <- matched_cols[1]


taxa_values <- trimws(as.character(es_p$data[[taxa_col]]))
all_taxa <- unique(taxa_values[!is.na(taxa_values)])

ordered_taxa <- c(
  taxa_order[taxa_order %in% all_taxa],
  setdiff(all_taxa, taxa_order)
)


es_p$data[[taxa_col]] <- factor(
  taxa_values,
  levels = rev(ordered_taxa)
)


es_p <- es_p +
  scale_color_manual(
    values = c("#FDDED7", "#A28CC2", "#CCD376", "#5CB0C3")
  ) +
  theme(
    panel.background = element_rect(
      fill = "grey96",
      color = NA
    ),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(
      color = "white",
      linetype = "solid"
    ),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

print(es_p)


# Save figure
ggsave(
  filename = file.path(path.plot, "11_biomarker/es_p-1.jpeg"),
  plot = es_p,
  width = 8,
  height = 10
)

ggsave(
  filename = file.path(path.plot, "11_biomarker/es_p-1.pdf"),
  plot = es_p,
  width = 8,
  height = 10
)



