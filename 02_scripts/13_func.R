# Fig5A
rm(list = ls())
getwd()
path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

# install file2meco package (>= 0.7.1)
library(file2meco)
library(microeco)
library(magrittr)
library(ggplot2)
library(vegan)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(microbiomeMarker)
library(dplyr)
library(tibble)
library(phyloseq)

ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))

ps <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 123, replace = FALSE)
sample_sums(ps)


?microeco
?microtable

mt <- phyloseq2meco(ps)

#remove OTUs with the taxonomic assignments “mitochondria” or “chloroplast” and are not assigned in the Kingdom “k__Archaea” or “k__Bacteria”
mt$filter_pollution(taxa = c("mitochondria", "chloroplast"))
mt$tax_table %<>% base::subset(Kingdom == "k__Archaea" | Kingdom == "k__Bacteria")
#To make the OTU and sample information consistent across all files in the object
mt$tidy_dataset()

mt$sample_sums() %>% range


# create object of trans_func
t1 <- trans_func$new(mt)
# for fungi example, see https://chiliubio.github.io/microeco_tutorial/other-dataset.html#fungi-data
# default database for prokaryotes is FAPROTAX database

t1$cal_spe_func(prok_database = "FAPROTAX")# NJC19 or fungi_database = c("FUNGuild", "FungalTraits")[1],FUNGuild_confidence = c("Highly Probable", "Probable", "Possible")
t1$res_spe_func[1:5, 1:2]

# assignment coverage
# 5074/6870=73.86%

t1$cal_spe_func_perc(abundance_weighted = T)
t1$res_spe_func_perc[1:5, 1:2]

t1$func_group_list
class(t1$func_group_list)


fun_tab<-t(t1$res_spe_func_perc)
colnames(fun_tab)<-c("SG1","SG2","SG3","SG4","SG5","SG6","PA1","PA2","PA3","PA4","PA5","PA6",
                     "SV1","SV2","SV3","SV4","SV5","SV6","CK1","CK2","CK3","CK4","CK5","CK6")

# annotation_col
annotation_col = as.matrix.data.frame(sample_data(ps))[,c("Part","Group")]
rownames(annotation_col)<-c("SG1","SG2","SG3","SG4","SG5","SG6","PA1","PA2","PA3","PA4","PA5","PA6",
                     "SV1","SV2","SV3","SV4","SV5","SV6","CK1","CK2","CK3","CK4","CK5","CK6")
annotation_col<-as.data.frame(annotation_col)

# annotation_row
annotation_row<-t1$func_group_list$FAPROTAX
annotation_row <- data.frame(
  Annotation = rep(names(annotation_row), sapply(annotation_row, length)),
  Function = unlist(annotation_row)
)

# annotation_row <- annotation_row %>%
#   mutate_all(~ gsub("Energy source", "Energy_source", .))

index<-annotation_row[, 2] %in% rownames(fun_tab)
annotation_row<-annotation_row[index,]

filter <- !(rownames(fun_tab) %in% annotation_row[, 2])  
other_func <- rownames(fun_tab)[filter]  

other_data<-data.frame(row.names = paste("Other",1:length(other_func),sep = ""),
                        Annotation = rep("Other",length(other_func)),
                       Function = other_func)
annotation_row_new<-rbind(annotation_row,other_data)
rownames(annotation_row_new) <- NULL
annotation_row_new1<-column_to_rownames(annotation_row_new,var = "Function")

annotation_row_new1$Annotation[rownames(annotation_row_new1) == "chemoheterotrophy"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "phototrophy"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "oxygenic_photoautotrophy"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "photosynthetic_cyanobacteria"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "aromatic_compound_degradation"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "dark_hydrogen_oxidation"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "aerobic_anoxygenic_phototrophy"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "anoxygenic_photoautotrophy_H2_oxidizing"] <- "C-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "dark_thiosulfate_oxidation"] <- "S-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "nitrate_ammonification"] <- "N-cycle"
annotation_row_new1$Annotation[rownames(annotation_row_new1) == "nitrite_ammonification"] <- "N-cycle"

annotation_row_new1 <- annotation_row_new1 %>%
  mutate(Annotation = factor(Annotation, levels = c("C-cycle", "N-cycle", "S-cycle","Energy source","Other"))) %>%
  arrange(Annotation)

fun_tab <- fun_tab[rownames(annotation_row_new1), ]

ann_colors = list(
  Part = c(rhizosphere = "#1B9E77", bulk = "firebrick"),
  Group = c(CK = "#1B9E77", PA = "#D95F02",SG = "#E7298A", SV = "#66A61E"),
  Annotation = c("Energy source"="#7570B3", "C-cycle"="#E7298A", "S-cycle"="#66A61E","N-cycle"="#F46F43",Other="blue")
)

jpeg(file.path(path.plot, "09_func/faprotax/09_func_faprotax-1.jpeg"), width = 10, height = 10, units = "in", res = 300)
ComplexHeatmap::pheatmap(fun_tab, annotation_col = annotation_col[,c(2,1)], 
                         annotation_row = annotation_row_new1,
                         annotation_colors = ann_colors,
                         border_color = "gray90",
                         cluster_rows = F,
                         cluster_cols = F,
                         color =colorRampPalette(c("white", "pink", "red"))(100),
                         heatmap_legend_param = list(title_position = "topcenter",
                                                     legend_direction = "horizontal"),
                         column_split = annotation_col$Group
                         )
dev.off()
# ------------------------------------------------------------------------------------------------------------------------------------

# Fig5B
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
colnames(tmp)<-c("SG1","SG2","SG3","SG4","SG5","SG6","PA1","PA2","PA3","PA4","PA5","PA6",
                            "SV1","SV2","SV3","SV4","SV5","SV6","CK1","CK2","CK3","CK4","CK5","CK6")
View(tmp)

metadata_func<-annotation_row_new
rownames(metadata_func)<-metadata_func$Function
faprotax_mt <- microtable$new(otu_table = tmp, tax_table = metadata_func, sample_table = annotation_col)
faprotax_mt$tidy_dataset()


t2 <- trans_diff$new(dataset = faprotax_mt, 
                     method = "KW_dunn", 
                     group = "Group", 
                     taxa_level = "Function",
                     alpha = 0.05,
                     p_adjust_method = "bh",
                     filter_thres = 0,
                     transformation = NULL)
View(t2$res_diff)
View(t2$res_abund)

# top10 function
p2<-t2$plot_diff_abund(use_number = 1:10,
                       coord_flip = T,add_sig = T,
                       keep_prefix = F,
                       xtext_size = 15,
                       ytitle_size = 15,
                       order_x_mean = T)+
  ggplot2::theme_bw()+
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title.x = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 15),  # 调整图例文本大小
        legend.title = element_text(size = 15, face = "bold"))+
  ggplot2::ylab("Abundance");p2

View(t2$res_diff)
View(t2$res_abund)

# calculate relative abundances
diff_fun<-c("chemoheterotrophy", "phototrophy","photoautotrophy","oxygenic_photoautotrophy","photosynthetic_cyanobacteria",
            "anaerobic_chemoheterotrophy","fermentation","nitrogen_fixation","ureolysis",
            "photoheterotrophy",
            "human_associated","animal_parasites_or_symbionts",
            "denitrification","nitrate_denitrification","nitrite_denitrification","aerobic_anoxygenic_phototrophy",
            "chitinolysis","human_pathogens_all","iron_respiration","manganese_oxidation")

t3 <- trans_diff$new(dataset = faprotax_mt, 
                     method = "KW_dunn", 
                     group = "Group", 
                     taxa_level = "Function",
                     alpha = 0.05,
                     p_adjust_method = "bh",
                     filter_thres = 0,
                     transformation = NULL)

p3<-t3$plot_diff_abund(use_number = 1:72,
                       coord_flip = T,add_sig = T,
                       select_taxa = diff_fun,
                       keep_prefix = F,
                       xtext_size = 15,
                       ytitle_size = 15,
                       order_x_mean = F)+
  ggplot2::theme_bw()+
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 15,face = "bold"),
        axis.title.x = element_text(size = 15,face = "bold"),
        legend.text = element_text(size = 15),  # 调整图例文本大小
        legend.title = element_text(size = 15, face = "bold"))+
  ggplot2::ylab("Abundance");p3

ggsave(file.path(path.plot,"09_func/faprotax/09_func_faprotax-3.jpeg"),plot=p3, width = 8, height = 12)

# http://www.ehbio.com/test/venn/#/
fapro_tab<-as.data.frame(t(t1$res_spe_func_perc))
fapro_mt <- microtable$new(otu_table = fapro_tab, tax_table = metadata_func, sample_table = t1$sample_table)
fapro_ps<-meco2phyloseq(fapro_mt)

for (i in c("CK", "SV", "SG", "PA")) {
  fapro_ps<-meco2phyloseq(fapro_mt)
  fapro_ps<-phyloseq::subset_samples(fapro_ps, Group == i)
  fapro_ps<- phyloseq::prune_taxa(taxa_sums(fapro_ps) > 0, fapro_ps)
  print(paste("Group:", i))
  print(ntaxa(fapro_ps))
}
# 结果
# [1] "Group: CK"
# [1] 57
# [1] "Group: SV"
# [1] 58
# [1] "Group: SG"
# [1] 49
# [1] "Group: PA"
# [1] 59

# -----------------------------------------------------------------------------

# Fig6
# bugbase
# ref:https://github.com/mestaki/qiime2-to-BugBase

conda env list
conda activate qiime2-amplicon-2024.5

wd=/mnt/e/b_myproject/kan_2024_full16s/03_data/analysis-individual/kan_2024/03_qiime2
classifier=/mnt/e/b_myproject/kan_2024_full16s/03_data/analysis-individual/CLUSTER/taxonomy/classifier_greengene
rawfq=/mnt/e/b_myproject/kan_2024_full16s/03_data/analysis-individual/kan_2024/00_raw_fastq
soft=/mnt/e/b_myproject/kan_2024_full16s/02_soft

qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $classifier/97_otus.fasta \
--output-path $classifier/97_otus.qza

qiime vsearch cluster-features-closed-reference \
--i-table $wd/otutab.qza \
--i-sequences $wd/repseq.qza \
--i-reference-sequences $classifier/97_otus.qza \
--p-perc-identity 0.97 \
--o-clustered-table $wd/gg_otutab.qza \
--o-clustered-sequences $wd/gg_seqs.qza \
--o-unmatched-sequences $wd/gg_unmatched_seqs

qiime feature-table summarize \
--i-table $wd/gg_otutab.qza \
--o-visualization $wd/gg_otutab.qzv


qiime feature-table rarefy \
--i-table $wd/gg_otutab.qza \
--p-sampling-depth 1094 \  
--o-rarefied-table $wd/gg_otutab_rarefy.qza

qiime tools export \
--input-path $wd/gg_otutab_rarefy.qza \
--output-path $wd/gg_otutab_exported

biom convert \
-i $wd/gg_otutab_exported/feature-table.biom \
-o $wd/gg_otutab.tsv \
--to-tsv
#in gg_otutab.tsv, remove # Constructed from biom file
# 在线运行，输入gg_otutab.tsv，ImageGP:https://www.bic.ac.cn/ImageGP/index.php/Home/Index/BugBase.html

# assignment coverage
# 420/1830=22.95%

library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggsci)

rm(list = ls())
getwd()

path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

bugbase_tab <- read.table(
  file.path(path.plot, "09_func/bugbase/predicted_phenotypes/predictions.txt"),
  row.names = 1, header = TRUE, sep = "\t", comment.char = ""
)

metadata <- read.table(
  file.path(path.data, "analysis-individual/kan_2024/01_Metadata/metadata-bugbase.txt"),
  row.names = 1, header = TRUE, sep = "\t", comment.char = ""
)


metadata <- metadata[rownames(bugbase_tab), , drop = FALSE]


metadata$Group <- factor(
  metadata$Group,
  levels = c("CK", "SG", "PA", "SV")
)


functions <- c(
  "Aerobic",
  "Anaerobic",
  "Contains_Mobile_Elements",
  "Facultatively_Anaerobic",
  "Forms_Biofilms",
  "Gram_Negative",
  "Gram_Positive",
  "Potentially_Pathogenic",
  "Stress_Tolerant"
)


symnum.args <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  symbols = c("****", "***", "**", "*", "ns")
)


path.out <- file.path(
  path.plot,
  "09_func/bugbase/plots"
)

dir.create(
  path.out,
  recursive = TRUE,
  showWarnings = FALSE
)


plot_list <- list()
dunn_list <- list()


for (fun in functions) {
  
  message("Processing: ", fun)
  
  
  if (!fun %in% colnames(bugbase_tab)) {
    warning("Column not found: ", fun)
    next
  }
  
  
  dat <- data.frame(
    Group = metadata$Group,
    Abundance = as.numeric(bugbase_tab[[fun]])
  )
  
  dat <- dat %>%
    filter(!is.na(Group), !is.na(Abundance))
  
  dat$Group <- droplevels(dat$Group)
  
  
  if (nlevels(dat$Group) < 2) {
    warning("Fewer than two groups for: ", fun)
    next
  }
  
  
  dunn_res <- dat %>%
    dunn_test(
      Abundance ~ Group,
      p.adjust.method = "BH"
    ) %>%
    mutate(
      p.adj.signif = as.character(
        symnum(
          p.adj,
          corr = FALSE,
          cutpoints = symnum.args$cutpoints,
          symbols = symnum.args$symbols
        )
      )
    )
  
  dunn_list[[fun]] <- dunn_res
  
  
  y_max <- max(dat$Abundance, na.rm = TRUE)
  y_min <- min(dat$Abundance, na.rm = TRUE)
  y_range <- y_max - y_min
  
  
  if (y_range == 0) {
    y_range <- max(abs(y_max), 0.1)
  }
  
  dunn_res$y.position <- seq(
    y_max + 0.10 * y_range,
    y_max + 0.40 * y_range,
    length.out = nrow(dunn_res)
  )
  
  
  p <- ggplot(
    dat,
    aes(x = Group, y = Abundance, color = Group)
  ) +
    geom_boxplot(
      aes(fill = Group),
      alpha = 0.6,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.12,
      size = 2
    ) +
    scale_color_manual(
      values = pal_npg("nrc")(9)
    ) +
    scale_fill_manual(
      values = pal_npg("nrc")(9)
    ) +
    theme_bw() +
    stat_pvalue_manual(
      dunn_res,
      label = "p.adj.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      hide.ns = FALSE
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.10))
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 1
      ),
      axis.title.x = element_blank()
    ) +
    labs(
      y = "Predicted phenotype abundance",
      title = fun
    )
  
  
  plot_list[[fun]] <- p
}


dunn_all <- bind_rows(
  lapply(names(dunn_list), function(fun) {
    dunn_list[[fun]] %>%
      mutate(Function = fun, .before = 1)
  })
)

library(patchwork)

combined_plot <- wrap_plots(
  plot_list[
    c(
      "Aerobic",
      "Anaerobic",
      "Contains_Mobile_Elements",
      "Facultatively_Anaerobic",
      "Forms_Biofilms",
      "Gram_Negative",
      "Gram_Positive",
      "Potentially_Pathogenic",
      "Stress_Tolerant"
    )
  ],
  ncol = 3
)


print(combined_plot)


ggsave(
  filename = file.path(
    path.plot,
    "09_func/bugbase/bugbase.pdf"
  ),
  plot = combined_plot,
  width = 12,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)


write.csv(
  dunn_all,
  file.path(path.out, "BugBase_Dunn_BH_all_functions.csv"),
  row.names = FALSE
)
# ============================================================

# Table3
library(dplyr)
library(tidyr)
library(rstatix)
library(multcompView)


functions <- c(
  "Aerobic",
  "Anaerobic",
  "Contains_Mobile_Elements",
  "Facultatively_Anaerobic",
  "Forms_Biofilms",
  "Gram_Negative",
  "Gram_Positive",
  "Potentially_Pathogenic",
  "Stress_Tolerant"
)

group_levels <- c("CK", "SG", "PA", "SV")

scale_factor <- 100


metadata <- metadata[rownames(bugbase_tab), , drop = FALSE]


bugbase_long <- bugbase_tab %>%
  tibble::rownames_to_column("SampleID") %>%
  left_join(
    data.frame(
      SampleID = rownames(metadata),
      Group = metadata$Group
    ),
    by = "SampleID"
  ) %>%
  mutate(
    Group = factor(Group, levels = group_levels)
  ) %>%
  pivot_longer(
    cols = all_of(functions),
    names_to = "Function",
    values_to = "Abundance"
  ) %>%
  filter(
    !is.na(Group),
    !is.na(Abundance)
  )



summary_tab <- bugbase_long %>%
  group_by(Function, Group) %>%
  summarise(
    Mean = mean(Abundance) * scale_factor,
    SD = sd(Abundance) * scale_factor,
    .groups = "drop"
  )


kw_tab <- bugbase_long %>%
  group_by(Function) %>%
  kruskal_test(Abundance ~ Group) %>%
  ungroup() %>%
  mutate(
    P_adj_overall = p.adjust(p, method = "BH")
  ) %>%
  select(Function, p, P_adj_overall)


dunn_tab <- bugbase_long %>%
  group_by(Function) %>%
  dunn_test(
    Abundance ~ Group,
    p.adjust.method = "BH"
  ) %>%
  ungroup()


letter_list <- lapply(functions, function(fun) {
  
  overall_p <- kw_tab$P_adj_overall[
    match(fun, kw_tab$Function)
  ]
  
  # Default: all groups share the same letter
  letter_result <- data.frame(
    Function = fun,
    Group = group_levels,
    Letter = "a"
  )
  

  if (length(overall_p) == 1 &&
      !is.na(overall_p) &&
      overall_p < 0.05) {
    
    current_dunn <- dunn_tab %>%
      filter(Function == fun)
    
    if (nrow(current_dunn) > 0) {
      
      pvals <- setNames(
        current_dunn$p.adj,
        paste(
          current_dunn$group1,
          current_dunn$group2,
          sep = "-"
        )
      )
      
      cld <- multcompView::multcompLetters(
        pvals,
        threshold = 0.05
      )$Letters
      
      letter_result$Letter <- unname(
        cld[letter_result$Group]
      )
      
      
      letter_result$Letter[
        is.na(letter_result$Letter)
      ] <- "a"
    }
  }
  
  letter_result
})

letter_tab <- bind_rows(letter_list)


table_values <- summary_tab %>%
  left_join(
    letter_tab,
    by = c("Function", "Group")
  ) %>%
  mutate(
    Value = sprintf(
      "%.2f \u00B1 %.2f%s",
      Mean,
      SD,
      Letter
    )
  ) %>%
  select(Function, Group, Value) %>%
  pivot_wider(
    names_from = Group,
    values_from = Value
  )

# Add adjusted overall P values
final_table <- table_values %>%
  left_join(
    kw_tab %>%
      select(Function, P_adj_overall),
    by = "Function"
  ) %>%
  mutate(
    P_adj_overall = sprintf("%.4f", P_adj_overall)
  ) %>%
  select(
    Function,
    CK,
    SG,
    PA,
    SV,
    P_adj_overall
  )


print(final_table, n = Inf)

write.csv(
  final_table,
  file.path(
    path.plot,
    "09_func/bugbase/BugBase_summary_table.csv"
  ),
  row.names = FALSE,
  na = ""
)

# -----------------------------------------------------------------------------------------------------------------------
