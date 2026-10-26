# compare differences the number of reads by group

library(dplyr)
library(tidyr)
library(rstatix)
library(multcompView)

dat<- read.table(file.path(path.data,"analysis-individual/kan_2024/04_other/read_count.txt"),
                 row.names= 1, header=T, sep="\t",comment.char="")
dat2 <- dat %>%
  mutate(
    Group = as.character(Group),
    Group = factor(Group, levels = c("CK", "SG", "PA", "SV")),
    count = as.numeric(count)
  ) %>%
  filter(!is.na(Group), !is.na(count))

print(table(dat2$Group))

desc <- dat2 %>%
  group_by(Group) %>%
  summarise(
    Mean = mean(count),
    SD = sd(count),
    .groups = "drop"
  ) %>%
  mutate(
    Result = sprintf("%.2f ± %.2f", Mean, SD)
  )

kw <- dat2 %>%
  kruskal_test(count ~ Group) %>%
  mutate(
    P_adj_overall = p.adjust(p, method = "BH")
  )

dunn <- dat2 %>%
  dunn_test(
    count ~ Group,
    p.adjust.method = "BH"
  )

print(dunn)


group_levels <- levels(droplevels(dat2$Group))

p_mat <- matrix(
  1,
  nrow = length(group_levels),
  ncol = length(group_levels),
  dimnames = list(group_levels, group_levels)
)

for (i in seq_len(nrow(dunn))) {
  g1 <- as.character(dunn$group1[i])
  g2 <- as.character(dunn$group2[i])
  
  p_mat[g1, g2] <- dunn$p.adj[i]
  p_mat[g2, g1] <- dunn$p.adj[i]
}

letters <- multcompView::multcompLetters(
  p_mat,
  compare = "<",
  threshold = 0.05
)$Letters

letter_df <- data.frame(
  Group = names(letters),
  Letter = unname(letters)
)


final_table <- desc %>%
  left_join(letter_df, by = "Group") %>%
  mutate(
    Result = paste0(Result, Letter)
  ) %>%
  select(Group, Result) %>%
  pivot_wider(
    names_from = Group,
    values_from = Result
  ) %>%
  mutate(
    P_adj_overall = kw$P_adj_overall[1]
  ) %>%
  select(any_of(c("CK", "SG", "PA", "SV")), P_adj_overall)

print(final_table)


write.csv(
  final_table,
  file.path(
    path.plot,
    "read_count_significance_letters.csv"
  ),
  row.names = FALSE,
  na = "",
  fileEncoding = "UTF-8"
)
#-----------------------------------------------------------------------------

# Table 2

rm(list = ls())

path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(phyloseq)
library(dplyr)
library(tidyr)
library(multcompView)

# read data
ps_net = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))

cat("Total reads:", sum(sample_sums(ps_net)), "\n")
cat("Total ASVs:", ntaxa(ps_net), "\n")

# Check group information
metadata <- data.frame(sample_data(ps_net))
print(table(metadata$Group))

group_order <- c("CK", "SG", "PA", "SV")

stopifnot(all(group_order %in% as.character(metadata$Group)))

sample_data(ps_net)$Group <- factor(
  as.character(sample_data(ps_net)$Group),
  levels = group_order
)

# calculate the relative abundance (%) of each sample and combine them to the phylum level.
ps_rel <- transform_sample_counts(
  ps_net,
  function(x) x / sum(x) * 100
)

ps_phylum <- tax_glom(
  ps_rel,
  taxrank = "Phylum",
  NArm = FALSE
)

# Organize one record for each sample and each phylum.
dat <- psmelt(ps_phylum) %>%
  mutate(
    Phylum = as.character(Phylum),
    Phylum = ifelse(
      is.na(Phylum) | Phylum == "",
      "Unclassified",
      Phylum
    ),
    Group = factor(as.character(Group), levels = group_order)
  ) %>%
  group_by(Sample, Group, Phylum) %>%
  summarise(
    RelAbundance = sum(Abundance),
    .groups = "drop"
  )

# The 11 fungal phyla with the highest average relative abundance were selected.
top_phyla <- dat %>%
  group_by(Phylum) %>%
  summarise(
    OverallMean = mean(RelAbundance),
    .groups = "drop"
  ) %>%
  arrange(desc(OverallMean)) %>%
  slice_head(n = 11) %>%
  pull(Phylum)

dat_top <- dat %>%
  filter(Phylum %in% top_phyla)

# Calculate the mean ± SD for each group.
summary_table <- dat_top %>%
  group_by(Phylum, Group) %>%
  summarise(
    Mean = mean(RelAbundance),
    SD = sd(RelAbundance),
    n = n(),
    .groups = "drop"
  )

# Perform a general Kruskal-Wallis test on each phylum.
kw_table <- dat_top %>%
  group_by(Phylum) %>%
  summarise(
    P_overall = kruskal.test(
      RelAbundance ~ Group
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    P_adj_overall = p.adjust(P_overall, method = "BH")
  )

# Pairwise Wilcoxon test
pairs <- combn(group_order, 2, simplify = FALSE)

pair_table <- bind_rows(lapply(top_phyla, function(ph) {
  d <- dat_top %>% filter(Phylum == ph)
  
  bind_rows(lapply(pairs, function(g) {
    x <- d$RelAbundance[d$Group == g[1]]
    y <- d$RelAbundance[d$Group == g[2]]
    
    data.frame(
      Phylum = ph,
      Comparison = paste(g, collapse = "-"),
      P_raw = wilcox.test(
        x, y,
        exact = FALSE
      )$p.value
    )
  }))
})) %>%
  mutate(
    P_adj = p.adjust(P_raw, method = "BH")
  )

# Significance letters are generated based on the corrected pairwise comparison results.
letter_table <- bind_rows(lapply(top_phyla, function(ph) {
  d <- pair_table %>% filter(Phylum == ph)
  
  p_values <- d$P_adj
  names(p_values) <- d$Comparison
  
  letters <- multcompView::multcompLetters(
    p_values,
    threshold = 0.05
  )$Letters
  
  data.frame(
    Phylum = ph,
    Group = names(letters),
    Letter = unname(letters)
  )
}))

final_table <- summary_table %>%
  left_join(letter_table, by = c("Phylum", "Group")) %>%
  mutate(
    Mean_SD = sprintf(
      "%.2f ± %.2f%s",
      Mean, SD, Letter
    )
  ) %>%
  select(Phylum, Group, Mean_SD) %>%
  pivot_wider(
    names_from = Group,
    values_from = Mean_SD
  ) %>%
  left_join(kw_table, by = "Phylum") %>%
  mutate(
    P_adj_overall = format.pval(
      P_adj_overall,
      digits = 3,
      eps = 0.001
    )
  ) %>%
  select(
    Phylum, CK, SG, PA, SV,
    P_adj_overall
  ) %>%
  arrange(
    match(Phylum, top_phyla)
  )

print(final_table, n = Inf)

write.csv(
  final_table,
  file.path(path.plot, "Table2_phylum_mean_SD_significance.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Save the pairwise comparison results
write.csv(
  pair_table,
  file.path(path.plot, "Table2_pairwise_tests.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("Table 2 and pairwise test results exported.\n")
# --------------------------------------------------------------------------------------------------

# statistics reads abundance
rm(list = ls())
getwd()
path.root <- "D:/test/kan_2024_full16s"
path.plot <- file.path(path.root, "05_plot")
path.data <- file.path(path.root, "03_data")

library(phyloseq)
# ps = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
# 
# # Statistical reads (excluding unassigned phylum)
# sum(sample_sums(ps))#109537
# sort(sample_sums(ps))
# ntaxa(ps) #View total ASV count


# statistics on reads
ps_net = base::readRDS(file.path(path.data,"phyloseq-objects/physeq-with-phylotree/phyloseq_network_kan2024.rds"))
sum(sample_sums(ps_net))#109537
sort(sample_sums(ps_net))
# ---------------------------------------------------------------------------------------------------------------------

ps_net_phylum <- tax_glom(ps_net,taxrank = rank_names(ps_net)[2]);ps_net_phylum
nm1<-as.data.frame(taxa_sums(ps_net_phylum));nm1
class(nm1)
nm2<-as.data.frame(phyloseq::tax_table(ps_net_phylum));nm2
class(nm2)
nm2$nm_ASV<-nm1[,1]
nm2

nm3<-as.data.frame(table(as.data.frame(phyloseq::tax_table(ps_net))$Phylum));nm3
nm3$Phylum<-nm3[,1];nm3

nm4<-merge.data.frame(nm2,nm3,by="Phylum")[,c(1,8,10)];nm4
class(nm4)

library(scales)
nm4[2:3] <- sapply(nm4[2:3], function (x)  percent(x / sum(x), accuracy = 0.01))
nm4

#______________________________________________________________________________________________________________________________
   
