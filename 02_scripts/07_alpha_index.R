# Table1
# Install once if needed
install.packages(c("dplyr", "tidyr", "FSA", "multcompView", "openxlsx"))

library(dplyr)
library(tidyr)
library(FSA)
library(multcompView)
library(openxlsx)


indices <- c(
  "Chao1", "ACE", "Observed", "Pielou",
  "Shannon", "Simpson", "InvSimpson", "Fisher"
)

group_order <- c("CK", "SG", "PA", "SV")

alpha$Group <- factor(
  as.character(alpha$Group),
  levels = group_order
)

# Ensure the output directory exists
dir.create(path.out, recursive = TRUE, showWarnings = FALSE)


get_letters <- function(p_adj, groups, alpha_level = 0.05) {

  # p_adj must be a named vector, e.g. "CK-SG", "CK-PA"
  # Values are adjusted pairwise P values.

  if (length(p_adj) == 0 || all(is.na(p_adj))) {
    return(setNames(rep("a", length(groups)), groups))
  }

  # multcompView requires comparisons named with group pairs
  cld <- multcompView::multcompLetters(
    p_adj,
    threshold = alpha_level
  )

  result <- setNames(rep(NA_character_, length(groups)), groups)
  result[names(cld$Letters)] <- cld$Letters

  # Fill any missing group labels conservatively
  result[is.na(result)] <- "a"

  result
}


desc_list <- list()
overall_list <- list()
pairwise_list <- list()
letter_list <- list()

for (idx in indices) {

  dat <- alpha %>%
    select(Group, Value = all_of(idx)) %>%
    filter(!is.na(Group), !is.na(Value)) %>%
    mutate(Group = factor(as.character(Group), levels = group_order))

  # Descriptive statistics
  desc <- dat %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      Mean = mean(Value),
      SD = sd(Value),
      Mean_SD = sprintf("%.2f ± %.2f", Mean, SD),
      .groups = "drop"
    ) %>%
    mutate(Index = idx)

  desc_list[[idx]] <- desc

  # Overall Kruskal-Wallis test
  kw <- kruskal.test(Value ~ Group, data = dat)

  overall_list[[idx]] <- data.frame(
    Index = idx,
    Statistic = unname(kw$statistic),
    df = unname(kw$parameter),
    P_value = kw$p.value
  )

  # Dunn's post-hoc pairwise comparisons
  dunn <- FSA::dunnTest(
    Value ~ Group,
    data = dat,
    method = "bh"
  )$res

  # FSA returns comparisons such as "CK - SG"
  # Normalize them to names accepted by multcompView
  dunn$Comparison <- gsub(
    "\\s*-\\s*", "-", dunn$Comparison
  )

  pvec <- dunn$P.adj
  names(pvec) <- dunn$Comparison

  # Store pairwise results
  pairwise_list[[idx]] <- data.frame(
    Index = idx,
    Comparison = dunn$Comparison,
    Z = dunn$Z,
    P_unadjusted = dunn$P.unadj,
    P_adjusted_BH = dunn$P.adj
  )

  # Compact letter display
  letters <- get_letters(
    p_adj = pvec,
    groups = group_order
  )

  letter_list[[idx]] <- data.frame(
    Index = idx,
    Group = names(letters),
    Letter = unname(letters)
  )
}


overall <- bind_rows(overall_list) %>%
  mutate(
    P_adj_BH = p.adjust(P_value, method = "BH")
  )

desc_all <- bind_rows(desc_list)
letters_all <- bind_rows(letter_list)
pairwise_all <- bind_rows(pairwise_list)

# Add letters to mean +/- SD
desc_final <- desc_all %>%
  left_join(letters_all, by = c("Index", "Group")) %>%
  mutate(
    Mean_SD_Letter = paste0(Mean_SD, Letter)
  )

# Wide publication table
table_mean <- desc_final %>%
  select(Index, Group, Mean_SD_Letter) %>%
  pivot_wider(
    names_from = Group,
    values_from = Mean_SD_Letter
  )

final_table <- table_mean %>%
  left_join(overall, by = "Index") %>%
  mutate(
    P_value = signif(P_value, 4),
    P_adj_BH = signif(P_adj_BH, 4)
  ) %>%
  select(
    Index, CK, SG, PA, SV,
    P_value, P_adj_BH
  )

print(final_table)


write.csv(
  final_table,
  file.path(path.out, "alpha_diversity_with_letters.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  desc_final,
  file.path(path.out, "alpha_diversity_descriptive_with_letters.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  pairwise_all,
  file.path(path.out, "alpha_diversity_Dunn_pairwise_BH.csv"),
  row.names = FALSE,
  na = ""
)

write.csv(
  overall,
  file.path(path.out, "alpha_diversity_overall_KW_BH.csv"),
  row.names = FALSE,
  na = ""
)

wb <- createWorkbook()

addWorksheet(wb, "Summary with letters")
writeData(wb, "Summary with letters", final_table)

addWorksheet(wb, "Descriptive statistics")
writeData(wb, "Descriptive statistics", desc_final)

addWorksheet(wb, "Dunn pairwise tests")
writeData(wb, "Dunn pairwise tests", pairwise_all)

addWorksheet(wb, "Overall tests")
writeData(wb, "Overall tests", overall)

saveWorkbook(
  wb,
  file.path(path.out, "alpha_diversity_with_letters.xlsx"),
  overwrite = TRUE
)

message("Finished. Output directory: ", path.out)

#------------------------------------------------------------

