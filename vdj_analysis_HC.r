###LOAD LIBRARIES
#%%
#load required libraries

library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggraph)
library(circlize)
library(tidyr)
library(readr)
library(vscDebugger)
library(purrr)
library(scales)
library(languageserver)
library(httpgd)

#read excel file into tibble
HC_master_df <- summary_tsv_all %>%
 filter(chain == "heavy")  

#define BR_Code list
BR_code <- c(
             "2405", "2443", "2778", "2929", "2989", "3094", "3851", "5873",
             "6527", "20045", "1079", "1215", "1299", "1400", "1455", "1468",
             "1470", "1507", "1515", "1533", "1551", "1602", "1763", "1767",
             "1783", "1792", "1823", "1839", "1924", "1956", "1963", "1982",
             "2035", "2133", "2140", "2149", "2161", "2162", "2177", "2180",
             "2227", "2265", "HD279", "UTSW001", "UTSW002", "UTSW003",
             "UTSW009", "UTSW013", "UTSW014", "UTSW015", "UTSW017",
             "UTSW018", "UTSW019", "UTSW020", "UTSW022", "UTSW023",
             "UTSW033", "227", "6082")

project <- "CYSLOOP"

new_plots_dir <- paste0("/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/ggplot_outputs_group/", project)
 if (!dir.exists(new_plots_dir)) {
   dir.create(new_plots_dir, recursive = TRUE)
 }
new_csv_dir <- paste0("/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/csv_outputs_group/", project)
 if (!dir.exists(new_csv_dir)) {
   dir.create(new_csv_dir, recursive = TRUE)
 }

plots_output_dir <- new_plots_dir
csv_output_dir <- new_csv_dir

experimental_group <- c("3094", "6527", "5873")
control_group <- c("1763", "2265", "3851")

exp_group <- "cysloop"
ctrl_group <- "control"

expected_genes_VHfam <- c(paste0("VH", 1:7))
expected_genes_JHfam <- c(paste0("JH", 1:6))
expected_genes_VJHpairs <- c(
  with(expand.grid(V = paste0("VH", 1:7), J = paste0("JH", 1:6)), paste(V, J, sep = ":")))

##FUNCTIONS
#filter for experimental and control groups
HC_master_df <- HC_master_df %>%
  mutate(group_ID = case_when(
    BR_code %in% experimental_group ~ exp_group,
    BR_code %in% control_group ~ ctrl_group,
    TRUE ~ NA_character_
  ), .after = 1)

HC_comparison_df <- HC_master_df %>%
  filter(group_ID %in% c(exp_group, ctrl_group)) %>%
  mutate(VJ_pair = paste0(v_call, ":", j_call)) %>%
  mutate(
    V_gene_mut = str_extract(v_call, "IGHV(\\d+)") %>%
      str_replace("IGHV", "VH"),
    J_gene_mut = str_extract(j_call, "IGHJ(\\d+)") %>%
      str_replace("IGHJ", "JH"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)) %>%
  relocate(V_gene_mut, J_gene_mut, VJ_only_gene, VJ_pair, .before = v_call) %>%
  select(-LC_isotype)


#summary_df is a grouped df so all df made from it need to be grouped 
#ryan said this could be condensed using case_when, can come back to this
HC_V_summary_df <- HC_comparison_df %>%
    group_by(group_ID, BR_code, V_gene_mut) %>%
    summarize(
      hit_count = n(),
      HC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE))%>%
    rename(gene = V_gene_mut) %>% 
   tidyr::complete(gene = expected_genes_VHfam,
      fill = list(
      hit_count = 0, HC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
      explicit = FALSE) %>%  
    mutate(type = "V_gene") %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits)

HC_J_summary_df <- HC_comparison_df %>%
    group_by(group_ID, BR_code, J_gene_mut) %>%
    summarize(
      hit_count = n(),
      HC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE))%>% 
    rename(gene = J_gene_mut) %>% 
    tidyr::complete(gene = expected_genes_JHfam,
      fill = list(
      hit_count = 0, HC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
      explicit = FALSE) %>% 
    mutate(type = "J_gene") %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits)

HC_VJ_summary_df <- HC_comparison_df %>%
    group_by(group_ID, BR_code, VJ_only_gene) %>%
    summarize(
      hit_count = n(),
      HC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    rename(gene = VJ_only_gene) %>%
    tidyr::complete(gene = expected_genes_VJHpairs,
      fill = list(
      hit_count = 0, HC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
      explicit = FALSE) %>%  
    mutate(type = "VJ_pair") %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits)

HC_summary_df <- bind_rows(
  HC_V_summary_df,
  HC_J_summary_df,
  HC_VJ_summary_df
)

HC_gene_means_df <- HC_summary_df %>%
    group_by(group_ID, gene, type) %>%
    summarize(
      hit_count_gene = sum(hit_count),
      HC_cdr3_aa_charge_gene = mean(HC_cdr3_aa_charge)) %>%
    ungroup() %>%
    group_by(group_ID, type) %>%
    mutate(
      total_hits_family = sum(hit_count_gene),
      percent_gene = if_else(total_hits_family == 0, 0, ((hit_count_gene / total_hits_family) * 100)),
      HC_cdr3_aa_charge_family = mean(HC_cdr3_aa_charge_gene)) 


###HEAVY CHAIN PLOTS

# Variable kappa and lambda Family Distribution Plots 

#counts plot is irrelevent - not normalized
HC_V_family_distribution_counts_plot <- HC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = hit_count_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = HC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = hit_count, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "# of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_counts_plot, width = 6, height = 4)

HC_V_family_distribution_percent_plot <- HC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = HC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "% of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_percent_plot, width = 6, height = 4)

    
#CDR3 charge distribution plot for V gene families 
HC_V_family_distribution_charge_plot <- HC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = HC_cdr3_aa_charge_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = HC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = HC_cdr3_aa_charge, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variable Gene Family CDR3 Charges",
       x = "Heavy Chain Gene Family", y = "Avg CDR3 Charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_charge_plot, width = 6, height = 4)


# Joint kappa and lambda Family Distribution Plots
#count plots are irrelevant
HC_J_family_distribution_counts_plot <- HC_summary_df %>%
  filter(type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = hit_count_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = HC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = hit_count, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "# of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_counts_plot, width = 6, height = 4)

HC_J_family_distribution_percent_plot <- HC_summary_df %>%
  filter(type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = HC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "% of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_percent_plot, width = 6, height = 4)



#CDR3 charge distribution plot for J gene families
HC_J_family_distribution_charge_plot <- HC_summary_df %>%
  filter(type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = HC_cdr3_aa_charge_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = HC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = HC_cdr3_aa_charge, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Gene Family CDR3 Charges",
       x = "Heavy Chain Gene Family", y = "Avg CDR3 Charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_charge_plot, width = 6, height = 4)


#V:J Gene Pairings counts and percentage plots

#V:J pairs top 10 counts 
VJ_pairs_HC_top10_count_plot <- HC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair") %>%
  pivot_longer(cols = c("hit_count_gene"),
               names_to = "measure", values_to = "hit_count_gene") %>%
    arrange(desc(hit_count_gene)) %>%
    slice_head(n = 10) %>%
    filter(hit_count_gene > 0) %>%
  ggplot(aes(x = gene, y = hit_count_gene, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top 10 V:J pairs Heavy Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("HC_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_HC_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage
VJ_pairs_HC_top10_percent_plot <- HC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair") %>%
  pivot_longer(cols = c("percent_gene"),
               names_to = "measure", values_to = "percent_gene") %>%
    arrange(desc(percent_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
ggplot(aes(x = gene, y = percent_gene, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 V:J pairs Heavy Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("HC_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_HC_top10_percent_plot, width = 8, height = 4)



#CDR3 Charge  
HC_charge_plot <- HC_gene_means_df %>%
  dplyr::filter(type == "V_gene" | type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "V_gene" | type == "J_gene") %>%
    group_by(group_ID, type) %>%
    distinct(HC_cdr3_aa_charge_family),
    mapping = aes(x = type, y = HC_cdr3_aa_charge_family, fill = group_ID), 
    position = "dodge") +
  geom_point(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "V_gene" | type == "J_gene") %>%
    group_by(group_ID, type),
    mapping = aes(x = type, y = HC_cdr3_aa_charge_gene, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #geom_text(data = summary_df %>%
              #dplyr::filter(type == "V_gene" | type == "J_gene") %>%
              #group_by(group_ID, ID, type) %>%
              #summarize(hit_count = sum(hit_count), .groups = "drop") %>%
              #mutate(measure = "HC_cdr3_aa_charge"),
            #aes(x = ID, y = 0, label = hit_count, fill = group_ID)) +
  #facet_grid(cols = vars(type), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heavy Chain CDR3 Charges",
       x = "Heavy Chain Gene Family", y = "Avg CDR3 charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_charge_plot, width = 6, height = 4)

#CDR3 lengths 

#V:J pairs set chord colors
row.col = adjustcolor(c(
  VH1 = "#0356BA", VH2 = "#008B90", VH3 = "#16B900", VH4 = "#905100", VH5 = "#b6033cff", VH6 = "#6f02b3ff",
  VH7 = "#4603a5cb"), alpha.f = 0.5)

grid.col = c(
  VH1 = "#0356BA", VH2 = "#008B90", VH3 = "#16B900", VH4 = "#905100", VH5 = "#F4004F", VH6 = "#9F00FF",
  VH7 = "#5800d4cb", 
  JH1 = "#0356BA", JH2 = "#008B90", JH3 = "#16B900", JH4 = "#905100", JH5 = "#F4004F", JH6 = "#9F00FF")
  

#V:J pairs chord diagram 
VJ_pairs_HC_chord_diagram_exp <- HC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == exp_group) %>%
  separate(gene, into = c("VH", "JH", sep = ":")) %>%
  group_by(VH, JH) %>%
  summarise(percent = percent_gene, .groups = 'drop') %>%
  pivot_wider(names_from = "JH", values_from = percent, values_fill = 0) 
  matrix_exp <- as.matrix(VJ_pairs_HC_chord_diagram_exp[, -1])
  matrix_exp[matrix_exp == 0] <- 1e-6
  rownames(matrix_exp) <- VJ_pairs_HC_chord_diagram_exp$VH
  circos.clear()
  circos.par(gap.after = c(rep(5, nrow(matrix_exp)-1), 15, rep(5, ncol(matrix_exp)-1), 15))
  png(filename = file.path(plots_output_dir, paste0("HC_VJ_pairs_chord_diagram_", exp_group, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    matrix_exp, link.sort = TRUE, row.col = row.col,
    grid.col = grid.col, reduce = 0, transparency = 0.2,
    annotationTrack = "grid", link.visible = TRUE, preAllocateTracks = list(track.height = 0.15))
  circos.trackPlotRegion(
    track.index = 1, panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(
        x = mean(get.cell.meta.data("xlim")),
        y = mean(get.cell.meta.data("ylim")),
        labels = sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0.5),
        cex = 1.3)
    },
    bg.border = NA)
  title(paste0("Heavy Chain V:J Gene Pairings in ", exp_group), cex.main = 2)
  dev.off()

VJ_pairs_HC_chord_diagram_control <- HC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == ctrl_group) %>%
  separate(gene, into = c("VH", "JH", sep = ":")) %>%
  group_by(VH, JH) %>%
  summarise(percent = percent_gene, .groups = 'drop') %>%
  pivot_wider(names_from = "JH", values_from = percent, values_fill = 0) 
  matrix_exp <- as.matrix(VJ_pairs_HC_chord_diagram_control[, -1])
    matrix_exp[matrix_exp == 0] <- 1e-6
  rownames(matrix_exp) <- VJ_pairs_HC_chord_diagram_control$VH
  circos.clear()
  circos.par(gap.after = c(rep(5, nrow(matrix_exp)-1), 15, rep(5, ncol(matrix_exp)-1), 15))
  png(filename = file.path(plots_output_dir, paste0("HC_VJ_pairs_chord_diagram_", ctrl_group, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    matrix_exp, link.sort = TRUE, row.col = row.col,
    grid.col = grid.col, reduce = 0, transparency = 0.2,
    annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))
  circos.trackPlotRegion(
    track.index = 1, panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(
        x = mean(get.cell.meta.data("xlim")),
        y = mean(get.cell.meta.data("ylim")),
        labels = sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0.5, 0.5),
        cex = 1.3)
    },
    bg.border = NA)
  title(paste0("Heavy Chain V:J Gene Pairings in ", ctrl_group), cex.main = 2)
  dev.off()

##CSV output

HC_export_individual_df <- HC_comparison_df %>%

HC_export_group_df <- HC_summary_df %>%


#write_csv(HC_export_individual_df, file.path(csv_output_dir, paste0(project, "_LightChain_VJ_individualdata_", exp_group, "_vs_", ctrl_group, ".csv")))
#write_csv(HC_export_group_df, file.path(csv_output_dir, paste0(project, "_LightChain_VJ_groupdata_", exp_group, "_vs_", ctrl_group, ".csv")))

write_csv(HC_summary_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_VJ_summarydata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_gene_means_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_VJ_meandata_", exp_group, "_vs_", ctrl_group, ".csv")))
#write_csv(HC_kappalambda_ratio_df, file.path(csv_output_dir, paste0(project, "_LightChain_kappalambda_ratios_", exp_group, "_vs_", ctrl_group, ".csv")))