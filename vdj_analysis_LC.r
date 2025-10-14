###LOAD LIBRARIES
#%%
#load required libraries
###LOAD LIBRARIES
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
LC_master_df <- summary_tsv_all %>%
 filter(chain == "Light")  

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
             "UTSW033", "227")

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


##FUNCTIONS
#filter for experimental and control groups
LC_master_df <- LC_master_df %>%
  mutate(group_ID = case_when(
    BR_code %in% experimental_group ~ exp_group,
    BR_code %in% control_group ~ ctrl_group,
    TRUE ~ NA_character_
  ), .after = 1)

LC_comparison_df <- LC_master_df %>%
  filter(group_ID %in% c(exp_group, ctrl_group)) 


#Count V:J gene pairings function
LC_comparison_df <- LC_comparison_df %>%
  mutate(VJ_pair = paste0(v_call, ":", j_call)) %>%
  mutate(
    V_gene_mut = str_extract(v_call, "IGKV(\\d+)|IGLV(\\d+)") %>%
      str_replace("IGKV", "VK") %>%
      str_replace("IGLV", "VL"),
    J_gene_mut = str_extract(j_call, "IGKJ(\\d+)|IGLJ(\\d+)") %>%
      str_replace("IGKJ", "JK") %>%
      str_replace("IGLJ", "JL"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)) %>%
  relocate(V_gene_mut, J_gene_mut, VJ_only_gene, VJ_pair, .before = v_call)


#summary_df is a grouped df so all df made from it need to be grouped 
#ryan said this could be condensed using case_when, can come back to this
#add group_by BR_code - should be a simple addition 
LC_summary_df <- bind_rows(
  LC_comparison_df %>%
    group_by(group_ID, BR_code, V_gene_mut) %>%
    summarize(
      hit_count = n(),
      LC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    mutate(
      percent_value = hit_count / sum(hit_count),
      percent = percent_value * 100) %>%  
    rename(gene = V_gene_mut) %>%
    mutate(
      type = "V_gene", 
      ID = if_else(str_detect(gene, "VK"), "kappa", "lambda")),
                          
  LC_comparison_df %>%
    group_by(group_ID, BR_code, J_gene_mut) %>%
    summarize(
      hit_count = n(),
      LC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE))%>%
    mutate(
      percent_value = hit_count / sum(hit_count),
      percent = percent_value * 100) %>%
    rename(gene = J_gene_mut) %>%
    mutate(
      type = "J_gene", 
      ID = if_else(str_detect(gene, "JK"), "kappa", "lambda")),

  LC_comparison_df %>%
    group_by(group_ID, BR_code, VJ_only_gene) %>%
    summarize(
      hit_count = n(),
      LC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    mutate(
      percent_value = hit_count / sum(hit_count),
      percent = percent_value * 100) %>%
    rename(gene = VJ_only_gene) %>%
    mutate(
      type = "VJ_pair", 
      ID = if_else(str_detect(gene, "VK"), "kappa", "lambda")))
    

###LIGHT CHAIN PLOTS

# Variable kappa and lambda Family Distribution Plots
LC_V_family_distribution_counts_plot <- LC_summary_df %>%
  dplyr::filter(type == "V_gene") %>%
  pivot_longer(cols = c("hit_count"),
               names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = gene, y = value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge2(0.8, preserve = "single")) +
  facet_grid(rows = vars(measure), cols = vars(ID), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Light Chain Gene Family", y = "# of Hits", fill = "Group")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_V_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_V_family_distribution_counts_plot, width = 6, height = 4)

LC_V_family_distribution_percent_plot <- LC_summary_df %>%
  dplyr::filter(type == "V_gene") %>%
  pivot_longer(cols = c("percent"),
               names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = gene, y = value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge(0.8, preserve = "single")) +
  facet_grid(rows = vars(measure), cols = vars(ID), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Light Chain Gene Family", y = "% of Hits", fill = "Group") 
    ggsave(filename = file.path(plots_output_dir, paste0("LC_V_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_V_family_distribution_percent_plot, width = 6, height = 4)
    
#CDR3 charge distribution plot for V gene families 
LC_V_family_distribution_charge_plot <- LC_summary_df %>%
  dplyr::filter(type == "V_gene") %>%
  pivot_longer(cols = c("LC_cdr3_aa_charge"),
               names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = gene, y = value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge2(0.8, preserve = "single")) +
  facet_grid(rows = vars(measure), cols = vars(ID), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variable Gene Family CDR3 Charges",
       x = "Light Chain Gene Family", y = "Avg CDR3 charge", fill = "Group")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_V_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_V_family_distribution_charge_plot, width = 6, height = 4)

# Joint kappa and lambda Family Distribution Plots
LC_J_family_distribution_counts_plot <- LC_summary_df %>%
  dplyr::filter(type == "J_gene") %>%
  pivot_longer(cols = c("hit_count"),
               names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = gene, y = value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge(0.8, preserve = "single")) +
  facet_grid(rows = vars(measure), cols = vars(ID), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Light Chain Gene Family", y = "# of Hits", fill = "Group") 
    ggsave(filename = file.path(plots_output_dir, paste0("LC_J_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_J_family_distribution_counts_plot, width = 6, height = 4)

LC_J_family_distribution_percent_plot <- LC_summary_df %>%
  dplyr::filter(type == "J_gene") %>%
  pivot_longer(cols = c("percent"),
               names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = gene, y = value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge(0.8, preserve = "single")) +
  facet_grid(rows = vars(measure), cols = vars(ID), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Light Chain Gene Family", y = "% of Hits", fill = "Group") 
    ggsave(filename = file.path(plots_output_dir, paste0("LC_J_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_J_family_distribution_percent_plot, width = 6, height = 4)

#CDR3 charge distribution plot for J gene families
LC_J_family_distribution_charge_plot <- LC_summary_df %>%
  dplyr::filter(type == "J_gene") %>%
  pivot_longer(cols = c("LC_cdr3_aa_charge"),
               names_to = "measure", values_to = "value") %>%
  ggplot(aes(x = gene, y = value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge2(0.8, preserve = "single")) +
  facet_grid(rows = vars(measure), cols = vars(ID), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Gene Family CDR3 Charges",
       x = "Light Chain Gene Family", y = "Avg CDR3 charge", fill = "Group")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_J_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_J_family_distribution_charge_plot, width = 6, height = 4)


#V:J Gene Pairings counts, heatmap, and percentage plots

#V:J pairs top 10 counts kappa and lambda
VJ_pairs_LC_top10_count_plot <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair") %>%
  pivot_longer(cols = c("hit_count"),
               names_to = "measure", values_to = "hit_count") %>%
    arrange(desc(hit_count)) %>%
    slice_head(n = 10) %>%
  ggplot(aes(x = gene, y = hit_count, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top 10 V:J pairs Light Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("LC_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_LC_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage kappa and lambda
VJ_pairs_LC_top10_percent_plot <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair") %>%
  pivot_longer(cols = c("percent"),
               names_to = "measure", values_to = "percent") %>%
    arrange(desc(percent)) %>%
    slice_head(n = 10) %>%
ggplot(aes(x = gene, y = percent, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 V:J pairs Light Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_LC_top10_percent_plot, width = 8, height = 4)

#V:J pairs top 10 counts kappa
VJ_pairs_LC_kappa_top10_count_plot <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair", ID == "kappa") %>%
  pivot_longer(cols = c("hit_count"),
               names_to = "measure", values_to = "hit_count") %>%
    arrange(desc(hit_count)) %>%
    slice_head(n = 10) %>%
  ggplot(aes(x = gene, y = hit_count, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top 10 V:J pairs Kappa Light Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("LC_kappa_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_LC_kappa_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage kappa
VJ_pairs_LC_kappa_top10_percent_plot <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair", ID == "kappa") %>%
  pivot_longer(cols = c("percent"),
               names_to = "measure", values_to = "percent") %>%
    arrange(desc(percent)) %>%
    slice_head(n = 10) %>%
ggplot(aes(x = gene, y = percent, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 V:J pairs Kappa Light Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_kappa_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_LC_kappa_top10_percent_plot, width = 8, height = 4)

#V:J pairs top 10 counts lambda
VJ_pairs_LC_lambda_top10_count_plot <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair", ID == "lambda") %>%
  pivot_longer(cols = c("hit_count"),
               names_to = "measure", values_to = "hit_count") %>%
    arrange(desc(hit_count)) %>%
    slice_head(n = 10) %>%
  ggplot(aes(x = gene, y = hit_count, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top 10 V:J pairs Lambda Light Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("LC_lambda_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_LC_lambda_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage lambda
VJ_pairs_LC_lambda_top10_percent_plot <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair", ID == "lambda") %>%
  pivot_longer(cols = c("percent"),
               names_to = "measure", values_to = "percent") %>%
    arrange(desc(percent)) %>%
    slice_head(n = 10) %>%
ggplot(aes(x = gene, y = percent, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 V:J pairs Lambda Light Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_lambda_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_LC_lambda_top10_percent_plot, width = 8, height = 4)


#V:J pairs chord diagram 
VJ_pairs_LC_chord_diagram_exp <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == exp_group) %>%
  separate(gene, into = c("VK", "JK", sep = ":")) %>%
  group_by(VK, JK) %>%
  summarise(percent = sum(percent), .groups = 'drop') %>%
  pivot_wider(names_from = "JK", values_from = percent, values_fill = 0) 
  matrix_exp <- as.matrix(VJ_pairs_LC_chord_diagram_exp[, -1])
  rownames(matrix_exp) <- VJ_pairs_LC_chord_diagram_exp$VK
  circos.clear()
  circos.par(gap.after = c(rep(5, nrow(matrix_exp)-1), 15, rep(5, ncol(matrix_exp)-1), 15))
  png(filename = file.path(plots_output_dir, paste0("LC_VJ_pairs_chord_diagram_", exp_group, ".png")), 
      width = 600, height = 600)
  chordDiagram(matrix_exp, link.sort = TRUE, row.col = 1:13, transparency = 0.5)
  title(paste0("Light Chain V:J Gene Pairings in ", exp_group))
  dev.off()

VJ_pairs_LC_chord_diagram_control <- LC_summary_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == ctrl_group) %>%
  separate(gene, into = c("VK", "JK", sep = ":")) %>%
  group_by(VK, JK) %>%
  summarise(percent = sum(percent), .groups = 'drop') %>%
  pivot_wider(names_from = "JK", values_from = percent, values_fill = 0) 
  matrix_control <- as.matrix(VJ_pairs_LC_chord_diagram_control[, -1])
  rownames(matrix_control) <- VJ_pairs_LC_chord_diagram_control$VK
  circos.clear()
  circos.par(gap.after = c(rep(5, nrow(matrix_control)-1), 15, rep(5, ncol(matrix_control)-1), 15))
  png(filename = file.path(plots_output_dir, paste0("LC_VJ_pairs_chord_diagram_", ctrl_group, ".png")), 
      width = 600, height = 600)
  chordDiagram(matrix_control, link.sort = TRUE, row.col = 1:13, transparency = 0.5)
  title(paste0("Light Chain V:J Gene Pairings in ", ctrl_group))
  dev.off()

#CDR3 Charge Kappa vs Lambda plot - need to add counts and values - this is fucked up 
LC_Kappa_and_Lambda_charge_plot <- LC_summary_df %>%
  dplyr::filter(type == "V_gene" | type == "J_gene") %>%
  pivot_longer(cols = c("LC_cdr3_aa_charge"),
               names_to = "measure", values_to = "value") %>%
  group_by(group_ID, ID, type, measure) %>%
  summarize(mean_value = mean(value, na.rm = TRUE)) %>%
  ggplot(aes(x = ID, y = mean_value, fill = group_ID)) +
  geom_col(width = 0.7, position = position_dodge2(0.8, preserve = "single")) +
  geom_text(data = summary_df %>%
              dplyr::filter(type == "V_gene" | type == "J_gene") %>%
              group_by(group_ID, ID, type) %>%
              summarize(hit_count = sum(hit_count), .groups = "drop") %>%
              mutate(measure = "LC_cdr3_aa_charge"),
            aes(x = ID, y = 0, label = hit_count, fill = group_ID)) +
  facet_grid(rows = vars(measure), cols = vars(type), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Light Chain CDR3 Charges",
       x = "Light Chain Gene Family", y = "Avg CDR3 charge", fill = "Group")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_Kappa_and_Lambda_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_Kappa_and_Lambda_charge_plot, width = 6, height = 4)

##CSV output

LC_export_individual_df <- LC_comparison_df %>%

LC_export_group_df <- LC_summary_df %>%


write_csv(LC_export_individual_df, file.path(csv_output_dir, paste0(project, "_LightChain_VJ_individualdata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_export_group_df, file.path(csv_output_dir, paste0(project, "_LightChain_VJ_groupdata_", exp_group, "_vs_", ctrl_group, ".csv")))
