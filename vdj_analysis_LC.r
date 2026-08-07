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
#library(vscDebugger)
library(purrr)
library(scales)
library(languageserver)
#library(httpgd)
library(forcats)
library(ggbeeswarm)
library(gtools)

#read excel file into df
#filter for light chain sequences
LC_master_df <- summary_tsv_all %>%
 filter(chain == "light")  

#define BR_Code list
BR_code <- c(
             "2405", "2443", "2778", "2929", "2989", "3094", "3851", "5873",
             "6527", "20045", "1079", "1215", "1299", "1400", "1455", "1468",
             "1470", "1507", "1515", "1533", "1551", "1602", "1763", "1767",
             "1783", "1792", "1823", "1839", "1924", "1956", "1963", "1982",
             "2035", "2133", "2140", "2149", "2161", "2162", "2177", "2180",
             "2227", "2265", "HD279", "UTSW0001", "UTSW0002", "UTSW0003",
             "UTSW0009", "UTSW0013", "UTSW0014", "UTSW0015", "UTSW0017",
             "UTSW0018", "UTSW0019", "UTSW0020", "UTSW0022", "UTSW0023",
             "UTSW0033", "6082",
             "1455", "2177",	"2180",	"991",	"1823",	"1123",	"1299",	"1490",	"1507",	
             "1844",	"2133",	"1574",	"3500", "J10",	"J130",	"J20",	"J201",	"J203",
             "J218",	"J220",	"J24",	"J26",	"J34",	"J42",	"J46",	"J8",	"J9",	"J91",
             "J93", "3895",	"1559",	"1561",	"1607",	"1618",	"1629",	"1647",	"1690",	
             "1702",	"1707",	"1713",	"1724",	"1790", "1831",	"1837",	"1848",	"1850",	
             "3843",	"1551",	"3889",	"32",	"3891",	"3880",	"3911",	"176",	"3864",	"3870",	
             "3873",	"1400",	"3933",	"3934", "3945",	"3961",	"3976",	"3892",	"3780",	"3782",	
             "3587",	"3796",	"3834",	"3814", "3836",	"3825",	"3813",	"3839",	"3811",	"3805",
             "3808",	"2899",	"3845",	"3846",	"3847",	"3848",	"3850",	"3853", "3061",	"3255",
             "3267",	"381",	"1136",	"1145",	"1226",	"1288",	"1292",	"1313",	"1397",	"1403",
             "1446",	"1524",	"1700", "2024",	"2218",	"1503",	"3892",	"1589",	"1606",	"1611",
             "1617",	"1636",	"1646",	"1650",	"1658",	"1692",	"1698",	"1717", "1722",	"1734",
             "1741",	"1747",	"1787",	"1796",	"1801",	"1825",	"1851",	"2666",	"2708",	"3300",	
             "2946",	"3108",	"547", "3287",	"3469",	"3806",	"3905",	"606",	"1173",	"1238",
             "1301",	"1468", "1579",	"1864",	"1922",	"1823",	"1574",	"USCHC017", "USCHC018",
             "USCHC019",	"USCHC020", "USCHC021",	"USCHC022",	"USCHC023",	"USCHC024",	"USCHC025",
             "USCHC027", "USCHC028",	"USCHC029",	"USCHC030",	"USCHC031",	"USCHC032",	"USCHC033",
             "USCHC034",	"USCHC035",	"USCHC036", "USCHC037",	"USCHC038",	"USCHC040",	"USCH001",
             "USCH002",	"USCH006",	"USCH004",	"USCH003",	"USCH005",	"USCH007",	"USCH010",	
             "USCH011",	"USCH013",	"USCH015",	"USCH014",	"USCH012",	"USCH016",	"USCH009",
             "USCH008",	"1503",	"USCHC001",	"USCHC002", "USCHC003",	"USCHC004",	"USCHC006",
             "USCHC007",	"USCHC008",	"USCHC009",	"USCHC010",	"USCHC011",	"USCHC012",	
             "USCHC013",	"USCHC014",	"USCHC015",	"USCHC016",	"USCH019",	"USCH018",	"USCH024",
             "USCH029",	"USCH030",	"USCH031", "USCH033",	"USCH037",	"USCH034",	"USCH017",
             "USCH040",	"USCH041",	"USCH044",	"USCH022",	"USCH023",	"USCH025",	"USCH026",
             	"USCH027",	"USCH035",	"USCH036",	"USCH039",	"USCH028",	"USCH042",
              "3895",	"1483",	"1616",	"1638",	"1697",	"3931", "USCHC026",	"USCHC039",
              "USCHC005",	"USCH043",	"USCH020",	"USCH021",	"USCH032"
)

project <- "CYSLOOP_AAG"

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


experimental_group <- c("6527")
control_group <- c("1763", "3851")

exp_group <- "cys:AAG"
ctrl_group <- "control"

expected_genes_VKfam <- c(paste0("VK", 1:6))
expected_genes_VLfam <- c(paste0("VL", 1:10))
expected_genes_JKfam <- c(paste0("JK", 1:5))
expected_genes_JLfam <- c(paste0("JL", 1:3))
expected_genes_VJKpairs <- c(
  with(expand.grid(V = paste0("VK", 1:6), J = paste0("JK", 1:5)), paste(V, J, sep = ":")))
expected_genes_VJLpairs <- c(  
  with(expand.grid(V = paste0("VL", 1:10), J = paste0("JL", 1:3)), paste(V, J, sep = ":")))


##FUNCTIONS
#filter for experimental and control groups
LC_master_df <- LC_master_df %>%
  mutate(group_ID = case_when(
    BR_code %in% experimental_group ~ exp_group,
    BR_code %in% control_group ~ ctrl_group,
    TRUE ~ NA_character_
  ), .after = 1)

LC_comparison_df <- LC_master_df %>%
  filter(group_ID %in% c(exp_group, ctrl_group)) %>%
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
LC_V_summary_df <- LC_comparison_df %>%
    group_by(group_ID, BR_code, V_gene_mut) %>%
    summarize(
      hit_count = n(),
      LC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE))%>%
    rename(gene = V_gene_mut) %>% 
    tidyr::complete(gene = expected_genes_VKfam,
      fill = list(
        hit_count = 0, LC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
        explicit = FALSE) %>% 
    tidyr::complete(gene = expected_genes_VLfam,
      fill = list(
        hit_count = 0, LC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
        explicit = FALSE) %>%   
    mutate(
     type = "V_gene",
     isotype = if_else(str_detect(gene, "VK"), "kappa", "lambda")) %>%  
    ungroup() %>%
    group_by(group_ID, BR_code, isotype) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene, isotype) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE),
          cdr3_stdev = sd(LC_cdr3_aa_charge))%>%
    ungroup() %>%
    group_by(gene) %>%
    mutate(
          percent_p_value = tryCatch(
            t.test(percent ~ group_ID)$p.value,
            error = function(e)NA_real_),
          cdr3_p_value = tryCatch(
            t.test(LC_cdr3_aa_charge ~ group_ID)$p.value,
            error = function(e)NA_real_))

LC_J_summary_df <- LC_comparison_df %>%
    group_by(group_ID, BR_code, J_gene_mut) %>%
    summarize(
      hit_count = n(),
      LC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE))%>% 
    rename(gene = J_gene_mut) %>% 
    tidyr::complete(gene = expected_genes_JKfam,
      fill = list(
        hit_count = 0, LC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
        explicit = FALSE) %>% 
    tidyr::complete(gene = expected_genes_JLfam,
      fill = list(
        hit_count = 0, LC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
        explicit = FALSE) %>%  
    mutate(
      type = "J_gene",
      isotype = if_else(str_detect(gene, "JK"), "kappa", "lambda")) %>%
    ungroup() %>%
    group_by(group_ID, BR_code, isotype) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene, isotype) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE),
          cdr3_stdev = sd(LC_cdr3_aa_charge))%>%
    ungroup() %>%
    group_by(gene) %>%
    mutate(  
          percent_p_value = tryCatch(
            t.test(percent ~ group_ID)$p.value,
            error = function(e)NA_real_),
          cdr3_p_value = tryCatch(
            t.test(LC_cdr3_aa_charge ~ group_ID)$p.value,
            error = function(e)NA_real_))

LC_VJ_summary_df <- LC_comparison_df %>%
    group_by(group_ID, BR_code, VJ_only_gene) %>%
    summarize(
      hit_count = n(),
      LC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    rename(gene = VJ_only_gene) %>%
    tidyr::complete(gene = expected_genes_VJKpairs,
      fill = list(
        hit_count = 0, LC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
        explicit = FALSE) %>% 
    tidyr::complete(gene = expected_genes_VJLpairs,
      fill = list(
        hit_count = 0, LC_cdr3_aa_charge = 0, percent_value = 0, percent = 0), 
        explicit = FALSE) %>%  
    mutate(
      type = "VJ_pair", 
      isotype = if_else(str_detect(gene, "VK"), "kappa", "lambda")) %>%
    ungroup() %>%
    group_by(group_ID, BR_code, isotype) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene, isotype) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE),
           cdr3_stdev = sd(LC_cdr3_aa_charge))%>%
    ungroup() %>%
    group_by(gene) %>%
    mutate(  
          percent_p_value = tryCatch(
            t.test(percent ~ group_ID)$p.value,
            error = function(e)NA_real_),
          cdr3_p_value = tryCatch(
            t.test(LC_cdr3_aa_charge ~ group_ID)$p.value,
            error = function(e)NA_real_))

LC_summary_df <- bind_rows(
  LC_V_summary_df,
  LC_J_summary_df,
  LC_VJ_summary_df
)

LC_gene_means_df <- LC_summary_df %>%
    group_by(group_ID, gene, type, isotype) %>%
    summarize(
      hit_count_gene = sum(hit_count),
      LC_cdr3_aa_charge_gene = mean(LC_cdr3_aa_charge),
      percent_stdev = unique(percent_stdev),
      percent_p_value = unique(percent_p_value),
      cdr3_stdev = unique(cdr3_stdev),
      cdr3_p_value = unique(cdr3_p_value)) %>%
    ungroup() %>%
    group_by(group_ID, type, isotype) %>%
    mutate(
      total_hits_isotype = sum(hit_count_gene),
      percent_gene = if_else(total_hits_isotype == 0, 0, ((hit_count_gene / total_hits_isotype) * 100)),
      LC_cdr3_aa_charge_isotype = mean(LC_cdr3_aa_charge_gene)) 

LC_cdr3_avgBR_charge <- LC_comparison_df %>%
    group_by(group_ID, BR_code, LC_isotype) %>%
    filter(chain == "light") %>%
    summarize(
      mean_br_cdr3charge = mean(cdr3_aa_charge)) %>%
    ungroup() %>%
    group_by(group_ID, LC_isotype) %>%
    mutate(group_cdr3charge = mean(mean_br_cdr3charge), 
           stdev = sd(mean_br_cdr3charge, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(LC_isotype) %>%
    mutate(
      p_value = tryCatch(
      t.test(mean_br_cdr3charge ~ group_ID)$p.value,
      error = function(e)NA_real_))

#kappa lambda ratio df - using v_gene, but the j_gene would be same number    
LC_kappalambda_ratio_df <- LC_summary_df %>%
    group_by(group_ID, BR_code) %>%
    filter(type == "V_gene", BR_code != 2180) %>%
    summarize(
      kappa_sum_br = sum(ifelse(isotype == "kappa", hit_count, 0)),
      lambda_sum_br = sum(ifelse(isotype == "lambda", hit_count, 0)),
      kappa_lambda_ratio_br = (kappa_sum_br / lambda_sum_br)) %>%
    ungroup() %>%
    group_by(group_ID) %>%
    mutate(
      kappa_sum_group = sum(kappa_sum_br),
      lambda_sum_group = sum(lambda_sum_br),
      kappa_lambda_ratio_group = (kappa_sum_group / lambda_sum_group))

#make wide df for presentations 
LC_summary_wide_df <- LC_summary_df %>%
    ungroup() %>%
    filter(type != "VJ_pair") %>%
    select(-isotype, -hit_count, -type, -LC_cdr3_aa_charge, -percent_stdev, 
            -percent_p_value, -cdr3_p_value, -cdr3_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent")

LC_means_wide_df <- LC_gene_means_df %>%
    ungroup() %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -LC_cdr3_aa_charge_gene, -total_hits_isotype, 
           -LC_cdr3_aa_charge_isotype, -cdr3_stdev, -percent_p_value, -cdr3_p_value, -percent_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent_gene")

LC_percent_stdev_wide_df <- LC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -LC_cdr3_aa_charge_gene, -total_hits_isotype,
           -LC_cdr3_aa_charge_isotype, -percent_gene, -percent_p_value, -cdr3_p_value, -cdr3_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent_stdev")

LC_percent_pvalue_wide_df <- LC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -LC_cdr3_aa_charge_gene, -total_hits_isotype,
           -LC_cdr3_aa_charge_isotype, -percent_gene, -percent_stdev, -cdr3_p_value, -cdr3_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent_p_value")

LC_cdr3_stdev_wide_df <- LC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -LC_cdr3_aa_charge_gene, -total_hits_isotype,
           -LC_cdr3_aa_charge_isotype, -percent_gene, -percent_p_value, -cdr3_p_value, -percent_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "cdr3_stdev")

LC_cdr3_pvalue_wide_df <- LC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -LC_cdr3_aa_charge_gene, -total_hits_isotype,
           -LC_cdr3_aa_charge_isotype, -percent_gene, -percent_stdev, -percent_stdev, 
           -cdr3_stdev, -percent_p_value) %>%
    pivot_wider(names_from = "gene", values_from = "cdr3_p_value")

LC_cell_counts_df <- LC_summary_df %>%
    group_by(group_ID, BR_code) %>%
    filter(type == "V_gene") %>%
    summarize(cells = sum(hit_count))


###LIGHT CHAIN PLOTS

# Variable kappa and lambda Family Distribution Plots 

LC_V_family_distribution_counts_plot <- LC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot() +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = hit_count_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = LC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = hit_count, shape = group_ID, group = group_ID), 
    size = 1.0, position = position_dodge(width = 0.8),
      shape = 21, color = "black", fill = "gray", alpha = 0.5) +
  facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Light Chain Gene Family", y = "# of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_V_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_V_family_distribution_counts_plot, width = 6, height = 4)

LC_V_family_distribution_percent_plot <- LC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot(
    data = LC_summary_df, mapping = aes(x = gene, y = percent, group = group_ID)
  ) +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    mutate(
      gene_numeric = as.numeric(str_extract(gene, "\\d+")),
      gene = fct_reorder(gene, gene_numeric)) %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), width = 0.8,
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_beeswarm(
    data = LC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    mutate(
      gene_numeric = as.numeric(str_extract(gene, "\\d+")),
      gene = fct_reorder(gene, gene_numeric)) %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
    dodge.width = 0.8, method = "center", 
    priority = "density", corral = "wrap", 
    corral.width = 0.15, alpha = 0.7
    #size = 1.3, position = position_jitterdodge(jitter.height = 0, jitter.width = 0.0)
    ) +
 #the text is too squished together 
  #geom_text(
   # data = LC_gene_means_df %>%
   ## dplyr::filter(type == "V_gene") %>%
    #group_by(group_ID) %>%
   # dplyr::distinct(gene, .keep_all = TRUE),
   # mapping = aes(x = gene, y = - 1, label = round(percent_gene, 2)),
   # vjust = 4, size = 1.8, colour = "black", 
   # position = position_dodge2(width = 1, preserve = "total")) +
  #coord_cartesian(ylim = c(-5, max(LC_gene_means_df$percent_gene))) +
  facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Light Chain Gene Family", y = "% of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_V_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_V_family_distribution_percent_plot, width = 15, height = 4)


# Joint kappa and lambda Family Distribution Plots
LC_J_family_distribution_counts_plot <- LC_summary_df %>%
  filter(type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = hit_count_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = LC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = hit_count, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Light Chain Gene Family", y = "# of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_J_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_J_family_distribution_counts_plot, width = 6, height = 4)

LC_J_family_distribution_percent_plot <- LC_summary_df %>%
  filter(type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_beeswarm(
    data = LC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
    dodge.width = 0.8, method = "center", 
    priority = "density", corral = "wrap", 
    corral.width = 0.10, alpha = 0.7
    #size = 1.3, position = position_jitterdodge(jitter.width = 0.2)
    ) +
  #geom_text(
  #  data = LC_gene_means_df %>%
  #  dplyr::filter(type == "J_gene") %>%
  #  group_by(group_ID) %>%
  #  dplyr::distinct(gene, .keep_all = TRUE),
  ##  mapping = aes(x = gene, y = - 1, label = round(percent_gene, 2)),
   # vjust = 4, size = 1.8, colour = "black", 
   # position = position_dodge2(width = 1, preserve = "total")) +
  #coord_cartesian(ylim = c(-5, max(LC_gene_means_df$percent_gene))) +
  facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Light Chain Gene Family", y = "% of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_J_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_J_family_distribution_percent_plot, width = 10, height = 4)

#V:J Gene Pairings counts, heatmap, and percentage plots

#V:J pairs top 10 counts kappa and lambda
VJ_pairs_LC_top10_count_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair") %>%
  pivot_longer(cols = c("hit_count_gene"),
               names_to = "measure", values_to = "hit_count_gene") %>%
    arrange(desc(hit_count_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
  ggplot(aes(x = gene, y = hit_count_gene, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    labs(title = "Top 10 V:J pairs Light Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("LC_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_LC_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage kappa and lambda
VJ_pairs_LC_top10_percent_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair") %>%
  group_by(group_ID) %>%
  pivot_longer(cols = c("percent_gene"),
               names_to = "measure", values_to = "percent_gene") %>%
    arrange(desc(percent_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
ggplot(aes(x = gene, y = percent_gene, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = -2, margin=margin(r=15)), 
        axis.title.x = element_text(hjust = 0.5, size = 14, margin=margin(t=20))) +
  labs(title = "Top 10 V:J pairs Light Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_LC_top10_percent_plot, width = 18, height = 4)

#V:J pairs top 10 counts kappa
VJ_pairs_LC_kappa_top10_count_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", isotype == "kappa") %>%
  pivot_longer(cols = c("hit_count_gene"),
               names_to = "measure", values_to = "hit_count_gene") %>%
    arrange(desc(hit_count_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
  ggplot(aes(x = gene, y = hit_count_gene, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    labs(title = "Top 10 V:J pairs Kappa Light Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("LC_kappa_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_LC_kappa_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage kappa
VJ_pairs_LC_kappa_top10_percent_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", isotype == "kappa") %>%
  pivot_longer(cols = c("percent_gene"),
               names_to = "measure", values_to = "percent_gene") %>%
    arrange(desc(percent_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
ggplot(aes(x = gene, y = percent_gene, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Top 10 V:J pairs Kappa Light Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_kappa_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_LC_kappa_top10_percent_plot, width = 8, height = 4)

#V:J pairs top 10 counts lambda
VJ_pairs_LC_lambda_top10_count_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", isotype == "lambda") %>%
  pivot_longer(cols = c("hit_count_gene"),
               names_to = "measure", values_to = "hit_count_gene") %>%
    arrange(desc(hit_count_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
  ggplot(aes(x = gene, y = hit_count_gene, fill = group_ID)) +
    geom_col(position = position_dodge(0.5)) +
    facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
    theme_grey(base_size = 14) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    labs(title = "Top 10 V:J pairs Lambda Light Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("LC_lambda_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_LC_lambda_top10_count_plot, width = 8, height = 4)

#V:J pairs percentage lambda
VJ_pairs_LC_lambda_top10_percent_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", isotype == "lambda") %>%
  pivot_longer(cols = c("percent_gene"),
               names_to = "measure", values_to = "percent_gene") %>%
    arrange(desc(percent_gene)) %>%
    slice_head(n = 10) %>%
    filter(percent_gene > 0) %>%
ggplot(aes(x = gene, y = percent_gene, fill = group_ID)) +
  geom_col(position = position_dodge(0.5)) +
  facet_grid(rows = vars(measure), cols = vars(group_ID), scales = "free") +
  theme_gray(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Top 10 V:J pairs Lambda Light Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_lambda_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_LC_lambda_top10_percent_plot, width = 8, height = 4)

#CDR3 Charge Kappa vs Lambda plot 
LC_Kappa_and_Lambda_charge_plot <- LC_gene_means_df %>%
  dplyr::filter(type == "V_gene" | type == "J_gene") %>%
  #group_by(group_ID, isotype, type) %>%
  ggplot() +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "V_gene" | type == "J_gene") %>%
    group_by(group_ID, isotype, type) %>%
    distinct(LC_cdr3_aa_charge_isotype),
    mapping = aes(x = isotype, y = LC_cdr3_aa_charge_isotype, fill = group_ID), 
    position = "dodge") +
  geom_point(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "V_gene" | type == "J_gene") %>%
    group_by(group_ID, isotype, type),
    mapping = aes(x = isotype, y = LC_cdr3_aa_charge_gene, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  #geom_text(data = summary_df %>%
              #dplyr::filter(type == "V_gene" | type == "J_gene") %>%
              #group_by(group_ID, ID, type) %>%
              #summarize(hit_count = sum(hit_count), .groups = "drop") %>%
              #mutate(measure = "LC_cdr3_aa_charge"),
            #aes(x = ID, y = 0, label = hit_count, fill = group_ID)) +
  facet_grid(cols = vars(type), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Light Chain CDR3 Charges",
       x = "Light Chain Gene Family", y = "Avg CDR3 charge")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_Kappa_and_Lambda_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_Kappa_and_Lambda_charge_plot, width = 6, height = 4)

#CDR3 CHARGE PER BR
LC_charge_plot_perBR <- LC_cdr3_avgBR_charge %>%
  ggplot() +
  geom_col(
    data = LC_cdr3_avgBR_charge %>%
    group_by(group_ID, LC_isotype) %>%
    distinct(group_cdr3charge),
    mapping = aes(x = group_ID, y = group_cdr3charge, fill = group_ID), 
    position = "dodge") +
  geom_beeswarm(
    data = LC_cdr3_avgBR_charge %>%
    group_by(group_ID, BR_code, LC_isotype),
    mapping = aes(x = group_ID, y = mean_br_cdr3charge, shape = group_ID, group = group_ID), 
        dodge.width = 0.8, 
        size = 1.2,
        method = "swarm",
        color = "black",
        corral = "gutter",
        corral.width = 0.8, 
        alpha = 0.7
    #size = 1.3, position = position_jitterdodge(jitter.width = 0.2
    ) +
  #geom_text(data = summary_df %>%
              #dplyr::filter(type == "V_gene" | type == "J_gene") %>%
              #group_by(group_ID, ID, type) %>%
              #summarize(hit_count = sum(hit_count), .groups = "drop") %>%
              #mutate(measure = "HC_cdr3_aa_charge"),
            #aes(x = ID, y = 0, label = hit_count, fill = group_ID)) +
  facet_grid(cols = vars(LC_isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Light Chain CDR3 Charges",
       x = "Group", y = "Avg CDR3 charge")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_CDR3_charge_perBR_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_charge_plot_perBR, width = 6, height = 4)

#CDR3 charge distribution plot for V gene families 
LC_V_family_distribution_charge_plot <- LC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot() +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    mutate(
      gene_numeric = as.numeric(str_extract(gene, "\\d+")),
      gene = fct_reorder(gene, gene_numeric)) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = LC_cdr3_aa_charge_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_beeswarm(
    data = LC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    mutate(
      gene_numeric = as.numeric(str_extract(gene, "\\d+")),
      gene = fct_reorder(gene, gene_numeric)) %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = LC_cdr3_aa_charge, shape = group_ID, group = group_ID), 
    dodge.width = 0.8, method = "center", priority = "density", corral = "wrap", 
    corral.width = 0.15, alpha = 0.7) +
    #size = 1.0, position = position_jitterdodge(jitter.width = 0.2), alpha = 0.7) +
  #too squished together
 # geom_text(
  #  data = LC_gene_means_df %>%
   # dplyr::filter(type == "V_gene") %>%
   # group_by(group_ID) %>%
   # dplyr::distinct(gene, .keep_all = TRUE),
   # mapping = aes(x = gene, y = - 1, label = round(LC_cdr3_aa_charge_gene, 2)),
   # vjust = 4, size = 1.8, colour = "black", 
   # position = position_dodge2(width = 1, preserve = "total")) +
  facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Variable Gene Family CDR3 Charges",
       x = "Light Chain Gene Family", y = "Avg CDR3 Charge")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_V_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_V_family_distribution_charge_plot, width = 12, height = 4)

#CDR3 charge distribution plot for J gene families
LC_J_family_distribution_charge_plot <- LC_summary_df %>%
  filter(type == "J_gene") %>%
  ggplot() +
  geom_col(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = LC_cdr3_aa_charge_gene, fill = group_ID), 
    position = position_dodge2(width = 0.8, preserve = "total")) +
  geom_point(
    data = LC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = LC_cdr3_aa_charge, shape = group_ID, group = group_ID), 
    size = 1.3, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_text(
    data = LC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = - 1, label = round(LC_cdr3_aa_charge_gene, 2)),
    vjust = 4, size = 1.8, colour = "black", 
    position = position_dodge2(width = 1, preserve = "total")) +
  facet_grid(cols = vars(isotype), scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Joint Gene Family CDR3 Charges",
       x = "Light Chain Gene Family", y = "Avg CDR3 Charge")
    ggsave(filename = file.path(plots_output_dir, paste0("LC_J_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = LC_J_family_distribution_charge_plot, width = 6, height = 4)

#kappa:lambda ratios plot 
kappa_lambda_ratio_plot <- LC_kappalambda_ratio_df %>%
  ggplot() +
  geom_col(
    data = LC_kappalambda_ratio_df %>%
    distinct(group_ID, kappa_lambda_ratio_group), 
    aes(x = group_ID, y = kappa_lambda_ratio_group, fill = group_ID)) +
  #geom_errorbar(
   # data = LC_kappalambda_ratio_df %>%
   ## distinct(group_ID, kappa_lambda_ratio_group),
   # aes(x = group_ID, y = kappa_lambda_ratio_group,
   #  ymin = mean - sd, ymax = mean + sd), width = 0.3) +
  geom_beeswarm(
    data = LC_kappalambda_ratio_df,
    mapping = aes(x = group_ID, y = kappa_lambda_ratio_br, shape = group_ID), 
                  method = "center", priority = "density", corral = "wrap", corral.width = 0.2,
                  alpha = 0.7) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Light Chain Kappa to Lambda Ratio")
  ggsave(filename = file.path(plots_output_dir, paste0("LC_kappa_lambda_ratio_", exp_group, "_vs_", ctrl_group, ".png")),
    plot = kappa_lambda_ratio_plot, width = 6, height = 4)



#V:J pairs set chord colors
row.col = adjustcolor(c(
  VK1 = "#0356BA", VK2 = "#00759A", VK3 = "#008B90", VK4 = "#00A07F", VK5 = "#16B900", VK6 = "#AAB200",
  VL1 = "#905100", VL2 = "#CD2C00", VL3 = "#F4004F", VL4 = "#FF00AC", VL5 = "#F900FF", VL6 = "#9F00FF", 
  VL7 = "#a200ffaf", VL8 = "#6A00FF", VL9 = "#5800d4cb", VL10 = "#2000d4cb"), alpha.f = 0.5)

grid.col = c(
  VK1 = "#0356BA", VK2 = "#00759A", VK3 = "#008B90", VK4 = "#00A07F", VK5 = "#16B900", VK6 = "#AAB200",
  VL1 = "#905100", VL2 = "#CD2C00", VL3 = "#F4004F", VL4 = "#FF00AC", VL5 = "#F900FF", VL6 = "#9F00FF", 
  VL7 = "#a200ffaf", VL8 = "#6A00FF", VL9 = "#5800d4cb", VL10 = "#2000d4cb",
  JK1 = "#2f6fbdff", JK2 = "#34849cff", JK3 = "#348f92ff", JK4 = "#3b9785ff", JK5 = "#53bd45ff", 
  JL1 = "#997240ff", JL2 = "#cc6f55ff", JL3 = "#f55e8eff")

#V:J pairs chord diagram 
VJ_pairs_LC_chord_diagram_exp <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == exp_group) %>%
  separate(gene, into = c("VK", "JK", sep = ":")) %>%
  group_by(VK, JK) %>%
  summarise(percent = percent_gene, .groups = 'drop')
  VJ_pairs_LC_chord_diagram_exp$percent[VJ_pairs_LC_chord_diagram_exp$percent == 0] <- 1e-6
  VJ_pairs_LC_chord_diagram_exp$VK <- factor(VJ_pairs_LC_chord_diagram_exp$VK, levels = 
                                              mixedsort(unique(VJ_pairs_LC_chord_diagram_exp$VK)))
  VJ_pairs_LC_chord_diagram_exp$JK <- factor(VJ_pairs_LC_chord_diagram_exp$JK, levels = 
                                              mixedsort(unique(VJ_pairs_LC_chord_diagram_exp$JK)))
  circos.clear()
  all_sectors <- c(levels(VJ_pairs_LC_chord_diagram_exp$VK),
                   levels(VJ_pairs_LC_chord_diagram_exp$JK))
  gap_vec[length(levels(VJ_pairs_LC_chord_diagram_exp$VK))] <- 15
  gap_vec <- rep(5, length(all_sectors))
  circos.par(gap.after = gap_vec)
  png(filename = file.path(plots_output_dir, paste0("LC_VJ_pairs_chord_diagram_", exp_group, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    VJ_pairs_LC_chord_diagram_exp, 
    order = c(levels(VJ_pairs_LC_chord_diagram_exp$VK),
              levels(VJ_pairs_LC_chord_diagram_exp$JK)),
    link.sort = TRUE, 
    grid.col = grid.col, 
    reduce = 0,
    col = adjustcolor(grid.col[VJ_pairs_LC_chord_diagram_exp$VK], alpha.f = 0.5),
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
  title(paste0("Light Chain V:J Gene Pairings in ", exp_group), cex.main = 2)
  dev.off()

VJ_pairs_LC_chord_diagram_control <- LC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == ctrl_group) %>%
  separate(gene, into = c("VK", "JK", sep = ":")) %>%
  group_by(VK, JK) %>%
  summarise(percent = percent_gene, .groups = 'drop')
  VJ_pairs_LC_chord_diagram_control$percent[VJ_pairs_LC_chord_diagram_control$percent == 0] <- 1e-6
  VJ_pairs_LC_chord_diagram_control$VK <- factor(VJ_pairs_LC_chord_diagram_control$VK, levels = 
                                              mixedsort(unique(VJ_pairs_LC_chord_diagram_control$VK)))
  VJ_pairs_LC_chord_diagram_control$JK <- factor(VJ_pairs_LC_chord_diagram_control$JK, levels = 
                                              mixedsort(unique(VJ_pairs_LC_chord_diagram_control$JK)))
  circos.clear()
  all_sectors <- c(levels(VJ_pairs_LC_chord_diagram_control$VK),
                   levels(VJ_pairs_LC_chord_diagram_control$JK))
  gap_vec[length(levels(VJ_pairs_LC_chord_diagram_control$VK))] <- 15
  gap_vec <- rep(5, length(all_sectors))
  circos.par(gap.after = gap_vec)
  png(filename = file.path(plots_output_dir, paste0("LC_VJ_pairs_chord_diagram_", ctrl_group, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    VJ_pairs_LC_chord_diagram_control, 
    order = c(levels(VJ_pairs_LC_chord_diagram_control$VK),
              levels(VJ_pairs_LC_chord_diagram_control$JK)),
    link.sort = TRUE, 
    grid.col = grid.col, 
    reduce = 0,
    col = adjustcolor(grid.col[VJ_pairs_LC_chord_diagram_control$VK], alpha.f = 0.5),
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
  title(paste0("Light Chain V:J Gene Pairings in ", ctrl_group), cex.main = 2)
  dev.off()



##CSV output

write_csv(LC_summary_df, file.path(csv_output_dir, paste0(project, "_LightChain_VJ_summarydata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_gene_means_df, file.path(csv_output_dir, paste0(project, "_LightChain_VJ_meandata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_kappalambda_ratio_df, file.path(csv_output_dir, paste0(project, "_LightChain_kappalambda_ratios_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_cdr3_avgBR_charge, file.path(csv_output_dir, paste0(project, "_LightChain_cdr3_avgBR_charge_ratios_", exp_group, "_vs_", ctrl_group, ".csv")))

write_csv(LC_summary_wide_df, file.path(csv_output_dir, paste0(project, "_LightChain_summary_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_means_wide_df, file.path(csv_output_dir, paste0(project, "_LightChain_means_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_percent_stdev_wide_df, file.path(csv_output_dir, paste0(project, "_LightChain_percent_stdev_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_percent_pvalue_wide_df, file.path(csv_output_dir, paste0(project, "_LightChain_percent_pvalue_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_cdr3_stdev_wide_df, file.path(csv_output_dir, paste0(project, "_LightChain_cdr3_stdev_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_cdr3_pvalue_wide_df, file.path(csv_output_dir, paste0(project, "_LightChain_cdr3_pvalue_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(LC_cell_counts_df, file.path(csv_output_dir, paste0(project, "_LightChain_cellcount_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))