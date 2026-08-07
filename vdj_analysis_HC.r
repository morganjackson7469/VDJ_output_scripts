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
library(purrr)
library(scales)
library(languageserver)
library(httpgd)
library(ggbeeswarm)
library(gtools)

#read summary file into tibble
HC_master_df <- summary_tsv_NGS %>%
 filter(chain == "heavy")

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
             "1455", "2177",	"2180",	"992",	"1823",	"1123",	"1299",	"1490",	"1507",	
             "1844",	"2133",	"1574",	"3500", "J10",	"J130",	"J20",	"J201",	"J203",
             "J218",	"J220",	"J24",	"J26",	"J34",	"J42",	"J46",	"J8",	"J9",	"J91",
             "J93", "6634",	"1559",	"1561",	"1607",	"1618",	"1629",	"1647",	"1690",	
             "1702",	"1707",	"1713",	"1724",	"1790", "1831",	"1837",	"1848",	"1850",	
             "6608",	"6611",	"6617",	"6618",	"6624",	"6658",	"6670",	"6691",	"6696",	"6714",	
             "6716",	"6741",	"6779",	"6781", "6794",	"6829",	"6854",	"6623",	"3780",	"3782",	
             "3587",	"3796",	"3834",	"3814", "3836",	"3825",	"3813",	"3839",	"3811",	"3805",
             "3808",	"2899",	"3845",	"3846",	"3847",	"3848",	"3850",	"3853", "3061",	"3255",
             "3267",	"381",	"1136",	"1145",	"1226",	"1288",	"1292",	"1313",	"1397",	"1403",
             "1446",	"1524",	"1700", "2024",	"2218",	"1503",	"1589",	"1606",	"1611",
             "1617",	"1636",	"1646",	"1650",	"1658",	"1692",	"1698",	"1717", "1722",	"1734",
             "1741",	"1747",	"1787",	"1796",	"1801",	"1825",	"1851",	"2726", "3197", "3369",
             "3936", "4284", "4643", "4720", "5139", "6565", "6655",	"606",	"1173",	"1238",
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
              "1483",	"1616",	"1638",	"1697",	"6774", "USCHC026",	"USCHC039",
              "USCHC005",	"USCH043",	"USCH020",	"USCH021",	"USCH032")

project <- "MAV_NGS"

new_plots_dir <- paste0("/mnt/md0/s440792/output/", project, "/ggplots")
 if (!dir.exists(new_plots_dir)) {
   dir.create(new_plots_dir, recursive = TRUE)
 }
new_csv_dir <- paste0("/mnt/md0/s440792/output/", project, "/csv_outputs")
 if (!dir.exists(new_csv_dir)) {
   dir.create(new_csv_dir, recursive = TRUE)
 }

plots_output_dir <- new_plots_dir
csv_output_dir <- new_csv_dir

experimental_group <- c("J10", "J130", "J20", "J201", "J203", "J218",
                      "J220", "J24", "J26", "J34", "J42", "J46", "J8",
                      "J9", "J91", "J93", "6634", "6608", "6611",
                      "6617", "6618", "6624", "6658", "6670", "6691",
                      "6696", "6714", "6716", "6741", "6779",
                      "6781", "6794", "6829", "6854", "6623", "6634",
                      "1455", "2177", "2180", "992", "1823", "1123",
                        "1299", "1490", "1507", "1844", "2133",
                        "1574", "2726", "3197", "3369", "3936", "4284",
                        "4643", "4720", "5139", "6565", "6655",
                        "606", "1173", "1238", "1301", "1468", "1579",
                        "1864", "1922")
control_group <- c("381", "1136", "1145", "1226", "1288", "1292", "1313",
                   "1397", "1403", "1446", "1503", "1524", "1700", "2024",
                   "2218")

exp_group <- "BAMS_all"
ctrl_group <- "LAMS"

BR_codes_BAMS_ty <- c("J10", "J130", "J20", "J201", "J203", "J218",
                      "J220", "J24", "J26", "J34", "J42", "J46", "J8",
                      "J9", "J91", "J93", "6634", "6608", "6611",
                      "6617", "6618", "6624", "6658", "6670", "6691",
                      "6696", "6714", "6716", "6741", "6779",
                      "6781", "6794", "6829", "6854", "6623", "6634")

BR_codes_BAMS <- c("1455", "2177", "2180", "992", "1823", "1123",
                        "1299", "1490", "1507", "1844", "2133",
                        "1574", "2726", "3197", "3369", "3936", "4284",
                        "4643", "4720", "5139", "6565", "6655",
                        "606", "1173", "1238", "1301", "1468", "1579",
                        "1864", "1922")
BR_codes_BAHC <- c("1559", "1561", "1607", "1618", "1629", "1647", "1690",
                    "1702", "1707", "1713", "1724", "1790", "1831", "1837",
                    "1848", "1850", "1589", "1606", "1611", "1617", "1636",
                    "1646", "1650", "1658", "1692", "1698", "1717", "1722",
                    "1734", "1741", "1747", "1787", "1796", "1801", "1825",
                    "1851")
BR_codes_LAMS <- c("381", "1136", "1145", "1226", "1288", "1292", "1313",
                   "1397", "1403", "1446", "1503", "1524", "1700", "2024",
                   "2218")


expected_genes_VHfam <- c(paste0("VH", 1:7))
expected_genes_JHfam <- c(paste0("JH", 1:6))
expected_genes_VJHpairs <- c(
  with(expand.grid(V = paste0("VH", 1:7), J = paste0("JH", 1:6)), paste(V, J, sep = ":")))

##FUNCTIONS
#filter for experimental and control groups
HC_master_df <- HC_master_df %>%
  mutate(group_ID = case_when(
    BR_code %in% experimental_group ~ exp_group,
    BR_code %in% control_group ~ ctrl_group)) %>%
  drop_na(cdr3_aa_charge, v_call, j_call)

HC_comparison_df <- HC_master_df %>%
  filter(group_ID %in% c(exp_group, ctrl_group)) %>%
  mutate(VJ_pair = paste0(v_call, ":", j_call)) %>%
  mutate(
    V_gene_mut = str_extract(v_call, "IGHV(\\d+)") %>%
      str_replace("IGHV", "VH"),
    J_gene_mut = str_extract(j_call, "IGHJ(\\d+)") %>%
      str_replace("IGHJ", "JH"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)) %>%
  relocate(V_gene_mut, J_gene_mut, VJ_only_gene, VJ_pair, .before = v_call)
  #select(-LC_isotype)


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
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE),
           cdr3_stdev = sd(HC_cdr3_aa_charge))%>%
    ungroup() %>%
    group_by(gene) %>%
    mutate(
          percent_p_value = tryCatch(
            t.test(percent ~ group_ID)$p.value,
            error = function(e)NA_real_),
          cdr3_p_value = tryCatch(
            t.test(HC_cdr3_aa_charge ~ group_ID)$p.value,
            error = function(e)NA_real_))
            


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
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene) %>%
    mutate(percent_stdev = sd(percent),
           cdr3_stdev = sd(HC_cdr3_aa_charge)) %>%
        ungroup() %>%
    group_by(gene) %>%
    mutate(
          percent_p_value = tryCatch(
            t.test(percent ~ group_ID)$p.value,
            error = function(e)NA_real_),
          cdr3_p_value = tryCatch(
            t.test(HC_cdr3_aa_charge ~ group_ID)$p.value,
            error = function(e)NA_real_))

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
    select(-total_hits) %>%
        ungroup() %>%
    group_by(group_ID, gene) %>%
    mutate(percent_stdev = sd(percent),
           cdr3_stdev = sd(HC_cdr3_aa_charge)) %>%
        ungroup() %>%
    group_by(gene) %>%
    mutate(percent_p_value = tryCatch(
            t.test(percent ~ group_ID)$p.value,
            error = function(e)NA_real_),
           cdr3_p_value = tryCatch(
           t.test(HC_cdr3_aa_charge ~ group_ID)$p.value,
           error = function(e)NA_real_))

HC_summary_df <- bind_rows(
  HC_V_summary_df,
  HC_J_summary_df,
  HC_VJ_summary_df
)

HC_gene_means_df <- HC_summary_df %>%
    group_by(group_ID, gene, type) %>%
    summarize(
      hit_count_gene = sum(hit_count),
      HC_cdr3_aa_charge_gene = mean(HC_cdr3_aa_charge),
      percent_stdev = unique(percent_stdev),
      percent_p_value = unique(percent_p_value),
      cdr3_stdev = unique(cdr3_stdev),
      cdr3_p_value = unique(cdr3_p_value)) %>%
    ungroup() %>%
    group_by(group_ID, type) %>%
    mutate(
      total_hits_family = sum(hit_count_gene),
      percent_gene = if_else(total_hits_family == 0, 0, ((hit_count_gene / total_hits_family) * 100)),
      HC_cdr3_aa_charge_family = mean(HC_cdr3_aa_charge_gene)) 

#to calculate cdr3 charge per BR code regardless of gene family - will need to add to csv export
HC_cdr3_avgBR_charge <- HC_comparison_df %>%
    group_by(group_ID, BR_code) %>%
    filter(chain == "heavy") %>%
    summarize(
      mean_br_cdr3charge = mean(cdr3_aa_charge)) %>%
    ungroup() %>%
    group_by(group_ID) %>%
    mutate(
      group_cdr3charge = mean(mean_br_cdr3charge), 
      stdev = sd(mean_br_cdr3charge, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      p_value = tryCatch(
      t.test(mean_br_cdr3charge ~ group_ID)$p.value,
      error = function(e)NA_real_))

#make wide df for presentations 
HC_summary_wide_df <- HC_summary_df %>%
    filter(type != "VJ_pair") %>%
    select(-hit_count, -type, -HC_cdr3_aa_charge, -percent_p_value, -cdr3_stdev,
           -cdr3_p_value, -percent_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent")

HC_means_wide_df <- HC_gene_means_df %>%
    ungroup() %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family, -HC_cdr3_aa_charge_family,
            -percent_p_value, -cdr3_stdev, -cdr3_p_value, -percent_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent_gene")

HC_percent_stdev_wide_df <- HC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family,
           -HC_cdr3_aa_charge_family, -percent_gene, -percent_p_value, -cdr3_stdev, -cdr3_p_value) %>%
    pivot_wider(names_from = "gene", values_from = "percent_stdev")

HC_percent_pvalue_wide_df <- HC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family,
           -HC_cdr3_aa_charge_family, -percent_gene, -percent_stdev, -cdr3_p_value, -cdr3_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "percent_p_value")

HC_cdr3gene_stdev_wide_df <- HC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family,
           -HC_cdr3_aa_charge_family, -percent_gene, -percent_p_value, -percent_stdev, -cdr3_p_value) %>%
    pivot_wider(names_from = "gene", values_from = "cdr3_stdev")

HC_cdr3gene_pvalue_wide_df <- HC_gene_means_df %>%
    group_by(group_ID) %>%
    filter(type != "VJ_pair") %>%
    select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family,
           -HC_cdr3_aa_charge_family, -percent_gene, -percent_stdev, -percent_p_value, -cdr3_stdev) %>%
    pivot_wider(names_from = "gene", values_from = "cdr3_p_value")

HC_cell_counts_df <- HC_summary_df %>%
    group_by(group_ID, BR_code) %>%
    filter(type == "V_gene") %>%
    summarize(cells = sum(hit_count))


###HEAVY CHAIN PLOTS

# Variable Family Distribution Plots 
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
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "# of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_counts_plot, width = 6, height = 4)

##!!! use this one 
HC_V_family_distribution_percent_plot <- HC_summary_df %>%
  filter(type == "V_gene") %>%
  ggplot() +
  geom_col(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), 
    position = position_dodge2(width = 1, preserve = "total")) +
  geom_beeswarm(
    data = HC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
    dodge.width = 0.8, method = "center",  
    priority = "density", corral = "wrap", 
    corral.width = 0.18, alpha = 0.7) +
  #geom_errorbar(
   # data = HC_gene_means_df %>%
   # dplyr::filter(type == "V_gene") %>%
   # group_by(group_ID) %>%
   # dplyr::distinct(gene, .keep_all = TRUE),
   # mapping = aes(x = gene, y= percent_gene, ymin = percent_gene - stdev, ymax = percent_gene + stdev,
   # group = group_ID),
   # width = 0.1, position = position_dodge2(width = 1, preserve = "total")
  #)
  #geom_text(
   # data = HC_gene_means_df %>%
    #dplyr::filter(type == "V_gene") %>%
    #group_by(group_ID) %>%
    #dplyr::distinct(gene, .keep_all = TRUE),
    #mapping = aes(x = gene, y = - 1, label = round(percent_gene, 2)),
    #vjust = 4, size = 1.8, colour = "black", 
    #position = position_dodge2(width = 1, preserve = "total")) +
  #coord_cartesian(ylim = c(-5, max(HC_gene_means_df$percent_gene))) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Variable Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "% of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_percent_plot, width = 6, height = 4)


# Joint Family Distribution Plots
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
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "# of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_counts_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_counts_plot, width = 6, height = 4)

##!! USE THIS ONE 
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
  geom_beeswarm(
    data = HC_summary_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
    dodge.width = 0.8, method = "center",  
    priority = "density", corral = "wrap", 
    corral.width = 0.18, alpha = 0.7
    #size = 1.3, position = position_jitterdodge(jitter.width = 0.2)
    ) +
  #geom_text(
   # data = HC_gene_means_df %>%
   # dplyr::filter(type == "J_gene") %>%
   # group_by(group_ID) %>%
   # dplyr::distinct(gene, .keep_all = TRUE),
   # mapping = aes(x = gene, y = - 1, label = round(percent_gene, 2)),
   # vjust = 4, size = 1.8, colour = "black", 
   # position = position_dodge2(width = 1, preserve = "total")) +
  #coord_cartesian(ylim = c(-5, max(HC_gene_means_df$percent_gene))) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Joint Gene Family Distribution",
       x = "Heavy Chain Gene Family", y = "% of Hits")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_percent_plot, width = 6, height = 4)



#V:J PAIRINGS counts and percentage plots

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
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    labs(title = "Top 10 V:J pairs Heavy Chain",
        x = "V:J pairs", y = "Count")
      ggsave(filename = file.path(plots_output_dir, paste0("HC_V:J_pairs_top10_count_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = VJ_pairs_HC_top10_count_plot, width = 12, height = 4)

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
  theme(axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = -2, margin=margin(r=15)), 
        axis.title.x = element_text(hjust = 0.5, size = 14, margin=margin(t=20))) +
  labs(title = "Top 10 V:J pairs Heavy Chain",
       x = "V:J pairs", y = "Percentage")
  ggsave(filename = file.path(plots_output_dir, paste0("HC_V:J_pairs_top10_percent_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = VJ_pairs_HC_top10_percent_plot, width = 18, height = 8)



#CDR3 CHARGES
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
  ylim(-3.0, 0.5) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heavy Chain CDR3 Charges",
       x = "Heavy Chain Gene Family", y = "Avg CDR3 charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_charge_plot, width = 6, height = 4)

#CDR3 charge averages per BR code (avg of total br)
#not average of families per br code 
HC_charge_plot_perBR <- HC_cdr3_avgBR_charge %>%
  ggplot() +
  geom_col(
    data = HC_cdr3_avgBR_charge %>%
    group_by(group_ID) %>%
    distinct(group_cdr3charge),
    mapping = aes(x = group_ID, y = group_cdr3charge, fill = group_ID), 
    position = "dodge") +
  geom_beeswarm(
    data = HC_cdr3_avgBR_charge %>%
    group_by(group_ID, BR_code),
    mapping = aes(x = group_ID, y = mean_br_cdr3charge, shape = group_ID, group = group_ID), 
      dodge.width = 0.8, method = "center",  
      priority = "density", corral = "wrap", 
      corral.width = 0.18, alpha = 0.7
    #size = 1.3, position = position_jitterdodge(jitter.width = 0.2)
    ) +
  #geom_text(data = summary_df %>%
              #dplyr::filter(type == "V_gene" | type == "J_gene") %>%
              #group_by(group_ID, ID, type) %>%
              #summarize(hit_count = sum(hit_count), .groups = "drop") %>%
              #mutate(measure = "HC_cdr3_aa_charge"),
            #aes(x = ID, y = 0, label = hit_count, fill = group_ID)) +
  ylim(-3.0, 0.6) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Heavy Chain CDR3 Charges",
       x = "Group", y = "Avg CDR3 charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_CDR3_charge_perBR_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_charge_plot_perBR, width = 6, height = 4)

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
  geom_beeswarm(
    data = HC_summary_df %>%
    dplyr::filter(type == "V_gene") %>%
    group_by(group_ID), 
    mapping = aes(x = gene, y = HC_cdr3_aa_charge, shape = group_ID, group = group_ID), 
      dodge.width = 0.8, method = "center",  
      priority = "density", corral = "wrap", 
      corral.width = 0.18, alpha = 0.7
    #size = 1.3, position = position_jitterdodge(jitter.width = 0.2)
    ) +
  #geom_text(
   # data = HC_gene_means_df %>%
   # dplyr::filter(type == "V_gene") %>%
   # group_by(group_ID) %>%
   # dplyr::distinct(gene, .keep_all = TRUE),
   # mapping = aes(x = gene, y = - 1, label = round(HC_cdr3_aa_charge_gene, 2)),
   # vjust = 4, size = 1.8, colour = "black", 
   # position = position_dodge2(width = 1, preserve = "total")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Variable Gene Family CDR3 Charges",
       x = "Heavy Chain Gene Family", y = "Avg CDR3 Charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_charge_plot, width = 6, height = 4)

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
  geom_text(
    data = HC_gene_means_df %>%
    dplyr::filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE),
    mapping = aes(x = gene, y = - 1, label = round(HC_cdr3_aa_charge_gene, 2)),
    vjust = 4, size = 1.8, colour = "black", 
    position = position_dodge2(width = 1, preserve = "total")) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(title = "Joint Gene Family CDR3 Charges",
       x = "Heavy Chain Gene Family", y = "Avg CDR3 Charge")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_CDR3_charge_", exp_group, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_charge_plot, width = 6, height = 4)


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
  summarise(percent = percent_gene, .groups = 'drop')
  VJ_pairs_HC_chord_diagram_exp$percent[VJ_pairs_HC_chord_diagram_exp$percent == 0] <- 1e-6
  VJ_pairs_HC_chord_diagram_exp$VH <- factor(VJ_pairs_HC_chord_diagram_exp$VH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_exp$VH)))
  VJ_pairs_HC_chord_diagram_exp$JH <- factor(VJ_pairs_HC_chord_diagram_exp$JH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_exp$JH)))
  circos.clear()
  all_sectors <- c(levels(VJ_pairs_HC_chord_diagram_exp$VH),
                   levels(VJ_pairs_HC_chord_diagram_exp$JH))
  gap_vec <- c(
                #big gap between VH and JH
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp$VH))-1), 15,
                #closing gap
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp$JH)) -1), 15)
  circos.par(gap.after = gap_vec)
  png(filename = file.path(plots_output_dir, paste0("HC_VJ_pairs_chord_diagram_", exp_group, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    VJ_pairs_HC_chord_diagram_exp, 
    order = c(levels(VJ_pairs_HC_chord_diagram_exp$VH),
              levels(VJ_pairs_HC_chord_diagram_exp$JH)),
    link.sort = TRUE, 
    grid.col = grid.col, 
    reduce = 0,
    col = adjustcolor(grid.col[VJ_pairs_HC_chord_diagram_exp$VH], alpha.f = 0.5),
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
  summarise(percent = percent_gene, .groups = 'drop')
  VJ_pairs_HC_chord_diagram_control$percent[VJ_pairs_HC_chord_diagram_control$percent == 0] <- 1e-6
  VJ_pairs_HC_chord_diagram_control$VH <- factor(VJ_pairs_HC_chord_diagram_control$VH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_control$VH)))
  VJ_pairs_HC_chord_diagram_control$JH <- factor(VJ_pairs_HC_chord_diagram_control$JH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_control$JH)))
  circos.clear()
  all_sectors <- c(levels(VJ_pairs_HC_chord_diagram_control$VH),
                   levels(VJ_pairs_HC_chord_diagram_control$JH))
  gap_vec <- c(
                #big gap between VH and JH
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp$VH))-1), 15,
                #closing gap
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp$JH)) -1), 15)
  circos.par(gap.after = gap_vec)
  png(filename = file.path(plots_output_dir, paste0("HC_VJ_pairs_chord_diagram_", ctrl_group, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    VJ_pairs_HC_chord_diagram_control, 
    order = c(levels(VJ_pairs_HC_chord_diagram_control$VH),
              levels(VJ_pairs_HC_chord_diagram_control$JH)),
    link.sort = TRUE, 
    grid.col = grid.col, 
    reduce = 0,
    col = adjustcolor(grid.col[VJ_pairs_HC_chord_diagram_control$VH], alpha.f = 0.5),
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
  title(paste0("Heavy Chain V:J Gene Pairings in ", ctrl_group), cex.main = 2)
  dev.off()


##CSV output

write_csv(HC_summary_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_VJ_summarydata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_gene_means_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_VJ_meandata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cdr3_avgBR_charge, file.path(csv_output_dir, paste0(project, "_HeavyChain_cdr3_avgBR_charge_", exp_group, "_vs_", ctrl_group, ".csv")))

write_csv(HC_summary_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_summary_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_means_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_means_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_percent_stdev_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_percent_stdev_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_percent_pvalue_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_percent_pvalue_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cdr3gene_stdev_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_cdr3_stdev_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cdr3gene_pvalue_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_cdr3_pvalue_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cell_counts_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_cellcount_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))