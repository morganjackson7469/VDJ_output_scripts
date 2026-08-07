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
library(ggprism)
library(forcats)
library(ggpubr)
library(lmtest)
library(broom)
library(rstatix)

load("NGS_summary_workspace.RData")

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

project <- "MAV_NGS_3"

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

experimental_group_A <- c("1455", "2177", "2180", "992", "1823", "1123",
                        "1299", "1490", "1507", "1844", "2133",
                        "1574", "2726", "3197", "3369", "3936", "4284",
                        "4643", "4720", "5139", "6565", "6655",
                        "606", "1173", "1238", "1301", "1468", "1579",
                        "1864", "1922")
experimental_group_B <- c("J10", "J130", "J20", "J201", "J203", "J218",
                      "J220", "J24", "J26", "J34", "J42", "J46", "J8",
                      "J9", "J91", "J93", "6634", "6608", "6611",
                      "6617", "6618", "6624", "6658", "6670", "6691",
                      "6696", "6714", "6716", "6741", "6779",
                      "6781", "6794", "6829", "6854", "6623")
control_group <- c("1559", "1561", "1607", "1618", "1629", "1647", "1690",
                    "1702", "1707", "1713", "1724", "1790", "1831", "1837",
                    "1848", "1850", "1589", "1606", "1611", "1617", "1636",
                    "1646", "1650", "1658", "1692", "1698", "1717", "1722",
                    "1734", "1741", "1747", "1787", "1796", "1801", "1825",
                    "1851")

## T TEST WILL NOT WORK WITH 3 GROUPS 
#how to set it up to automatically detect if there are 2 groups or 3 or more? 
exp_group_A <- "BAMS_nonNAT"
exp_group_B <- "BAMS_NAT"
ctrl_group <- "BAHC"

BR_codes_BAMS_NAT <- c("J10", "J130", "J20", "J201", "J203", "J218",
                      "J220", "J24", "J26", "J34", "J42", "J46", "J8",
                      "J9", "J91", "J93", "6634", "6608", "6611",
                      "6617", "6618", "6624", "6658", "6670", "6691",
                      "6696", "6714", "6716", "6741", "6779",
                      "6781", "6794", "6829", "6854", "6623", "6634")
BR_codes_BAMS_nonNAT <- c("1455", "2177", "2180", "992", "1823", "1123",
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

#colors and shapes for the points
#change to exp_group_A, exp_group_B, cntl_group 
group_shapes <- c(
  "BAHC" = 15, #square
  "BAMS_nonNAT" = 16, #circle
  "BAMS_NAT" = 17 #triangle
)

group_colors <- c(
  "BAHC" = "#73ACE0", 
  "BAMS_nonNAT" = "#6EDF72", 
  "BAMS_NAT" = "#CEA57A")


#scott and sams scripts do not account for zeros if a gene or family is missing for a BR code
expected_genes_VHfam <- c(paste0("VH", 1:7))
expected_genes_JHfam <- c(paste0("JH", 1:6))
expected_genes_VJHpairs <- c(
 with(expand.grid(V = paste0("VH", 1:7), J = paste0("JH", 1:6)), paste(V, J, sep = ":")))

##FUNCTIONS
#filter for experimental and control groups
HC_master_df <- HC_master_df %>%
  mutate(group_ID = case_when(
    BR_code %in% experimental_group_A ~ exp_group_A,
    BR_code %in% experimental_group_B ~ exp_group_B,
    BR_code %in% control_group ~ ctrl_group)) %>%
  drop_na(cdr3_aa_charge, v_call, j_call)

HC_comparison_df <- HC_master_df %>%
  filter(group_ID %in% c(exp_group_A, exp_group_B, ctrl_group)) %>%
  mutate(VJ_pair = paste0(v_call, ":", j_call)) %>%
  mutate(
    V_gene_mut = str_extract(v_call, "IGHV(\\d+)") %>%
      str_replace("IGHV", "VH"),
    J_gene_mut = str_extract(j_call, "IGHJ(\\d+)") %>%
      str_replace("IGHJ", "JH"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)) %>%
  relocate(V_gene_mut, J_gene_mut, VJ_only_gene, VJ_pair, .before = v_call)
  #select(-LC_isotype)

# Individual V genes 
## THIS ANALYSIS NEEDS TO BE CORRECTED FOR THE MEAN CALCULATION percent_gene
HC_V_call_df <- HC_comparison_df %>%
    #remove allele number and extra allele calls 
    group_by(group_ID, BR_code, v_call) %>%
    mutate(
      v_call = map_chr(v_call, function(x) {
        v_genes <- str_split(x, ",")[[1]]
        v_genes <- sub("\\*.*$", "", v_genes)
        v_genes <- v_genes[!grepl("D$", v_genes)]
        v_genes <- sub("^IGHV", "VH", v_genes)
        #v_genes <- sub("^VH4-4$", "VH4-04", v_genes)
        v_genes <- unique(v_genes)
        paste(v_genes, collapse = ",")
      })) %>%
    summarize(
      hit_count = n()) %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, v_call) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE))
#this means calculation is incorrect 
HC_V_call_means_df <- HC_V_call_df %>%
    group_by(group_ID, v_call) %>%
      summarize(
        hit_count_gene = sum(hit_count),
        percent_stdev = unique(percent_stdev)) %>%
      mutate(
        total_hits_gene = sum(hit_count_gene),
        percent_gene = if_else(total_hits_gene == 0, 0, ((hit_count_gene / total_hits_gene) * 100)))

# individual V gene and J gene pairings 
HC_VJ_call_df <- HC_comparison_df %>%
    #remove allele number and extra allele calls 
    group_by(group_ID, BR_code, VJ_pair) %>%
    mutate(
    V_gene_mut = map_chr(v_call, function(x) {
        v_genes <- str_split(x, ",")[[1]]
        v_genes <- sub("\\*.*$", "", v_genes)
        v_genes <- sub("^IGHV", "VH", v_genes)
        v_genes <- sub("^VH4-4$", "VH4-04", v_genes)
        v_genes <- unique(v_genes)
        paste(v_genes, collapse = ",")}),
    J_gene_mut = str_extract(VJ_pair, "(?<=:)IGHJ.+") %>%
      str_remove("\\*.*$") %>%
      str_replace("^IGHJ", "JH"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)) %>%
    ungroup() %>%
    group_by(group_ID, BR_code, VJ_only_gene) %>%
    summarize(
      V_gene_mut = first(V_gene_mut),
      J_gene_mut = first(J_gene_mut),
      VJ_only_gene = first(VJ_only_gene),
      hit_count = n(),
      cdr3_aa_length = mean(cdr3_aa_length, na.rm = TRUE),
      cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    #select(total_hits) %>%
    ungroup() %>%
    group_by(group_ID, VJ_only_gene) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE)) %>%
    mutate(
      mean_cdr3_length = mean(cdr3_aa_length, na.rm = TRUE),
      mean_cdr3_charge = mean(cdr3_aa_charge, na.rm = TRUE))
#this means calculation is incorrect 
HC_VJ_call_means_df <- HC_VJ_call_df %>%
    group_by(group_ID, VJ_only_gene) %>%
    summarize(
      hit_count_gene = sum(hit_count),
      percent_stdev = unique(percent_stdev),
      mean_cdr3_length = mean(cdr3_aa_length, na.rm = TRUE),
      mean_cdr3_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    mutate(
      total_hits_gene = sum(hit_count_gene),
      percent_gene = if_else(total_hits_gene == 0, 0, ((hit_count_gene / total_hits_gene) * 100)))

#Individual V genes to J gene pairings - VH4 ONLY
#need to make VH4-4 to VH4-04 and move it to the front 
HC_VH4_J_call_df <- HC_comparison_df %>%
    #remove allele number and extra allele calls
    #this grouping and filter will not work here before V_gene_mut is made 
    filter(V_gene_mut == "VH4") %>%
    group_by(group_ID, BR_code, VJ_pair) %>%
    mutate(
    V_gene_mut = map_chr(v_call, function(x) {
      #the [1] removes everything after the first comma for ones with multiple gene calls
      # this will need to actually be split into equal amounts for each gene represented 
      # scott said the ones with multiple calls are equal homology 
      # if theres two, then 0.5 should go into each call, not sure how to do that with percentages
      # maybe 0.5 for each hit count? 
        v_genes <- str_split(x, ",")[[1]][1]
        v_genes <- sub("\\*.*$", "", v_genes)
        v_genes <- sub("^IGHV", "VH", v_genes)
        v_genes <- sub("^VH4-4$", "VH4-04", v_genes)
        v_genes <- unique(v_genes)
        paste(v_genes, collapse = ",")}),
    J_gene_mut = str_extract(VJ_pair, "(?<=:)IGHJ.+") %>%
      str_remove("\\*.*$") %>%
      str_replace("^IGHJ", "JH"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)) %>%
    ungroup() %>%
    group_by(group_ID, BR_code, VJ_only_gene) %>%
    summarize(
      V_gene_mut = first(V_gene_mut),
      J_gene_mut = first(J_gene_mut),
      VJ_only_gene = first(VJ_only_gene),
      hit_count = n(),
      cdr3_aa_length = mean(cdr3_aa_length, na.rm = TRUE),
      cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    #select(total_hits) %>%
    ungroup() %>%
    group_by(group_ID, VJ_only_gene) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE)) %>%
    mutate(
      mean_cdr3_length = mean(cdr3_aa_length, na.rm = TRUE),
      mean_cdr3_charge = mean(cdr3_aa_charge, na.rm = TRUE))
#this mean calculation is incorrect 
HC_VH4_J_call_means_df <- HC_VH4_J_call_df %>%
    group_by(group_ID, VJ_only_gene) %>%
    summarize(
      hit_count_gene = sum(hit_count),
      percent_stdev = unique(percent_stdev)) %>%
    mutate(
      total_hits_gene = sum(hit_count_gene),
      percent_gene = if_else(total_hits_gene == 0, 0, ((hit_count_gene / total_hits_gene) * 100))
    )

#Individual VH4 genes only NO J GENE 
HC_VH4_call_df <- HC_comparison_df %>%
    #remove allele number and extra allele calls
    filter(V_gene_mut == "VH4") %>%
    select(-J_gene_mut) %>%
    #group_by(group_ID, BR_code, V_gene_mut) %>%
    mutate(
    V_gene_mut = map_chr(v_call, function(x) {
      #the [1] removes everything after the first comma for ones with multiple gene calls
        v_genes <- str_split(x, ",")[[1]][1]
        v_genes <- sub("\\*.*$", "", v_genes)
        v_genes <- sub("^IGHV", "VH", v_genes)
        v_genes <- sub("^VH4-4$", "VH4-04", v_genes)
        v_genes <- unique(v_genes)
        paste(v_genes, collapse = ",")})) %>%
    filter(V_gene_mut != "VH4/OR15-8") %>%
    #ungroup() %>%
    group_by(group_ID, BR_code, V_gene_mut) %>%
    summarize(
      V_gene_mut = first(V_gene_mut),
      hit_count = n(),
      cdr3_aa_length = mean(cdr3_aa_length, na.rm = TRUE),
      cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    mutate(
      total_hits = sum(hit_count),
      percent = if_else(total_hits == 0, 0, ((hit_count / total_hits) * 100))) %>%
    #select(total_hits) %>%
    ungroup() %>%
    group_by(group_ID, V_gene_mut) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE)) %>%
    mutate(
      mean_cdr3_length = mean(cdr3_aa_length, na.rm = TRUE),
      mean_cdr3_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    ungroup()

HC_VH4_call_outliers <- HC_VH4_call_df %>%
    group_by(group_ID, V_gene_mut) %>%
    identify_outliers(percent)

#percent_gene needs to be recalculated this is all wrong, format it like VH4_bar
HC_VH4_call_means_df <- HC_VH4_call_df %>%
    group_by(group_ID, V_gene_mut) %>%
    summarize(
      hit_count_gene = sum(hit_count),
      percent_stdev = unique(percent_stdev),
      n_gene = n()) %>%
    mutate(
      total_hits_gene = sum(hit_count_gene),
      #FIX THIS HERE 
      percent_gene = if_else(total_hits_gene == 0, 0, ((hit_count_gene / total_hits_gene) * 100))
        )

      
write_csv(VH4_bar_df, file.path(csv_output_dir, paste0(project, "_HC_VH4_call_means_df.csv")))


HC_VH4_4_J6_family_distribution_percent_plot <- HC_VJ_call_df %>%
  filter(V_gene_mut == "VH4" &
          J_gene_mut == "JH6" &
          V_gene_mut != "OR15-8") %>%
  ggplot() +
  geom_col(
    data = HC_VJ_call_df %>%
    dplyr::filter(V_gene_mut == "VH4-4") %>%
    group_by(group_ID, J_gene_mut) %>%
    dplyr::distinct(V_gene_mut, .keep_all = TRUE),
    mapping = aes(x = VJ_only_gene, y = percent, fill = group_ID), 
    position = position_dodge2(width = 1, preserve = "total")) +
  geom_beeswarm(
    data = HC_VJ_call_df %>%
    dplyr::filter(V_gene_mut == "VH4-4") %>%
    group_by(group_ID, J_gene_mut), 
    mapping = aes(x = VJ_only_gene, y = percent, shape = group_ID, group = group_ID), 
    dodge.width = 0.8, method = "center", preserve.data.axis = TRUE,
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
  labs(title = "VH4-4 Gene Pairing",
       x = " ", y = "% of All Gene Pairings") 
  ggsave(filename = file.path(plots_output_dir, paste0("3HC_VH4:JH6_dis_percent_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".png")),
       plot = HC_VH4_family_distribution_percent_plot, width = 10, height = 4)

#make group ID factor to be reordered in the graph
HC_VH4_call_df$group_ID <- factor(HC_VH4_call_df$group_ID, levels=c('BAHC','BAMS_nonNAT','BAMS_NAT')) 


#is it possible to build somekind of template at the top to fill these dataframes in? (group_by(x, y) then define them here? idk)
#oneway.test is for welch's anova
 VH4_bar_df <- HC_VH4_call_df %>%
    group_by(group_ID, V_gene_mut) %>%
    summarize(
      avg_percent = mean(percent, na.rm = TRUE),
      percent_sd = sd(percent, na.rm = TRUE),
    .groups = "drop")
  VH4_point_df <- HC_VH4_call_df %>%
    group_by(group_ID, V_gene_mut) %>%
    ungroup()

#stats for just VH4-04
HC_VH4_04_filtered <- HC_VH4_call_df %>%
                      filter(V_gene_mut == "VH4-04") 

HC_VH4_04_model <- oneway.test(percent ~ group_ID, data = HC_VH4_04_filtered, var.equal = FALSE)
                print(HC_VH4_04_model)

HC_VH4_04_ttest <- HC_VH4_04_filtered %>%
  pairwise_t_test(
    percent ~ group_ID,
    p.adjust.method = "BH",
    pool.sd = FALSE
  ) %>%
  add_xy_position(x = "V_gene_mut", 
                #step.increase = 0.25)
  )
print(HC_VH4_04_ttest)

ymax <- VH4_bar_df %>%
    filter(V_gene_mut == "VH4-04") %>%
    summarise(y = max(avg_percent)) %>%
    pull(y)

#stats for all VH4 genes 
#HC_VH4_04_filtered <- HC_VH4_call_df %>%
 #                     filter(V_gene_mut == "VH4-04") 

HC_VH4_04_model <- oneway.test(percent ~ group_ID, data = HC_VH4_04_filtered, var.equal = FALSE)
                print(HC_VH4_04_model)

HC_VH4_04_ttest <- HC_VH4_04_filtered %>%
  pairwise_t_test(
    percent ~ group_ID,
    p.adjust.method = "BH",
    pool.sd = FALSE
  ) %>%
  add_xy_position(x = "V_gene_mut", 
                #step.increase = 0.25)
  )

HC_VH4_only_family_distribution_percent_plot <- HC_VH4_call_df %>%
  #filter(V_gene_mut == "VH4") %>%
  #filter(V_gene_mut != "VH4/OR15-8") %>%
  ggplot() +
  geom_col(
    data = VH4_bar_df,
    mapping = aes(x = V_gene_mut, y = avg_percent, fill = group_ID), 
    position = position_dodge2(width = 0.9),
    colour = "black", linewidth = 1.0) +
  #geom_beeswarm(
  #data = VH4_point_df,
   #mapping = aes(x = V_gene_mut, y = percent, shape = group_ID, group = group_ID), 
    #              dodge.width = 0.9, 
     #             method = "center", 
      #            preserve.data.axis = TRUE,
       #           priority = "density", 
        #          corral = "wrap", 
         #         corral.width = 0.05, 
          #        alpha = 0.9
           #       ) +
  #geom_errorbar(
   #   data = VH4_bar_df,
    #mapping = aes(x = V_gene_mut, y = avg_percent, ymin = avg_percent, ymax = avg_percent + percent_sd,
     #group = group_ID),
   #width = 0.5, 
   #position = position_dodge2(width = 0.9, padding = 0.6)
    # ) +
  #geom_text(
   # data = HC_gene_means_df %>%
    #dplyr::filter(type == "V_gene") %>%
    #group_by(group_ID) %>%
    #dplyr::distinct(gene, .keep_all = TRUE),
    #mapping = aes(x = gene, y = - 1, label = round(percent_gene, 2)),
    #vjust = 4, size = 1.8, colour = "black", 
    #position = position_dodge2(width = 1, preserve = "total")) +
  #coord_cartesian(ylim = c(-5, max(HC_gene_means_df$percent_gene))) + 
  stat_pvalue_manual(
    HC_VH4_04_ttest,
    size = 1, 
    p.format.style = "graphpad",
    hide.ns = TRUE) +
  theme_prism(
      palette = "black_and_white",
      base_size = 14, base_family = "arial") +
  theme(
        axis.title.y = element_text(size = 18),
        plot.title = element_text(
          size = 24, 
          face = "bold")
        ) +
  #theme(axis.title.y = element_text(size = 18)) +
  scale_fill_manual(values = group_colors) +
  #scale_shape_manual(values = group_shapes) +
  scale_y_continuous(expand = c(0,0.2), 
                    limits = c(0, 100)
                    ) +
  labs(title = "VH4 Genes",
       x = " ", y = "% VH4 Gene Frequency") 
  ggsave(filename = file.path(plots_output_dir, paste0("HC_VH4_only_bars_dis_3color_percent_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".png")),
       plot = HC_VH4_only_family_distribution_percent_plot, width = 16, height = 8, dpi = 600)

# V AND J GENE FAMILIES
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
       hit_count = NA_integer_, 
       HC_cdr3_aa_charge = NA_real_
       )) %>%  
    mutate(type = "V_gene") %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count, na.rm = TRUE),
      percent = if_else(total_hits == 0, NA_real_, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE),
           cdr3_stdev = sd(HC_cdr3_aa_charge, na.rm = TRUE))
           #%>%
    #ungroup() %>%
    #group_by(gene) %>%
    #mutate(
     #     percent_p_value = tryCatch(
      #      t.test(percent ~ group_ID)$p.value,
       #     error = function(e)NA_real_),
        #  cdr3_p_value = tryCatch(
         #   t.test(HC_cdr3_aa_charge ~ group_ID)$p.value,
          #  error = function(e)NA_real_))
            


HC_J_summary_df <- HC_comparison_df %>%
    group_by(group_ID, BR_code, J_gene_mut) %>%
    summarize(
      hit_count = n(),
      HC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE))%>% 
    rename(gene = J_gene_mut) %>%
    tidyr::complete(gene = expected_genes_JHfam,
     fill = list(
     hit_count = NA_integer_, 
     HC_cdr3_aa_charge = NA_real_
     #percent_value = 0, percent = 0), 
     # explicit = FALSE
    )) %>% 
    mutate(type = "J_gene") %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count, na.rm = TRUE),
      percent = if_else(total_hits == 0, NA_real_, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
    ungroup() %>%
    group_by(group_ID, gene) %>%
    mutate(percent_stdev = sd(percent, na.rm = TRUE),
           cdr3_stdev = sd(HC_cdr3_aa_charge, na.rm = TRUE)) #%>%
        #ungroup() %>%
    #group_by(gene) %>%
    #mutate(
     #     percent_p_value = tryCatch(
      #      t.test(percent ~ group_ID)$p.value,
       #     error = function(e)NA_real_),
        #  cdr3_p_value = tryCatch(
         #   t.test(HC_cdr3_aa_charge ~ group_ID)$p.value,
          #  error = function(e)NA_real_))

HC_VJ_summary_df <- HC_comparison_df %>%
    group_by(group_ID, BR_code, VJ_only_gene) %>%
      summarize(
        hit_count = n(),
        HC_cdr3_aa_charge = mean(cdr3_aa_charge, na.rm = TRUE)) %>%
    rename(gene = VJ_only_gene) %>%
     tidyr::complete(gene = expected_genes_VJHpairs,
      fill = list(
       hit_count = NA_integer_, 
       HC_cdr3_aa_charge = NA_real_ 
        #percent_value = 0, percent = 0), 
        #explicit = FALSE)
        ))%>%  
    mutate(type = "VJ_pair") %>%
    ungroup() %>%
    group_by(group_ID, BR_code) %>%
    mutate(
      total_hits = sum(hit_count, na.rm = TRUE),
      percent = if_else(total_hits == 0, NA_real_, ((hit_count / total_hits) * 100))) %>%
    select(-total_hits) %>%
        ungroup() %>%
    group_by(group_ID, gene) %>%
      mutate(percent_stdev = sd(percent, na.rm = TRUE),
            cdr3_stdev = sd(HC_cdr3_aa_charge, na.rm = TRUE)) #%>%
        #ungroup() %>%
    #group_by(gene) %>%
    #mutate(percent_p_value = tryCatch(
     #       t.test(percent ~ group_ID)$p.value,
      #      error = function(e)NA_real_),
       #    cdr3_p_value = tryCatch(
        #   t.test(HC_cdr3_aa_charge ~ group_ID)$p.value,
         #  error = function(e)NA_real_))

HC_summary_df <- bind_rows(
  HC_V_summary_df,
  HC_J_summary_df,
  HC_VJ_summary_df
)


HC_gene_means_df <- HC_summary_df %>%
    group_by(group_ID, gene, type) %>%
    summarize(
      hit_count_gene = sum(hit_count, na.rm = TRUE),
      HC_cdr3_aa_charge_gene = mean(HC_cdr3_aa_charge, na.rm = TRUE),
      percent_gene = mean(percent, na.rm = TRUE),
      percent_stdev = unique(percent_stdev),
      #percent_p_value = unique(percent_p_value),
      cdr3_stdev = unique(cdr3_stdev)) %>%
      #cdr3_p_value = unique(cdr3_p_value)) %>%
    ungroup() %>%
    group_by(group_ID, type) %>%
    mutate(
      total_hits_family = sum(hit_count_gene, na.rm = TRUE),
      #percent_gene = if_else(total_hits_family == 0, 0, ((hit_count_gene / total_hits_family) * 100)),
      HC_cdr3_aa_charge_family = mean(HC_cdr3_aa_charge_gene, na.rm = TRUE)) 

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
      select(-hit_count, -type, -HC_cdr3_aa_charge) %>%
              #-percent_p_value, -cdr3_stdev,
            #-cdr3_p_value, -percent_stdev) %>%
        pivot_wider(names_from = "gene", values_from = "percent")

  HC_means_wide_df <- HC_gene_means_df %>%
      ungroup() %>%
      filter(type != "VJ_pair") %>%
      select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family, -HC_cdr3_aa_charge_family) %>%
              #-percent_p_value, -cdr3_stdev, -cdr3_p_value, -percent_stdev) %>%
      pivot_wider(names_from = "gene", values_from = "percent_gene")

  HC_percent_stdev_wide_df <- HC_gene_means_df %>%
      group_by(group_ID) %>%
      filter(type != "VJ_pair") %>%
      select(-type, -hit_count_gene, -HC_cdr3_aa_charge_gene, -total_hits_family, -HC_cdr3_aa_charge_family, -percent_gene) %>% 
              #-percent_p_value, -cdr3_stdev, -cdr3_p_value) %>%
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

#make factor for correct order on graph 
HC_summary_df$group_ID <- factor(HC_summary_df$group_ID, levels=c('BAHC','BAMS_nonNAT','BAMS_NAT')) 
HC_gene_means_df$group_ID <- factor(HC_gene_means_df$group_ID, levels=c('BAHC','BAMS_nonNAT','BAMS_NAT')) 

#filter data for plot
 VH_bar_df <- HC_gene_means_df %>%
    filter(type == "V_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE) 
  VH_point_df <- HC_summary_df %>%
    filter(type == "V_gene") %>%
    group_by(group_ID) 

#plot data 
HC_V_family_distribution_percent_plot <- HC_summary_df %>%
  ggplot() +
  geom_col(
    data = VH_bar_df,
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), 
    position = position_dodge2(width = 0.9),
    colour = "black", linewidth = 0.5) +
  geom_beeswarm(
    data = VH_point_df, 
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
                  dodge.width = 0.9, 
                  method = "center",  
                  preserve.data.axis = TRUE,
                  priority = "density", 
                  corral = "wrap", 
                  corral.width = 0.05, 
                  alpha = 0.9) +
  geom_errorbar(
   data = VH_bar_df,
   mapping = aes(x = gene, y= percent_gene, 
                ymin = pmax(0, percent_gene - percent_stdev), 
                ymax = percent_gene + percent_stdev,
                group = group_ID),
                #width = 0.1, 
                position = position_dodge2(width = 0.9, padding = 0.6)
              ) +
  #geom_text(
    #mapping = aes(x = gene, y = - 1, label = round(percent_gene, 2)),
    #vjust = 4, size = 1.8, colour = "black", 
    #position = position_dodge2(width = 1, preserve = "total")) +
  #coord_cartesian(ylim = c(-5, max(HC_gene_means_df$percent_gene))) +
  theme_prism(
    palette = "black_and_white",
    base_size = 12,
    base_family = "arial") +
  theme(
    plot.title = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    size = 24,
    face = "bold") +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  scale_y_continuous(expand = mult(c(0, 0.02))) +
  labs(title = "VH Families",
       x = "Heavy Chain Gene Family", y = "% VH Family Usage")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_V_famdis_percent_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".png")),
       plot = HC_V_family_distribution_percent_plot, width = 14, height = 8)
### does plot() work here?? 

# Joint Family Distribution Plots.
#filter data for plot
 JH_bar_df <- HC_gene_means_df %>%
    filter(type == "J_gene") %>%
    group_by(group_ID) %>%
    dplyr::distinct(gene, .keep_all = TRUE) 
  JH_point_df <- HC_summary_df %>%
    filter(type == "J_gene") %>%
    group_by(group_ID) 

#plot data 
HC_J_family_distribution_percent_plot <- HC_summary_df %>%
  ggplot() +
  geom_col(
    data = JH_bar_df,
    mapping = aes(x = gene, y = percent_gene, fill = group_ID), 
    position = position_dodge2(width = 0.9),
    colour = "black", linewidth = 0.5) +
  geom_beeswarm(
    data = JH_point_df,
    mapping = aes(x = gene, y = percent, shape = group_ID, group = group_ID), 
                  dodge.width = 0.9, 
                  method = "center",  
                  preserve.data.axis = TRUE,
                  priority = "density", 
                  corral = "wrap", 
                  corral.width = 0.05, 
                  alpha = 0.9
                  #size = 1.3, position = position_jitterdodge(jitter.width = 0.2)
                  ) +
  geom_errorbar(
   data = JH_bar_df,
   mapping = aes(x = gene, y= percent_gene, 
                ymin = pmax(0, percent_gene - percent_stdev), 
                ymax = percent_gene + percent_stdev,
                group = group_ID),
                #width = 0.1, 
                position = position_dodge2(width = 0.9, padding = 0.6)
              ) +
  theme_prism(
    palette = "black_and_white",
    base_size = 12,
    base_family = "arial") +
  theme(
    plot.title = element_text(size = 18),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    size = 24,
    face = "bold") +
  scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  scale_y_continuous(expand = c(0, 0.25)) +
  labs(title = "JH Families",
       x = "Heavy Chain Gene Family", y = "% JH Family Usage")
    ggsave(filename = file.path(plots_output_dir, paste0("HC_J_famdis_percent_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".png")),
       plot = HC_J_family_distribution_percent_plot, width = 14, height = 8)



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
VJ_pairs_HC_chord_diagram_exp_group_A <- HC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == exp_group_A) %>%
  separate(gene, into = c("VH", "JH", sep = ":")) %>%
  group_by(VH, JH) %>%
  summarise(percent = percent_gene, .groups = 'drop')
  VJ_pairs_HC_chord_diagram_exp_group_A$percent[VJ_pairs_HC_chord_diagram_exp_group_A$percent == 0] <- 1e-6
  VJ_pairs_HC_chord_diagram_exp_group_A$VH <- factor(VJ_pairs_HC_chord_diagram_exp_group_A$VH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_exp_group_A$VH)))
  VJ_pairs_HC_chord_diagram_exp_group_A$JH <- factor(VJ_pairs_HC_chord_diagram_exp_group_A$JH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_exp_group_A$JH)))
  circos.clear()
  all_sectors <- c(levels(VJ_pairs_HC_chord_diagram_exp_group_A$VH),
                   levels(VJ_pairs_HC_chord_diagram_exp_group_A$JH))
  gap_vec <- c(
                #big gap between VH and JH
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp_group_A$VH))-1), 15,
                #closing gap
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp_group_A$JH)) -1), 15)
  circos.par(gap.after = gap_vec)
  png(filename = file.path(plots_output_dir, paste0("2HC_VJ_pairs_chord_diagram_", exp_group_A, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    VJ_pairs_HC_chord_diagram_exp_group_A, 
    order = c(levels(VJ_pairs_HC_chord_diagram_exp_group_A$VH),
              levels(VJ_pairs_HC_chord_diagram_exp_group_A$JH)),
    link.sort = TRUE, 
    grid.col = grid.col, 
    reduce = 0,
    col = adjustcolor(grid.col[VJ_pairs_HC_chord_diagram_exp_group_A$VH], alpha.f = 0.5),
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
        cex = 2.0)
    },
    bg.border = NA)
  title(paste0("Heavy Chain V:J Gene Pairings in ", exp_group_A), cex.main = 2)
  dev.off()

VJ_pairs_HC_chord_diagram_exp_group_B <- HC_gene_means_df %>%
  dplyr::filter(type == "VJ_pair", group_ID == exp_group_B) %>%
  separate(gene, into = c("VH", "JH", sep = ":")) %>%
  group_by(VH, JH) %>%
  summarise(percent = percent_gene, .groups = 'drop')
  VJ_pairs_HC_chord_diagram_exp_group_B$percent[VJ_pairs_HC_chord_diagram_exp_group_B$percent == 0] <- 1e-6
  VJ_pairs_HC_chord_diagram_exp_group_B$VH <- factor(VJ_pairs_HC_chord_diagram_exp_group_B$VH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_exp_group_B$VH)))
  VJ_pairs_HC_chord_diagram_exp_group_B$JH <- factor(VJ_pairs_HC_chord_diagram_exp_group_B$JH, levels = 
                                              mixedsort(unique(VJ_pairs_HC_chord_diagram_exp_group_B$JH)))
  circos.clear()
  all_sectors <- c(levels(VJ_pairs_HC_chord_diagram_exp_group_B$VH),
                   levels(VJ_pairs_HC_chord_diagram_exp_group_B$JH))
  gap_vec <- c(
                #big gap between VH and JH
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp_group_B$VH))-1), 15,
                #closing gap
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp_group_B$JH)) -1), 15)
  circos.par(gap.after = gap_vec)
  png(filename = file.path(plots_output_dir, paste0("HC_VJ_pairs_chord_diagram_", exp_group_B, ".png")), 
      width = 800, height = 800)
  par(mar = c(0, 0, 3, 0))
  chordDiagram(
    VJ_pairs_HC_chord_diagram_exp_group_B, 
    order = c(levels(VJ_pairs_HC_chord_diagram_exp_group_B$VH),
              levels(VJ_pairs_HC_chord_diagram_exp_group_B$JH)),
    link.sort = TRUE, 
    grid.col = grid.col, 
    reduce = 0,
    col = adjustcolor(grid.col[VJ_pairs_HC_chord_diagram_exp_group_B$VH], alpha.f = 0.5),
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
        cex = 2.0)
    },
    bg.border = NA)
  title(paste0("Heavy Chain V:J Gene Pairings in ", exp_group_B), cex.main = 2)
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
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp_group_A$VH))-1), 15,
                #closing gap
                rep(5, length(levels(VJ_pairs_HC_chord_diagram_exp_group_A$JH)) -1), 15)
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
        cex = 2.0)
    },
    bg.border = NA)
  title(paste0("Heavy Chain V:J Gene Pairings in ", ctrl_group), cex.main = 2)
  dev.off()


##CSV output
write_csv(summary_tsv_NGS, file.path(csv_output_dir, paste0(project, "summary_tsv.csv")))
write_csv(HC_master_df, file.path(csv_output_dir, paste0(project, "HC_df.csv")))
write_csv(HC_V_call_df, file.path(csv_output_dir, paste0(project, "HC_v_call_df.csv")))
write_csv(HC_V_call_means_df, file.path(csv_output_dir, paste0(project, "HC_v_call_means_df.csv")))
write_csv(HC_VJ_call_df, file.path(csv_output_dir, paste0(project, "_HC_VJ_call_df.csv")))
write_csv(HC_VJ_call_means_df, file.path(csv_output_dir, paste0(project, "_HC_VJ_call_means_df.csv")))

write_csv(HC_VH4_J_call_df, file.path(csv_output_dir, paste0(project, "_HC_VH4_J_call_df.csv")))
write_csv(HC_VH4_J_call_means_df, file.path(csv_output_dir, paste0(project, "_HC_VH4_J_call_means_df.csv")))

write_csv(HC_summary_df, file.path(csv_output_dir, 
            paste0(project, "_HeavyChain_VJ_summarydata_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".csv")))
write_csv(HC_gene_means_df, file.path(csv_output_dir, 
            paste0(project, "_HeavyChain_VJ_meandata_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cdr3_avgBR_charge, file.path(csv_output_dir, paste0(project, "_HeavyChain_cdr3_avgBR_charge_", exp_group, "_vs_", ctrl_group, ".csv")))

write_csv(HC_summary_wide_df, file.path(csv_output_dir, 
          paste0(project, "_HeavyChain_summary_widedata_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".csv")))
write_csv(HC_means_wide_df, file.path(csv_output_dir, 
          paste0(project, "_HeavyChain_means_widedata_", exp_group_A, "_vs_", exp_group_B, "_vs_", ctrl_group, ".csv")))
write_csv(HC_percent_stdev_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_percent_stdev_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_percent_pvalue_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_percent_pvalue_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cdr3gene_stdev_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_cdr3_stdev_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cdr3gene_pvalue_wide_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_cdr3_pvalue_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))
write_csv(HC_cell_counts_df, file.path(csv_output_dir, paste0(project, "_HeavyChain_cellcount_widedata_", exp_group, "_vs_", ctrl_group, ".csv")))