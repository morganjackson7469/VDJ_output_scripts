###LOAD LIBRARIES
#%%
#load required libraries
library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(tidyr)
library(readr)
library(vscDebugger)

#read excel file into tibble
vdj_pairmaster_jcvi <- read_excel("/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/vdj_pairmaster.xlsx")

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

csv_output_dir <- "/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/csv_outputs"
plots_output_dir <- "/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/ggplot_outputs_individual"

HC_master_df <- tibble()


##FUNCTIONS
#%%
#filter for BR 1763, 2265, 3851
cysloop_controls <- vdj_pairmaster_jcvi %>%
  filter(BR_code %in% c("3094", "6527", "5873"))

#function to count gene family usage
count_gene_family <- function(df, column, patterns) {
  sapply(patterns, function(p) sum(grepl(p, df[[column]])))
}

#Top 10 gene distribution
top_genes <- function(df, column_name = "HC_v_call", top_n = 10,
                      filter_pattern = NULL) {
  if (!is.null(filter_pattern)) {
    df <- df %>% filter(str_detect(.data[[column_name]], filter_pattern))
  }
  df %>%
    count(.data[[column_name]], sort = TRUE) %>%
    slice_head(n = top_n)
}

#Count V:J gene pairings function
cysloop_controls <- cysloop_controls %>%
  mutate(VJ_pair = paste0(HC_v_call, ":", HC_j_call)) %>%
  mutate(
    V_gene_mut = str_extract(HC_v_call, "IGHV(\\d+)") %>% str_replace("IGHV", "VH"),
    J_gene_mut = str_extract(HC_j_call, "IGHJ(\\d+)") %>% str_replace("IGHJ", "JH"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)
  )

count_vj_pairs <- function(df, v_pattern, j_pattern, v_column = "HC_v_call", j_column = "HC_j_call") {
  df %>%
    filter(str_detect(.data[[v_column]], v_pattern) & str_detect(.data[[j_column]], j_pattern)) %>%
    nrow()
}

percent_vj_pairs <- function(df, vj_pattern, vj_column = "VJ_only_gene") {
  total <- nrow(df)
  count <- df %>%
    filter(str_detect(.data[[vj_column]], vj_pattern)) %>%
    nrow()
  percent <- (count / total) * 100
  return(percent)
}


##HEAVY CHAIN ANALYSIS
#%%

#IGHV and IGHJ families
IGHV_families <- paste0("IGHV", 1:7)
IGHJ_families <- paste0("IGHJ", 1:7)

#count V and J gene families
IGHV_counts <- count_gene_family(cysloop_controls, "HC_v_call", IGHV_families)
IGHJ_counts <- count_gene_family(cysloop_controls, "HC_j_call", IGHJ_families)

#create distribution tibbles
V_fam_distribution <- tibble(HC_v_fam = IGHV_families, count = IGHV_counts) %>%
mutate(distribution = 100 * count / sum(count))

J_fam_distribution <- tibble(HC_j_fam = IGHJ_families, count = IGHJ_counts) %>%
mutate(distribution = 100 * count / sum(count))

#HC gene distribution
HC_v_gene_distribution <- cysloop_controls %>%
filter(str_detect(HC_v_call, "IGHV"))

HC_j_gene_distibution <- cysloop_controls %>%
filter(str_detect(HC_j_call, "IGHJ"))

#get top 10 V and J gene hits 
top10_v_genes <- top_genes(cysloop_controls, column_name = "HC_v_call")
top10_j_genes <- top_genes(cysloop_controls, column_name = "HC_j_call")

#V:J pairs
top_vj_pairings <- cysloop_controls %>%
  count(VJ_only_gene, sort = TRUE) %>%
  slice_max(n, n = 10)

heatmap_vj_pairs <- cysloop_controls %>%
  count(V_gene_mut, J_gene_mut)

vj_pair_percentage_table <- cysloop_controls %>%
  count(VJ_only_gene) %>%
  mutate(percent = (n / sum(n)) * 100) %>%
  slice_max(percent, n = 10)

#this can call the specific percentage 
percent_vj_pairs(cysloop_controls, "VH3:JH3")
#if want to show only top 10 percentage use slice_max(percent, n = 10)



##PLOT AND SAVE V & J
#%%
# IGHV Family Distribution
ggplot(V_fam_distribution, aes(x = HC_v_fam, y = distribution)) +
geom_bar(stat = "identity", fill = "steelblue") +
theme_grey(base_size = 14) +
labs(title = "IGHV Family Distribution in Cysloop",
x = "V Gene Family", y = "Percentage (%)")
ggsave(filename = file.path(plots_output_dir, "V_fam_distribution_cysloop.png"),
       width = 6, height = 4)

#IGHV top 10 genes 
ggplot(top10_v_genes, aes(x = reorder(HC_v_call, -n), y = n)) +
geom_bar(stat = "identity", fill = "steelblue") +
theme_grey(base_size = 14) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "IGHV Top 10 Genes in Cysloop",
x = "V Top 10 Genes", y = "Count")
ggsave(filename = file.path(plots_output_dir, "V_genes_top10_cysloop.png"),
       width = 6, height = 4)

# IGHJ Family Distribution
ggplot(J_fam_distribution, aes(x = HC_j_fam, y = distribution)) +
geom_bar(stat = "identity", fill = "steelblue") +
theme_grey(base_size = 14) +
labs(title = "IGHJ Family Distribution in Cysloop",
x = "J Gene Family", y = "Percentage (%)")
ggsave(filename = file.path(plots_output_dir, "J_fam_distribution_cysloop.png"),
       width = 6, height = 4)

#IGHV top 10 genes 
ggplot(top10_j_genes, aes(x = reorder(HC_j_call, -n), y = n)) +
geom_bar(stat = "identity", fill = "steelblue") +
theme_grey(base_size = 14) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "IGHJ Top 10 Genes in Cysloop",
x = "J Top 10 Genes", y = "Count")
ggsave(filename = file.path(plots_output_dir, "J_genes_top10_cysloop.png"),
       width = 6, height = 4)

#V:J pairs
ggplot(top_vj_pairings, aes(x = reorder(VJ_only_gene, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 V:J pairs in Cysloop",
       x = "V:J pairs", y = "Count")
ggsave(filename = file.path(plots_output_dir, "V:J_pairs_top10_cysloop.png"),
       width = 6, height = 4)

ggplot(heatmap_vj_pairs, aes(x = J_gene_mut, y = V_gene_mut, fill = n)) +
  geom_tile(color = "#ffffff60") +
  scale_fill_gradient(low = "#edcdcde4", high = "darkred") +
  labs(title = "VH:JH Gene Pairings in Cysloop",
       x = "JH Gene", y = "VH Gene", fill = "Count") +
  theme_grey()
ggsave(filename = file.path(plots_output_dir, "VH:JH_pairs_cysloop.png"),
        width = 6, height = 4)

#V:J pairs percentage
ggplot(vj_pair_percentage_table, aes(x = reorder(VJ_only_gene, -percent), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "V:J pairs percentage in Cysloop",
       x = "V:J pairs", y = "Percentage")
ggsave(filename = file.path(plots_output_dir, "V:J_pairs_percentage_cysloop.png"),
       width = 6, height = 4)




#%%
##CSV output
selected_codes <- c("3094", "6527", "5873")

for (code in selected_codes) {
  csv_df <- cysloop_controls %>%
    filter(BR_code == code)


IGHV_counts <- count_gene_family(csv_df, "HC_v_call", IGHV_families)
IGHJ_counts <- count_gene_family(csv_df, "HC_j_call", IGHJ_families)

V_fam_df <- tibble(
    BR_code = code,
    Gene_family = IGHV_families,
    Type = "IGHV",
    count = IGHV_counts,
    distribution = 100 * count / sum(count)
  )

J_fam_df <- tibble(
  BR_code = code,
  Gene_family = IGHJ_families,
  Type = "IGHJ",
  count = IGHJ_counts,
  distribution = 100 * count / sum(count)
)

write_csv(V_fam_df, file.path(csv_output_dir, paste0("IGHV_family_counts_", code, ".csv")))
write_csv(J_fam_df, file.path(csv_output_dir, paste0("IGHJ_family_counts_", code, ".csv")))
HC_master_df <- bind_rows(HC_master_df, V_fam_df, J_fam_df)
}

