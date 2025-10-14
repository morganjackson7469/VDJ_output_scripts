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
LC_master_df <- tibble()


##FUNCTIONS
#%%
#filter for BR 1763, 2265, 3851
cysloop <- vdj_pairmaster_jcvi %>%
  filter(BR_code %in% c("3094", "6527", "5873"))

#function to count gene family usage
count_gene_family <- function(df, column, patterns) {
  sapply(patterns, function(p) sum(grepl(p, df[[column]])))
}

#Top 10 gene distribution
top_genes <- function(df, column_name = "LC_v_call", top_n = 10,
                      filter_pattern = NULL) {
  if (!is.null(filter_pattern)) {
    df <- df %>% filter(str_detect(.data[[column_name]], filter_pattern))
  }
  df %>%
    count(.data[[column_name]], sort = TRUE) %>%
    slice_head(n = top_n)
}

#Count V:J gene pairings function
cysloop <- cysloop %>%
  mutate(VJ_pair = paste0(LC_v_call, ":", LC_j_call)) %>%
  mutate(
    V_gene_mut = str_extract(LC_v_call, "IGKV(\\d+)|IGLV(\\d+)") %>%
      str_replace("IGKV", "VK") %>%
      str_replace("IGLV", "VL"),
    J_gene_mut = str_extract(LC_j_call, "IGKJ(\\d+)|IGLJ") %>%
      str_replace("IGKJ", "JK") %>%
      str_replace("IGLJ", "LJ"),
    VJ_only_gene = paste0(V_gene_mut, ":", J_gene_mut)
  )

count_vj_pairs <- function(df, v_pattern, j_pattern, v_column = "LC_v_call", j_column = "LC_j_call") {
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


###LIGHT CHAIN ANALYSIS
#%%

IGKV_families <- paste0("IGKV", 1:7)
IGLV_families <- paste0("IGLV", 1:7)

#count V and J gene families
IGKV_counts <- count_gene_family(cysloop_controls, "LC_v_call", IGKV_families)
IGLV_counts <- count_gene_family(cysloop_controls, "LC_v_call", IGLV_families)

#create distribution tibbles
K_fam_distribution <- tibble(KC_v_fam = IGKV_families, count = IGKV_counts) %>%
  mutate(distribution = 100 * count / sum(count))

L_fam_distribution <- tibble(LC_v_fam = IGLV_families, count = IGLV_counts) %>%
  mutate(distribution = 100 * count / sum(count))

#light chain gene distribution
KC_v_gene_distribution <- cysloop_controls %>%
  filter(str_detect(LC_v_call, "IGKV"))

LC_v_gene_distribution <- cysloop_controls %>%
  filter(str_detect(LC_v_call, "IGLV"))

#get top 10 K and L gene hits, remove allele name  
top10_k_genes <- top_genes(KC_v_gene_distribution, column_name = "LC_v_call") 
top10_k_genes <- top10_k_genes %>%
  mutate(gene_only = str_remove(LC_v_call, "\\*.*$"))

top10_l_genes <- top_genes(LC_v_gene_distribution, column_name = "LC_v_call") 
top10_l_genes <- top10_l_genes %>%
  mutate(gene_only = str_remove(LC_v_call, "\\*.*$"))

#get top 10 V and J gene hits 
top10_v_genes <- top_genes(cysloop_controls, column_name = "LC_v_call")
top10_j_genes <- top_genes(cysloop_controls, column_name = "LC_j_call")

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


###LIGHT CHAIN PLOTS
#%%
# IGKV Family Distribution
ggplot(K_fam_distribution, aes(x = KC_v_fam, y = distribution)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  labs(title = "IGKV Family Distribution in Cysloop",
       x = "LC Kappa Gene Family", y = "Percentage (%)")
ggsave(filename = file.path(plots_output_dir, "K_fam_distribution_cysloop.png"),
       width = 6, height = 4)

#IGKV top 10 genes 
ggplot(top10_k_genes, aes(x = reorder(gene_only, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IGKV Top 10 Genes in Cysloop",
       x = "KV Top 10 Genes", y = "Count")
ggsave(filename = file.path(plots_output_dir, "KV_genes_top10_cysloop.png"),
       width = 6, height = 4)

# IGLV Family Distribution
ggplot(L_fam_distribution, aes(x = LC_v_fam, y = distribution)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  labs(title = "IGLC Family Distribution in Cysloop",
       x = "LC Lambda Gene Family", y = "Percentage (%)")
ggsave(filename = file.path(plots_output_dir, "L_fam_distribution_cysloop.png"),
       width = 6, height = 4)

#IGLV top 10 genes 
ggplot(top10_l_genes, aes(x = reorder(gene_only, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IGLV Top 10 Genes in Cysloop",
       x = "LV Top 10 Genes", y = "Count")
ggsave(filename = file.path(plots_output_dir, "LV_genes_top10_cysloop.png"),
       width = 6, height = 4)

#V:J pairs
ggplot(top_vj_pairings, aes(x = reorder(VJ_only_gene, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Top 10 V:J pairs in Cysloop",
       x = "V:J pairs", y = "Count")
ggsave(filename = file.path(plots_output_dir, "LC_V:J_pairs_top10_cysloop.png"),
       width = 6, height = 4)

ggplot(heatmap_vj_pairs, aes(x = J_gene_mut, y = V_gene_mut, fill = n)) +
  geom_tile(color = "#ffffff60") +
  scale_fill_gradient(low = "#edcdcde4", high = "darkred") +
  labs(title = "VH:JH Gene Pairings in Cysloop",
       x = "JH Gene", y = "VH Gene", fill = "Count") +
  theme_grey()
ggsave(filename = file.path(plots_output_dir, "LC_VH:JH_pairs_cysloop.png"),
        width = 6, height = 4)

#V:J pairs percentage
ggplot(vj_pair_percentage_table, aes(x = reorder(VJ_only_gene, -percent), y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_grey(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "V:J pairs percentage in Cysloop",
       x = "V:J pairs", y = "Percentage")
ggsave(filename = file.path(plots_output_dir, "LC_V:J_pairs_percentage_cysloop.png"),
       width = 6, height = 4)

#%%
##CSV output
selected_codes <- c("3094", "6527", "5873")

for (code in selected_codes) {
  csv_df <- cysloop_controls %>%
    filter(BR_code == code)


IGHV_counts <- count_gene_family(csv_df, "LC_v_call", IGHV_families)
IGHJ_counts <- count_gene_family(csv_df, "LC_j_call", IGHJ_families)

V_fam_df <- tibble(
    BR_code = code,
    Gene_family = IGHV_families,
    Type = "IGV",
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

# Write master CSV
write_csv(HC_master_df, file.path(csv_output_dir, "HC_gene_family_counts.csv"))

##LIGHT CHAIN CSV
selected_codes <- c("3094", "6527", "5873")
for (code in selected_codes) {
  csv_df <- cysloop_controls %>%
    filter(BR_code == code)

IGKV_counts <- count_gene_family(csv_df, "LC_v_call", IGKV_families)
IGLV_counts <- count_gene_family(csv_df, "LC_v_call", IGLV_families)

KV_fam_df <- tibble(
    BR_code = code,
    Gene_family = IGKV_families,
    Type = "IGKV",
    count = IGKV_counts,
    distribution = 100 * count / sum(count)
  )

LV_fam_df <- tibble(
  BR_code = code,
  Gene_family = IGLV_families,
  Type = "IGLV",
  count = IGLV_counts,
  distribution = 100 * count / sum(count)
)

write_csv(KV_fam_df, file.path(csv_output_dir, paste0("IGKV_family_counts_", code, ".csv")))
write_csv(LV_fam_df, file.path(csv_output_dir, paste0("IGLV_family_counts_", code, ".csv")))
LC_master_df <- bind_rows(LC_master_df, KV_fam_df, LV_fam_df)
}

# Write master CSV
write_csv(LC_master_df, file.path(csv_output_dir, "LC_gene_family_counts.csv"))