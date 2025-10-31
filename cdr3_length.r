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

individual_outputs_dir <- ("/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/individual_outputs")

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

cdr3_length_files <- list.files(path = individual_outputs_dir, pattern = "\\.igblast\\.junction_aa_wo_duplicates_length\\.tsv$", full.names = TRUE)

file_match <- summary_tsv_all %>%
  mutate(group_ID = case_when(
    BR_code %in% experimental_group ~ exp_group,
    BR_code %in% control_group ~ ctrl_group),
    .after = 1) %>%
   filter(group_ID %in% c(exp_group, ctrl_group)) %>%
   group_by(BR_code) %>%
   slice_head(n = 1) %>%
   select(BR_code, group_ID, source_file) %>%
   ungroup() %>%
   mutate(
    match_key = str_match(source_file, "^\\d+\\.([^.]+)\\.")[,2])

cdr3_matched_files <- cdr3_length_files[
    grep(paste(file_match$match_key, collapse = "|"), cdr3_length_files)]

cdr3_length_df <- purrr::map_df(
    cdr3_matched_files, ~{
        readr::read_tsv(.x) %>%
    mutate(source_file = basename(.x)) %>%
    mutate(
        match_key = str_match(source_file, "^[^.]+")) %>%
    left_join(file_match %>%
        select(match_key, BR_code, group_ID), by = "match_key") %>%
    select(-source_file, -match_key) %>%
    mutate(CDR3_percent = CDR3_RELATIVE * 100)
    }
)

cdr3_length_avg_df <- cdr3_length_df %>%
    group_by(group_ID, CDR3_LENGTH) %>%
    summarize(
        cdr3_avg_count = mean(CDR3_COUNT)) %>%
    ungroup() %>%
    group_by(group_ID) %>%
    mutate(
        cdr3_total_count = sum(cdr3_avg_count),
        cdr3_avg_percent = ifelse(cdr3_avg_count == 0, 0, ((cdr3_avg_count / cdr3_total_count) * 100)
    )
)



#PLOTS

cdr3_length_counts_plot <- cdr3_length_df %>%
    ggplot() +
    geom_line(
        aes(x = CDR3_LENGTH, y = CDR3_COUNT, colour = BR_code, linetype = group_ID)) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "CDR3 Lengths",
       x = "CDR3 Length", y = "CDR3 Count")
    ggsave(filename = file.path(plots_output_dir, paste0("CDR3_length_counts_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = cdr3_length_counts_plot, width = 8, height = 4)

cdr3_length_percent_plot <- cdr3_length_df %>%
    ggplot() +
    geom_line(
        aes(x = CDR3_LENGTH, y = CDR3_percent, colour = BR_code, linetype = group_ID)) +
    theme_classic(base_size = 14) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "CDR3 Lengths",
       x = "CDR3 Length", y = "CDR3 percent")
    ggsave(filename = file.path(plots_output_dir, paste0("CDR3_length_percent_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = cdr3_length_percent_plot, width = 8, height = 4)

cdr3_avg_length_counts_plot <- cdr3_length_avg_df %>%
    ggplot() +
    geom_line(
        aes(x = CDR3_LENGTH, y = cdr3_avg_count, group = group_ID, colour = group_ID)) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average CDR3 Lengths",
       x = "CDR3 Length", y = "Avg CDR3 Count")
    ggsave(filename = file.path(plots_output_dir, paste0("CDR3_avg_length_counts_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = cdr3_avg_length_counts_plot, width = 8, height = 4)

cdr3_avg_length_percent_plot <- cdr3_length_avg_df %>%
    ggplot() +
    geom_line(
        aes(x = CDR3_LENGTH, y = cdr3_avg_percent, group = group_ID, colour = group_ID)) +
    theme_classic(base_size = 14) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average CDR3 Lengths",
       x = "CDR3 Length", y = "Avg CDR3 Percent")
    ggsave(filename = file.path(plots_output_dir, paste0("CDR3_avg_length_percent_", exp_group, "_vs_", ctrl_group, ".png")),
        plot = cdr3_avg_length_percent_plot, width = 8, height = 4)