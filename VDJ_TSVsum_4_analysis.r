##SCRIPT to convert VDJserver tsv output files to readable format for R analysis of heavy and light chains 
library(readr)
library(dplyr)
library(stringr)

BR_code <- c(
             "2405", "2443", "2778", "2929", "2989", "3094", "3851", "5873",
             "6527", "20045", "1079", "1215", "1299", "1400", "1455", "1468",
             "1470", "1507", "1515", "1533", "1551", "1602", "1763", "1767",
             "1783", "1792", "1823", "1839", "1924", "1956", "1963", "1982",
             "2035", "2133", "2140", "2149", "2161", "2162", "2177", "2180",
             "2227", "2265", "HD279", "UTSW0001", "UTSW0002", "UTSW0003",
             "UTSW0009", "UTSW0013", "UTSW0014", "UTSW0015", "UTSW0017",
             "UTSW0018", "UTSW0019", "UTSW0020", "UTSW0022", "UTSW0023",
             "UTSW0033")

tsv_dir <- "/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/summary_files/"
output_dir <- "/home/morganjackson/bioinformatics/data/vdj_outputs/VDJserver_JCVIsamples_20250519/output/"

#list all tsv files in the directory 
tsv_files <- list.files(path = tsv_dir, pattern = "\\.tsv$", full.names = TRUE)

#read and combine all tsv files into df, add source_file column to identify original file 
summary_tsv_all <- bind_rows(
    lapply(tsv_files, function(file) {
      read_tsv(file, col_names = TRUE) %>%
    mutate(source_file = basename(file))
  }))

#convert productive column from logical value to character string
#filter only productive sequences
#extract BR code from sequence_id column, add as a new column 
#
summary_tsv_all <- summary_tsv_all %>%
    mutate(productive = as.character(productive)) %>%
    filter(productive == "TRUE") %>%
    mutate(BR_code = str_extract(sequence_id, paste0(BR_code, collapse = "|"))) %>%
    mutate(sample_type = case_when(
      str_detect(sequence_id, "PBMC") ~ "PBMC",
      str_detect(sequence_id, "CSF") ~ "CSF",
      str_detect(sequence_id, "Bld") ~ "Blood",
      TRUE ~ "NA")) %>%
    mutate(chain = case_when(
      str_detect(v_call, "IGH") ~ "heavy",
      str_detect(v_call, "IGK|IGL") ~ "light",
      TRUE ~ "NA")) %>%
    mutate(LC_isotype = case_when(
      str_detect(v_call, "IGKV") ~ "kappa",
      str_detect(v_call, "IGLV") ~ "lambda",
      TRUE ~ "NA")) %>%
    relocate(BR_code, chain, sample_type, .before = v_call) %>%
    relocate(source_file, .after = sequence_id) %>%
    relocate(LC_isotype, .after = chain)
