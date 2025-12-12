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
             "UTSW0033", "6082",
             "1455", "2177",	"2180",	"992",	"1823",	"1123",	"1299",	"1490",	"1507",	
             "1844",	"2133",	"1574",	"3500", "J10",	"J130",	"J20",	"J201",	"J203",
             "J218",	"J220",	"J24",	"J26",	"J34",	"J42",	"J46",	"J8",	"J9",	"J91",
             "J93", "6634",	"1559",	"1561",	"1607",	"1618",	"1629",	"1647",	"1690",	
             "1702",	"1707",	"1713",	"1724",	"1790", "1831",	"1837",	"1848",	"1850",	
             "6608",	"6611",	"6617",	"32",	"6624",	"6658",	"6670",	"6691",	"6696",	"6714",	
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
              "USCHC005",	"USCH043",	"USCH020",	"USCH021",	"USCH032"
)

tsv_dir <- "/mnt/md0/Projects/MonsonLab/MAV/summary"

#list all tsv files in the directory 
tsv_files <- list.files(path = tsv_dir, pattern = "\\.tsv$", full.names = TRUE)

#read and combine all tsv files into df, add source_file column to identify original file 
summary_tsv_NGS <- bind_rows(
    lapply(tsv_files, function(file) {
      read_tsv(file, col_names = TRUE) %>%
    mutate(source_file = basename(file))
  }))

#convert productive column from logical value to character string
#filter only productive sequences
#extract BR code from sequence_id column, add as a new column 
#prepare safe ordered regex from BR_code vector
codes <- unique(BR_code)

#sort by descending length to prefer longer matched (J91 over J9)
codes <- codes[order(nchar(codes), decreasing = TRUE)]

#escape regex metacharacters in codes
codes_esc <- str_replace_all(codes, "([\\^$.|?*+()\\[\\]{}\\\\-])", "\\\\\\1")

#case-insensitive pattern: start code after non-alpha num char
#does not require seperator (2177ASC)
pattern <- paste0("(?i)(?:^|[^A-Za-z0-9])(", paste(codes_esc, collapse = "|"), ")")


summary_tsv_NGS <- summary_tsv_NGS %>%
    mutate(productive = as.character(productive)) %>%
    filter(productive == "TRUE") %>%
    mutate(BR_code = {
        m <- str_match(source_file, pattern)
        m[,2]}) %>%
    mutate(BR_code_matched = !is.na(BR_code)) %>%
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

print("finished creating summary dataframe")