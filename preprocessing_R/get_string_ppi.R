suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(STRINGdb))
suppressMessages(library(igraph))

#### Load the STRING database including physical interactions only and prep for downstream PCSF

#First get the STRING network
string_db <- STRINGdb$new(version="11.5", species=9606, network_type = 'physical', score_threshold = 1)
full_string_db <- string_db$load_all()

string_db_df <- as_data_frame(full_string_db)

#Format the PPI for downstream PCSF analysis
ppiPrepped <- data.frame(
  "protein1" = string_db_df$from,
  "protein2" = string_db_df$to,
  "cost" = 1 - (string_db_df$combined_score / 1000)
)

setwd("~/Bioinformatics/Network Analysis of HFpEF/Hahn_FGN_Analysis/ref")
write.table(ppiPrepped, 'prepped_STRING_9606.protein.links.v11.5.txt', row.names = FALSE, sep = '\t',  quote = FALSE)