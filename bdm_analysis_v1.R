# Packages ----
library(pacman) ; library(tidyverse)
p_load("BiocManager","dada2","Biostrings")

# below doesnt work, memory allocation fails

# # Input data ----
# midori_db <- "midori_datasets/MIDORI2_UNIQ_NUC_SP_GB269_CO1_DADA2.fasta"
# all_otu_fasta <- "bdm_otus/all_otus.fasta"
#
# # Reading in sequences ----
# dna <- readDNAStringSet(all_otu_fasta)
#
# # Taxonomy assignment ----
# taxa <- assignTaxonomy(seqs = as.character(dna),
#                        refFasta = midori_db, # reference database for matching
#                        multithread = T, # 20 threads for my machine, can set this as TRUE for automatic multithreading
#                        tryRC = T,
#                        minBoot = 50, # 50 standard for COI
#                        outputBootstraps = T,
#                        verbose = T)
#
# # Convert to a clean data frame ----
# taxa_df <- as.data.frame(taxa) %>%
#   rownames_to_column("Sequence")
#
# # memory allocation issue on default run, spitting into chunks first
# n_chunks <- ceiling(length(dna)/25)
# dna_chunks <- split(dna, rep(1:n_chunks, each = 25, length.out = length(dna)))
#
# # Run assignment in a loop and combine results
# taxa_list <- lapply(dna_chunks, function(x) {
#   assignTaxonomy(as.character(x),
#                  midori_db,
#                  multithread = T,
#                  tryRC = T,
#                  minBoot = 50)
# })
#
# # Merge chunks back into one matrix
# taxa <- do.call(rbind, taxa_list)


# Blast query ran through terminal using midori reference database----

# Data formatting

# input file and col names
blast_cols <- c("OTU_ID", "Subject_ID", "Perc_Identity", "Alignment_Length",
                "Mismatches", "Gaps", "Q_Start", "Q_End", "S_Start", "S_End",
                "E_Value", "Bit_Score")

blast_res <- read_csv("midori_datasets/bat_diet_blast_res_v3.csv", col_names = blast_cols)

# Best hit selection for each OTU (Highest Bit_Score and Perc_Identity)
best_hits <- blast_res %>%
  group_by(OTU_ID) %>%
  slice_max(order_by = tibble(Bit_Score, Perc_Identity), n = 1, with_ties = FALSE) %>%
  ungroup()

# Reformat taxonomy file for readability
taxa <- best_hits %>%
  separate(Subject_ID,
           into = c("TaxID", "Root", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right")

View(taxa)

# clear out non essential cols
taxa_table <- taxa %>%
  select(OTU_ID, TaxID, Order, Family, Perc_Identity)

View(taxa_table)

# filter for Arthropoda, remove otus matched to bats, fungi, a mollusc as well

diet_orders <- taxa_table %>%
  filter(TaxID == "Arthropoda_6656") %>%
  group_by(Order) %>%
  mutate(count = n()) %>%
  select(Order,count) %>% unique()

View(diet_orders)





