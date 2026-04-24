# packages
library(tidyverse)

# load the .tsv taxa assigned file paths
f_paths <- list.files(path = "fastq_trim_midori_kraken_conf_0-1_excel/",
                      pattern = ".xlsx$",
                      full.names = T)

# apply read_xlsx to all files at their paths, column types are made to default to character type.
# Also adds a column for sample id to the df, pulls the number from the input file name itself as it applies to each file

data_list <- lapply(f_paths, function(f) {
  df <- read_xlsx(f,
                  skip = 1,
                  range = cell_cols(1:5),
                  col_names = T,
                  col_types = NULL,
                  trim_ws = T)
  df$sample_id <- gsub("\\.xlsx","",basename(f)) # Extracts just the number
  return(df)
})

# bind all files read in as one df
bound_df <- do.call(rbind, data_list)
View(bound_df)

# check all possible rank codes
unique(bound_df$rank)

# diet spp. kraken df

d_kr_df <- bound_df %>%
  filter(rank =="S") %>%
  group_by(taxID) %>%
  mutate(taxon_count = n())

# running a taxid match with taxonomizr package

library(taxonomizr)

# increase timeout cutoff to an hour and use https protocol, connection keeps failing
options(timeout = 3600)

#prepare local database
prepareDatabase(sqlFile = "testdb.sqlite",
                getAccessions = F,
                protocol = 'https',
                vocal = T)

# isolate taxID

d_kr_taxid <- unique(d_kr_df$taxID)

# get taxa names for the tax ids
d_kr_taxa <- getTaxonomy(ids = d_kr_taxid,
                        sqlFile = "testdb.sqlite",
                        desiredTaxa = c("phylum","order","genus","species"),
                        getNames = T)

# convert the output to a df
d_kr_taxa_matched <- d_kr_taxa %>%
  as.data.frame() %>%
  rownames_to_column(var = "taxid") %>%
  mutate(taxid = str_trim(taxid))

# change taxID to numeric
d_kr_taxa_matched$taxid <- as.numeric(d_kr_taxa_matched$taxid)

unique(d_kr_taxa_matched$phylum)
unique(d_kr_taxa_matched$species)

# checks for NA entries
sum(is.na(d_kr_taxa_matched$taxid))

sum(is.na(d_kr_taxa_matched$species))
sum(is.na(d_kr_taxa_matched$genus))

sum(is.na(d_kr_taxa_matched$order))
sum(is.na(d_kr_taxa_matched$phylum))


# join assigned spp names with taxid key to df with sample id and other details

d_kr_tax_join <- full_join(d_kr_df,
                           d_kr_taxa_matched,
                           by = join_by(taxID == taxid),
                           relationship = "many-to-one")


# filter taxa name assigned df for arthropods
d_kr_taxa_named <- d_kr_tax_join %>%
  filter(phylum == "Arthropoda")

# total numner of unique sample ids
length(unique(d_kr_taxa_named$sample_id))

# create concise table for the species and its freq of occ data
d_kr_taxa_foo <- d_kr_taxa_named %>%
  ungroup() %>%
  select(species,taxon_count) %>%
  unique() %>%
  arrange(desc(taxon_count))

write_csv(d_kr_taxa_foo,"diet_kraken_spp_freq.csv")
