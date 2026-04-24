# packages
library(tidyverse)

# load the .tsv taxa assigned file paths
f_paths <- list.files(path = "bdm_asv_blast_tax/",
                      pattern = ".tsv$",
                      full.names = T)

# apply read_tsv to all files at their paths, column types are made to default to character type.
# Also adds a column for sample id to the df, pulls the number from the input file name itself as it applies to each file

data_list <- lapply(f_paths, function(f) {
  df <- read_tsv(f,
                 col_names = F,
                 col_types = cols(.default = "c"),
                 show_col_types = F)
  df$sample_id <- gsub("\\.tsv", "", basename(f)) # Extracts just the number
  return(df)
})

# bind all files read in as one df
bound_df <- do.call(rbind, data_list)
View(bound_df)

# replace | to ; separator
bound_df$X2 <- gsub("|",";",bound_df$X2,fixed = T)

# separate columns at ";" as a separator
diet_df <- bound_df %>%
  separate_wider_delim(
    X2,
    delim = ";",
    names_sep = "_",
    too_few = "align_start") %>% # fill empty right-hand columns with NA
  separate_wider_delim(
    X1,
    delim = ";",
    names_sep = "_",
    too_few = "align_start") %>%
  separate_wider_delim(
    X9,
    delim = ";",
    names_sep = "_",
    too_few = "align_start")

View(diet_df)

# remove size label from the size column and rename column to size
diet_df$X1_2 <- gsub("size=","",diet_df$X1_2,fixed = T)

# rename size and taxid column
diet_df <- diet_df %>%
  rename(size = X1_2,
         taxid = X2_3)

# drop column that has taxid label
diet_df <- diet_df %>% select(-X2_2)

# check if certain columns have only NAs
unique(diet_df$X7)
unique(diet_df$X8)

# drop the two uninformative cols
diet_df <- diet_df %>%
  select(-X7,-X8)

# check more cols
unique(diet_df$X9_31)
unique(diet_df$X9_32)
unique(diet_df$X9_33)

# few species names show up, retained

# move up sample id column
diet_df <- diet_df %>% relocate(sample_id,taxid)

# rename more cols
diet_df <- diet_df %>%
  rename(seq = X2_1,
         p_ident = X3,
         length = X4,
         evalue = X5,
         bitscore = X6)

diet_df <- diet_df %>% relocate(taxid, .after = seq)

# seq id?
diet_df <- diet_df %>%
  rename(seq_id = X1_1)

diet_df <- diet_df %>%
  rename(qseq_id = seq_id)

# pick out the best of the two hits for each query sequence
# arrange by least evalue, highest bitscore and highest p-ident and select the first row

d_hits <- diet_df %>%
  group_by(qseq_id) %>%
  arrange(evalue, desc(bitscore), desc(p_ident)) %>%
  slice(1) %>%
  ungroup()

View(d_hits)

# check that for every query, one out of two records are retained in the output df - below two lines should match in number
nrow(diet_df)/2
nrow(d_hits)

d_hits_filter <- d_hits %>%
  select(-(X9_1:X9_9))

# check types of taxa at order column - some formats are different than the most, need filtration here
unique(d_hits_filter$X9_10)

d_hits_clean <- d_hits_filter %>%
  filter(X9_10 == "Vertebrata" | X9_10 == "Arthropoda" | X9_10 == "Mollusca")

# running a taxid match with taxonomizr

library(taxonomizr)

# increase timeout cutoff to an hour and use https protocol, connection keeps failing
options(timeout = 3600)

#prepare local database
prepareDatabase(sqlFile = "testdb.sqlite",
                getAccessions = F,
                protocol = 'https',
                vocal = T)

# check what all arthropods are present in each sample and how many hits each got

arth_df <- d_hits_clean %>%
  filter(X9_10 == "Arthropoda") %>%
  group_by(taxid,sample_id) %>%
  mutate(count=n()) %>%
  select(sample_id,taxid,count,size) %>%
  arrange(sample_id,taxid) %>%
  ungroup() %>%
  group_by(sample_id,taxid) %>%
  mutate(size = as.numeric(size)) %>%
  mutate(total_seq_count = sum(size)) %>%
  select(-size) %>%
  unique()

# list out unique taxid's of arthropods
arth_ids <- unique(arth_df$taxid)

str(arth_ids)

# get taxa names for the arth tax ids
arth_spp <- getTaxonomy(ids = arth_ids,
            sqlFile = "testdb.sqlite",
            desiredTaxa = "species",
            getNames = T)

# write spp file
# write_csv(as.data.frame(arth_spp), file = "arth_sp_list.csv")

# convert the output to a df

arth_spp_df <- arth_spp %>%
  as.data.frame() %>%
  rownames_to_column(var = "taxid") %>%
  mutate(taxid = str_trim(taxid))

# join assigned spp names to with taxid key to df with sample id and seq hit details

arth_tax_join <- full_join(arth_df,
                           arth_spp_df,
                           by = "taxid",
                           relationship = "many-to-one")

# total numner of samples each species has been detected in, freq of occurrence, sp_foo
arth_taxa <- arth_tax_join %>%
  group_by(taxid) %>%
  mutate(sp_foo = n())

length(unique(arth_taxa$sample_id))

length(unique(d_hits_clean$sample_id))

length(unique(arth_df$sample_id))

length(unique(d_hits_filter$sample_id))

length(unique(d_hits$sample_id))

arth_sample_count_test <- d_hits_clean %>%
  filter(X9_10 == "Arthropoda")

length(unique(arth_sample_count_test$sample_id))
