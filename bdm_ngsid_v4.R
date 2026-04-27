# packages
library(tidyverse)

# load the .tsv taxa assigned file paths
f_paths <- list.files(path = "bdm_ngsid_midori-blast_asv/",
                      pattern = ".tsv$",
                      full.names = T)

# Apply read_tsv to all files at their paths, column types are guessed
# Adds a column for sample id to the df, pulls the number from the input file name itself as it applies to each file

data_list <- lapply(f_paths, function(f) {
  df <- read_tsv(f,
                  col_select = 1:6,
                  col_names = F,
                  col_types = NULL,
                  trim_ws = T)
  df$sample_id <- gsub("\\.tsv","",basename(f)) # Extracts just the number
  return(df)
})

# bind all files read in as one df
bound_df <- do.call(rbind, data_list)
View(bound_df)

# Reformat columns
# Extract the read number from the end of X1
# Extract seq and taxid from X2
# Extract consensus cl id from X1
# Remove redundant cols, Rename cols, Reorder cols

ngsid_df_v1 <- bound_df %>%
  extract(X1,
          into = "sup_reads",
          regex = "_(\\d+)$",
          remove = F, convert = T) %>%
  extract(X2,
          into = c("seq_id", "taxid"),
          regex = "seq(\\d+)\\|taxid\\|(\\d+)",
          remove = F, convert = T) %>%
  extract(X1,
          into = "cl_id",
          regex = "consensus_cl_id_(\\d+)_",
          remove = F, convert = T) %>%
  select(-X1,-X2) %>%
  rename(p_ident = X3,
         length = X4,
         evalue = X5,
         bitscore = X6) %>%
 select(sample_id, taxid, seq_id, cl_id, sup_reads, p_ident, length, evalue, bitscore)


# check how many entries per cl ID
table(table(ngsid_df_v1$cl_id))

# arrange by increasing evalue, desc bitscore and desc percent identity
# This will organise best hit on top, then pick this top, first hit

ngsid_df_v2 <- ngsid_df_v1 %>%
  group_by(sample_id,cl_id) %>%
  arrange(evalue, desc(bitscore), desc(p_ident)) %>%
  slice(1) %>%
  ungroup()

# check number of rows for each cl ID
id_counts <- ngsid_df_v2 %>%
  count(cl_id) %>%
  arrange(desc(n))

print(id_counts)

# total number of unique cl ids
length(unique(ngsid_df_v2$cl_id))

# total number of sample ids
length(unique(ngsid_df_v2$sample_id))

# check if there are duplicate cl ids within each sample id
duplicates <- ngsid_df_v2 %>%
  group_by(sample_id, cl_id) %>%
  filter(n() > 1) %>%
  ungroup()

print(duplicates) # this should print 0 rows. There should be only one row of each cl id inside one sample id

# check how many rows are reduced
nrow(ngsid_df_v2)/nrow(ngsid_df_v1)

# get taxon names from taxID

# running a taxid match with taxonomizr package

library(taxonomizr)

# increase timeout cutoff to an hour and use https protocol, connection keeps failing
options(timeout = 3600)

#prepare local database
prepareDatabase(sqlFile = "testdb.sqlite",
                getAccessions = F,
                protocol = 'https',
                vocal = T)

# pull unique taxIDs from the df
taxids <- as.numeric(unique(ngsid_df_v2$taxid))

# get taxa names for the tax ids
ngsid_df_taxa <- getTaxonomy(ids = taxids,
                         sqlFile = "testdb.sqlite",
                         desiredTaxa = c("phylum","order","family","genus","species"),
                         getNames = T)

# update the output to a df
ngsid_taxa_reshape <- ngsid_df_taxa %>%
  as.data.frame() %>%
  rownames_to_column(var = "taxid") %>%
  mutate(taxid = str_trim(taxid))

# change taxid col type in taxa name df
ngsid_taxa_reshape$taxid <- as.numeric(ngsid_taxa_reshape$taxid)

# join dfs
ngsid_join <- full_join(ngsid_df_v2,
                        ngsid_taxa_reshape,
                           by = join_by(taxid == taxid),
                           relationship = "many-to-one")


# Some checks

# total arthropod records
sum(ngsid_join$phylum == "Arthropoda",na.rm = T)

# check number of arthropod hits with less than 5 reads
sum(ngsid_join$sup_reads < 5 & ngsid_join$phylum == "Arthropoda") # not too many, can discard

# unique phyla
unique(ngsid_join$phylum)

# Arthropod counts with greater than 90, 85 and 80% match
sum(ngsid_join$p_ident >= 90 & ngsid_join$phylum == "Arthropoda", na.rm = T)
sum(ngsid_join$p_ident >= 85 & ngsid_join$phylum == "Arthropoda", na.rm = T)
sum(ngsid_join$p_ident >= 80 & ngsid_join$phylum == "Arthropoda", na.rm = T)

# Very few records lost at 80% p ident
ngsid_join %>% filter(phylum == "Arthropoda") %>% pull(p_ident) %>% IQR(na.rm = T)


# Clean up non essential cols from df and count total supporting reads for each taxid
ngsid_v3 <- ngsid_join %>%
  filter(phylum == "Arthropoda", p_ident >= 80) %>% # tighten filter on taxon matches here
  group_by(taxid) %>%
  mutate(taxid_total_reads = sum(sup_reads)) %>%
  select(-(phylum:genus)) %>%
  ungroup()

#Within one sample id, there can be more than one cl_id pointing to the same taxid
# To correctly count sp_foo then, have to remove other variables, such as cl_id, then group by taxid and sample id
# This avoids duplicate rows for the same taxid/species within one sample

# This is done in the below chunk - we pick the one taxon hit with the best evalue, bitscore, p ident and support reads, in that order of priority
# sample IDs are retained here. When those are mapped to season, seasonal patterns can be analysed

ngsid_v4 <- ngsid_v3 %>%
  group_by(sample_id, taxid) %>%
  arrange(evalue, desc(bitscore), desc(p_ident), desc(sup_reads)) %>%
  slice(1) %>%
  ungroup()

# check if there are duplicate tax ids within each sample id
duplicates <- ngsid_v4 %>%
  group_by(sample_id, tax_id) %>%
  filter(n() > 1) %>%
  ungroup()

print(duplicates) # this should print 0 rows. There should be only one row of each tax id inside one sample id

# Now each sample id has one instance for every tax id. We can count total occurrences of each taxon
ngsid_v5 <- ngsid_v4 %>%
  select(sample_id,taxid,species,taxid_total_reads) %>%
  group_by(taxid) %>%
  mutate(sp_foo = n()) %>%
  ungroup()

# Simpler table with freq of occurrence for each arthropod diet species, sample_ids removed
sp_foo_no_sid <- ngsid_v5 %>%
  select(species,sp_foo, taxid_total_reads) %>%
  unique() %>%
  arrange(desc(sp_foo), desc(taxid_total_reads)) %>%
  mutate(sp_group = fct_lump_n(species, n = 10, w = sp_foo, other_level = "Others")) %>%
  filter(sp_group != "Others")

# Pest classification of diet spp.

sp_foo_pest_class <- sp_foo_no_sid %>%
  select(-sp_group) %>%
  mutate(sp_type = case_when(species == "Eysarcoris ventralis" ~ "P",
                             species == "Pyrilla perpusilla" ~ "P",
                             species == "Elasmolomus pallens" ~ "P",
                             species == "Bagrada hilaris" ~ "P",
                             species == "Odontotermes obesus" ~ "P",
                             species == "Nilaparvata lugens" ~ "P",
                             .default = "NP"))


# Plot prepping ====

# font setup
library(showtext)
showtext_auto()
font_add(family = "Bahnschrift", regular = "Bahnschrift.ttf")

# plot features
point_size = 4
pointtext_size = 5
pointtext_family = "Bahnschrift"

linecolor = "darkgrey"
linetype = "dashed"

labeltext_size = 5
labeltext_family = "Bahnschrift"
# graph axis labels

# check total number of samples
length(unique(ngsid_v5$sample_id))

x_lab = "Frequency of occurrence (FOO) [n = 125]"
y_lab = "Species (Top 10 by FOO)"

# plot element adjustments
plot_theme <- theme_classic(base_family = "Bahnschrift") +
  theme(axis.title = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.text.y = element_text(margin = margin(r = 5), ),
        legend.title = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 14, color = "black"))


# Plotting ----

ggplot(sp_foo_pest_class,
       aes(x = sp_foo, y = reorder(species,sp_foo))) +
  geom_segment(aes(xend = 0, yend = species),
               color = linecolor,
               linetype = linetype) +
  geom_point(aes(shape = sp_type, color = sp_type),
             size = point_size,
             show.legend = T) +
  scale_color_manual(name = "Species Type",
                     labels = c("Not Pest","Pest"),
                     values = c("lightgreen","red3")) +
  scale_shape_manual(name = "Species Type",
                     labels = c("Not Pest","Pest"),
                     values = c(16,17)) +
  geom_text(aes(label = sp_foo),
            nudge_x = 3.5,
            size = labeltext_size,
            family = labeltext_family) +
  expand_limits(x = max(sp_foo_pest_class$sp_foo) * 1.2) +
  labs(x = x_lab, y = y_lab) +
  plot_theme

# Export plot out as a png file to use. This will take last created plot in the running session
# Saved in the source project directory
ggsave(filename = "diet_foo_plot.png",plot = last_plot(),
       dpi = 300, units = "px", width = 1600, height = 1080)


