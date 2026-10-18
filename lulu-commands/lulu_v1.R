# packages ----

library(lulu) ; library(tidyverse) ; library(data.table)

# Data input ----

# otu table boldigger 95
otu_tbl <- read_tsv("lulu-commands/BOLDigger_OTU_table_95.tsv")

# matchlist made with blast+
match_ls <- data.table::fread("lulu-commands/match_list.txt")


# Format OTU table for lulu ----

otu_tbl_clean <- otu_tbl |>
  select(OTU_ID,bat_01:bat_99)

otu_tbl_clean$OTU_ID <-   sub("_0+","_", otu_tbl_clean$OTU_ID)


# check if all the OTU IDs in match_ls are present in otu_tbl
# sum(...) should be equal to length(unique(match_ls$V1))
length(unique(match_ls$V1))
sum(unique(match_ls$V1) %in% otu_tbl_clean$OTU_ID)

# set otu ids as rownames in the otu table since lulu wants only numeric cols
otu_tbl_clean <- as.data.frame(otu_tbl_clean)
otu_tbl_num <- otu_tbl_clean[, -1]
rownames(otu_tbl_num) <- otu_tbl_clean[, 1]

# lulu run ---- this step takes about 4 to 5 hours on the latest dataset
lulu_out <- lulu(otutable = otu_tbl_num, matchlist = match_ls)

# Save LULU output object
saveRDS(lulu_out, file = "lulu_output_v1_18_oct_2026.rds")

# lulu outputs ----
lulu_out$curated_count ; lulu_out$discarded_count
lulu_out$runtime
lulu_out$curated_table

# lulu output has OTU IDs as rownames - covert these to a col and export

lulu_otu_tbl <- rownames_to_column(lulu_out$curated_table, "otu_id")
data.table::fwrite(lulu_otu_tbl,
                   file = "lulu_otu_tbl_boldigger95_v1.csv",
                   col.names = T)

# data-table method, that modifies the object, does not create a copy
# setDT(lulu_out$curated_table,keep.rownames = "otu_id")

# OTUs' taxon matches need to be matched using the OTU IDs as the key and joined with the lulu curated OTU table
# Match the OTU ID format betwee nthe two tables, then join the tables, i.e. bring the taxon names to the lulu otu table







