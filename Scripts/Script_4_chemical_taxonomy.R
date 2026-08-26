# ============================================================
# Script 4
# Get ClassyFire chemical classifications
# ============================================================
#
# Purpose:
#
# This script uses the final physicochemical information table created in
# Data Wrangling pt. 3 and adds chemical classification information from
# ClassyFire.
#
# The script first reads in the final antibiotic and pollutant physicochemical
# information table. It then selects the compound names and PubChem CIDs.
#
# PubChem is searched using the webchem package to retrieve the InChIKey for
# each compound from its CID.
#
# The InChIKeys are then used with the classyfireR package to retrieve chemical
# classifications from ClassyFire. The classification levels extracted are:
#
# - kingdom
# - superclass
# - class
# - subclass
#
# The ClassyFire classifications are then joined back to the compound table and
# saved as an Excel file.
#
# ============================================================


## ---------- 0) Setup ----------

# Clear the global environment.
# Use with care because this removes all existing objects.
rm(list = ls())

# Load required packages.
library(webchem)      # retrieve InChIKeys from PubChem using CIDs
library(classyfireR)  # retrieve ClassyFire chemical classifications
library(dplyr)        # data wrangling: select, rename, joins and pipes
library(purrr)        # safely run ClassyFire queries and combine outputs
library(tibble)       # create tidy tables
library(writexl)      # save Excel files
library(readxl)       # read Excel files
library(openxlsx)     # create Excel workbooks with multiple sheets
library(here)         # make file paths portable across different computers
# by using paths relative to the project root


# These can be used in the methods section or supplementary information if needed.
citation("webchem")
citation("classyfireR")
citation("purrr")
citation("readxl")
citation("openxlsx")
citation("here")


#Check R version and session information

# Check R version.
R.version.string

# Record session information so the package versions used are documented.
sessionInfo()


#Set working directory and read in final data

# Set working directory to where the final physicochemical information file is saved.
setwd(here("Files", "Files_for_coding"))


# List files to check that the final workbook is there.
list.files()

# Read in the final antibiotic and pollutant physicochemical information table.
df <- read.csv("Final_AB_Pol_Physchem_info.csv")


## ---------- 1) Select compound names and PubChem CIDs ----------

# Make a smaller data frame containing only the columns needed to get InChIKeys.
df2 <- df %>%
  dplyr::select(
    Compound_Name,
    CID
  )

# Check the data frame.
head(df2)
str(df2)


## ---------- 2) Get InChIKeys from PubChem ----------

# Use PubChem CIDs to retrieve InChIKeys.
inchikey_tbl <- pc_prop(
  df2$CID,
  from       = "cid",
  properties = c("InChIKey")
)

# Check the InChIKey table.
head(inchikey_tbl)

# Join the InChIKeys back to the compound name and CID table.
df3 <- df2 %>%
  left_join(
    inchikey_tbl %>%
      select(
        CID = CID,
        InChIKey
      ),
    by = "CID"
  )

# Check for any missing InChIKeys.
df3 %>%
  filter(is.na(InChIKey))
#[1] Compound_Name CID           InChIKey     
#<0 rows> (or 0-length row.names)

#none are missing 

## ---------- 3) Get ClassyFire classifications ----------

# Get a ClassyFire classification object for each InChIKey.
#
# safely() is used so that if one compound fails, the whole script does not stop.
# Failed classifications are returned as NULL and handled in the next section.
class_list <- map(
  df3$InChIKey,
  safely(get_classification)
)


## ---------- 4) Extract ClassyFire classification levels ----------

# Extract kingdom, superclass, class and subclass from each ClassyFire result.
class_df <- map_dfr(seq_along(class_list), function(i) {
  
  # Get the ClassyFire result for the current compound.
  res <- class_list[[i]]$result
  
  # Keep the InChIKey from the original table so everything stays in the same order.
  ik <- df3$InChIKey[i]
  
  # If the InChIKey failed classification, return NA values.
  if (is.null(res)) {
    return(tibble(
      InChIKey   = ik,
      kingdom    = NA_character_,
      superclass = NA_character_,
      class      = NA_character_,
      subclass   = NA_character_
    ))
  }
  
  # res is a ClassyFire object.
  # The classification slot contains a table with:
  # - Level
  # - Classification
  # - CHEMONT
  cf_tbl <- res@classification
  
  # Small helper function to extract each classification level.
  get_level <- function(level_name) {
    x <- cf_tbl$Classification[cf_tbl$Level == level_name]
    if (length(x) == 0) NA_character_ else x[1]
  }
  
  # Return one row per compound.
  tibble(
    InChIKey   = ik,
    kingdom    = get_level("kingdom"),
    superclass = get_level("superclass"),
    class      = get_level("class"),
    subclass   = get_level("subclass")
  )
})

# Check the classification table.
head(class_df)

# Check for any missing ClassyFire classifications.
class_df %>%
  filter(
    is.na(kingdom) |
      is.na(superclass) |
      is.na(class) |
      is.na(subclass)
  )


## ---------- 5) Join ClassyFire classifications back to compound table ----------

# Join the ClassyFire classes back to the compound table.
df4 <- df3 %>%
  left_join(
    class_df,
    by = "InChIKey"
  )

# Check the final classification table.
head(df4)
str(df4)


## ---------- 6) Save ClassyFire classification output ----------

# set working directory
setwd(here("Files", "Tables"))

# Save the final table containing compound names, CIDs, InChIKeys and
# ClassyFire classifications.
# These will be used in Supplementary Table S1.1 and to create Table_S2_3_Chemical_Taxonomy
write_xlsx(
  df4,
  "Table_S2_3_Chemical_Taxonomy.xlsx"
)
