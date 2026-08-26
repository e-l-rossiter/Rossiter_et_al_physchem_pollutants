# ============================================================
# Script 2
# Prepare Chemicalize and eNTRYway inputs and join outputs
# ============================================================
#
# Purpose:
#
# This script uses the PubChem information generated in Script 1
# to prepare input files for Chemicalize.com and eNTRYway.
#
# First, the PubChem workbook is read in and the compounds with MW <= 600 g/mol
# are selected. These are the compounds that we want to use for physicochemical
# property calculations.
#
# A text file containing SMILES and compound names is then created for bulk
# upload to Chemicalize.com. The Chemicalize output is downloaded as an SDF file
# and read back into R using ChemmineR.
#
# The Chemicalize SDF data block is converted into a table containing the
# calculated physicochemical properties. The Chemicalize SMILES are then used to
# create a second text file for bulk upload to eNTRYway.
#
# The eNTRYway output is read back into R as a CSV file. Compound names are then
# checked and renamed where needed, because spaces were removed from compound
# names during the upload process.
#
# Finally, the Chemicalize output, eNTRYway output and PubChem information are
# joined together by compound name. The final joined table is saved as an Excel
# workbook for later analysis.
#
# ============================================================


## ---------- 0) Setup ----------

# Clear the R environment. Use with care because it removes all existing objects.
rm(list = ls())

# Load required packages.
library(ChemmineR)  # read SDF files and extract SDF IDs/data blocks
library(dplyr)      # data wrangling: filter, select, mutate, rename, joins and pipes
library(readxl)     # read Excel files and check Excel sheet names
library(openxlsx)   # save Excel workbooks
library(here)      # make file paths portable across different computers
# by using paths relative to the project root

##Get package citations

# Get the citations for the R packages used in this script.
# These can be used in the methods section or supplementary information if needed.

citation("ChemmineR")
citation("dplyr")
citation("readxl")
citation("openxlsx")
citation("here")

#Check R version and session information

# Check R version.
R.version.string
# "R version 4.4.2 (2024-10-31 ucrt)"

# For this part, we are using an older version of R as some packages may not
# work with newer versions.

# Record session information so the package versions used are documented.
sessionInfo()

# Set working directory to where the PubChem information workbook is saved from Script 1.

setwd(here("Files", "Files_for_coding"))


# List files in the working directory to check that the PubChem workbook is there.
list.files()



# Read in the CSV version.
pubchem_info_small <- read.csv(
  "Pubchem_Downloads.csv",
  stringsAsFactors = FALSE
)

# name it pubchem_info_small again as it is the same as in Script 1

## ---------- 1) Check PubChem data  ----------

# check it
str(pubchem_info_small)
head(pubchem_info_small)
colnames(pubchem_info_small)


# check 
pubchem_info_small %>%
  count(Less_equal_600_da)

#  Less_equal_600_da   n
#1              TRUE 219

# this is the correct file we don't want any false


## ---------- 2) Create Chemicalize upload file ----------

# Make a text file for Chemicalize with just SMILES and compound name.

# We need this so as we can do a bulk prediction on chemicalize # using the SMILEs

upload_file_1 <- pubchem_info_small %>%
  select(
    SMILES,
    Compound_Name
  )

# Remove spaces from compound names.
# This is done because the Chemicalize/eNTRYway upload process can be easier
# when compound names do not contain spaces.

upload_file_1 <- upload_file_1 %>%
  mutate(
    Compound_Name = gsub(" ", "", Compound_Name)
  )

# Check the upload file.
head(upload_file_1)




## ---------- 3) Save Chemicalize upload file ----------

# Set working directory to Files_for_uploading directory 

setwd(here("Files", "Files_for_uploading"))


# Save as a tab-delimited text file so it can be uploaded to Chemicalize.
write.table(
  upload_file_1,
  file      = "Pubchem_SMILES_for_chemicalise_upload.txt",
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

# The text file was then uploaded to Chemicalize.com as a bulk upload.
# The Chemicalize output was downloaded as an SDF file and read back into R below.


## ---------- 4) Read in Chemicalize SDF output ----------

# Set working directory to where the Chemicalize SDF file is saved.
setwd(here("Files", "Downloaded_files"))

# Check files in the working directory.
list.files()

# Read in the SDF file downloaded from Chemicalize.
sdfset <- read.SDFset("Chemicalise_sdf_batch_calculation.sdf")


# Check how many compounds are in the SDF file.
length(sdfset)

# Look at IDs from the SDF header.
sdfid(sdfset)


## ---------- 5) Convert Chemicalize SDF data block into a table ----------

# Convert the SDF data fields into a matrix.
# Rows are molecules and columns are the Chemicalize output fields.
block_ma <- datablock2ma(
  datablocklist = datablock(sdfset)
)

# Make a data frame containing the SDF ID and all Chemicalize data fields.
chemicalize_tbl <- data.frame(
  SDFID = sdfid(sdfset),
  block_ma,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Check the column names in the Chemicalize output table.
colnames(chemicalize_tbl)

## ---------- 6) Tidy the Chemicalize batch ----------


# Make a copy
chemicalize_tbl_clean <- chemicalize_tbl

# Remove anything like [MProp:scalar:double], [MProp:scalar:integer], etc.
colnames(chemicalize_tbl_clean) <- gsub(
  "\\s*\\[MProp:scalar:[^]]+\\]",
  "",
  colnames(chemicalize_tbl_clean)
)

# Make duplicate column names unique after cleaning
colnames(chemicalize_tbl_clean) <- make.unique(colnames(chemicalize_tbl_clean), sep = "_")

# Check cleaned names
colnames(chemicalize_tbl_clean)

# we have some duplicate rows so lets check what is happening here 

duplicate_cols_check <- chemicalize_tbl_clean %>%
  select(
    SDFID,
    `Strongest acidic pKa`,
    `Strongest acidic pKa_1`,
    `Strongest basic pKa`,
    `Strongest basic pKa_1`,
    `Isoelectric point`,
    `Isoelectric point_1`,
    `Intrinsic solubility [logS]`,
    `Intrinsic solubility [logS]_1`
  )

View(duplicate_cols_check)

# Check unique values in each duplicate column pair

unique(chemicalize_tbl_clean$`Strongest acidic pKa`)
unique(chemicalize_tbl_clean$`Strongest acidic pKa_1`)

#`Strongest acidic pKa` has the actual values while  `Strongest acidic pKa_1`  has N/As se we need to get rid of Strongest acidic pKa_1

unique(chemicalize_tbl_clean$`Strongest basic pKa`)
unique(chemicalize_tbl_clean$`Strongest basic pKa_1`)


#`Strongest basic pKa_1` has the actual values while  `Strongest basic pKa` has N/As so we need to remove `Strongest basic pKa`


unique(chemicalize_tbl_clean$`Isoelectric point`)
unique(chemicalize_tbl_clean$`Isoelectric point_1`)


#`Isoelectric point_1` has the actual values while `Isoelectric point`, just has N/As so we need to remove `Isoelectric point`,

unique(chemicalize_tbl_clean$`Intrinsic solubility [logS]`)
unique(chemicalize_tbl_clean$`Intrinsic solubility [logS]_1`)


#`Intrinsic solubility [logS]`, has the actual values while  `Intrinsic solubility [logS]_1  has N/As so we need to remove  `Intrinsic solubility [logS]_1

# remove these 

chemicalize_tbl_clean_final <- chemicalize_tbl_clean %>%
  select(
    -`Strongest acidic pKa_1`,
    -`Strongest basic pKa`,
    -`Isoelectric point`,
    -`Intrinsic solubility [logS]_1`
  )

# rename 
chemicalize_tbl_clean_final <- chemicalize_tbl_clean_final %>%
  rename(
    `Strongest basic pKa` = `Strongest basic pKa_1`,
    `Isoelectric point` = `Isoelectric point_1`
  )

# check colnames
colnames(chemicalize_tbl_clean_final)

# reorder 

chemicalize_tbl_clean_final <- chemicalize_tbl_clean_final %>%
  relocate(
    `Strongest basic pKa`,
    `Isoelectric point`,
    .after = `Strongest acidic pKa`
  )


# check colnames
colnames(chemicalize_tbl_clean_final)



## ---------- 7) Create eNTRYway upload file ----------

# reneame chemicalize_tbl_clean_final
chemicalize_tbl_2 <- chemicalize_tbl_clean_final



# Take the Chemicalize SMILES and SDF IDs to upload to eNTRYway.
chemicalize_smiles <- chemicalize_tbl_2 %>%
  select(
    SMILES,
    SDFID
  )

# set working directory 

setwd(here("Files", "Files_for_uploading"))


# Save as a tab-delimited text file so it can be uploaded to eNTRyway.
write.table(
  chemicalize_smiles,
  file      = "Chemicalise_SMILES_for_eNTRy_upload.txt",
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

# This file is then uploaded to:
# https://entryway.igb.illinois.edu/
#
# The eNTRYway output is downloaded as a CSV file and read back into R below.



## ---------- 8) Read in eNTRYway output ----------
# set working directory to where eNTRy way result was downloaded

setwd(here("Files", "Downloaded_files"))

# Check files in the working directory.
list.files()

#The downloaded eNTRy way file is eNTRy_way_Batch

# Read in the eNTRYway results CSV file.
results_entry_way <- read.csv(
  "eNTRy_way_Batch.csv"
)


## ---------- 9) Check compound-name matching between PubChem and eNTRYway ----------


# Because spaces were removed in the upload files, we need to add them back
# before joining the tables.

# Make another copy of the Chemicalize batch table.
chemicalize_tbl_3 <- chemicalize_tbl_2

# Check column names before renaming.
colnames(chemicalize_tbl_3)
colnames(results_entry_way)

# Rename eNTRYway columns so that they match later
results_entry_way <- results_entry_way %>%
  rename(
    Compound_Name = name,
    SMILES        = smiles
  )


# Check names in PubChem but NOT in eNTRYway.
pub_chem_only <- setdiff(
  pubchem_info_small$Compound_Name,
  results_entry_way$Compound_Name
)

pub_chem_only

#[1] "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonic acid"
#[2] "(2,4-Dichlorophenoxy)Acetic acid"                            
#[3] "2-methyl-4-chlorophenoxyacetic acid"                         
#[4] "4-Phenoxybutyric acid"                                       
#[5] "Acetaminophen sulfate"                                       
#[6] "Acetylsalicylic acid"                                        
#[7] "Clofibric acid"                                              
#[8] "Ethylenediaminetetraacetic acid"                             
#[9] "Fipronil sulfone"                                            
#[10] "Metazachlor ESA"                                           
#[11] "Methylchlorophenoxypropionic acid"                         
#[12] "Salicylic acid"                                            
#[13] "Kanamycin a"                                               
#[14] "Nalidixic acid"   
#total 14

# Check names in eNTRYway but NOT in PubChem.
results_only <- setdiff(
  results_entry_way$Compound_Name,
  pubchem_info_small$Compound_Name
)

results_only

#[1] "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonicacid"
#[2] "(2,4-Dichlorophenoxy)Aceticacid"                            
#[3] "2-methyl-4-chlorophenoxyaceticacid"                         
#[4] "4-Phenoxybutyricacid"                                       
#[5] "Acetylsalicylicacid"                                        
#[6] "Clofibricacid"                                              
#[7] "Ethylenediaminetetraaceticacid"                             
#[8] "Fipronilsulfone"                                            
#[9] "MetazachlorESA"                                             
#[10] "Methylchlorophenoxypropionicacid"                           
#[11] "Acetaminophensulfate"                                       
#[12] "Salicylicacid"                                              
#[13] "Kanamycina"                                                 
#[14] "Nalidixicacid"


# These are the same compounds, but the eNTRYway names have no spaces.
# We rename them so they match the original PubChem compound names.


## ---------- 10) Rename eNTRYway compound names ----------

results_entry_way <- results_entry_way %>%
  mutate(
    Compound_Name = case_when(
      Compound_Name == "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonicacid" ~ "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonic acid",
      Compound_Name == "(2,4-Dichlorophenoxy)Aceticacid"                             ~ "(2,4-Dichlorophenoxy)Acetic acid",
      Compound_Name == "2-methyl-4-chlorophenoxyaceticacid"                         ~ "2-methyl-4-chlorophenoxyacetic acid",
      Compound_Name == "4-Phenoxybutyricacid"                                       ~ "4-Phenoxybutyric acid",
      Compound_Name == "Acetylsalicylicacid"                                        ~ "Acetylsalicylic acid",
      Compound_Name == "Clofibricacid"                                              ~ "Clofibric acid",
      Compound_Name == "Ethylenediaminetetraaceticacid"                             ~ "Ethylenediaminetetraacetic acid",
      Compound_Name == "Fipronilsulfone"                                            ~ "Fipronil sulfone",
      Compound_Name == "MetazachlorESA"                                             ~ "Metazachlor ESA",
      Compound_Name == "Methylchlorophenoxypropionicacid"                           ~ "Methylchlorophenoxypropionic acid",
      Compound_Name == "Acetaminophensulfate"                                       ~ "Acetaminophen sulfate",
      Compound_Name == "Salicylicacid"                                              ~ "Salicylic acid",
      Compound_Name == "Kanamycina"                                                 ~ "Kanamycin a",
      Compound_Name == "Nalidixicacid"                                              ~ "Nalidixic acid",
      TRUE                                                                          ~ Compound_Name
    )
  )

# Check again for names in eNTRYway but NOT in PubChem.
results_only <- setdiff(
  results_entry_way$Compound_Name,
  pubchem_info_small$Compound_Name
)

results_only
# character(0)

# Check again for names in PubChem but NOT in eNTRYway.
pub_chem_only <- setdiff(
  pubchem_info_small$Compound_Name,
  results_entry_way$Compound_Name
)

pub_chem_only
# character(0)

# eNTRy way compound names renamed sucessfully 



## ---------- 11) Rename Chemicalize compound name column ----------

# Rename the Chemicalize SDF ID column so that it matches the other tables.
chemicalize_tbl_3 <- chemicalize_tbl_3 %>%
  rename(
    Compound_Name = SDFID
  )

# Check names in Chemicalize but NOT in eNTRYway.
chemicalise_batch_only <- setdiff(
  chemicalize_tbl_3$Compound_Name,
  results_entry_way$Compound_Name
)

chemicalise_batch_only
#[1] "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonicacid"
#[2] "(2,4-Dichlorophenoxy)Aceticacid"                            
#[3] "2-methyl-4-chlorophenoxyaceticacid"                         
#[4] "4-Phenoxybutyricacid"                                       
#[5] "Acetylsalicylicacid"                                        
#[6] "Clofibricacid"                                              
#[7] "Ethylenediaminetetraaceticacid"                             
#[8] "Fipronilsulfone"                                            
#[9] "MetazachlorESA"                                             
#[10] "Methylchlorophenoxypropionicacid"                           
#[11] "Acetaminophensulfate"                                       
#[12] "Salicylicacid"                                              
#[13] "Kanamycina"                                                 
#[14] "Nalidixicacid"


## ---------- 12) Rename Chemicalize compound names ----------

# Rename the Chemicalize compound names in the same way as the eNTRYway names.
chemicalize_tbl_3 <- chemicalize_tbl_3 %>%
  mutate(
    Compound_Name = case_when(
      Compound_Name == "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonicacid" ~ "((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonic acid",
      Compound_Name == "(2,4-Dichlorophenoxy)Aceticacid"                             ~ "(2,4-Dichlorophenoxy)Acetic acid",
      Compound_Name == "2-methyl-4-chlorophenoxyaceticacid"                         ~ "2-methyl-4-chlorophenoxyacetic acid",
      Compound_Name == "4-Phenoxybutyricacid"                                       ~ "4-Phenoxybutyric acid",
      Compound_Name == "Acetylsalicylicacid"                                        ~ "Acetylsalicylic acid",
      Compound_Name == "Clofibricacid"                                              ~ "Clofibric acid",
      Compound_Name == "Ethylenediaminetetraaceticacid"                             ~ "Ethylenediaminetetraacetic acid",
      Compound_Name == "Fipronilsulfone"                                            ~ "Fipronil sulfone",
      Compound_Name == "MetazachlorESA"                                             ~ "Metazachlor ESA",
      Compound_Name == "Methylchlorophenoxypropionicacid"                           ~ "Methylchlorophenoxypropionic acid",
      Compound_Name == "Acetaminophensulfate"                                       ~ "Acetaminophen sulfate",
      Compound_Name == "Salicylicacid"                                              ~ "Salicylic acid",
      Compound_Name == "Kanamycina"                                                 ~ "Kanamycin a",
      Compound_Name == "Nalidixicacid"                                              ~ "Nalidixic acid",
      TRUE                                                                          ~ Compound_Name
    )
  )

# Check again for names in Chemicalize but NOT in eNTRYway.
chemicalise_batch_only <- setdiff(
  chemicalize_tbl_3$Compound_Name,
  results_entry_way$Compound_Name
)

chemicalise_batch_only
# character(0)

#success

## ---------- 13) Save the Chemicalize output and eNTRy way outputs as Table_S2_2  ----------

# set working directory 

setwd(here("Files", "Tables"))

# Save Table_S2_2_Chemicalize_batch as Excel file
write.xlsx(
  chemicalize_tbl_3,
  file = "Table_S2_2_Chemicalize_batch.xlsx",
  overwrite = TRUE
)

# Save Table_S2_5_eNTRy_way_Batch as Excel file
write.xlsx(
  results_entry_way,
  file = "Table_S2_5_eNTRy_way_Batch.xlsx",
  overwrite = TRUE
)


## ---------- 14) Rename Chemicalize output columns ----------

# Check the Chemicalize batch table before renaming columns.
head(chemicalize_tbl_3)
str(chemicalize_tbl_3)

# Check number of columns and column names.
ncol(chemicalize_tbl_3)
colnames(chemicalize_tbl_3)


# Create new names for the Chemicalize batch output.
new_names <- c(
  "Compound_Name",
  "molar_mass_g_mol",
  "exact_mass_da",
  "formula",
  "composition",
  "lipinski_ro5",
  "atom_count",
  "heavy_atom_count",
  "asymmetric_atom_count",
  "rotatable_bond_count",
  "ring_count",
  "aromatic_ring_count",
  "hetero_ring_count",
  "fsp3",
  "hbd_count",
  "hba_count",
  "formal_charge",
  "tpsa_a2",
  "polarizability_a3",
  "molar_refractivity_cm3_mol",
  "vdw_volume_a3",
  "min_projection_area_a2",
  "max_projection_area_a2",
  "min_projection_radius_a",
  "max_projection_radius_a",
  "vdw_surface_area_a2",
  "solvent_accessible_surface_area_a2",
  "pka_acid_strongest",
  "pka_base_strongest",
  "isoelectric_point",
  "intrinsic_solubility_logS",
  "solubility_category",
  "logp",
  "logd_ph0",
  "logd_ph1",
  "logd_ph2",
  "logd_ph3",
  "logd_ph4",
  "logd_ph5",
  "logd_ph6",
  "logd_ph7",
  "logd_ph8",
  "logd_ph9",
  "logd_ph10",
  "logd_ph11",
  "logd_ph12",
  "logd_ph13",
  "logd_ph14",
  "hlb",
  "iupac_name",
  "smiles",
  "inchi",
  "inchikey",
  "traditional_name",
  "common_names",
  "cas_registry_numbers"
)

# Check that the number of new names matches the number of columns.
stopifnot(ncol(chemicalize_tbl_3) == length(new_names))

# Apply the new column names.
names(chemicalize_tbl_3) <- new_names

# Check renamed columns.
colnames(chemicalize_tbl_3)

# Make a copy of the Chemicalize batch table.
chemicalize_tbl_4 <- chemicalize_tbl_3


## ---------- 15) Rename overlapping columns before joining ----------

# Check column names before joining.
colnames(chemicalize_tbl_4)
colnames(results_entry_way)
colnames(pubchem_info_small)

# Rename PubChem columns so that it is clear they came from PubChem.
pubchem_info_small <- pubchem_info_small %>%
  rename(
    smiles_PC       = SMILES,
    MW_Da_PC           = MW_Da,
    CAS_PC          = CAS  )%>%
  select(-Less_equal_600_da)

# Rename eNTRYway SMILES column so that it is clear it came from eNTRYway.
results_entry_way <- results_entry_way %>%
  rename(
    smiles_ew = SMILES
  )

# Rename eNTRYway columns so that it is clear they came from eNTRYway.
results_entry_way <- results_entry_way %>%
  rename(
    molwt_EW      = molwt,
    formula_EW    = formula,
    rb_EW         = rb,
    glob_EW       = glob,
    pbf_EW        = pbf,
    func_group_EW = func_group
  )

# Check column names after renaming.
colnames(chemicalize_tbl_4)
colnames(results_entry_way)
colnames(pubchem_info_small)



## ---------- 16) Join Chemicalize, eNTRYway and PubChem outputs ----------

# Join Chemicalize and eNTRYway outputs by compound name.
joined_chemicalize_eNTRy_way <- chemicalize_tbl_4 %>%
  left_join(
    results_entry_way,
    by = "Compound_Name"
  )

# Check column names after first join.
colnames(joined_chemicalize_eNTRy_way)

# Join the PubChem information to the Chemicalize/eNTRYway table.
joined_chemicalize_eNTRy_way_pub_chem <- joined_chemicalize_eNTRy_way %>%
  left_join(
    pubchem_info_small,
    by = "Compound_Name"
  )

# Check final column names.
colnames(joined_chemicalize_eNTRy_way_pub_chem)


## ---------- 17) Save final joined table ----------
# Set working directory

setwd(here("Files", "Files_for_coding"))

# Save the final joined Chemicalize, eNTRYway and PubChem table for use in Script 3.
write.xlsx(
  joined_chemicalize_eNTRy_way_pub_chem,
  "joined_chemicalize_eNTRy_way_pub_chem.xlsx"
)
