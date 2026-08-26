# ============================================================
# Script 1
# Get PubChem identifiers and SMILES for compounds of interest
# ============================================================
#
# Purpose:
#
# This script creates a PubChem reference table for the pollutants and
# antibiotics of interest, so that their SMILES strings can be used later to
# calculate physicochemical properties on Chemicalize.com.
#
# The script also selects only compounds with MW <= 600 Da, as we did not want
# to compute physicochemical properties for compounds that were not of interest.
#
#
# The script first defines the compound lists used in the study, including:
#
# - pollutants of interest the sources of which are referenced in the paper
# - antibiotics from O'Shea, Rosemarie, and Heinz E. Moser:
# "Physicochemical properties of antibacterial compounds: implications for
# drug discovery." Journal of Medicinal Chemistry,
# doi:10.1021/jm700967e
# - additional antibiotics from other sources which are referenced in the paper
#
#

# PubChem is searched using the webchem package to retrieve PubChem compound
# identifiers, known as CIDs, for each compound name. Where PubChem does not
# automatically return a CID, we manually add the correct CID for Polymyxin B1.
#
#
# Using the PubChem CIDs, the webchem package retrieves compound information
# from PubChem, including:
#
# - SMILES
# - molecular weight
# - CAS numbers
#
# The webchem package is also used to retrieve PubChem synonyms for each
# compound. We then extract the first valid CAS Registry Number where available.
#
# The final PubChem information table, containing data extracted from PubChem,
# includes:
#
# - Compound_Name
# - PubChem CID
# - SMILES
# - molecular weight (g/mol)
# - CAS number
# - CAS number with dashes removed
#
#
# The compounds are split into two groups based on molecular weight:
#
# - compounds with molecular weight less than or equal to 600 g/mol
# - compounds with molecular weight greater than 600 g/mol
#
# This is because we are only interested in compounds that are <= 600 Da.
# 600 Da is exactly equivalent to 600 g/mol.
#
#
# Finally, the PubChem information is saved into one Excel workbook with three
# sheets:
#
# - All_comps: all compounds with PubChem information
# - Comps_less_or_equal_600_g_mol: compounds with MW <= 600 g/mol
# - Comps_greater_than_600_g_mol: compounds with MW > 600 g/mol
#
#
# These outputs are used later for bulk upload to Chemicalize.com, to calculate
# physicochemical properties using the SMILES strings in the web-based
# application.
#
#
# NOTE for anyone wanting to use Chemicalize:
# Unfortunately, ChemAxon is retiring Chemicalize on 30 June 2027.

## ---------- 0) Setup ----------

# Clear the environment so that the script starts from a clean workspace.
rm(list = ls())

# Load required packages.
library(webchem)   # search PubChem and retrieve CIDs, SMILES, MWs, synonyms and CAS numbers
library(dplyr)     # data wrangling: filter, mutate, rename, joins and pipes
library(purrr)     # loop over PubChem synonym lists and combine outputs
library(tibble)    # create tidy tables
library(writexl)   # save simple Excel files
library(openxlsx)  # create Excel workbooks with multiple sheets
library(here)      # make file paths portable across different computers
                   # by using paths relative to the project root

# Get the citations for the R packages used in this script.

citation("webchem")
citation("dplyr")
citation("purrr")
citation("tibble")
citation("writexl")
citation("openxlsx")
citation("here")

# Check R version.
R.version.string
# "R version 4.4.2 (2024-10-31 ucrt)"
# For this part, we are using an older version of R as some packages may not
# work with newer versions.

# Record session information so the package versions used are documented.
sessionInfo()

#
## ---------- 1) Define compound lists used in the study ----------

# These are the pollutants by name.
pollutants <- c("((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonic acid",
                "(2,4-Dichlorophenoxy)Acetic acid",
                "10,11-dihydroxycarbamazepine",
                "1-Decanol",
                "2(3H)-Benzothiazolone",
                "2,6-Dichlorobenzamide",
                "2-Hydroxyibuprofen",
                "2-methyl-4-chlorophenoxyacetic acid",
                "4-Phenoxybutyric acid",
                "5-Amino-2-methylphenol",
                "5-Methyl-1H-benzotriazole",
                "Acenaphthene",
                "Acetaminophen",
                "Acetaminophen sulfate",
                "Acetamiprid",
                "Acetylsalicylic acid",
                "Aclonifen",
                "Amitriptyline",
                "Amphetamine",
                "Aniline",
                "Atenolol",
                "Atrazine",
                "Bentazone",
                "Benzophenone",
                "1H-Benzotriazole",
                "Benzoylecgonine",
                "Caffeine",
                "Carbamazepine",
                "Carbamazepine-10,11-epoxide",
                "Cetirizine",
                "Chlorantraniliprole",
                "Chloridazon-desphenyl",
                "Citalopram",
                "Clofibric acid",
                "Clothianidin",
                "Cocaine",
                "Codeine",
                "Cotinine",
                "Cyantraniliprole",
                "Cyclohexanone",
                "Dicamba",
                "Dichlobenil",
                "Diclofenac",
                "Diphenylamine",
                "Diuron",
                "Ethylenediaminetetraacetic acid",
                "Estradiol",
                "Estriol",
                "Ethylbenzene",
                "Fipronil",
                "Fipronil sulfone",
                "Flufenacet",
                "Fluoranthene",
                "Fluoxetine",
                "Flupyradifurone",
                "Gabapentin",
                "Gliclazide",
                "Guanylurea",
                "Hydrocodone",
                "Ibuprofen",
                "Imidacloprid",
                "Isoprocarb",
                "Ketamine",
                "Lamotrigine",
                "Lidocaine",
                "Linalool",
                "Sotalol",
                "Melamine",
                "Metaflumizone",
                "Metazachlor",
                "Metazachlor ESA",
                "Metformin",
                "Methylamine",
                "Methylchlorophenoxypropionic acid",
                "Modafinil",
                "Morphine",
                "N,N,N',N'-Tetraacetylethylenediamine",
                "N,N-diethyl-m-toluamide",
                "N,N-Dimethylaniline",
                "Naproxen",
                "Nitrobenzene",
                "p-Benzoquinone",
                "Propyzamide",
                "Ranitidine",
                "Salicylic acid",
                "Schradan",
                "Sertraline",
                "Simazine",
                "Sitagliptin",
                "Spirotetramat",
                "Sulisobenzone",
                "Terbutryn",
                "Terpineol",
                "Tonalide",
                "Tramadol",
                "N-Phenyl-1-naphthylamine",
                "2-Naphthylamine",
                "1-Naphthylamine", 
                "N-nitrosodiphenylamine",
                "Venlafaxine"
)

#see Supplementary Table S1.1 for references 


## these are the antibiotics by name from (O’Shea and Moser 2008)

GN_antibiotics_oshea <- c(
  "Kanamycin a",
  "Paromomycin",
  "Streptomycin",
  "Amikacin",
  "Arbekacin",
  "Dibekacin",
  "Gentamicin",
  "Isepamicin",
  "Neomycin",
  "Netilmicin",
  "Sisomicin",
  "Tobramycin",
  "Chloramphenicol",
  "Loracarbef",
  "Ertapenem",
  "Imipenem",
  "Meropenem",
  "Tomopenem", #*
  "Cefetamet",
  "Ceftibuten",
  "Cefaclor",
  "Cefadroxil",
  "Cefamandole",
  "Cefazolin",
  "Cefdinir",
  "Cefditoren",
  "Cefixime",
  "Cefmetazole",
  "Cefoperazone",
  "Cefotaxime",
  "Cefotetan",
  "Cefoxitin",
  "Cefpodoxime",
  "Cefprozil",
  "Ceftizoxime",
  "Ceftriaxone",
  "Cefuroxime",
  "Cephalexin",
  "Cephalothin",
  "Cephradine",
  "Ceftaroline", #*
  "Cefepime",
  "Cefpirome",
  "Ceftazidime",
  "Ceftobiprole",
  "Iclaprim",
  "Trimethoprim",
  "Nalidixic acid",
  "Delafloxacin", #*
  "Ciprofloxacin",
  "Clinafloxacin",
  "Danofloxacin",
  "Difloxacin",
  "Dx-619", #*
  "Enoxacin",
  "Fleroxacin",
  "Garenoxacin",
  "Gatifloxacin",
  "Gemifloxacin",
  "Grepafloxacin",
  "Levofloxacin",
  "Lomefloxacin",
  "Moxifloxacin",
  "Nadifloxacin",
  "Norfloxacin",
  "Pefloxacin",
  "Rufloxacin",
  "Sitafloxacin",
  "Sparfloxacin",
  "Temafloxacin",
  "Trovafloxacin",
  "Azithromycin",
  "Aztreonam",
  "Fosfomycin",
  "Faropenem",
  "Doripenem",
  "Amoxicillin",
  "Ampicillin",
  "Carbenicillin",
  "Mezlocillin",
  "Ticarcillin",
  "Azlocillin",
  "Piperacillin",
  "Sulfabenzamide",
  "Sulfacetamide",
  "Sulfachlorpyridazine",
  "Sulfadiazine",
  "Sulfadimethoxine",
  "Sulfaguanidine",
  "Sulfamerazine",
  "Sulfameter",
  "Sulfamethazine",
  "Sulfamethizole",
  "Sulfamethoxazole",
  "Sulfamethoxypyridazine",
  "Sulfamonomethoxine",
  "Sulfanitran",
  "Sulfaphenazole",
  "Sulfapyridine",
  "Sulfaquinoxaline",
  "Sulfathiazole",
  "Sulfisoxazole",
  "Polymyxin b1",
  "Demeclocycline",
  "Doxycycline",
  "Meclocycline",
  "Methacycline",
  "Minocycline",
  "Oxytetracycline",
  "Omadacycline", #*
  "Tetracycline",
  "Tigecycline",
  "Chlortetracycline"
)

#* note in O'shea and Mosser 2008 

# Some antibiotics are listed in O'Shea and Moser 2008 using development
# or research code names 
#
# In the database https://antibioticdb.com/
#
# R-115685 is listed as Tomopenem: 
# https://antibioticdb.com/AdbCompoundDisplayForward?id=321 
# Which reached a highest development stage of Phase 2
# It's Development status is Inactive
# and its "Reason dropped" is given as:
# "Discontinued for unknown reasons"
#
# T-91825 is listed as Ceftaroline: 
# https://antibioticdb.com/AdbCompoundDisplayForward?id=1444  
# It's Highest development stage is listed as:
# "Approved by FDA in 2010; Approved for paediatric use in 2016
# It's Development status is Approved
#
# PTK-0796 is listed as Omadacycline: 
# https://antibioticdb.com/AdbCompoundDisplayForward?id=246
# It's Highest development stage is listed as: 
# Approved by FDA in 2018
# It's Development status is Approved
#
# Dx-619 is listed as Dx-619
# https://antibioticdb.com/AdbCompoundDisplayForward?id=498
# It's Highest development stage is listed as:
# Phase 1
# Development status is listed as	Inactive
#
# it is listed in antibioticdb.com as a 
# gram positive antibiotic
# 
# In O'shea and Mosser 2008 SI
# Dx-619 is listed as having:
# E.coli       MIC50:       0.03 µg/mL
# P aeruginosa MIC50:       1 µg/mL 
# S aureus     MIC50:       0.008 µg/mL 
#
# In O'shea and Mosser 2008: 
# Compounds with >100-fold difference between Gram-positive 
# and -negative MIC values were declared inactive against 
# Gram-negative bacteria, even if their MIC value is 8 
# or lower
#
# Hence why it is included in our study as a gram-negative 
# compound
#
# ABT-492 is not found in the database:
# https://antibioticdb.com/ 
# but it is listed on pubchem as Delafloxacin
# https://pubchem.ncbi.nlm.nih.gov/compound/487101
# According to pubchem:
# "It was approved in June 2017 under the trade name 
# Baxdela for use in the treatment of acute bacterial skin 
# and skin structure infections.


## these are the antibiotics by name from other sources
GN_antibiotics_other <- c("Ofloxacin", 
                       "Tigemonam",
                       "Lenapenem",
                       "Sulopenem",
                       "Biapenem",
                       "Tebipenem",
                       "Rolitetracycline",
                       "Eravacycline",
                       "Tetroxoprim",  
                       "Aditoprim",
                       "Florfenicol",
                       "Thiamphenicol") 

#see "Table 2: Gram-negative antibiotics evaluated in this study" for references 

# Combine all compounds into one list 


compounds <- c(pollutants, GN_antibiotics_oshea, GN_antibiotics_other)

## ---------- 2) Get PubChem CIDs for compound names ----------

# Get PubChem CIDs for names.
# This will take a few min minute to load.

cids <- webchem::get_cid(
  compounds,
  from = "name",
  match = "ask"
)

# match = "ask" is useful where PubChem returns more than one possible match.

# In this case, we just had more than one match for Ceftaroline
#
# Select result for 'Ceftaroline':
# 1: 9938701
# 2: 9852981

# NOTE the order maybe be flipped so double check that option 1 is the same for you before running this code

# Here is how it was resolved:

# These CIDs 9852981 and 9938701 were checked on PubChem 

# Option 1- 9938701 is for Ceftaroline. According to PubChem Ceftaroline. is a cephalosporin that is the active metabolite of the prodrug ceftaroline fosamil.

# Option 2- 9852981 is for Ceftaroline Fosamil. According to PubChem is an N-phosphono prodrug of the fifth-generation cephalosporin derivative ceftaroline with antibacterial activity. Ceftaroline fosamil is hydrolyzed to the active form ceftaroline in vivo. 


# We want Ceftaroline and not the prodrug therefore we select Option 1
#
# Selection: 1



# Check the first few results.

head(cids)

# A tibble: 6 × 2
#query                                                        cid     
#<chr>                                                        <chr>   
# 1 ((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonic acid 16212225
#2 (2,4-Dichlorophenoxy)Acetic acid                             1486    
#3 10,11-dihydroxycarbamazepine                                 83852   
#4 1-Decanol                                                    8174    
#5 2(3H)-Benzothiazolone                                        13625   
#6 2,6-Dichlorobenzamide                                        16183  

# Look for any N/As.
cids %>% 
  filter(is.na(cid))


# A tibble: 2 × 2
#query                         cid  
#<chr>                         <chr>
##  1 Polymyxin b1                  NA   

## ---------- 3) Manually fix any missing PubChem CIDs ----------

# Polymyxin B1 was not automatically returned by PubChem in this search.
# The correct PubChem CID is added manually.
#
# Polymyxin B1:
# https://pubchem.ncbi.nlm.nih.gov/compound/11228650
#CID = 11228650

#fix
cids <- cids %>%
  mutate(
    cid = if_else(
      query == "Polymyxin b1" & is.na(cid),
      "11228650",   
      cid
    )
  )

# Check that the manual correction worked.
cids %>% 
  dplyr::filter(query == "Polymyxin b1")


# A tibble: 1 × 2
#query        cid     
#<chr>        <chr>   
#  1 Polymyxin b1 11228650

#Success it worked 

# Look for any N/As to double check
cids %>% 
  filter(is.na(cid))

# A tibble: 0 × 2
# ℹ 2 variables: query <chr>, cid <chr>

#Success no more N/As 


## ---------- 4) Get SMILES and molecular weight from PubChem ----------

# Use the PubChem CIDs to retrieve SMILES and molecular weight.

props <- pc_prop(
  cids$cid, 
  from = "cid", 
  properties = c("SMILES",# SMILES: canonical with stereochemical and isotopic information.
                 "MolecularWeight"
  ) 
)
# Look at the first few results.
head(props)

# A tibble: 6 × 3
# CID      MolecularWeight SMILES                                   
# <chr>    <chr>           <chr>                                    
# 1 16212225 275.30          CC(C)N(C1=CC=C(C=C1)F)C(=O)CS(=O)(=O)O
# 2 1486     221.03          C1=CC(=C(C=C1Cl)Cl)OCC(=O)O      
# 3 83852    270.28          C1=CC=C2C(=C1)C(C(C3=CC=CC=C3N2C(=O)N)O)O
# 4 8174     158.28          CCCCCCCCCCO
# 5 13625    151.19          C1=CC=C2C(=C1)NC(=O)S2    
# 6 16183    190.02          C1=CC(=C(C(=C1)Cl)C(=O)N)Cl

# Check for missing SMILES.
props %>% dplyr::filter(is.na(SMILES))

# A tibble: 0 × 3
# ℹ 3 variables: CID <chr>, MolecularWeight <chr>, SMILES <chr>

#no N/As, no missing SMILES.

# Check for missing molecular weights
props %>% dplyr::filter(is.na(MolecularWeight))
# A tibble: 0 × 3
# ℹ 3 variables: CID <chr>, MolecularWeight <chr>, SMILES <chr>

#no N/As, missing molecular weights


# Rename MolecularWeight to MW.
props <- props %>%
  rename(
    MW = MolecularWeight
  )

# Check the renamed table.

head(props)

# A tibble: 6 × 3
#CID      MW     SMILES                                   
#<chr>    <chr>  <chr>                                    
#  1 16212225 275.30 CC(C)N(C1=CC=C(C=C1)F)C(=O)CS(=O)(=O)O
#2 1486     221.03 C1=CC(=C(C=C1Cl)Cl)OCC(=O)O            
#3 83852    270.28 C1=CC=C2C(=C1)C(C(C3=CC=CC=C3N2C(=O)N)O)O
#4 8174     158.28 CCCCCCCCCCO                    
#5 13625    151.19 C1=CC=C2C(=C1)NC(=O)S2                
#6 16183    190.02 C1=CC(=C(C(=C1)Cl)C(=O)N)Cl


## ---------- 5) Get PubChem synonyms and extract CAS numbers ----------

# Get synonyms for each CID from PubChem.
# CAS numbers are stored as synonyms in PubChem, so we retrieve all synonyms
# and then extract the ones that look like valid CAS Registry Numbers.

syns <- pc_synonyms(cids$cid, from = "cid", match = "all")

# The code above will take a little while to run.


# Extract the first valid CAS number for each CID, where available.
cas_tbl <- map2_dfr(
  syns,
  cids$cid,
  ~ {
    cas_vec <- .x[is.cas(.x)]         # keep only things that look like CAS RNs
    tibble(
      CID = .y,
      CAS = if (length(cas_vec)) cas_vec[1] else NA_character_
    )
  }
) %>%
  dplyr::distinct(CID, .keep_all = TRUE)   # <- ensure one row per CID

##lets break the above down
#map2_dfr is a purr function
#map2: iterate over two vectors/lists at the same time
#dfr: return each result and combine them row-wise into a data frame
#in the df syns we have a big list of all of the synonms stored in pubchem for each of our compound IDs
#so map2_dfr(syns, cids$cid, ~ { ... })Loops over two things in parallel:
#syns:a list where each element contains synonyms for one compound
#cids$cid: the PubChem CID for that compound
#For each pair:.x is the current synonym vector and .y is the current CID
#_dfr means each result is returned as a row-like object and then all results are combined into one data frame.
#cas_vec <- .x[is.cas(.x)]Checks which synonyms in .x look like valid CAS Registry Numbers

#is.cas is a function from webchem designed to extract whatever looks like a cas number
#Keeps only those
#So cas_vec becomes a vector of candidate CAS numbers

#tibble(
#CID = .y,
#CAS = if (length(cas_vec)) cas_vec[1] else NA_character_
#)
#Creates a one-row tibble for that compound
#CID = .y stores the current PubChem CID
#CAS = ... stores:the first CAS-looking synonym, if any exist, otherwise NA
#So after map2_dfr(...), you get a data frame with one row per input pair, containing: CID
#first detected CAS number or NA
#dplyr::distinct(CID, .keep_all = TRUE)Removes duplicate rows with the same CID,Keeps only the first row for each CID .keep_all = TRUE means it keeps the other columns too, not just CID


# View the CAS table.

cas_tbl

# check for missing CAS numbers
cas_tbl %>% dplyr::filter(is.na(CAS))

# CID      CAS  
# <chr>    <chr>
# 1 22722247 NA 

#lets check what this compound without as CAS number is 

cids %>% 
  dplyr::filter(cid == "22722247")

# A tibble: 1 × 2
# query                             cid     
# <chr>                             <chr>   
# 1 Methylchlorophenoxypropionic acid 22722247

#lets search the PubChem CID for Methylchlorophenoxypropionic acid and see why there is no CAS number

# https://pubchem.ncbi.nlm.nih.gov/compound/22722247

# There is no given CAS number for this compound


## ---------- 6) Build the final PubChem information table ----------
# Check the column names. 

colnames(cids)
#[1] "query" "cid"  
colnames(props)
#[1] "CID"    "MW"     "SMILES"
colnames(cas_tbl)
#[1] "CID" "CAS"


# Rename cid in cids to CID so that it matches the other tables.
cids <- cids %>%
  rename(
    CID = cid
  )

# Build a table containing compound name, CID, SMILES and MW.
# CAS numbers are then added by joining with cas_tbl.
smiles_pubchem <- tibble(
  Compound_Name = cids$query,
  CID           = cids$CID,
  SMILES        = props$SMILES,
  MW            = props$MW
) %>%
  left_join(cas_tbl, by = "CID")   

# Check the table.
head(smiles_pubchem)

# Rename this as the final PubChem information table.
pubchem_info <- smiles_pubchem

## ---------- 7) Split compounds by TYPE ----------

# Make new column in the pubchem_info df called TYPE identiying if the compound is in the lists made earlier:

# pollutants
# GN_antibiotics_oshea
# GN_antibiotics_other 

pubchem_info <- pubchem_info %>%
  mutate(
    TYPE = case_when(
      Compound_Name %in% pollutants ~ "pollutant",
      Compound_Name %in% GN_antibiotics_oshea ~ "GN_antibiotics_oshea",
      Compound_Name %in% GN_antibiotics_other ~ "antibiotics_other",
      TRUE ~ NA_character_
    )
  )

# check it worked:
pubchem_info %>%
  count(TYPE)

# A tibble: 3 × 2
#TYPE                  n
#<chr>             <int>
#  1 GN_antibiotics_oshea   113
#2 antibiotics_other    12
#3 pollutant           100

#it worked

#Check N/As
pubchem_info %>%
  filter(is.na(TYPE)) %>%
  select(Compound_Name)

# A tibble: 0 × 1
# ℹ 1 variable: Compound_Name <chr>
# no N/As

## ---------- 8) Split compounds by molecular weight ----------

# Make sure MW is numeric and add TRUE/FALSE column for MW <= 600 Da.
pubchem_info <- pubchem_info %>%
  mutate(
    MW = as.numeric(MW),
    Less_equal_600_da = MW <= 600
  )


# Compounds with MW greater than 600 Da (same as 600 g/mol which are the units used by pubchem).
pubchem_info_big <- pubchem_info %>%
  filter(!Less_equal_600_da)

# Compounds in pubchem_info_big will be excluded due to being over the porin mediated up-take cut off weight of600 Da (same as 600 g/mol which are the units used by pubchem).

# Check the compounds that are greater than 600 Da.

# Convert to a data frame for display, as tibbles can print rounded values.
as.data.frame(
  pubchem_info_big %>%
    select(
      Compound_Name,
      MW
    )
)

# Compound_Name     MW
# 1   Paromomycin  615.6
# 2      Neomycin  614.6
# 3  Cefoperazone  645.7
# 4   Ceftaroline  604.7
# 5  Azithromycin  749.0
# 6  Polymyxin b1 1203.5

# These were the six antibiotics from O'Shea et al. that were > 600 Da.


# Compounds with MW less than or equal to 600 Da (same as 600 g/mol which are the units used by pubchem).
pubchem_info_small <- pubchem_info %>%
  filter(Less_equal_600_da)

# Check the compounds that are >=600 Da.

# Convert to a data frame for display, as tibbles can print rounded values.

as.data.frame(
  pubchem_info_small %>%
    select(
      Compound_Name,
      MW
    )
)

# count how many will be included in the final dataset 

pubchem_info_small %>%
  summarise(
    total_antibiotics_count = sum(TYPE %in% c("GN_antibiotics_oshea", "antibiotics_other")),
    total_pollutants_count = sum(TYPE == "pollutant")
  )



# A tibble: 1 × 2
#total_antibiotics_count total_pollutants_count
#<int>                  <int>
#1                     119                    100


## ---------- 9) Save PubChem information as an Excel workbook to make Table_S2_1_Pubchem_Downloads ----------

# Set working directory before saving outputs.
setwd(here("Files", "Tables"))

# now lets order the df the wat we want it 
# check column names 
colnames(pubchem_info)

# change

pubchem_info_small <- pubchem_info_small%>%
  select(Compound_Name,
         TYPE,
         CID, 
         SMILES,
         CAS, 
         MW_Da = MW,
         CAS,               
         Less_equal_600_da)

# check 
head(pubchem_info_small)


# save the table Table_S2_1_Pubchem_Downloads
write.xlsx(
  pubchem_info_small,
  "Table_S2_1_Pubchem_Downloads.xlsx"
)

# Also save as a CSV file to use in later scripts.

setwd(here("Files", "Files_for_coding"))

write.csv(
  pubchem_info_small,
  "Pubchem_Downloads.csv",
  row.names = FALSE
)



