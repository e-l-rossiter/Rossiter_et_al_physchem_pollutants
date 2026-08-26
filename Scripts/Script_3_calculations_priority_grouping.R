
# ============================================================
# Script 3
# Linear interpolation, relative PSA, eNTRY rules and grouping
# ============================================================
#
# Purpose:
#
# This script uses the final joined PubChem, Chemicalize and eNTRYway output
# created in Script 2
#
# The script calculates additional physicochemical variables needed for the
# analysis, including:
#
# - interpolated logD values at  7.4, 7.5 and 8.5
# - relative polar surface area
# - hydrophilic compound status
# - low globularity status
# - low flexibility status
# - eNTRY compound status
#
# Chemicalize bulk predictions only returned cLogD values at integer pH values,
# so logD at pH 7.4 was calculated by linear interpolation between logD at pH 7
# and logD at pH 8. 
#
# Ionisation at pH 7.4 was described using compound net charge, and the major
# microspecies was defined as the species present at the highest fractional
# abundance at pH 7.4. As these values were not included in the Chemicalize bulk
# upload, they were manually transcribed from the Chemicalize web-based
# interface and joined into the main data frame.
#
# For compounds ionised at pH 7.4 but lacking an amine, the presence of an
# ionisable nitrogen was also assessed manually. These manually checked
# functional-group assignments were then joined into the main data frame.
#
# The script also joins compound use/class information so that pollutants and
# antibiotics can be separated and summarised by compound class.
#
# The final data frame is then used to assign compounds into:
#
# - Group 1
# - Group 2
# - Group 3
# - No group
#
# These groups are based on hydrophilicity, eNTRY rule status and charge state
# at pH 7.4.
#
# Finally, the script creates separate tables for pollutants and antibiotics in
# each group, checks the compounds assigned to each group, and saves the final
# physicochemical-property table as an Excel file.

### ---------- 0) Setup ----------

# Clear the environment so that the script starts from a clean workspace.
rm(list = ls())

# Load required packages.
library(dplyr)     # data wrangling: filter, mutate, rename, joins and pipes
library(readxl)    # read Excel files
library(writexl)   # save simple Excel files
library(openxlsx)  # create Excel workbooks with multiple sheets
library(here)      # make file paths portable across different computers
# by using paths relative to the project root
library(stringr)      # work with text strings, including detecting text patterns and changing case

#Get package citations

# Get the citations for the R packages used in this script.
# These can be used in the methods section or supplementary information if needed.

citation("dplyr")
citation("readxl")
citation("writexl")
citation("openxlsx")
citation("here")
citation("stringr")



# Check R version and session information

# Check R version.
R.version.string
# "R version 4.4.2 (2024-10-31 ucrt)"

# For this part, we are using an older version of R as some packages may not
# work with newer versions.

# Record session information so the package versions used are documented.
sessionInfo()


### ---------- 1) Read in joined PubChem, Chemicalize and eNTRYway file ----------


# Set working directory
setwd(here("Files", "Files_for_coding"))

# List files

list.files()

# Read in the joined table created in Script 2
joined_df <- read_xlsx("joined_chemicalize_eNTRy_way_pub_chem.xlsx")


# Make a working copy of the joined data frame.
joined_df_2 <- joined_df

### ---------- 2) Calculate interpolated logD values ----------

# Chemicalize bulk predictions returned cLogD values only at integer pH values.
# Therefore, logD at pH 7.4 was calculated by linear interpolation between
# logD at pH 7 and logD at pH 8.
#
# Additional interpolated values were also calculated for pH 6.5, 7.5 and 8.5.

# Convert the relevant logD columns to numeric before interpolation.
logD6 <- as.numeric(joined_df_2[["logd_ph6"]])
logD7 <- as.numeric(joined_df_2[["logd_ph7"]])
logD8 <- as.numeric(joined_df_2[["logd_ph8"]])


# Linear interpolation between logD7 and logD8 to get pH 7.4
joined_df_2$logD_pH7_4_interp <- logD7 + 0.4 * (logD8 - logD7)


# ---------- 3.4) Reorder columns and calculate relative PSA ----------

# Check column names after adding interpolated logD values.
colnames(joined_df_2)


# Reorder the joined data frame so that the main physicochemical variables are 
# grouped together.
joined_df_2 <- joined_df_2 %>%
  dplyr::select(
    Compound_Name,
    TYPE,
    molar_mass_g_mol, 
    exact_mass_da,
    formula,
    lipinski_ro5,
    atom_count, 
    heavy_atom_count,
    asymmetric_atom_count, 
    rotatable_bond_count,
    ring_count,
    aromatic_ring_count,
    hetero_ring_count,
    fsp3,
    hbd_count,
    hba_count,
    formal_charge,
    tpsa_a2, 
    polarizability_a3,
    molar_refractivity_cm3_mol, 
    vdw_volume_a3,
    min_projection_area_a2, 
    max_projection_area_a2,
    min_projection_radius_a, 
    max_projection_radius_a,  
    vdw_surface_area_a2, 
    solvent_accessible_surface_area_a2,
    pka_acid_strongest, 
    pka_base_strongest,
    isoelectric_point,
    intrinsic_solubility_logS, 
    solubility_category, 
    logp,
    logd_ph0,
    logd_ph1,
    logd_ph2, 
    logd_ph3,
    logd_ph4,
    logd_ph5,
    logd_ph6,
    logd_ph7, 
    logD_pH7_4_interp,
    logd_ph8,
    logd_ph9,
    logd_ph10,
    logd_ph11, 
    logd_ph12,
    logd_ph13,
    logd_ph14,
    hlb, 
    iupac_name,
    smiles, 
    inchi,  
    inchikey,
    traditional_name,
    common_names,
    cas_registry_numbers,
    smiles_ew,   
    formula_EW,
    molwt_EW,
    rb_EW,
    glob_EW, 
    pbf_EW,
    func_group_EW,
    CID,
    smiles_PC,
    CAS_PC, 
    MW_Da_PC)

# Calculate relative polar surface area.
# This is calculated as topological polar surface area divided by exact mass.
tpsa <- as.numeric(joined_df_2[["tpsa_a2"]])
vdw_surface_area_a2 <- as.numeric(joined_df_2[["vdw_surface_area_a2"]])

joined_df_2$rel_polar_sa <- tpsa / vdw_surface_area_a2


# Reorder the joined data frame so that the main physicochemical variables are
# grouped together.

joined_df_2 <- joined_df_2 %>%
  dplyr::select(
    Compound_Name,
    molar_mass_g_mol, 
    exact_mass_da,
    formula,
    lipinski_ro5,
    atom_count, 
    heavy_atom_count,
    asymmetric_atom_count, 
    rotatable_bond_count,
    ring_count,
    aromatic_ring_count,
    hetero_ring_count,
    fsp3,
    hbd_count,
    hba_count,
    formal_charge,
    tpsa_a2,
    rel_polar_sa,
    polarizability_a3,
    molar_refractivity_cm3_mol, 
    vdw_volume_a3,
    min_projection_area_a2, 
    max_projection_area_a2,
    min_projection_radius_a, 
    max_projection_radius_a,  
    vdw_surface_area_a2, 
    solvent_accessible_surface_area_a2,
    pka_acid_strongest, 
    pka_base_strongest,
    isoelectric_point,
    intrinsic_solubility_logS, 
    solubility_category, 
    logp,
    logd_ph0,
    logd_ph1,
    logd_ph2, 
    logd_ph3,
    logd_ph4,
    logd_ph5,
    logd_ph6,
    logd_ph7, 
    logD_pH7_4_interp,
    logd_ph8,
    logd_ph9,
    logd_ph10,
    logd_ph11, 
    logd_ph12,
    logd_ph13,
    logd_ph14,
    hlb, 
    iupac_name,
    smiles, 
    inchi,  
    inchikey,
    traditional_name,
    common_names,
    cas_registry_numbers,
    smiles_ew,   
    formula_EW,
    molwt_EW,
    rb_EW,
    glob_EW, 
    pbf_EW,
    func_group_EW,
    CID,
    smiles_PC,
    CAS_PC, 
    MW_Da_PC)




### ---------- 5) Read in manually transcribed charge information ----------

# Ionisation at pH 7.4 was described using compound net charge.
# The major microspecies was defined as the species present at the highest
# fractional abundance at pH 7.4.
#
# These values were not included in the Chemicalize bulk upload, so they were
# manually transcribed from the Chemicalize web-based interface.

# Set working directory to where the manually transcribed files are saved.

setwd(here("Files", "Tables"))

# List files to check that the charge workbook is present.
list.files()

# Read in manually transcribed charge information.
# Sheet 1 contains pollutants and sheet 2 contains antibiotics.
df_pols_charges <- read_excel("Table_S2_4_Predominant_charge.xlsx", sheet = 1)
df_abs_charges  <- read_excel("Table_S2_4_Predominant_charge.xlsx", sheet = 2)

# Check column names.
colnames(df_pols_charges)
colnames(df_abs_charges)

# Combine pollutant and antibiotic charge information into one table.
df_charges <- dplyr::bind_rows(
  df_pols_charges,
  df_abs_charges
)



# Join manually transcribed charge information to the main data frame.

df_joined_final <- joined_df_2 %>%
  left_join(df_charges, by = "Compound_Name")

# Check whether any compounds are missing dominant charge information.
missing_df <- df_joined_final %>%
  filter(is.na(`per_dominant_charge_pH_7_4`)) %>%
  distinct(Compound_Name) %>%
  pull(Compound_Name)

missing_df
# character(0)

#no missing compounds


### ---------- 6) Read in manual functional group checks and compound use groups etc. ----------
# List files to check that the ionisable N  workbook is present.
list.files()
# For compounds ionised at pH 7.4 but lacking an amine, the presence of an
# ionisable nitrogen was assessed manually.

# Set working directory to where the manually transcribed files are saved.

setwd(here("Files", "Tables"))

# This file is Table S2.6 in the paper.
functional_groups_other <- read_excel("Table_S2_6_ionisable_N.xlsx")

# This file contains the uses/classes of compounds.
# It is the same as Tables 1 and 2 in the paper, but modified to help with
# plotting and downstream data wrangling.

setwd(here("Files", "Files_for_coding"))

Uses_of_compounds <- read_excel("Uses_of_compounds.xlsx")

# Check column names.
colnames(functional_groups_other)
colnames(Uses_of_compounds)

# rename Name to Compound_Name in functional_groups_other

functional_groups_other<-functional_groups_other%>%
  mutate(Compound_Name = Name)
# check that worked

colnames(Uses_of_compounds)


# Join the manually checked functional group information.
df_joined_final_2 <- df_joined_final %>%
  left_join(
    functional_groups_other,
    by = "Compound_Name"
  )


# Join compound use/class information.
df_joined_final_3 <- df_joined_final_2 %>%
  left_join(
    Uses_of_compounds,
    by = "Compound_Name"
  )


### ---------- 7) Assign eNTRY rules ----------

# eNTRY compounds are assigned using the following criteria:
#
# - rb_EW <= 5
# - glob_EW <= 0.25
# - func_group_other is one of:
#   "Primary Amine", "Secondary Amine", "Tertiary Amine", or "Other"
#
# Additional TRUE/FALSE are also created for hydrophilicity,
# low globularity and low flexibility.
# this is for when we assign groups later

df_joined_final_4 <- df_joined_final_3 %>%
  mutate(
    # TRUE if logD at pH 7.4 is less than or equal to 0.
    hydrophilic_compound = logD_pH7_4_interp <= 0,
    
    # TRUE if globularity from eNTRYway is less than or equal to 0.25.
    lowglob = glob_EW <= 0.25,
    
    # TRUE if rotatable bonds from eNTRYway are less than or equal to 5.
    lowflex = rb_EW <= 5,
    
    # TRUE if the compound meets the eNTRY rule criteria.
    ENTRY_compound = !is.na(rb_EW) & !is.na(glob_EW) &
      rb_EW <= 5 &
      glob_EW <= 0.25 &
      func_group_other %in% c(
        "Primary Amine",
        "Secondary Amine",
        "Tertiary Amine", 
        "Other"
      )
  )

# Check the eNTRY and hydrophilic flags.
df_joined_final_4 %>%
  select(
    ENTRY_compound,
    hydrophilic_compound
  )

### ---------- 8) Assign charge TRUE/FALSE variables ----------

# Assign negative charge TRUE/FALSE.
# This includes zwitterions because they have both a positive and negative charge.
df_joined_final_5 <- df_joined_final_4 %>%
  mutate(
    Neg_charge_T_F = Dominant_charge_pH_7_4_code %in% c("-", "--", "+-", "++--")
  )

# Assign positive charge TRUE/FALSE.
# This also includes zwitterions because they have both a positive and negative charge.
df_joined_final_5 <- df_joined_final_5 %>%
  mutate(
    Positive_charge_T_F = Dominant_charge_pH_7_4_code %in% c(
      "+", "++", "+++", "++++", "+++++", "+-", "++--"
    )
  )

# Assign zwitterionic charge TRUE/FALSE.
df_joined_final_5 <- df_joined_final_5 %>%
  mutate(
    Zwitt_charge_T_F = Dominant_charge_pH_7_4_code %in% c("+-", "++--")
  )

# Assign uncharged TRUE/FALSE.
df_joined_final_5 <- df_joined_final_5 %>%
  mutate(
    Uncharged_T_F = Dominant_charge_pH_7_4_code %in% c("0")
  )

# Check the charge variables.
df_joined_final_5 %>%
  select(
    ENTRY_compound,
    hydrophilic_compound,
    Neg_charge_T_F,
    Positive_charge_T_F,
    Zwitt_charge_T_F,
    Uncharged_T_F
  )

### ---------- 9) Assign Group 1, Group 2, Group 3 and No group ----------

# Assign compounds into groups based on hydrophilicity, eNTRY rule status and
# charge state at pH 7.4.
#
# Group 1:
# Hydrophilic eNTRY compounds, including compounds that are positively charged
# or zwitterionic and meet the eNTRY criteria.
#
# Group 2:
# Hydrophilic compounds that are positively charged or zwitterionic, but do not
# meet the eNTRY criteria.
#
# Group 3:
# Non-hydrophilic compounds that meet the eNTRY criteria.
#
# No group:
# Compounds that do not meet the above criteria.

df_joined_final_6 <- df_joined_final_5 %>%
  mutate(
    Group = case_when(
      hydrophilic_compound == TRUE & ENTRY_compound == TRUE &
        !(Positive_charge_T_F == TRUE | Zwitt_charge_T_F == TRUE) ~ "Group 1",
      
      (Positive_charge_T_F == TRUE | Zwitt_charge_T_F == TRUE) &
        ENTRY_compound == TRUE & hydrophilic_compound == TRUE ~ "Group 1",
      
      (Positive_charge_T_F == TRUE | Zwitt_charge_T_F == TRUE) &
        ENTRY_compound == FALSE & hydrophilic_compound == TRUE ~ "Group 2",
      
      hydrophilic_compound == FALSE & ENTRY_compound == TRUE ~ "Group 3",
      
      TRUE ~ "No group"
    )
  )



# Check column names.
colnames(df_joined_final_6)

unique(df_joined_final_6$Dominant_charge_pH_7_4_code)

# make a column which gives the chemical species at pH 7.4 in words

df_joined_final_6 <- df_joined_final_6 %>%
  mutate(
    spp_7_4 = case_when(
      Dominant_charge_pH_7_4_code %in% c("-", "--") ~ "Anion",
      Dominant_charge_pH_7_4_code == "0" ~ "Neutral",
      Dominant_charge_pH_7_4_code %in% c("+", "++", "+++", "++++", "+++++") ~ "Cation",
      Dominant_charge_pH_7_4_code == "+-" ~ "Zwitterion",
      TRUE ~ NA_character_
    )
  )

# check it 
df_joined_final_6 %>%
  count(Dominant_charge_pH_7_4_code, spp_7_4)


# A tibble: 9 × 3
# Dominant_charge_pH_7_4_code spp_7_4        n
# <chr>                       <chr>      <int>
# 1 +                           cation        20
# 2 ++                          cation         1
# 3 +++                         cation         2
# 4 ++++                        cation         4
# 5 +++++                       cation         4
# 6 +-                          zwitterion    36
# 7 -                           anion         72
# 8 --                          anion         16
# 9 0                           neutral       64

colnames(df_joined_final_6)

# make Cephradine and Cephalexin match the rest of the cephms 

df_joined_final_6 <- df_joined_final_6 %>%
  mutate(
    Compound_Name = case_when(
      Compound_Name == "Cephradine" ~ "Cefradine",
      Compound_Name == "Cephalexin" ~ "Cefalexin",
      TRUE ~ Compound_Name
    )
  )

# check it worked 
df_joined_final_6 %>%
  filter(Compound_Name %in% c("Cefradine", "Cefalexin", "Cephradine", "Cephalexin")) %>%
  select(Compound_Name)

# # A tibble: 2 × 1
# Compound_Name     
# <chr>    
# 1 Cefalexin
# 2 Cefradine
#
# it worked


### ---------- 10) Make smaller grouped frame for checking and making tables for paper ----------

# Make a smaller data frame containing the main variables needed to check group
# assignment.
df_joined_final_6_small <- df_joined_final_6 %>%
  select(
    Group,
    CODE_2,
    Class_abbrv,
    Compound_Name,
    exact_mass_da,
    rel_polar_sa,
    hbd_count,
    hba_count,
    logD_pH7_4_interp,
    logp,
    rb_EW,
    rotatable_bond_count,
    glob_EW,
    per_dominant_charge_pH_7_4,
    Charge_pH_7_4_num,
    func_group_updated,
    spp_7_4
  ) %>%
  mutate(
    exact_mass_da = round(as.numeric(exact_mass_da), 0),
    rel_polar_sa = round(as.numeric(rel_polar_sa), 2),
    logD_pH7_4_interp = round(as.numeric(logD_pH7_4_interp), 2),
    logp = round(as.numeric(logp), 2),
    glob_EW = round(as.numeric(glob_EW), 2),
    per_dominant_charge_pH_7_4 = round(as.numeric(per_dominant_charge_pH_7_4), 2),
    Charge_pH_7_4_num = round(as.numeric(Charge_pH_7_4_num), 2),
    hbd_count = as.integer(hbd_count),
    hba_count = as.integer(hba_count),
    rb_EW = as.integer(rb_EW),
    rotatable_bond_count = as.integer(rotatable_bond_count)
  )

df_joined_final_6_small



### ---------- 11) Split compounds by group ----------

# Create separate data frames for compounds assigned to Group 1, Group 2,
# Group 3 and No group. These are useful for checking the compounds assigned
# to each priority group.

group_1_compounds <- df_joined_final_6_small %>%
  filter(Group == "Group 1")

group_2_compounds <- df_joined_final_6_small %>%
  filter(Group == "Group 2")

group_3_compounds <- df_joined_final_6_small %>%
  filter(Group == "Group 3")

no_group_compounds <- df_joined_final_6_small %>%
  filter(Group == "No group")

#Check pollutant compounds by group

# Group 1 pollutants.
group_1_pols <- group_1_compounds %>%
  filter(CODE_2 == "Pollutant") %>%
  select(-CODE_2)

group_1_pols %>%
  select(
    Class_abbrv,
    Compound_Name
  ) %>%
  print(n = Inf)

# A tibble: 7 × 2
#Class_abbrv                    Compound_Name  
#<chr>                          <chr>          
#1 Non-antibiotic drug            Amphetamine    
#2 Non-antibiotic drug metabolite Benzoylecgonine
#3 Non-antibiotic drug            Codeine        
#4 Non-antibiotic drug metabolite Guanylurea     
#5 Non-antibiotic drug            Metformin      
#6 Industrial chemical            Methylamine    
#7 Non-antibiotic drug            Sitagliptin    

# Group 2 pollutants.
group_2_pols <- group_2_compounds %>%
  filter(CODE_2 == "Pollutant")%>%
  select(-CODE_2)

group_2_pols%>%
  select(Class_abbrv,
         Compound_Name)

# A tibble: 4 × 2
#Class_abbrv         Compound_Name
#<chr>               <chr>        
#1 Non-antibiotic drug Atenolol     
#2 Non-antibiotic drug Gabapentin   
#3 Non-antibiotic drug Sotalol      
#4 Non-antibiotic drug Morphine

# Group 3 pollutants.

group_3_pols <- group_3_compounds %>%
  filter(CODE_2 == "Pollutant")%>%
  select(-CODE_2)

group_3_pols%>%
  select(Class_abbrv,
         Compound_Name)

# A tibble: 7 × 2
#Class_abbrv         Compound_Name
#<chr>               <chr>        
#1 Non-antibiotic drug Amitriptyline
#2 Non-antibiotic drug Citalopram   
#3 Non-antibiotic drug Cocaine      
#4 Non-antibiotic drug Ketamine     
#5 Non-antibiotic drug Lidocaine    
#6 Non-antibiotic drug Sertraline   
#7 Non-antibiotic drug Venlafaxine  

# No group pollutants.

no_group_pols <- no_group_compounds %>%
  filter(CODE_2 == "Pollutant")%>%
  select(-CODE_2)

no_group_pols %>%
  select(
    Class_abbrv,
    Compound_Name
  ) %>%
  print(n = Inf)

#A tibble: 81 × 2
#Class_abbrv                    Compound_Name                                               
#<chr>                          <chr>                                                       
#1 Pesticide metabolite           ((4-Fluorophenyl)(propan-2-yl)carbamoyl)methanesulfonic acid
#2 Pesticide                      (2,4-Dichlorophenoxy)Acetic acid                            
#3 Non-antibiotic drug metabolite 10,11-dihydroxycarbamazepine                                
#4 Industrial chemical            1-Decanol                                                   
#5 Industrial chemical            2(3H)-Benzothiazolone                                       
#6 Pesticide metabolite           2,6-Dichlorobenzamide                                       
#7 Non-antibiotic drug metabolite 2-Hydroxyibuprofen                                          
#8 Pesticide                      2-methyl-4-chlorophenoxyacetic acid                         
#9 Industrial chemical            4-Phenoxybutyric acid                                       
#10 Industrial chemical            5-Amino-2-methylphenol                                      
#11 Industrial chemical            5-Methyl-1H-benzotriazole                                   
#12 Industrial chemical            Acenaphthene                                                
#13 Non-antibiotic drug            Acetaminophen                                               
#14 Pesticide                      Acetamiprid                                                 
#15 Non-antibiotic drug            Acetylsalicylic acid                                        
#16 Pesticide                      Aclonifen                                                   
#17 Industrial chemical            Aniline                                                     
#18 Pesticide                      Atrazine                                                    
#19 Pesticide                      Bentazone                                                   
#20 Industrial chemical            Benzophenone                                                
#21 Industrial chemical            1H-Benzotriazole                                            
#22 Lifestyle product              Caffeine                                                    
#23 Non-antibiotic drug            Carbamazepine                                               
#24 Non-antibiotic drug metabolite Carbamazepine-10,11-epoxide                                 
#25 Non-antibiotic drug            Cetirizine                                                  
#26 Pesticide                      Chlorantraniliprole                                         
#27 Pesticide metabolite           Chloridazon-desphenyl                                       
#28 Non-antibiotic drug metabolite Clofibric acid                                              
#29 Pesticide                      Clothianidin                                                
#30 Lifestyle product metabolite   Cotinine                                                    
#31 Pesticide                      Cyantraniliprole                                            
#32 Industrial chemical            Cyclohexanone                                               
#33 Pesticide                      Dicamba                                                     
#34 Pesticide                      Dichlobenil                                                 
#35 Non-antibiotic drug            Diclofenac                                                  
#36 Industrial chemical            Diphenylamine                                               
#37 Pesticide                      Diuron                                                      
#38 Industrial chemical            Ethylenediaminetetraacetic acid                             
#39 Non-antibiotic drug            Estradiol                                                   
#40 Non-antibiotic drug            Estriol                                                     
#41 Industrial chemical            Ethylbenzene                                                
#42 Pesticide                      Fipronil                                                    
#43 Pesticide metabolite           Fipronil sulfone                                            
#44 Pesticide                      Flufenacet                                                  
#45 Industrial chemical            Fluoranthene                                                
#46 Non-antibiotic drug            Fluoxetine                                                  
#47 Pesticide                      Flupyradifurone                                             
#48 Non-antibiotic drug            Gliclazide                                                  
#49 Non-antibiotic drug            Hydrocodone                                                 
#50 Non-antibiotic drug            Ibuprofen                                                   
#51 Pesticide                      Imidacloprid                                                
#52 Pesticide                      Isoprocarb                                                  
#53 Non-antibiotic drug            Lamotrigine                                                 
#54 Industrial chemical            Linalool                                                    
#55 Industrial chemical            Melamine                                                    
#56 Pesticide                      Metaflumizone                                               
#57 Pesticide                      Metazachlor                                                 
#58 Pesticide metabolite           Metazachlor ESA                                             
#59 Pesticide                      Methylchlorophenoxypropionic acid                           
#60 Non-antibiotic drug            Modafinil                                                   
#61 Industrial chemical            N,N,N',N'-Tetraacetylethylenediamine                        
#62 Pesticide                      N,N-diethyl-m-toluamide                                     
#63 Industrial chemical            N,N-Dimethylaniline                                         
#64 Non-antibiotic drug            Naproxen                                                    
#65 Industrial chemical            Nitrobenzene                                                
#66 Industrial chemical            p-Benzoquinone                                              
#67 Pesticide                      Propyzamide                                                 
#68 Non-antibiotic drug            Ranitidine                                                  
#69 Personal care product          Salicylic acid                                              
#70 Pesticide                      Schradan                                                    
#71 Pesticide                      Simazine                                                    
#72 Pesticide                      Spirotetramat                                               
#73 Personal care product          Sulisobenzone                                               
#74 Pesticide                      Terbutryn                                                   
#75 Industrial chemical            Terpineol                                                   
#76 Industrial chemical            Tonalide                                                    
#77 Non-antibiotic drug            Tramadol                                                    
#78 Industrial chemical            N-Phenyl-1-naphthylamine                                    
#79 Industrial chemical            2-Naphthylamine                                             
#80 Industrial chemical            1-Naphthylamine                                             
#81 Industrial chemical            N-nitrosodiphenylamine                             



#Check antibiotic compounds by group

# Group 1 antibiotics.

group_1_abs <- group_1_compounds %>%
  filter(CODE_2 == "Antibiotic")%>%
  select(-CODE_2)

group_1_abs%>%
  select(Class_abbrv,
         Compound_Name) %>%
  print(n = Inf)

# A tibble: 35 × 2
#Class_abbrv     Compound_Name    
#<chr>           <chr>            
#1 Cephm           Loracarbef       
#2 Cephm           Cefaclor         
#3 Cephm           Cefadroxil       
#4 Cephm           Cefprozil        
#5 Cephm           Cefalexin       
#6 Cephm           Cefradine       
#7 Fluoroquinolone Ciprofloxacin    
#8 Fluoroquinolone Clinafloxacin    
#9 Fluoroquinolone Danofloxacin     
#10 Fluoroquinolone Dx-619           
#11 Fluoroquinolone Enoxacin         
#12 Fluoroquinolone Fleroxacin       
#13 Fluoroquinolone Gatifloxacin     
#14 Fluoroquinolone Gemifloxacin     
#15 Fluoroquinolone Levofloxacin     
#16 Fluoroquinolone Lomefloxacin     
#17 Fluoroquinolone Moxifloxacin     
#18 Fluoroquinolone Norfloxacin      
#19 Fluoroquinolone Pefloxacin       
#20 Fluoroquinolone Rufloxacin       
#21 Fluoroquinolone Sitafloxacin     
#22 Penicillin      Amoxicillin      
#23 Penicillin      Ampicillin       
#24 Tetracycline    Demeclocycline   
#25 Tetracycline    Doxycycline      
#26 Tetracycline    Meclocycline     
#27 Tetracycline    Methacycline     
#28 Tetracycline    Minocycline      
#29 Tetracycline    Oxytetracycline  
#30 Tetracycline    Tetracycline     
#31 Tetracycline    Chlortetracycline
#32 Fluoroquinolone Ofloxacin        
#33 Penem           Biapenem         
#34 Tetracycline    Rolitetracycline 
#35 Tetracycline    Eravacycline


#35 antibiotics in group 1 

group_1_abs %>%
  count(Class_abbrv, sort = TRUE)

# A tibble: 5 × 2
#Class_abbrv         n
#<chr>           <int>
#1 Fluoroquinolone    16
#2 Tetracycline       10
#3 Cephm               6
#4 Penicillin          2
#5 Penem               1

#AB classes in Group 1 

# 16 Fluoroquinolones
#10 Tetracyclines
#9 Beta-lactams (#6 Cephms, #2 Penicillins, #1 Penem)



# Group 2 antibiotics.

group_2_abs <- group_2_compounds %>%
  filter(CODE_2 == "Antibiotic")%>%
  select(-CODE_2)


group_2_abs%>%
  select(Class_abbrv,
         Compound_Name) %>%
  print(n = Inf)

# A tibble: 17 × 2
#Class_abbrv    Compound_Name
#<chr>          <chr>        
#1 Aminoglycoside Kanamycin a  
#2 Aminoglycoside Streptomycin 
#3 Aminoglycoside Amikacin     
#4 Aminoglycoside Arbekacin    
#5 Aminoglycoside Dibekacin    
#6 Aminoglycoside Gentamicin   
#7 Aminoglycoside Isepamicin   
#8 Aminoglycoside Netilmicin   
#9 Aminoglycoside Sisomicin    
#10 Aminoglycoside Tobramycin   
#11 Penem          Imipenem     
#12 Penem          Meropenem    
#13 Penem          Tomopenem    
#14 Cephm          Cefepime     
#15 Cephm          Cefpirome    
#16 Penem          Doripenem    
#17 Penem          Lenapenem 

#17 ABs total in Group 2

group_2_abs %>%
  count(Class_abbrv, sort = TRUE)

# A tibble: 3 × 2
#Class_abbrv        n
#<chr>          <int>
# 1 Aminoglycoside    10
#2 Penem              5
#3 Cephm              2

#10 AGs, 7 BLs


# Group 3 antibiotics.

group_3_abs <- group_3_compounds %>%
  filter(CODE_2 == "Antibiotic")%>%
  select(-CODE_2)

group_3_abs
group_3_abs%>%
  select(Class_abbrv,
         Compound_Name) %>%
  print(n = Inf)

# A tibble: 6 × 2
#Class_abbrv     Compound_Name
#<chr>           <chr>        
#1 Fluoroquinolone Difloxacin   
#2 Fluoroquinolone Garenoxacin  
#3 Fluoroquinolone Grepafloxacin
#4 Fluoroquinolone Sparfloxacin 
#5 Fluoroquinolone Temafloxacin 
#6 Fluoroquinolone Trovafloxacin

#6 ABs in Group 3 all FQs

# No group antibiotics.

no_group_abs <- no_group_compounds %>%
  filter(CODE_2 == "Antibiotic")%>%
  select(-CODE_2)


no_group_abs%>%
  select(Class_abbrv,
         Compound_Name) %>%
  print(n = Inf)

# A tibble: 61 × 2
#Class_abbrv       Compound_Name         
#<chr>             <chr>                 
#1 Amphenicol        Chloramphenicol       
#2 Penem             Ertapenem             
#3 Cephm             Cefetamet             
#4 Cephm             Ceftibuten            
#5 Cephm             Cefamandole           
#6 Cephm             Cefazolin             
#7 Cephm             Cefdinir              
#8 Cephm             Cefditoren            
#9 Cephm             Cefixime              
#10 Cephm             Cefmetazole           
#11 Cephm             Cefotaxime            
#12 Cephm             Cefotetan             
#13 Cephm             Cefoxitin             
#14 Cephm             Cefpodoxime           
#15 Cephm             Ceftizoxime           
#16 Cephm             Ceftriaxone           
#17 Cephm             Cefuroxime            
#18 Cephm             Cephalothin           
#19 Cephm             Ceftazidime           
#20 Cephm             Ceftobiprole          
#21 Diaminopyrimidine Iclaprim              
#22 Diaminopyrimidine Trimethoprim          
#23 Quinolone         Nalidixic acid        
#24 Fluoroquinolone   Delafloxacin          
#25 Fluoroquinolone   Nadifloxacin          
#26 Monobactam        Aztreonam             
#27 Phosphonic acid   Fosfomycin            
#28 Penem             Faropenem             
#29 Penicillin        Carbenicillin         
#30 Penicillin        Mezlocillin           
#31 Penicillin        Ticarcillin           
#32 Penicillin        Azlocillin            
#33 Penicillin        Piperacillin          
#34 Sulphonamide      Sulfabenzamide        
#35 Sulphonamide      Sulfacetamide         
#36 Sulphonamide      Sulfachlorpyridazine  
#37 Sulphonamide      Sulfadiazine          
#38 Sulphonamide      Sulfadimethoxine      
#39 Sulphonamide      Sulfaguanidine        
#40 Sulphonamide      Sulfamerazine         
#41 Sulphonamide      Sulfameter            
#42 Sulphonamide      Sulfamethazine        
#43 Sulphonamide      Sulfamethizole        
#44 Sulphonamide      Sulfamethoxazole      
#45 Sulphonamide      Sulfamethoxypyridazine
#46 Sulphonamide      Sulfamonomethoxine    
#47 Sulphonamide      Sulfanitran           
#48 Sulphonamide      Sulfaphenazole        
#49 Sulphonamide      Sulfapyridine         
#50 Sulphonamide      Sulfaquinoxaline      
#51 Sulphonamide      Sulfathiazole         
#52 Sulphonamide      Sulfisoxazole         
#53 Tetracycline      Omadacycline          
#54 Tetracycline      Tigecycline           
#55 Monobactam        Tigemonam             
#56 Penem             Sulopenem             
#57 Penem             Tebipenem             
#58 Diaminopyrimidine Tetroxoprim           
#59 Diaminopyrimidine Aditoprim             
#60 Amphenicol        Florfenicol           
#61 Amphenicol        Thiamphenicol

no_group_abs %>%
  count(Class_abbrv, sort = TRUE)
# A tibble: 11 × 2
#Class_abbrv           n
#<chr>             <int>
# 1 Sulphonamide         19
#2 Cephm                18
#3 Penicillin            5
#4 Diaminopyrimidine     4
#5 Penem                 4
#6 Amphenicol            3
#7 Fluoroquinolone       2
#8 Monobactam            2
#9 Tetracycline          2
#10 Phosphonic acid       1
#11 Quinolone             1

#create antibiotic class counts 
antibiotic_class_counts <- bind_rows(
  group_1_abs %>%
    count(Class_abbrv, name = "n") %>%
    mutate(Group = "Group 1"),
  
  group_2_abs %>%
    count(Class_abbrv, name = "n") %>%
    mutate(Group = "Group 2"),
  
  group_3_abs %>%
    count(Class_abbrv, name = "n") %>%
    mutate(Group = "Group 3"),
  
  no_group_abs %>%
    count(Class_abbrv, name = "n") %>%
    mutate(Group = "No group")
) %>%
  select(Group, Class_abbrv, n) %>%
  arrange(
    factor(Group, levels = c("Group 1", "Group 2", "Group 3", "No group")),
    desc(n),
    Class_abbrv
  )

antibiotic_class_counts

# A tibble: 20 × 3
#Group    Class_abbrv           n
#<chr>    <chr>             <int>
#  1 Group 1  Fluoroquinolone      16
#2 Group 1  Tetracycline         10
#3 Group 1  Cephm                 6
#4 Group 1  Penicillin            2
#5 Group 1  Penem                 1
#6 Group 2  Aminoglycoside       10
#7 Group 2  Penem                 5
#8 Group 2  Cephm                 2
#9 Group 3  Fluoroquinolone       6
#10 No group Sulphonamide         19
#11 No group Cephm                18
#12 No group Penicillin            5
#13 No group Diaminopyrimidine     4
#14 No group Penem                 4
#15 No group Amphenicol            3
#16 No group Fluoroquinolone       2
#17 No group Monobactam            2
#18 No group Tetracycline          2
#19 No group Phosphonic acid       1
#20 No group Quinolone             1


### ---------- 12) Create prioritised pollutant table ----------

# Create a table of prioritised pollutants to save as Table 5.
#
# This table keeps only compounds which:
# - are pollutants
# - are assigned to Group 1, Group 2 or Group 3
#
# Compounds in "No group" are excluded.

Prioritised_pollutants <- df_joined_final_6_small %>%
  filter(
    CODE_2 == "Pollutant",
    Group %in% c("Group 1", "Group 2", "Group 3")
  ) %>%
  
  # Order the priority groups so that Group 1 appears first, followed by
  # Group 2 and Group 3.
  #
  # Add a blank Description column so descriptions can be filled in manually
  # for Table 5.
  mutate(
    Group = factor(Group, levels = c("Group 1", "Group 2", "Group 3")),
    Description = ""
  ) %>%
  
  # Sort the table by priority group.
  arrange(Group)


# Format percentage columns.
#
# rel_polar_sa and per_dominant_charge_pH_7_4 are stored as proportions, so they
# are multiplied by 100 and rounded to one decimal place.

Prioritised_pollutants <- Prioritised_pollutants %>%
  mutate(
    rel_polar_sa = round(rel_polar_sa * 100, 1),
    per_dominant_charge_pH_7_4 = round(per_dominant_charge_pH_7_4 * 100, 1)
  ) %>%
  select(-CODE_2)



# check col names
colnames(Prioritised_pollutants)

# Select and rename columns for Table 5.

Prioritised_pollutants <- Prioritised_pollutants %>%
  select(
    Priority_Group     = Group,
    Name               = Compound_Name,
    Description,
    MW_Da              = exact_mass_da,
    Rel_PSA_percent    = rel_polar_sa,
    HBD                = hbd_count,
    HBA                = hba_count,
    RB                 = rb_EW, 
    Glob               = glob_EW,
    clogP              = logp,
    clogD_7_4          = logD_pH7_4_interp,
    spp_7_4,
    spp_7_4_percent    = per_dominant_charge_pH_7_4,
    Net_charge_7_4     = Charge_pH_7_4_num,
    Ionisable_Nitrogen = func_group_updated
  )


# Manually add pollutant descriptions for Table 5.

Prioritised_pollutants <- Prioritised_pollutants %>%
  mutate(
    Description = case_when(
      Name == "Methylamine" ~ "Industrial chemical; industrial feedstock; natural product",
      Name == "Amphetamine" ~ "NAD; CNS stimulant; psychoactive compound",
      Name == "Sitagliptin" ~ "NAD; antidiabetic drug",
      Name == "Guanylurea" ~ "NAD metabolite; metformin metabolite",
      Name == "Metformin" ~ "NAD; antidiabetic drug",
      Name == "Codeine" ~ "NAD; analgesic; opioid; psychoactive compound",
      Name == "Benzoylecgonine" ~ "NAD metabolite; cocaine metabolite",
      Name == "Gabapentin" ~ "NAD; antiepileptic drug; psychoactive compound",
      Name == "Atenolol" ~ "NAD; antihypertensive drug; β-blocker",
      Name == "Sotalol" ~ "NAD; antihypertensive drug; β-blocker",
      Name == "Morphine" ~ "NAD; analgesic; opioid; psychoactive compound",
      Name == "Ketamine" ~ "NAD; anaesthetic; antidepressant; psychoactive compound",
      Name == "Sertraline" ~ "NAD; antidepressant; SSRI; psychoactive compound",
      Name == "Amitriptyline" ~ "NAD; antidepressant; TCA; psychoactive compound",
      Name == "Lidocaine" ~ "NAD; anaesthetic",
      Name == "Citalopram" ~ "NAD; antidepressant; SSRI; psychoactive compound",
      Name == "Venlafaxine" ~ "NAD; antidepressant; SNRI; psychoactive compound",
      Name == "Cocaine" ~ "NAD; anaesthetic; CNS stimulant; psychoactive compound",
      TRUE ~ Description
    )
  )


# Check that descriptions were added correctly.

Prioritised_pollutants %>%
  select(
    Priority_Group,
    Name,
    Description
  ) %>%
  print(n = Inf)


# Save prioritised pollutants as Table 5.

setwd(here("Files", "Tables"))

write.xlsx(
  Prioritised_pollutants,
  "Table_5_Prioritised_pollutants.xlsx"
)



### ---------- 13) Create prioritised antibiotic tables ----------

# Create a table of prioritised antibiotics to save as supplementary tables.
#
# This table keeps only compounds which:
# - are antibiotics
# - are assigned to Group 1, Group 2 or Group 3
#
# Compounds in "No group" are excluded.

Prioritised_antibiotics <- df_joined_final_6_small %>%
  filter(
    CODE_2 == "Antibiotic",
    Group %in% c("Group 1", "Group 2", "Group 3")
  ) %>%
  
  # Order the priority groups so that Group 1 appears first, followed by
  # Group 2 and Group 3.
  #
  # Add blank MOA and Class columns so these can be filled in manually.
  mutate(
    Group = factor(Group, levels = c("Group 1", "Group 2", "Group 3")),
    MOA = "",
    Class = ""
  ) %>%
  
  # Sort by priority group, then by antibiotic class and compound name.
  arrange(Group, Class_abbrv, Compound_Name)


# Convert percentage columns from proportions to percentages and round to one
# decimal place.

Prioritised_antibiotics <- Prioritised_antibiotics %>%
  mutate(
    rel_polar_sa = round(rel_polar_sa * 100, 1),
    per_dominant_charge_pH_7_4 = round(per_dominant_charge_pH_7_4 * 100, 1)
  ) %>%
  select(-CODE_2)


# Check column names and antibiotic classes before formatting.
colnames(Prioritised_antibiotics)
unique(Prioritised_antibiotics$Class_abbrv)

# Select and rename columns for the prioritised antibiotic supplementary tables.

Prioritised_antibiotics <- Prioritised_antibiotics %>%
  select(
    Priority_Group     = Group,
    Name               = Compound_Name,
    Class_abbrv,
    Class,
    MOA,
    MW_Da              = exact_mass_da,
    Rel_PSA_percent    = rel_polar_sa,
    HBD                = hbd_count,
    HBA                = hba_count,
    RB                 = rb_EW, # this is the same as rb_EW
    Glob               = glob_EW,
    clogP              = logp,
    clogD_7_4          = logD_pH7_4_interp,
    spp_7_4,
    spp_7_4_percent    = per_dominant_charge_pH_7_4,
    Net_charge_7_4     = Charge_pH_7_4_num,
    Ionisable_Nitrogen = func_group_updated
  )


# Split prioritised antibiotics into separate data frames by priority group.
# These will become supplementary Tables S3.1, S3.2 and S3.3.

# group 1 

Prioritised_antibiotics_g1 <- Prioritised_antibiotics %>%
  filter(Priority_Group == "Group 1")

Prioritised_antibiotics_g2 <- Prioritised_antibiotics %>%
  filter(Priority_Group == "Group 2")

Prioritised_antibiotics_g3 <- Prioritised_antibiotics %>%
  filter(Priority_Group == "Group 3")

# Check compound names before manually adding Class and MOA.
unique(Prioritised_antibiotics_g1$Name)
unique(Prioritised_antibiotics_g2$Name)
unique(Prioritised_antibiotics_g3$Name)

### ---------- 13.1) Add Class and MOA for Group 1 antibiotics ----------

# Manually fill in antibiotic class and mechanism of action for Group 1
# antibiotics.
#
# These values come from Table 2 and S1.2: Gram-negative antibiotics evaluated in this
# study.

ab_g1_lookup <- tribble(
  ~Name, ~Class_new, ~MOA_new,
  "Cefadroxil",        "Cephm - Cephalosporin", "Cell wall biosynthesis",
  "Cefaclor",          "Cephm - Cephalosporin", "Cell wall biosynthesis",
  "Loracarbef",        "Cephm - Carbacephem",   "Cell wall biosynthesis",
  "Cefradine",         "Cephm - Cephalosporin", "Cell wall biosynthesis",
  "Cefalexin",         "Cephm - Cephalosporin", "Cell wall biosynthesis",
  "Amoxicillin",       "Penicillin",            "Cell wall biosynthesis",
  "Ampicillin",        "Penicillin",            "Cell wall biosynthesis",
  "Cefprozil",         "Cephm - Cephalosporin", "Cell wall biosynthesis",
  "Biapenem",          "Penem - Carbapenem",    "Cell wall biosynthesis",
  "Clinafloxacin",     "Fluoroquinolone",       "DNA synth",
  "Sitafloxacin",      "Fluoroquinolone",       "DNA synth",
  "Gemifloxacin",      "Fluoroquinolone",       "DNA synth",
  "Dx-619",            "Fluoroquinolone",       "DNA synth",
  "Ciprofloxacin",     "Fluoroquinolone",       "DNA synth",
  "Enoxacin",          "Fluoroquinolone",       "DNA synth",
  "Norfloxacin",       "Fluoroquinolone",       "DNA synth",
  "Lomefloxacin",      "Fluoroquinolone",       "DNA synth",
  "Gatifloxacin",      "Fluoroquinolone",       "DNA synth",
  "Moxifloxacin",      "Fluoroquinolone",       "DNA synth",
  "Levofloxacin",      "Fluoroquinolone",       "DNA synth",
  "Ofloxacin",         "Fluoroquinolone",       "DNA synth",
  "Rufloxacin",        "Fluoroquinolone",       "DNA synth",
  "Pefloxacin",        "Fluoroquinolone",       "DNA synth",
  "Danofloxacin",      "Fluoroquinolone",       "DNA synth",
  "Fleroxacin",        "Fluoroquinolone",       "DNA synth",
  "Doxycycline",       "Tetracycline",          "30S Ribosome",
  "Demeclocycline",    "Tetracycline",          "30S Ribosome",
  "Meclocycline",      "Tetracycline",          "30S Ribosome",
  "Chlortetracycline", "Tetracycline",          "30S Ribosome",
  "Tetracycline",      "Tetracycline",          "30S Ribosome",
  "Methacycline",      "Tetracycline",          "30S Ribosome",
  "Oxytetracycline",   "Tetracycline",          "30S Ribosome",
  "Minocycline",       "Tetracycline",          "30S Ribosome",
  "Rolitetracycline",  "Tetracycline",          "30S Ribosome",
  "Eravacycline",      "Tetracycline",          "30S Ribosome"
)


# Join the lookup table and fill Class and MOA.

Prioritised_antibiotics_g1 <- Prioritised_antibiotics_g1 %>%
  left_join(
    ab_g1_lookup,
    by = "Name"
  ) %>%
  mutate(
    Class = coalesce(Class_new, Class),
    MOA = coalesce(MOA_new, MOA)
  ) %>%
  select(
    -Class_new,
    -MOA_new
  )


# Check for any compounds where Class or MOA did not join correctly.

Prioritised_antibiotics_g1 %>%
  filter(is.na(MOA) | MOA == "") %>%
  select(Name, Class, MOA)

# A tibble: 0 × 3
# This means all Group 1 antibiotics matched successfully.

### ---------- 13.2) Add Class and MOA for Group 2 antibiotics ----------

# Manually fill in antibiotic class and mechanism of action for Group 2
# antibiotics.
#
# These values come from Table 2 and S1.2: Gram-negative antibiotics evaluated in this
# study.

# manually fill in MOA and class

# Create lookup table for Group 2 antibiotic class and MOA.
ab_g2_lookup <- tribble(
  ~Name, ~Class_new, ~MOA_new,
  "Arbekacin",      "Aminoglycoside",              "30S Ribosome",
  "Dibekacin",      "Aminoglycoside",              "30S Ribosome",
  "Sisomicin",      "Aminoglycoside",              "30S Ribosome",
  "Tobramycin",     "Aminoglycoside",              "30S Ribosome",
  "Gentamicin",     "Aminoglycoside",              "30S Ribosome",
  "Netilmicin",     "Aminoglycoside",              "30S Ribosome",
  "Amikacin",       "Aminoglycoside",              "30S Ribosome",
  "Isepamicin",     "Aminoglycoside",              "30S Ribosome",
  "Kanamycin a",    "Aminoglycoside",              "30S Ribosome",
  "Streptomycin",   "Aminoglycoside",              "30S Ribosome",
  "Cefpirome",      "Cephm - Cephalosporin",       "Cell wall biosynthesis",
  "Cefepime",       "Cephm - Cephalosporin",       "Cell wall biosynthesis",
  "Doripenem",      "Penem - Carbapenem",          "Cell wall biosynthesis",
  "Meropenem",      "Penem - Carbapenem",          "Cell wall biosynthesis",
  "Imipenem",       "Penem - Carbapenem",          "Cell wall biosynthesis",
  "Lenapenem",      "Penem - Carbapenem",          "Cell wall biosynthesis",
  "Tomopenem",      "Penem - Carbapenem",          "Cell wall biosynthesis"
)

# Join the lookup table and fill Class and MOA.

Prioritised_antibiotics_g2 <- Prioritised_antibiotics_g2 %>%
  left_join(
    ab_g2_lookup,
    by = "Name"
  ) %>%
  mutate(
    Class = coalesce(Class_new, Class),
    MOA = coalesce(MOA_new, MOA)
  ) %>%
  select(
    -Class_new,
    -MOA_new
  )

# Check for any compounds where Class or MOA did not join corre

Prioritised_antibiotics_g2 %>%
  filter(is.na(MOA) | MOA == "") %>%
  select(Name, Class, MOA)

# A tibble: 0 × 3
# ℹ 3 variables: Name <chr>, Class <chr>, MOA <chr>

# it worked


### ---------- 13.3) Add Class and MOA for Group 2 antibiotics ----------

# Manually fill in antibiotic class and mechanism of action for Group 3
# antibiotics.
#
# These values come from Table 2 and S1.2: Gram-negative antibiotics evaluated in this
# study.


ab_g3_lookup <- tribble(
  ~Name, ~Class_new, ~MOA_new,
  "Trovafloxacin",  "Fluoroquinolone", "DNA synth",
  "Grepafloxacin",  "Fluoroquinolone", "DNA synth",
  "Sparfloxacin",   "Fluoroquinolone", "DNA synth",
  "Temafloxacin",   "Fluoroquinolone", "DNA synth",
  "Garenoxacin",    "Fluoroquinolone", "DNA synth",
  "Difloxacin",     "Fluoroquinolone", "DNA synth"
)

# Join the lookup table and fill Class and MOA.
Prioritised_antibiotics_g3 <- Prioritised_antibiotics_g3 %>%
  left_join(
    ab_g3_lookup,
    by = "Name"
  ) %>%
  mutate(
    Class = coalesce(Class_new, Class),
    MOA = coalesce(MOA_new, MOA)
  ) %>%
  select(
    -Class_new,
    -MOA_new
  )

# Check for any compounds where Class or MOA did not join corre

Prioritised_antibiotics_g3 %>%
  filter(is.na(MOA) | MOA == "") %>%
  select(Name, Class, MOA)

# A tibble: 0 × 3
# ℹ 3 variables: Name <chr>, Class <chr>, MOA <chr>

### ---------- 14) Save prioritised antibiotic tables ----------

# Save prioritised antibiotic tables for the supplementary information.
#
# Prioritised_antibiotics_g1 is Table S3.1.
# Prioritised_antibiotics_g2 is Table S3.2.
# Prioritised_antibiotics_g3 is Table S3.3.

setwd(here("Files", "Tables"))


# Remove the original abbreviated class column 

Prioritised_antibiotics_g1 <- Prioritised_antibiotics_g1 %>%
  select(-Class_abbrv)

Prioritised_antibiotics_g2 <- Prioritised_antibiotics_g2 %>%
  select(-Class_abbrv)

Prioritised_antibiotics_g3 <- Prioritised_antibiotics_g3 %>%
  select(-Class_abbrv)


# Save Table S3.1.
write.xlsx(
  Prioritised_antibiotics_g1,
  "Table_S3_1_prioritised_ABs.xlsx"
)

# Save Table S3.2.
write.xlsx(
  Prioritised_antibiotics_g2,
  "Table_S3_2_prioritised_ABs.xlsx"
)

# Save Table S3.3.
write.xlsx(
  Prioritised_antibiotics_g3,
  "Table_S3_3_prioritised_ABs.xlsx"
)


### ---------- 15) Save final physicochemical-property table ----------

# Save the final data frame containing:
#
# - PubChem information
# - Chemicalize physicochemical properties
# - eNTRYway descriptors
# - manually transcribed charge information
# - manually checked functional group information
# - compound use/class information
# - eNTRY rule flags
# - Group 1, Group 2, Group 3 or No group assignment
#
# This CSV file is used later in the PCA and PERMANOVA scripts.

setwd(here("Files", "Files_for_coding"))

write.csv(
  df_joined_final_6,
  "Final_AB_Pol_Physchem_info.csv",
  row.names = FALSE
)
### ---------- 16) Create supplementary physicochemical-property tables ----------

# Create supplementary tables for all antibiotics and all pollutants.
#
# These tables contain the selected physicochemical properties needed for the
# supplementary information.
#

# for antibiotics & pollutants we want the following cols for SI

# Name
# Priority Group
# Moleular Weight (Da)	
# Topological Polar Surface Area  [A^2]	
# Van der Waals surface area [A^2]	
# Relative Polar Surface Area 	
# Relative Polar Surface Area (%)	
# Hydrogen Bond Donor Count	#Hydrogen Bond Acceptor Count	
# Rotatable Bond Count	
# Globularity Score
# PFB
# cLogP	
# cLogD7.4	
# pKa Acid Strongest	
# pKa Base Strongest	
# Predominant spp7.4	 
# % spp7.4	
# Net Charge pH 7.4	
# Ionisable Nitrogen?	
# Hydrophilic compound?	#Low globularity?	
# Low Flexability?	
                                                                        # these columns in df_joined_final_6 are 
# Compound_Name = Name 
# Group = Priority Group
# exact_mass_da = Molecular Weight (Da)	
# tpsa_a2 = Topological Polar Surface Area  [A^2]	
# vdw_surface_area_a2 = Van der Waals surface area [A^2]	
# rel_polar_sa = Relative Polar Surface Area 	
# rel_polar_sa * 100 = Relative Polar Surface Area (%)	
# hbd_count = Hydrogen Bond Donor Count	
# hba_count = Hydrogen Bond Acceptor Count	
# we can use rb_EW or rotatable_bond_count = Rotatable Bond Count	as they are the same
# glob_EW = Globularity Score
# pbf_EW
# logp = cLogP	
# logD_pH7_4_interp = cLogD7.4	
# pka_acid_strongest = pKa Acid Strongest	
# pka_base_strongest = pKa Base Strongest	
# spp_7_4 = Predominant spp7.4	 
# per_dominant_charge_pH_7_4 * 100 = % spp7.4	
# Charge_pH_7_4_num = Net Charge pH 7.4	
# func_group_other = Ionisable Nitrogen?	
# hydrophilic_compound = Hydrophilic compound?	
# lowglob = Low globularity?	
# lowflex = Low Flexability?	                                                 

# Create SI table with selected and renamed columns.

df_joined_final_6_SI <- df_joined_final_6 %>%
  mutate(
    `Relative Polar Surface Area (%)` = rel_polar_sa * 100,
    `% spp7.4` = per_dominant_charge_pH_7_4 * 100
  ) %>%
  select(
    Name                                      = Compound_Name,
    `Priority Group`                         = Group,
    `Molecular Weight (Da)`                  = exact_mass_da,
    `Topological Polar Surface Area [A^2]`   = tpsa_a2,
    `Van der Waals surface area [A^2]`        = vdw_surface_area_a2,
    `Relative Polar Surface Area`            = rel_polar_sa,
    `Relative Polar Surface Area (%)`,
    `Hydrogen Bond Donor Count`              = hbd_count,
    `Hydrogen Bond Acceptor Count`           = hba_count,
    `Rotatable Bond Count`                   = rotatable_bond_count,
    `Globularity Score`                      = glob_EW,
    PBF                                      = pbf_EW,
    cLogP                                    = logp,
    cLogD7.4                                 = logD_pH7_4_interp,
    `pKa Acid Strongest`                     = pka_acid_strongest,
    `pKa Base Strongest`                     = pka_base_strongest,
    `Predominant spp7.4`                     = spp_7_4,
    `% spp7.4`,
    `Net Charge pH 7.4`                      = Charge_pH_7_4_num,
    `Ionisable Nitrogen?`                    = func_group_other,
    `Hydrophilic compound?`                  = hydrophilic_compound,
    `Low globularity?`                       = lowglob,
    `Low flexibility?`                       = lowflex,
    CODE_2
  )


# Create Table S2.7: all antibiotics.

Table_S2_7_all_abs <- df_joined_final_6_SI %>%
  filter(CODE_2 == "Antibiotic") %>%
  select(-CODE_2)


# Create Table S2.8: all pollutants.

Table_S2_8_all_pols <- df_joined_final_6_SI %>%
  filter(CODE_2 == "Pollutant") %>%
  select(-CODE_2)

### ---------- 16.1) Add detailed antibiotic class and MOA to Table S2.7 ----------
# For Table S2.7, add more detailed antibiotic class and mechanism of action
# information from the manually created workbook Detailed_uses_and_MOA.xlsx.

setwd(here("Files", "Files_for_coding"))

Detailed_uses_and_MOA_abs <- read.xlsx(
  "Detailed_uses_and_MOA.xlsx",
  sheet = 1
)

# Join detailed antibiotic information by compound name.

Table_S2_7_all_abs <- Table_S2_7_all_abs %>%
  left_join(
    Detailed_uses_and_MOA_abs,
    by = "Name"
  )


# Move Antibiotic_Class and Mechanism_of_Action next to Name.

Table_S2_7_all_abs <- Table_S2_7_all_abs %>%
  relocate(
    Antibiotic_Class,
    Mechanism_of_Action,
    .after = Name
  )

# Check column order.
colnames(Table_S2_7_all_abs)


# sucess
### ---------- 16.2) Add detailed pollutant use information to Table S2.8 ----------

# For Table S2.8, add detailed pollutant use type and description from
# sheet 2 of Detailed_uses_and_MOA.xlsx.

Detailed_uses_and_MOA_pols <- read.xlsx(
  "Detailed_uses_and_MOA.xlsx",
  sheet = 2
)

# Join detailed pollutant information by compound name.

Table_S2_8_all_pols <- Table_S2_8_all_pols %>%
  left_join(
    Detailed_uses_and_MOA_pols,
    by = "Name"
  )


# Move Pollutant_use_type and Description next to Name.

Table_S2_8_all_pols <- Table_S2_8_all_pols %>%
  relocate(
    Pollutant_use_type,
    Description,
    .after = Name
  )

# Check column order.
colnames(Table_S2_8_all_pols)

# sucess

### ---------- 17) Save supplementary physicochemical-property tables ----------

# Save Table S2.7 and Table S2.8 as Excel files.

setwd(here("Files", "Tables"))

# Save Table S2.7: all antibiotics.
write_xlsx(
  Table_S2_7_all_abs,
  "Table_S2_7_all_abs.xlsx"
)

# Save Table S2.8: all pollutants.
write_xlsx(
  Table_S2_8_all_pols,
  "Table_S2_8_all_pols.xlsx"
)



### ---------- 18) Counts of pollutant & antibiotic types ----------
### ---------- 18) Counts of pollutant and antibiotic types ----------

# This section creates summary counts for:
#
# - pollutant use types in Table S2.8
# - descriptions assigned to prioritised pollutants
# - antibiotic classes in Table S2.7
#
# These counts are used to help describe the composition of the pollutant and
# antibiotic datasets in the manuscript/supplementary information.


# Count pollutant use types 

# Check column names in the all-pollutants supplementary table.
colnames(Table_S2_8_all_pols)

# Check the unique pollutant use types.
unique(Table_S2_8_all_pols$Pollutant_use_type)

# Count the number of pollutants in each use type.
# sort = TRUE orders the output from the highest count to the lowest count.

pollutant_use_type_counts <- Table_S2_8_all_pols %>%
  count(Pollutant_use_type, sort = TRUE)

pollutant_use_type_counts

# Pollutant_use_type                 n
# <chr>                          <int>
# 1 Non-antibiotic drug               31
# 2 Industrial chemical               27
# 3 Pesticide                         26
# 4 Non-antibiotic drug metabolite     7
# 5 Pesticide metabolite               5
# 6 Personal care product              2
# 7 Lifestyle product                  1
# 8 Lifestyle product metabolite       1


##Count prioritised pollutant description types

# Check column names in the prioritised pollutant table.
colnames(Prioritised_pollutants)

# Check the unique descriptions assigned to prioritised pollutants.
unique(Prioritised_pollutants$Description)

# View prioritised pollutant names and descriptions.
Prioritised_pollutants %>%
  select(Name, Description)

# There are 18 prioritised pollutants in total.

# Count the number of prioritised pollutants described as non-antibiotic drugs.
#
# This searches for "NAD;" so that NADs are counted separately from
# NAD metabolites.

Total_NADs <- Prioritised_pollutants %>%
  summarise(
    n_NADs = sum(str_detect(Description, "NAD;"))
  )

Total_NADs
# 15


# Count the number of prioritised pollutants described as NAD metabolites.

Total_NAD_metabolites <- Prioritised_pollutants %>%
  summarise(
    n_NAD_metabolite = sum(str_detect(Description, "NAD metabolite"))
  )

Total_NAD_metabolites
# 2


# Calculate total NADs plus NAD metabolites.

Total_NADs_and_metabolites <- Total_NADs$n_NADs +
  Total_NAD_metabolites$n_NAD_metabolite

Total_NADs_and_metabolites
# 17


# Count the number of prioritised pollutants described as psychoactive compounds.

Total_psychoactives <- Prioritised_pollutants %>%
  summarise(
    n_psychoactive = sum(str_detect(Description, "psychoactive compound"))
  )

Total_psychoactives
# 10


# Count antibiotic classes 

# Count the number of antibiotics in each antibiotic class.
# The classes come from the detailed antibiotic class information added to
# Table S2.7.

antibiotic_class_counts <- Table_S2_7_all_abs %>%
  count(Antibiotic_Class) %>%
  arrange(Antibiotic_Class)

antibiotic_class_counts

# A tibble: 14 × 2
# Antibiotic_Class                      n
# <chr>                             <int>
# 1 "Aminoglycoside"                     10
# 2 "Amphenicol"                          3
# 3 "Diaminopyrimidine"                   4
# 4 "Phosphonic acid"                     1
# 5 "Quinolone"                           1
# 6 "Quinolone; Fluoroquinolone"         24
# 7 "Sulphonamide"                       19
# 8 "Tetracycline"                       12
# 9 "β-lactam; Cephm; Carbacephem"        1
# 10 "β-lactam; Cephm; Cephalosporin "   25
# 11 "β-lactam; Monobactam"               2
# 12 "β-lactam; Penem"                    1
# 13 "β-lactam; Penem; Carbapenem"        9
# 14 "β-lactam; Penicillin"               7



# Calculate broader antibiotic class totals used in the text.

total_quinolone <- Table_S2_7_all_abs %>%
  filter(str_detect(Antibiotic_Class, "Quinolone|Fluoroquinolone")) %>%
  summarise(n = n())

total_quinolone
# 25


# Calculate total cephm antibiotics.
# This includes carbacephems and cephalosporins.

total_cephm <- Table_S2_7_all_abs %>%
  filter(str_detect(Antibiotic_Class, "Cephm")) %>%
  summarise(n = n())

total_cephm
# 26


# Calculate total penem antibiotics.
# This includes penems and carbapenems.

total_penem <- Table_S2_7_all_abs %>%
  filter(str_detect(Antibiotic_Class, "Penem")) %>%
  summarise(n = n())

total_penem
# 10


# Calculate total β-lactam antibiotics.
# This includes cephms, monobactams, penems and penicillins.

total_beta_lactam <- Table_S2_7_all_abs %>%
  filter(str_detect(Antibiotic_Class, "β-lactam")) %>%
  summarise(n = n())

total_beta_lactam
# 45