# ============================================================
# Script 6
# pt 1. PCA of 8 physicochemical properties excluding molecular weight
# ============================================================
# Purpose:

#
# This script uses the final antibiotic and pollutant physicochemical
# information table created during data wrangling.
#
# The aim of this script is to run a PCA using eight physicochemical properties,
# excluding molecular weight. 
# The eight physicochemical properties used in the PCA are:
#
# - cLogP
# - cLogD at pH 7.4
# - relative polar surface area
# - hydrogen bond donor count
# - hydrogen bond acceptor count
# - net charge at pH 7.4
# - rotatable bond count
# - globularity score
#
# The script first reads in the final physicochemical information table and
# selects the columns needed for PCA and plotting.
#
# The PCA is then run using scaled and centred variables so that all
# physicochemical properties contribute on comparable scales.
#
# The script creates:
#
# - a scree plot showing the proportion of variance explained by the first
# eight principal components
# - PCA biplots for PC1 vs PC2, PC1 vs PC3 and PC2 vs PC3
# - a heatmap of PCA species scores for PC1-PC8
# - a PCA summary table containing eigenvalues, proportion explained and
# cumulative proportion explained
#
# The script also uses the PCA site scores to create more detailed PCA plots for
# visualising compound separation in physicochemical-property space.
#
# For these plots, the PCA site scores are combined with compound metadata,
# including compound name, compound type, class and priority group.
#
# Labelled plots are created for PC1 vs PC2 and PC1 vs PC3. These plots
# are not used directly in the paper, but are used to help
# with manual labelling of the final figures.
#
# The script creates supplementary PCA plots where:
#
# - antibiotics are coloured and shaped by antibiotic class
# - pollutants are coloured and shaped by pollutant use category
#
# Finally, the script creates the main PCA figures where compounds are coloured
# by priority group. Antibiotics are shaped by antibiotic class, while pollutants
# are shown using a single pollutant shape. 
# These outputs are saved as TIFF figures and Excel tables for use in the
# manuscript and supplementary information.
#

# ============================================================

## ---------- 1.0) Setup ----------

# Clear the global environment.
# Use with care because this removes all existing objects.
rm(list = ls())

# Load required packages.
library(dplyr)      # data wrangling: select, mutate, filter, rename and pipes
library(tidyr)      # reshape data for the PCA heatmap
library(tibble)     # convert row names to columns
library(ggplot2)    # make PCA plots and heatmaps
library(ggrepel)    # add non-overlapping labels to PCA plots
library(ggtext)     # allow markdown-style axis text in ggplot themes
library(ggnewscale) # allow multiple colour, fill and shape scales in one ggplot
library(vegan)      # run PCA using rda(), extract scores and eigenvalues
library(writexl)    # save Excel files
library(readxl)     # read Excel files, if needed
library(stringr)    # tidy character strings and remove extra spaces
library(scales)     # plotting scale handling, including squish for heatmaps
library(here)         # make file paths portable across different computers
# Get package citations  

# Get the citations for the R packages used in this script

# These can be used in the methods section or supplementary information if needed.

citation("dplyr")
citation("tidyr")
citation("tibble")
citation("ggplot2")
citation("ggrepel")
citation("ggtext")
citation("ggnewscale")
citation("vegan")
citation("writexl")
citation("readxl")
citation("stringr")
citation("scales")
citation("here")


# Check R version and session information

# Check R version.
R.version.string
# R version 4.5.0 (2025-04-11 ucrt)

# Record session information so the package versions used are documented.
sessionInfo()


# Set working directory and read in final data 

# Set working directory to where the final physicochemical information file is saved.
setwd(here("Files", "Files_for_coding"))

# List files to check that the final CSV file is there.
list.files()

# Read in the final antibiotic and pollutant physicochemical information table.
df <- read.csv("Final_AB_Pol_Physchem_info.csv")



## ---------- 1.2) Tidy the data frame ----------

# Check the data frame.
head(df)
str(df)
colnames(df)

# Select the columns needed for PCA, plotting and later interpretation.
# Rename columns so that the names are shorter and easier to use in the script.

df2 <- df %>%
  dplyr::select(
    Name = Compound_Name,
    TYPE = CODE_2,
    Class = Class_abbrv,
    Class_2 = Class_abbrv_2,
    MW = exact_mass_da,
    logp,
    logd7_4 = logD_pH7_4_interp,
    Rel_PSA_vdw = rel_polar_sa,
    HBD = hbd_count,
    HBA = hba_count,
    Dominant_charge_pH_7_4_code,
    Numeric_charge = Charge_pH_7_4_num,
    hydrophilic_compound,
    ENTRY_compound,
    func_group_updated,
    rb_EW,
    glob_EW,
    Group,
    Positive_charge_T_F,
    Zwitt_charge_T_F,
    Neg_charge_T_F
  )

# Check selected column names.
colnames(df2)


## ---------- 1.3) Check expected pollutant group assignments ----------

# Double check that the expected number of pollutants are assigned to each
# priority group.

df2 %>%
  filter(TYPE == "Pollutant" & Group == "Group 1") %>%
  select(Name)

# 7 pollutants, as expected.

df2 %>%
  filter(TYPE == "Pollutant" & Group == "Group 2") %>%
  select(Name)

# 4 pollutants, as expected.

df2 %>%
  filter(TYPE == "Pollutant" & Group == "Group 3") %>%
  select(Name)

# 7 pollutants, as expected.

## ---------- 1.4) Create PCA data frame ----------

# Make a data frame containing only the physicochemical properties used in PCA.

props_cols <- df2 %>%
  dplyr::select(
    clogp              = logp,
    clogD7_4           = logd7_4,
    Rel_PSA_vdw,
    HBD,
    HBA,
    Numeric_charge_7_4 = Numeric_charge,
    RBs                = rb_EW,
    Glob               = glob_EW
  )

# Make sure all PCA variables are numeric/intergers before running PCA.
props_cols <- props_cols %>%
  mutate(
    clogp              = as.numeric(clogp),
    clogD7_4           = as.numeric(clogD7_4),
    Rel_PSA_vdw        = as.numeric(Rel_PSA_vdw),
    HBD                = as.integer(HBD),
    HBA                = as.integer(HBA),
    Numeric_charge_7_4 = as.numeric(Numeric_charge_7_4),
    RBs                = as.integer(RBs),
    Glob               = as.numeric(Glob)
  )

# Make a data frame containing compound names and grouping variables.
names_cols <- df2 %>%
  dplyr::select(
    Name,
    TYPE,
    Class,
    Group
  )

# Combine names/grouping columns with the PCA property columns.
properties <- cbind(names_cols, props_cols)



## ---------- 1.5) Run PCA ----------

# Remove TYPE because it is not needed for the PCA object.
pca_df_props <- subset(properties, select = -c(TYPE))

# Make row names the compound names.
pca_props <- pca_df_props[, -1]
rownames(pca_props) <- pca_df_props[, 1]

# Remove non-numeric columns before PCA.
pca_props_num <- pca_props %>%
  dplyr::select(
    -Class,
    -Group
  )

# Run PCA using vegan::rda().
# The variables are scaled and centred so that all properties contribute on
# comparable scales.
PCA <- rda(
  pca_props_num,
  scale = TRUE,
  center = TRUE
)

## ---------- 1.6) Make PCA scree plot ----------

# Extract eigenvalues and calculate the proportion of variance explained by
# each principal component.

eig_prop <- PCA$CA$eig / sum(PCA$CA$eig)
eig_all <- PCA$CA$eig

# Make Figure S5.5:
# Scree plot showing the proportion of variance explained by the first eight
# principal components and the cumulative variance explained.

make_scree_plot <- function(eig_all, N = 8) {
  
  eig <- eig_all[1:N]
  
  eig_prop <- eig / sum(eig_all)
  cum_prop <- cumsum(eig_prop)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(
    mar = c(6, 4.5, 1.5, 1) + 0.1,
    xpd = NA,
    cex.axis = 0.95,
    cex.lab = 1.05,
    lwd = 1.2
  )
  
  y_top <- max(c(eig_prop, cum_prop)) * 1.20
  
  bp <- barplot(
    eig_prop,
    names.arg = paste0("PC", 1:N),
    las = 2,
    ylab = "Proportion of variance explained",
    ylim = c(0, y_top),
    border = "black"
  )
  
  lines(bp, cum_prop, type = "b", pch = 16, lwd = 1.2)
  
  text(
    bp,
    cum_prop,
    labels = paste0(round(100 * cum_prop, 1), "%"),
    pos = 3,
    cex = 0.75
  )
  
  legend(
    "topright",
    legend = c("Cumulative"),
    lty = c(NA, 1),
    pch = c(NA, 16),
    bty = "n",
    cex = 0.8
  )
  
  # Add Kaiser threshold line.
  kaiser <- 1 / sum(eig_all)
  
  abline(
    h = kaiser,
    col = "red",
    lty = 2,
    lwd = 1.2
  )
  
  text(
    x = bp[1],
    y = kaiser,
    labels = "Kaiser (eig = 1)",
    pos = 3,
    cex = 0.8
  )
}

# Set working directory for PCA/PERMANOVA without molecular weight 
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Figs"))
# Save Figure S5.5 as a TIFF.

tiff(
  "Figure_S5_5_scree_plot_no_MW.tiff",
  width = 7,
  height = 5,
  units = "in",
  res = 600,
  compression = "lzw"
)

make_scree_plot(eig_all, N = 8)

dev.off()

## ---------- 1.7) Check variance explained ----------

# Calculate the percentage of variance explained by the first two axes.
sum((as.vector(PCA$CA$eig) / sum(PCA$CA$eig))[1:2])
# [1] 0.6900039
# ~69%

# Calculate the percentage of variance explained by the first three axes.
sum((as.vector(PCA$CA$eig) / sum(PCA$CA$eig))[1:3])
# [1] 0.8077188
# ~81%


## ---------- 1.8) Extract PCA scores and save PCA summary table Table S5.1. ----------

# Extract site scores for further analyses.
sitePCA <- PCA$CA$u
sitePCA

# Extract species scores for further analyses.
speciesPCA <- PCA$CA$v
speciesPCA

# Convert site and species scores to data frames.
sitePCA_df <- as.data.frame(PCA$CA$u)

speciesPCA_df <- as.data.frame(PCA$CA$v) %>%
  tibble::rownames_to_column(var = "Parameter")

# Extract species scores for PC1-PC8.
load <- scores(PCA, display = "species", choices = 1:8)

# Get a summary of the PCA.
summary(PCA)

# Make Table S5.1:
# PCA eigenvalues and variance explained for PC1-PC8.

eig <- vegan::eigenvals(PCA)
prop <- eig / sum(eig)
cumprop <- cumsum(prop)

pca_summary_df <- data.frame(
  PC = paste0("PC", seq_along(eig)),
  Eigenvalue = as.numeric(eig),
  Proportion_Explained = as.numeric(prop),
  Cumulative_Proportion = as.numeric(cumprop)
)

pca_summary_df

# Save Table S5.1.

setwd(here("Files", "Tables"))


write_xlsx(
  pca_summary_df,
  "Table_S5_1_pca_summary_no_mw.xlsx"
)


## ---------- 1.9) Make PCA biplots ----------

# Make biplots for:
#
# - Figure S5.2: PC1 vs PC2
# - Figure S5.3: PC1 vs PC3
# - Figure S5.4: PC2 vs PC3

# Rename variable labels for plotting.
new_names <- c(
  clogD7_4 = "cLogD7.4",
  clogp = "cLogP",
  RBs = "#RB",
  HBD = "#HBD",
  HBA = "#HBA",
  Glob = "Globularity score",
  Rel_PSA_vdw = "Rel. PSA",
  Numeric_charge_7_4 = "Net Charge pH7.4"
)

# Use pca_summary_df to check the proportion of variance explained by PC1-PC3.
pca_summary_df %>%
  filter(PC %in% c("PC1", "PC2", "PC3")) %>%
  select(
    PC,
    Proportion_Explained
  )

# Make proportion explained labels for the biplot axes.
pc_var <- c(
  `1` = "55.2%",
  `2` = "13.8%",
  `3` = "11.8%"
)

# Function to create PCA biplots.

plot_pca_biplot <- function(PCA, axes = c(1, 2), scaling = 0,
                            label_cex = 0.6,
                            label_col = "#0072B2") {
  
  plot(
    PCA,
    type = "n",
    scaling = scaling,
    choices = axes,
    xlab = "",
    ylab = ""
  )
  
  points(
    PCA,
    display = "sites",
    scaling = scaling,
    choices = axes,
    pch = 16,
    cex = 0.4
  )
  
  sp <- scores(
    PCA,
    display = "species",
    scaling = scaling,
    choices = axes
  )
  
  mul <- ordiArrowMul(sp)
  
  x <- sp[, 1] * mul
  y <- sp[, 2] * mul
  
  labs <- rownames(sp)
  labs[labs %in% names(new_names)] <- new_names[labs[labs %in% names(new_names)]]
  
  arrows(0, 0, x, y, length = 0.08)
  
  pos <- ifelse(
    abs(x) > abs(y),
    ifelse(x > 0, 4, 2),
    ifelse(y > 0, 3, 1)
  )
  
  text(
    x,
    y,
    labels = labs,
    pos = pos,
    offset = 0.4,
    cex = label_cex,
    col = label_col
  )
  
  x_pc <- paste0("PC", axes[1])
  y_pc <- paste0("PC", axes[2])
  
  x_var <- pc_var[as.character(axes[1])]
  y_var <- pc_var[as.character(axes[2])]
  
  mtext(x_pc, side = 1, line = 2.2, cex = 1.1)
  mtext(y_pc, side = 2, line = 2.2, cex = 1.1)
  
  mtext(x_var, side = 1, line = 3.2, cex = 0.9)
  mtext(y_var, side = 2, line = 3.2, cex = 0.9)
}

# View biplots.
plot_pca_biplot(PCA, axes = c(1, 2), label_cex = 0.6)
plot_pca_biplot(PCA, axes = c(1, 3), label_cex = 0.6)
plot_pca_biplot(PCA, axes = c(2, 3), label_cex = 0.6)


## ---------- 1.10) Save PCA biplots ----------

setwd(here("Files", "PCA_PERMANOVA_NO_MW_Figs"))

# Save Figure S5.2: PC1 vs PC2.
tiff(
  "Figure_S5_2_PCA_PC1_PC2_no_MW.tiff",
  width = 7,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

plot_pca_biplot(PCA, axes = c(1, 2), label_cex = 0.6)

dev.off()

# Save Figure S5.3: PC1 vs PC3.
tiff(
  "Figure_S5_3_PCA_PC1_PC3_no_MW.tiff",
  width = 7,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

plot_pca_biplot(PCA, axes = c(1, 3), label_cex = 0.6)

dev.off()

# Save Figure S5.4: PC2 vs PC3.
tiff(
  "Figure_S5_4_PCA_PC2_PC3_no_MW.tiff",
  width = 7,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

plot_pca_biplot(PCA, axes = c(2, 3), label_cex = 0.6)

dev.off()


## ---------- 1.11) Make PCA species score heatmap ----------

# Make Figure S5.1:
# Heatmap of species scores for PC1-PC8.

# Extract species scores for PC1-PC8.
# from the fitted PCA object.
sp <- scores(
  PCA,
  display = "species",
  choices = 1:8,
  scaling = 0
)

# Order rows by clustering variables with similar PCA score patterns together.
row_order <- hclust(dist(sp))$order

sp_df <- as.data.frame(sp) %>%
  rownames_to_column("Parameter") %>%
  mutate(
    Parameter = factor(Parameter, levels = Parameter[row_order])
  )

# Convert to long format for ggplot.
sp_long <- sp_df %>%
  pivot_longer(
    -Parameter,
    names_to = "PC",
    values_to = "Loading"
  ) %>%
  mutate(
    label = sprintf("%.2f", Loading)
  )

# Map old variable names to cleaner plot labels.
param_map <- c(
  clogD7_4 = "cLogD7.4",
  clogp = "cLogP",
  Rel_PSA_vdw = "Relative PSA",
  HBD = "# H-Bond Donors",
  HBA = "# H-Bond Acceptors",
  RBs = "# Rotatable Bonds",
  Glob = "Globularity Score",
  Numeric_charge_7_4 = "Net Charge pH 7.4"
)

# Desired order from top to bottom.
param_order <- c(
  "cLogP",
  "cLogD7.4",
  "Relative PSA",
  "# H-Bond Donors",
  "# H-Bond Acceptors",
  "# Rotatable Bonds",
  "Globularity Score",
  "Net Charge pH 7.4"
)

sp_long <- sp_long %>%
  mutate(
    Parameter = recode(Parameter, !!!param_map),
    
    # ggplot puts the first factor level at the bottom of the y-axis,
    # so reverse the order to get the desired top-to-bottom order.
    Parameter = factor(Parameter, levels = rev(param_order))
  )

# Make heatmap using colour-blind friendly colours.
PCA_species_scores_heatmap_PC18 <- 
  ggplot(sp_long, aes(x = PC, y = Parameter, fill = Loading)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    midpoint = 0,
    low = "#15E635",
    mid = "white",
    high = "#FF096E",
    limits = c(-0.35, 0.36),
    oob = scales::squish
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    panel.grid = element_blank()
  )

# View heatmap.
PCA_species_scores_heatmap_PC18


## ---------- 1.12) Save PCA species score heatmap ----------

# Save Figure S5.1 as a TIFF.
# 300-600 dpi is typical for publication/Word outputs.
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Figs"))

ggsave(
  filename = "Figure_S5_1_PCA_species_scores_heatmap_no_MW.tiff",
  plot = PCA_species_scores_heatmap_PC18,
  device = "tiff",
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


## ---------- 1.13) Make Detailed PCA Plots with ggplot ----------

# Create a plotting data frame by combining compound information with PCA site scores.
# names_cols contains compound name, type, class and group information.
# sitePCA contains the PCA site scores for each compound.

pca_plot_df <- cbind(names_cols, sitePCA)

# Calculate the proportion of variance explained by each principal component.
pve <- PCA$CA$eig / sum(PCA$CA$eig)

# Create axis labels showing the percentage of variance explained.
lab1 <- paste0("PC1: ", round(100 * pve[1], 1), "%")
lab2 <- paste0("PC2: ", round(100 * pve[2], 1), "%")
lab3 <- paste0("PC3: ", round(100 * pve[3], 1), "%")

# Clean class names and create plotting variables.
# TYPE2 is used for clearer labelling of antibiotics and pollutants if needed.
# Group_plot keeps Group 1, Group 2 and Group 3 as priority groups, while all
# other compounds are labelled as "No Group".

pca_plot_df <- pca_plot_df %>%
  mutate(
    Class = str_squish(Class),
    TYPE2 = case_when(
      TYPE == "Antibiotic" ~ "GN Antibiotic",
      TYPE == "Pollutant"  ~ "Pollutant",
      TRUE                 ~ as.character(TYPE)
    ),
    Group_plot = ifelse(
      Group %in% c("Group 1", "Group 2", "Group 3"),
      Group,
      "No Group"
    )
  )


## ---------- 1.14) Make labelled PCA plots for manual checking ----------

# These plots are not used directly in the paper.
# They are used as checking plots and to help identify which pollutant labels
# should be added manually to the final figures.

# Define which pollutant groups should be labelled.
pol_label_groups <- c("Group 1", "Group 2", "Group 3")

# Define the column used for labels.
label_col <- "Name"


## ---------- 1.16) Plotting theme ----------

# Create a consistent plotting theme for PCA figures.

theme_journal <- theme_bw(base_size = 10, base_family = "Arial") + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), 
    aspect.ratio = 1, 
    axis.title.x = ggtext::element_markdown(size = 12), 
    axis.title.y = ggtext::element_markdown(size = 12), 
    axis.text.x = element_text(size = 9.5, colour = "black"), 
    axis.text.y = element_text(size = 9.5, colour = "black"), 
    legend.title = element_text(size = 9.5), 
    legend.text = element_text(size = 9), 
    legend.position = "right", 
    legend.direction = "vertical", 
    strip.background = element_blank(), 
    strip.text = element_text(size = 9), 
    plot.title = element_text(size = 11, face = "plain"), 
    plot.margin = margin(8, 8, 8, 8)
  )


## ---------- 1.15) Define priority group colours and transparency ----------

# Define colours for priority groups.
# No Group compounds are shown in grey and with lower transparency 

group_cols <- c(
  "Group 1" = "#EF1665",
  "Group 2" = "#FFC107",
  "Group 3" = "#1EB3E5",
  "No Group" = "grey70"
)

group_alpha <- c(
  "Group 1" = 1,
  "Group 2" = 1,
  "Group 3" = 1,
  "No Group" = 0.3
)


## ---------- 1.16) Plot priority groups with pollutant labels ----------

# Make a basic PCA plot for PC1 vs PC2.
# Points are coloured by priority group, shaped by compound type and labelled
# for prioritised pollutants.

Basic_group_PC1_2_labs <- ggplot(pca_plot_df, aes(PC1, PC2)) +
  geom_point(
    aes(
      color = Group_plot,
      shape = TYPE,
      alpha = Group_plot
    ),
    size = 3.5
  ) +
  scale_colour_manual(values = group_cols) +
  scale_alpha_manual(values = group_alpha) +
  ggnewscale::new_scale_colour() + 
  geom_text_repel(
    data = pca_plot_df %>%
      filter(TYPE == "Pollutant", Group %in% pol_label_groups),
    aes(label = .data[[label_col]]),
    size = 3,
    show.legend = FALSE,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.alpha = 0.8,
    arrow = arrow(
      length = unit(0.006, "npc"),
      type = "closed",
      ends = "last"
    )
  ) +
  labs(
    colour = "Group",
    shape = "TYPE",
    x = lab1,
    y = lab2
  ) +
  guides(alpha = "none", size = "none") +
  theme_journal


# Make a basic PCA plot for PC1 vs PC3.
Basic_group_PC1_3_labs <- ggplot(pca_plot_df, aes(PC1, PC3)) +
  geom_point(
    aes(
      color = Group_plot,
      shape = TYPE,
      alpha = Group_plot
    ),
    size = 3.5
  ) +
  scale_colour_manual(values = group_cols) +
  scale_alpha_manual(values = group_alpha) +
  ggnewscale::new_scale_colour() + 
  geom_text_repel(
    data = pca_plot_df %>%
      filter(TYPE == "Pollutant", Group %in% pol_label_groups),
    aes(label = .data[[label_col]]),
    size = 3,
    show.legend = FALSE,
    min.segment.length = 0,
    segment.size = 0.25,
    segment.alpha = 0.8,
    arrow = arrow(
      length = unit(0.006, "npc"),
      type = "closed",
      ends = "last"
    )
  ) +
  labs(
    colour = "Group",
    shape = "TYPE",
    x = lab1,
    y = lab3
  ) +
  guides(alpha = "none", size = "none") +
  theme_journal

# View labled plots.
Basic_group_PC1_2_labs
Basic_group_PC1_3_labs


## ---------- 1.17) Save labelled PCA plots ----------
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Figs"))

# Save PC1 vs PC2 labelled plot.
ggsave(
  filename = "Basic_group_PC1_2_labs_no_MW.tiff",
  plot = Basic_group_PC1_2_labs,
  device = "tiff",
  width = 7,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

# Save PC1 vs PC3 labelled plot.
ggsave(
  filename = "Basic_group_PC1_3_labs_no_MW.tiff",
  plot = Basic_group_PC1_3_labs,
  device = "tiff",
  width = 7,
  height = 7,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)
## ---------- 1.18) Make PCA plots coloured and shaped by compound class ----------

# Make plots for:
#
# Figure S5.7:
# PCA plot of physicochemical properties, PC1 vs PC2.
# Points represent compounds. Antibiotics are coloured and shaped by antibiotic
# class. Pollutants are coloured and shaped by use category.
#
# Figure S5.8:
# PCA plot of physicochemical properties, PC1 vs PC3.
# Points represent compounds. Antibiotics are coloured and shaped by antibiotic
# class. Pollutants are coloured and shaped by use category.

# Define antibiotic and pollutant plotting colours, 

# Define antibiotic colours by class.
ab_cols <- c(
  "Aminoglycoside" = "#a00e00",
  "Amphenicol" = "#FF4630",
  "Cephm" = "#A629CC",
  "Monobactam" = "#CC294F",
  "Penicillin" = "#0093A0",
  "Penem" = "#5429CC",
  "Diaminopyrimidine" = "#0F07F0",
  "Fluoroquinolone" = "#B800D0",
  "Quinolone" = "#3089FF",
  "Phosphonic acid" = "#881B6B",
  "Sulphonamide" = "#132b69",
  "Tetracycline" = "#135669"
)

# Check what shapes are available to us

# Create data frame of pch values
shape_df <- data.frame(
  shape = 0:25,
  x = rep(1:6, length.out = 26),
  y = rep(5:1, each = 6, length.out = 26)
)

# Plot all shapes
ggplot(shape_df, aes(x = x, y = y)) +
  geom_point(
    aes(shape = factor(shape)),
    size = 3,
    stroke = 1
  ) +
  geom_text(
    aes(label = shape),
    vjust = 2.2,
    size = 3.5
  ) +
  scale_shape_manual(values = shape_df$shape) +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "none"
  ) +
  labs(
    title = "ggplot2 point shapes 0-25"
  )


# Define antibiotic shapes by class.
ab_shapes <- c(
  "Aminoglycoside"      = 0,   # square, open
  "Amphenicol"          = 1,   # circle, open
  "Cephm"               = 7,   # square cross
  "Monobactam"          = 3,   # plus
  "Penicillin"          = 4,   # cross / x
  "Penem"               = 5,   # diamond, open
  "Diaminopyrimidine"   = 2,   # triangle, open
  "Fluoroquinolone"     = 6,   # inverted triangle, open
  "Quinolone"           = 11,  # star / superimposed triangles
  "Phosphonic acid"     = 13,  # circle cross
  "Sulphonamide"        = 9,   # diamond plus
  "Tetracycline"        = 8    # star / asterisk
)

# View the shapes for antibiotics 

plot.new()

legend(
  "center",
  legend = names(ab_shapes),
  pch = ab_shapes,
  pt.cex = 1.5,
  bty = "n",
  title = "Antibiotic class"
)


# Define pollutant colours by use category.
pol_cols <- c(
  "Industrial chemical" = "#f6c200",
  "Lifestyle product" = "#B0F600",
  "Lifestyle product metabolite" = "#3D49B8",
  "Non-antibiotic drug" = "#3DB86F",
  "Non-antibiotic drug metabolite" = "#3083FF",
  "Personal care product" = "#6E3DB8",
  "Pesticide" = "#d04e00",
  "Pesticide metabolite" = "#692613"
)

# Define pollutant shapes by use category.
pol_shapes <- c(
  "Industrial chemical" = 18,              # filled diamond
  "Lifestyle product" = 10,                # circle plus
  "Lifestyle product metabolite" = 12,     # square plus
  "Non-antibiotic drug" = 25,              # fillable inverted triangle
  "Non-antibiotic drug metabolite" = 17,   # filled triangle
  "Personal care product" = 20,            # small filled circle
  "Pesticide" = 16,                        # filled circle
  "Pesticide metabolite" = 15              # filled square
)

# View the shapes for pollutants 

plot.new()

legend(
  "center",
  legend = names(pol_shapes),
  pch = pol_shapes,
  pt.cex = 1.5,
  bty = "n",
  title = "Antibiotic class"
)

# Prepare antibiotic and pollutant plotting data 

# Order colour and shape vectors alphabetically so that legends are consistent.
ab_cols <- ab_cols[order(names(ab_cols))]
ab_shapes <- ab_shapes[order(names(ab_shapes))]
ab_breaks <- names(ab_cols)

pol_cols <- pol_cols[order(names(pol_cols))]
pol_shapes <- pol_shapes[order(names(pol_shapes))]
pol_breaks <- names(pol_cols)

# Clean class names and create separate class columns for antibiotics and pollutants.
# This helps ggplot apply separate colour and shape scales to antibiotics and pollutants.
pca_plot_df <- pca_plot_df %>%
  mutate(
    Class = str_squish(str_replace_all(Class, "\u00A0", " ")),
    Class_AB = if_else(
      TYPE == "Antibiotic" & Class %in% ab_breaks,
      Class,
      NA_character_
    ),
    Class_POL = if_else(
      TYPE == "Pollutant" & Class %in% pol_breaks,
      Class,
      NA_character_
    )
  )

# Make point sizes uniform.
plot_point_size <- 2
legend_point_size <- 2

# Define fillable ggplot shapes.
# Shapes 21-25 use both colour and fill, so they need to be plotted separately
# from non-fillable shapes.
fillable_shapes <- c(21, 22, 23, 24, 25)

# Split antibiotics into fillable and non-fillable shapes.
ab_fillable_classes <- names(ab_shapes)[ab_shapes %in% fillable_shapes]
ab_nonfillable_classes <- names(ab_shapes)[!ab_shapes %in% fillable_shapes]

ab_fill_dat <- pca_plot_df %>%
  filter(!is.na(Class_AB), Class_AB %in% ab_fillable_classes)

ab_nonfill_dat <- pca_plot_df %>%
  filter(!is.na(Class_AB), Class_AB %in% ab_nonfillable_classes)

# Split pollutants into fillable and non-fillable shapes.
pol_fillable_classes <- names(pol_shapes)[pol_shapes %in% fillable_shapes]
pol_nonfillable_classes <- names(pol_shapes)[!pol_shapes %in% fillable_shapes]

pol_fill_dat <- pca_plot_df %>%
  filter(!is.na(Class_POL), Class_POL %in% pol_fillable_classes)

pol_nonfill_dat <- pca_plot_df %>%
  filter(!is.na(Class_POL), Class_POL %in% pol_nonfillable_classes)

# legend appearance for mixed fillable/non-fillable antibiotic shapes
ab_legend_colour <- ifelse(
  ab_shapes[ab_breaks] %in% fillable_shapes,
  "black",
  unname(ab_cols[ab_breaks])
)

ab_legend_fill <- ifelse(
  ab_shapes[ab_breaks] %in% fillable_shapes,
  unname(ab_cols[ab_breaks]),
  "white"
)


# legend appearance for mixed fillable/non-fillable pollutant shapes
pol_legend_colour <- ifelse(
  pol_shapes[pol_breaks] %in% fillable_shapes,
  "black",
  unname(pol_cols[pol_breaks])
)

pol_legend_fill <- ifelse(
  pol_shapes[pol_breaks] %in% fillable_shapes,
  unname(pol_cols[pol_breaks]),
  "white"
)

#Figure S5.7: PC1 vs PC2, coloured and shaped by class 

Complex_shapes_PC1_2 <- ggplot(pca_plot_df, aes(PC1, PC2)) +
  
  # Add antibiotics.
  # These are non-fillable shapes, so they use colour rather than fill.
  geom_point(
    data = ab_nonfill_dat,
    aes(
      colour = Class_AB,
      shape = Class_AB
    ),
    size = plot_point_size,
    alpha = 0.85,
    stroke = 0.8,
    show.legend = TRUE
  ) +
  
  # Set antibiotic colours.
  scale_colour_manual(
    name = "Antibiotic",
    values = ab_cols,
    breaks = ab_breaks
  ) +
  
  # Set antibiotic shapes.
  scale_shape_manual(
    name = "Antibiotic",
    values = ab_shapes,
    breaks = ab_breaks,
    drop = FALSE
  ) +
  
  # Make sure antibiotic colour and shape are shown together.
  guides(
    colour = guide_legend(order = 1),
    shape  = guide_legend(order = 1)
  ) +
  
  # Start new colour, fill and shape scales for pollutants.
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale("shape") +
  # Add pollutants with non-fillable shapes.
  geom_point(
    data = pol_nonfill_dat,
    aes(
      colour = Class_POL,
      shape = Class_POL
    ),
    size = plot_point_size,
    alpha = 0.85,
    stroke = 0.5,
    show.legend = c(
      shape = TRUE,
      colour = FALSE
    )
  ) +
  
  # Add pollutants with fillable shapes.
  geom_point(
    data = pol_fill_dat,
    aes(
      fill = Class_POL,
      shape = Class_POL
    ),
    size = plot_point_size,
    alpha = 0.85,
    colour = "black",
    stroke = 0.3,
    show.legend = c(
      shape = TRUE,
      fill = FALSE
    )
  ) +
  
  scale_colour_manual(
    name = "Pollutants",
    values = pol_cols,
    breaks = pol_breaks,
    guide = "none"
  ) +
  
  scale_fill_manual(
    name = "Pollutants",
    values = pol_cols,
    breaks = pol_breaks,
    guide = "none"
  ) +
  
  scale_shape_manual(
    name = "Pollutants",
    values = pol_shapes,
    breaks = pol_breaks,
    drop = FALSE,
    guide = guide_legend(
      order = 2,
      override.aes = list(
        colour = pol_legend_colour,
        fill = pol_legend_fill,
        alpha = 1,
        size = legend_point_size,
        stroke = 0.4
      )
    )
  ) +
  
  labs(
    x = lab1,
    y = lab2
  ) +
  
  theme_journal +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.box.spacing = unit(0.2, "cm")
  )
# View Figure S5.7.
Complex_shapes_PC1_2

# Make Figure S5.8:
Complex_shapes_PC1_3 <- ggplot(pca_plot_df, aes(PC1, PC3)) +
  
  # Add antibiotics.
  # These are non-fillable shapes, so they use colour rather than fill.
  geom_point(
    data = ab_nonfill_dat,
    aes(
      colour = Class_AB,
      shape = Class_AB
    ),
    size = plot_point_size,
    alpha = 0.85,
    stroke = 0.8,
    show.legend = TRUE
  ) +
  
  # Set antibiotic colours.
  scale_colour_manual(
    name = "Antibiotic",
    values = ab_cols,
    breaks = ab_breaks
  ) +
  
  # Set antibiotic shapes.
  scale_shape_manual(
    name = "Antibiotic",
    values = ab_shapes,
    breaks = ab_breaks,
    drop = FALSE
  ) +
  
  # Make sure antibiotic colour and shape are shown together.
  guides(
    colour = guide_legend(order = 1),
    shape  = guide_legend(order = 1)
  ) +
  
  # Start new colour, fill and shape scales for pollutants.
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  ggnewscale::new_scale("shape") +
  # Add pollutants with non-fillable shapes.
  geom_point(
    data = pol_nonfill_dat,
    aes(
      colour = Class_POL,
      shape = Class_POL
    ),
    size = plot_point_size,
    alpha = 0.85,
    stroke = 0.5,
    show.legend = c(
      shape = TRUE,
      colour = FALSE
    )
  ) +
  
  # Add pollutants with fillable shapes.
  geom_point(
    data = pol_fill_dat,
    aes(
      fill = Class_POL,
      shape = Class_POL
    ),
    size = plot_point_size,
    alpha = 0.85,
    colour = "black",
    stroke = 0.3,
    show.legend = c(
      shape = TRUE,
      fill = FALSE
    )
  ) +
  
  scale_colour_manual(
    name = "Pollutants",
    values = pol_cols,
    breaks = pol_breaks,
    guide = "none"
  ) +
  
  scale_fill_manual(
    name = "Pollutants",
    values = pol_cols,
    breaks = pol_breaks,
    guide = "none"
  ) +
  
  scale_shape_manual(
    name = "Pollutants",
    values = pol_shapes,
    breaks = pol_breaks,
    drop = FALSE,
    guide = guide_legend(
      order = 2,
      override.aes = list(
        colour = pol_legend_colour,
        fill = pol_legend_fill,
        alpha = 1,
        size = legend_point_size,
        stroke = 0.4
      )
    )
  ) +
  
  labs(
    x = lab1,
    y = lab3
  ) +
  
  theme_journal +
  
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.box.spacing = unit(0.2, "cm")
  )
# View Figure S5.8.
Complex_shapes_PC1_3


## ----------  1.19) Save Figure S5.7 & S5.8 ----------
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Figs"))

# Save Figure S5.7 as a TIFF.
ggsave(
  filename = "Figure_S5_7_PC1_PC2_no_mw.tiff",
  plot = Complex_shapes_PC1_2,
  device = "tiff",
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# Save Figure S5.8 as a TIFF.
ggsave(
  filename = "Figure_S5_8_PC1_PC3_no_mw.tiff",
  plot = Complex_shapes_PC1_3,
  device = "tiff",
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


## ---------- 1.20) Make final PCA plots coloured by priority group ----------

# Make plots for:
#
# Figure 1:
# PCA plot, PC1 vs PC2, using eight physicochemical properties:
# hydrogen bond donor count, hydrogen bond acceptor count, cLogD7.4, cLogP,
# relative PSA, rotatable bond count, globularity score and net charge at pH 7.4.
# Points represent compounds, coloured by prioritisation group. Antibiotics are
# shaped by antibiotic class and pollutants are shown as triangles.
#
# Figure S5.6:
# PCA plot, PC1 vs PC3, using the same aesthetics.
#
# Labels were added manually later using Affinity Designer.

# Define antibiotic shapes for final group-coloured figures.
ab_shapes <- c(
  "Aminoglycoside"      = 11,  # superimposed triangles
  "Amphenicol"          = 10,  # circle plus
  "Cephm"               = 23,  # fillable diamond
  "Monobactam"          = 3,   # plus
  "Penicillin"          = 25,  # fillable inverted triangle
  "Penem"               = 22,  # fillable square
  "Diaminopyrimidine"   = 4,   # cross / x
  "Fluoroquinolone"     = 21,  # fillable circle
  "Quinolone"           = 13,  # circle cross
  "Phosphonic acid"     = 14,  # square and triangle
  "Sulphonamide"        = 9,   # diamond plus
  "Tetracycline"        = 8    # star / asterisk
)

# Order antibiotic and pollutant plotting vectors.
ab_cols <- ab_cols[order(names(ab_cols))]
ab_shapes <- ab_shapes[order(names(ab_shapes))]
ab_breaks <- names(ab_cols)

pol_cols <- pol_cols[order(names(pol_cols))]
pol_breaks <- names(pol_cols)

# Define colour-blind friendly group colours.
# Reference used when selecting colours:
# https://davidmathlogic.com/colorblind/
group_cols <- c(
  "Group 1" = "#EF1665",
  "Group 2" = "#FFC107",
  "Group 3" = "#1EB3E5",
  "No Group" = "grey70"
)

group_alpha <- c(
  "Group 1" = 1,
  "Group 2" = 1,
  "Group 3" = 1,
  "No Group" = 0.3
)

## ---------- 1.21) Prepare data for final group-coloured PCA plots ----------

# Clean class and group names and create plotting columns.
# Antibiotics are shaped by antibiotic class.
# Pollutants are all given one pollutant shape.
# Points are coloured by priority group.

pca_plot_df <- pca_plot_df %>%
  mutate(
    Class = str_squish(str_replace_all(Class, "\u00A0", " ")),
    Group = str_squish(str_replace_all(Group, "\u00A0", " ")),
    
    Group_plot = case_when(
      Group %in% c("Group 1", "Group 2", "Group 3") ~ Group,
      TRUE ~ "No Group"
    ),
    Group_plot = factor(Group_plot, levels = names(group_cols)),
    
    Class_AB = if_else(
      TYPE == "Antibiotic" & Class %in% ab_breaks,
      Class,
      NA_character_
    ),
    Shape_POL = if_else(
      TYPE == "Pollutant",
      "Pollutant",
      NA_character_
    )
  )

# Split antibiotics into fillable and non-fillable shapes.
fillable_shapes <- c(21, 22, 23, 24, 25)

ab_fillable_classes <- names(ab_shapes)[ab_shapes %in% fillable_shapes]
ab_nonfillable_classes <- names(ab_shapes)[!ab_shapes %in% fillable_shapes]

ab_fill_dat <- pca_plot_df %>%
  filter(!is.na(Class_AB), Class_AB %in% ab_fillable_classes)

ab_nonfill_dat <- pca_plot_df %>%
  filter(!is.na(Class_AB), Class_AB %in% ab_nonfillable_classes) %>%
  mutate(
    plot_size_nonfill = if_else(Group_plot == "No Group", 1.5, plot_point_size)
  )

# Create pollutant plotting data.
pol_dat <- pca_plot_df %>%
  filter(!is.na(Shape_POL))

# Use a darker grey for No Group non-fillable shapes so they can still be seen.
group_cols_nonfill <- group_cols
group_cols_nonfill["No Group"] <- "grey30"

# Define point and legend sizes.
plot_point_size <- 3
legend_point_size <- 2
group_legend_size <- 1.2

## ---------- 1.22) Plot  PC1 vs PC2 / PC1 vs PC 3 coloured by priority group ----------

## ploy Figure 1: PC1 vs PC2 coloured by priority group 

# PCA plot of eight physicochemical properties, PC1 vs PC2.
# Points represent compounds and are coloured by prioritisation group.
# Antibiotics are shaped by antibiotic class, while pollutants are shown as
# triangles.

PC1_2_colour_group <- ggplot(pca_plot_df, aes(PC1, PC2)) +
  
  # Add antibiotics with non-fillable shapes.
  # These shapes use colour for the point symbol.
  # No Group points are smaller using plot_size_nonfill.
  geom_point(
    data = ab_nonfill_dat,
    aes(
      shape = Class_AB,
      colour = Group_plot,
      alpha = Group_plot,
      size = plot_size_nonfill
    ),
    stroke = 0.8,
    show.legend = c(
      shape = TRUE,
      colour = TRUE,
      alpha = FALSE,
      size = FALSE
    )
  ) +
  
  # Add antibiotics with fillable shapes.
  # These shapes use fill for the inside of the point and black for the outline.
  geom_point(
    data = ab_fill_dat,
    aes(
      shape = Class_AB,
      fill = Group_plot,
      alpha = Group_plot
    ),
    size = plot_point_size,
    colour = "black",
    stroke = 0.3,
    show.legend = c(
      shape = TRUE,
      fill = FALSE,
      alpha = FALSE
    )
  ) +
  
  # Set fill colours for fillable antibiotic shapes.
  scale_fill_manual(
    values = group_cols,
    guide = "none"
  ) +
  
  # Set colours for non-fillable antibiotic shapes and create the Group legend.
  scale_colour_manual(
    name = "Group",
    values = group_cols_nonfill,
    breaks = names(group_cols),
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(
        shape = 16,
        size = group_legend_size,
        alpha = 1
      )
    )
  ) +
  
  # Set transparency by group.
  scale_alpha_manual(
    values = group_alpha,
    guide = "none"
  ) +
  
  # Use the point sizes already defined in the data.
  scale_size_identity(
    guide = "none"
  ) +
  
  # Set antibiotic shapes and create the antibiotic-class legend.
  scale_shape_manual(
    name = "Antibiotic",
    values = ab_shapes,
    breaks = ab_breaks,
    guide = guide_legend(
      override.aes = list(
        fill = "white",
        colour = "black",
        alpha = 1,
        size = legend_point_size,
        stroke = 0.4
      )
    )
  ) +
  
  # Start a new shape scale for pollutants.
  ggnewscale::new_scale("shape") +
  
  # Add pollutants.
  # All pollutants are shown using one pollutant shape and filled by group.
  geom_point(
    data = pol_dat,
    aes(
      shape = Shape_POL,
      fill = Group_plot,
      alpha = Group_plot
    ),
    size = plot_point_size,
    colour = "black",
    stroke = 0.3,
    show.legend = c(
      shape = TRUE,
      fill = FALSE,
      alpha = FALSE
    )
  ) +
  
  # Set pollutant shape and create the pollutant legend.
  scale_shape_manual(
    name = "Pollutant",
    values = c("Pollutant" = 24),
    breaks = "Pollutant",
    guide = guide_legend(
      override.aes = list(
        fill = "white",
        colour = "black",
        alpha = 1,
        size = legend_point_size,
        stroke = 0.3
      )
    )
  ) +
  
  # Add PCA axis labels.
  labs(
    x = lab1,
    y = lab2
  ) +
  
  # Apply journal-style theme.
  theme_journal +
  
  # Format legend size.
  theme(
    legend.key.size = unit(0.35, "cm")
  )

# View Figure 1.
PC1_2_colour_group


# Plot Figure S5.6
# PCA plot of eight physicochemical properties, PC1 vs PC3.
# Points represent compounds and are coloured by prioritisation group.
# Antibiotics are shaped by antibiotic class, while pollutants are shown as
# triangles.

PC1_3_colour_group <- ggplot(pca_plot_df, aes(PC1, PC3)) +
  
  # Add antibiotics with non-fillable shapes.
  # These shapes use colour for the point symbol.
  # No Group points are smaller using plot_size_nonfill.
  geom_point(
    data = ab_nonfill_dat,
    aes(
      shape = Class_AB,
      colour = Group_plot,
      alpha = Group_plot,
      size = plot_size_nonfill
    ),
    stroke = 0.8,
    show.legend = c(
      shape = TRUE,
      colour = TRUE,
      alpha = FALSE,
      size = FALSE
    )
  ) +
  
  # Add antibiotics with fillable shapes.
  # These shapes use fill for the inside of the point and black for the outline.
  geom_point(
    data = ab_fill_dat,
    aes(
      shape = Class_AB,
      fill = Group_plot,
      alpha = Group_plot
    ),
    size = plot_point_size,
    colour = "black",
    stroke = 0.3,
    show.legend = c(
      shape = TRUE,
      fill = FALSE,
      alpha = FALSE
    )
  ) +
  
  # Set fill colours for fillable antibiotic shapes.
  scale_fill_manual(
    values = group_cols,
    guide = "none"
  ) +
  
  # Set colours for non-fillable antibiotic shapes and create the Group legend.
  scale_colour_manual(
    name = "Group",
    values = group_cols_nonfill,
    breaks = names(group_cols),
    drop = FALSE,
    guide = guide_legend(
      override.aes = list(
        shape = 16,
        size = group_legend_size,
        alpha = 1
      )
    )
  ) +
  
  # Set transparency by group.
  scale_alpha_manual(
    values = group_alpha,
    guide = "none"
  ) +
  
  # Use the point sizes already defined in the data.
  scale_size_identity(
    guide = "none"
  ) +
  
  # Set antibiotic shapes and create the antibiotic-class legend.
  scale_shape_manual(
    name = "Antibiotic",
    values = ab_shapes,
    breaks = ab_breaks,
    guide = guide_legend(
      override.aes = list(
        fill = "white",
        colour = "black",
        alpha = 1,
        size = legend_point_size,
        stroke = 0.4
      )
    )
  ) +
  
  # Start a new shape scale for pollutants.
  ggnewscale::new_scale("shape") +
  
  # Add pollutants.
  # All pollutants are shown using one pollutant shape and filled by group.
  geom_point(
    data = pol_dat,
    aes(
      shape = Shape_POL,
      fill = Group_plot,
      alpha = Group_plot
    ),
    size = plot_point_size,
    colour = "black",
    stroke = 0.3,
    show.legend = c(
      shape = TRUE,
      fill = FALSE,
      alpha = FALSE
    )
  ) +
  
  # Set pollutant shape and create the pollutant legend.
  scale_shape_manual(
    name = "Pollutant",
    values = c("Pollutant" = 24),
    breaks = "Pollutant",
    guide = guide_legend(
      override.aes = list(
        fill = "white",
        colour = "black",
        alpha = 1,
        size = legend_point_size,
        stroke = 0.3
      )
    )
  ) +
  
  # Add PCA axis labels.
  labs(
    x = lab1,
    y = lab3
  ) +
  
  # Apply journal-style theme.
  theme_journal +
  
  # Format legend size.
  theme(
    legend.key.size = unit(0.35, "cm")
  )

# View Figure S5.6.
PC1_3_colour_group


## ---------- 1.23) Save Figure 1 & Figure S5.6 ----------

# Save Figure 1: PC1 vs PC2, coloured by priority group.
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Figs"))


ggsave(
  filename = "Figure_1_no_mw.tiff",
  plot = PC1_2_colour_group,
  device = "tiff",
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)



# Save Figure S5.6: PC1 vs PC3, coloured by priority group.
ggsave(
  filename = "Figure_S5_6_PC1_PC3_no_mw.tiff",
  plot = PC1_3_colour_group,
  device = "tiff",
  width = 7,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ============================================================
# pt. 2 PERMANOVA of 8 physicochemical properties excluding molecular weight
# ============================================================
#
# Purpose:
#
# This script uses the final antibiotic and pollutant physicochemical
# information table created during data wrangling.
#
# The aim of this script is to test whether antibiotics and pollutants differ
# in multivariate physicochemical-property space using PERMANOVA.
#
# Molecular weight is excluded from the analysis because compounds were already
# filtered to MW <= 600 Da earlier in the workflow.
#
# The eight physicochemical properties used in the analysis are:
#
# - cLogP
# - cLogD at pH 7.4
# - relative polar surface area
# - hydrogen bond donor count
# - hydrogen bond acceptor count
# - net charge at pH 7.4
# - rotatable bond count
# - globularity score
#
# The script first reads in the final physicochemical information table and
# selects the columns needed for PERMANOVA.
#
# A combined grouping variable is then created, called Group_type, which combines
# compound priority group and compound type. For example:
#
# - Group_1_pollutant
# - Group_1_antibiotic
# - Group_2_pollutant
# - Group_2_antibiotic
# - Group_3_pollutant
# - Group_3_antibiotic
# - NG_pollutant
# - NG_antibiotic
#
# The physicochemical variables are scaled and centred so that all properties
# contribute on comparable scales.
#
# Euclidean distances are then calculated using the scaled physicochemical
# variables.
#
# Before interpreting PERMANOVA, permutation-based tests of multivariate
# dispersion are carried out using betadisper and permutest. This checks whether
# significant PERMANOVA results could be influenced by differences in group
# dispersion rather than differences in group centroids.
#
# Pairwise betadisper tests are then run for selected group comparisons, and
# adjusted p-values are calculated using the Benjamini-Hochberg method.
#
# PERMANOVA is then run using adonis2, followed by pairwise PERMANOVA using
# pairwise.adonis2.
#
# Finally, the script extracts and formats the PERMANOVA and betadisper results,
# adds significance codes, creates smaller result tables for the main priority
# group comparisons, and saves the outputs as Excel files.
#
# ============================================================

## ---------- 2.0) Setup ----------

# Clear the global environment.
# Use with care because this removes all existing objects.
rm(list = ls())

# Load required packages.
library(dplyr)          # data wrangling: select, mutate, filter, case_when and pipes
library(tidyr)          # drop missing rows and separate contrast names
library(tibble)         # create tibbles and convert row names to columns
library(purrr)          # extract pairwise PERMANOVA outputs into one data frame
library(vegan)          # run betadisper, permutest and adonis2
library(writexl)        # save Excel files
library(pairwiseAdonis) # run pairwise PERMANOVA
library(here)         # make file paths portable across different computers


#Get package citations 

# Get the citations for the R packages used in this script.
# These can be used in the methods section or supplementary information if needed.

citation("dplyr")
citation("tidyr")
citation("tibble")
citation("purrr")
citation("vegan")
citation("writexl")
citation("pairwiseAdonis")
citation("here")

#Check R version and session information 

# Check R version.
R.version.string
# R version 4.5.0 (2025-04-11 ucrt)

# Record session information so the package versions used are documented.
sessionInfo()


# Set working directory and read in final data

# Set working directory to where the final physicochemical information file is saved.
setwd(here("Files", "Files_for_coding"))
# List files to check that the final CSV file is there.
list.files()

# Read in the final antibiotic and pollutant physicochemical information table.
df1 <- read.csv("Final_AB_Pol_Physchem_info.csv")

# Set working directory to where the PERMANOVA outputs will be saved.
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Results"))


## ---------- 2.1) Select and rename columns for PERMANOVA ----------

# Select the compound information and physicochemical variables needed for
# PERMANOVA.
#
# Molecular weight is not included because this analysis excludes MW.

df2 <- df1 %>%
  dplyr::select(
    Name = Compound_Name,
    TYPE = CODE_2,
    Class = Class_abbrv,
    Group,
    clogp = logp,
    clogD7_4 = logD_pH7_4_interp,
    Rel_PSA = rel_polar_sa,
    HBD_count = hbd_count,
    HBA_count = hba_count,
    Numeric_charge_7_4 = Charge_pH_7_4_num,
    RB_count = rb_EW,
    Glob_score = glob_EW
  )

# Check structure.
str(df2)


## ---------- 2.2) Prepare physicochemical property columns ----------

# Make a data frame containing only the physicochemical variables used in
# PERMANOVA.

props_cols <- df2 %>%
  dplyr::select(
    clogp,
    clogD7_4,
    Rel_PSA,
    HBD_count,
    HBA_count,
    Numeric_charge_7_4,
    RB_count,
    Glob_score
  )

# Make sure all physicochemical variables are numeric before scaling and
# calculating distances.

props_cols <- props_cols %>%
  mutate(
    clogp = as.numeric(clogp),
    clogD7_4 = as.numeric(clogD7_4),
    Rel_PSA = as.numeric(Rel_PSA),
    HBD_count = as.integer(HBD_count),
    HBA_count = as.integer(HBA_count),
    Numeric_charge_7_4 = as.numeric(Numeric_charge_7_4),
    RB_count = as.integer(RB_count),
    Glob_score = as.numeric(Glob_score)
  )


## ---------- 2.3) Prepare compound metadata columns ----------

# Make a data frame containing compound name, type, class and priority group.

names_cols <- df2 %>%
  dplyr::select(
    Name,
    Type = TYPE,
    Class,
    Group
  )

# Combine compound metadata with physicochemical properties.
properties <- cbind(names_cols, props_cols)



## ---------- 2.4) Create combined group/type variable ----------

# Make a new grouping variable called Group_type.
#
# This combines the compound priority group and compound type, so that
# antibiotics and pollutants can be compared within each priority group.

properties_2 <- properties %>%
  mutate(
    Group_type = case_when(
      Group == "Group 1" & Type == "Antibiotic" ~ "Group_1_antibiotic",
      Group == "Group 2" & Type == "Antibiotic" ~ "Group_2_antibiotic",
      Group == "Group 3" & Type == "Antibiotic" ~ "Group_3_antibiotic",
      Group == "No group" & Type == "Antibiotic" ~ "NG_antibiotic",
      
      Group == "Group 1" & Type == "Pollutant" ~ "Group_1_pollutant",
      Group == "Group 2" & Type == "Pollutant" ~ "Group_2_pollutant",
      Group == "Group 3" & Type == "Pollutant" ~ "Group_3_pollutant",
      Group == "No group" & Type == "Pollutant" ~ "NG_pollutant"
    )
  )

## ---------- 2.5) Prepare PERMANOVA data ----------

# Create a list of the physicochemical variables used in the PERMANOVA.
perm_vars <- c(
  "clogp",
  "clogD7_4",
  "Rel_PSA",
  "HBD_count",
  "HBA_count",
  "Numeric_charge_7_4",
  "RB_count",
  "Glob_score"
)

# Format grouping variables as factors and remove rows with missing
# physicochemical-property values.

perm_dat <- properties_2 %>%
  mutate(
    Type = factor(Type, levels = c("Pollutant", "Antibiotic")),
    Group = factor(Group, levels = c("Group 1", "Group 2", "Group 3", "No group")),
    Group_type = factor(
      Group_type,
      levels = c(
        "Group_1_pollutant",
        "Group_1_antibiotic",
        "Group_2_pollutant",
        "Group_2_antibiotic",
        "Group_3_pollutant",
        "Group_3_antibiotic",
        "NG_pollutant",
        "NG_antibiotic"
      )
    )
  ) %>%
  drop_na(all_of(perm_vars))

# Scale and centre the physicochemical variables so that all properties are on
# comparable scales.
X <- scale(perm_dat[, perm_vars], center = TRUE, scale = TRUE)

# Calculate Euclidean distance matrix using the scaled physicochemical variables.
d <- dist(X, method = "euclidean")

# Set number of permutations.
n_perm <- 10000


## ---------- 2.6) Test multivariate dispersion ----------

# Run betadisper to test whether groups differ in multivariate dispersion.
# This is important because PERMANOVA can be influenced by differences in
# dispersion as well as differences in group centroids.

bd <- betadisper(d, group = perm_dat$Group_type)

# Run ANOVA on the betadisper object.
anova(bd)

# Run permutation test on betadisper.
set.seed(123)
permutest <- permutest(bd, permutations = n_perm)

# Run pairwise permutation tests for multivariate dispersion.
pairwise_permutest <- permutest(
  bd,
  permutations = n_perm,
  pairwise = TRUE
)

pairwise_permutest

# Plot distances to group median.
boxplot(
  bd,
  ylab = "Distance to group median"
)


## ---------- 2.7) Run selected pairwise betadisper tests ----------

# Add the scaled variables back into the PERMANOVA data frame.
# This means pairwise betadisper tests use the same scaling as the full analysis.

perm_dat_scaled <- perm_dat
perm_dat_scaled[, perm_vars] <- X

# Function to run pairwise betadisper between two selected Group_type levels.
run_pairwise_betadisper <- function(dat, group1, group2, vars, n_perm = 10000) {
  
  sub <- dat %>%
    filter(Group_type %in% c(group1, group2)) %>%
    droplevels()
  
  if (nrow(sub) == 0) {
    stop("No rows found for one or both groups.")
  }
  
  if (length(unique(sub$Group_type)) != 2) {
    stop("Subset does not contain exactly two Group_type levels.")
  }
  
  # Use already-scaled variables.
  X_sub <- sub[, vars]
  
  # Calculate Euclidean distance matrix for the two-group subset.
  d_sub <- dist(X_sub, method = "euclidean")
  
  # Run betadisper for the two-group subset.
  bd <- betadisper(d_sub, group = sub$Group_type)
  
  # Run permutation test.
  bd_perm <- permutest(bd, permutations = n_perm)
  
  # Return results as one row.
  tibble(
    contrast = paste(group1, "vs", group2),
    group1 = group1,
    group2 = group2,
    n_group1 = sum(sub$Group_type == group1),
    n_group2 = sum(sub$Group_type == group2),
    F = bd_perm$tab[1, "F"],
    p = bd_perm$tab[1, "Pr(>F)"],
    mean_dist_group1 = unname(bd$group.distances[group1]),
    mean_dist_group2 = unname(bd$group.distances[group2])
  )
}


## ---------- 1.8) Define pairwise betadisper contrasts ----------

# Define selected contrasts for pairwise multivariate dispersion testing.
# The first three contrasts compare antibiotics and pollutants within each
# priority group.

contrast_df <- tribble(
  ~group1, ~group2,
  "Group_1_pollutant", "Group_1_antibiotic",
  "Group_2_pollutant", "Group_2_antibiotic",
  "Group_3_pollutant", "Group_3_antibiotic",
  "NG_pollutant",      "NG_antibiotic",
  "Group_1_pollutant", "NG_antibiotic",
  "Group_2_pollutant", "NG_antibiotic",
  "Group_3_pollutant", "NG_antibiotic",
  "Group_2_pollutant", "Group_1_antibiotic",
  "Group_3_pollutant", "Group_1_antibiotic",
  "Group_1_pollutant", "Group_2_antibiotic",
  "Group_3_pollutant", "Group_2_antibiotic",
  "Group_1_pollutant", "Group_3_antibiotic",
  "Group_2_pollutant", "Group_3_antibiotic",
  "NG_pollutant",      "Group_1_antibiotic",
  "NG_pollutant",      "Group_2_antibiotic",
  "NG_pollutant",      "Group_3_antibiotic"
)


## ---------- 2.9) Format pairwise betadisper results ----------

# Run pairwise betadisper for all selected contrasts.
set.seed(123)

betadisper_results <- bind_rows(
  lapply(seq_len(nrow(contrast_df)), function(i) {
    run_pairwise_betadisper(
      dat = perm_dat_scaled,
      group1 = contrast_df$group1[i],
      group2 = contrast_df$group2[i],
      vars = perm_vars,
      n_perm = n_perm
    )
  })
) %>%
  mutate(
    p_BH = p.adjust(p, method = "BH")
  ) %>%
  select(-group1, -group2)

# Rename columns so that it is clear these results are from betadisper.
betadisper_results <- betadisper_results %>%
  select(
    contrast,
    n_group1,
    n_group2,
    mean_dist_group1,
    mean_dist_group2,
    F_bd = F,
    p_bd = p,
    p_BH_bd = p_BH
  )

# Add significance codes based on adjusted p-values.
betadisper_results <- betadisper_results %>%
  mutate(
    sig_code = case_when(
      p_BH_bd <= 0.001 ~ "***",
      p_BH_bd <= 0.01  ~ "**",
      p_BH_bd <= 0.05  ~ "*",
      p_BH_bd <= 0.1   ~ ".",
      TRUE             ~ ""
    )
  )

## ---------- 2.10) Create betadisper results table for priority groups ----------

# Create a smaller betadisper results table showing only the antibiotic vs
# pollutant comparisons within Priority Groups 1-3.

betadisper_results_view <- betadisper_results %>%
  filter(
    contrast %in% c(
      "Group_1_pollutant vs Group_1_antibiotic",
      "Group_2_pollutant vs Group_2_antibiotic",
      "Group_3_pollutant vs Group_3_antibiotic"
    )
  ) %>%
  select(
    contrast,
    F_bd,
    p_BH_bd,
    sig_code
  )

# Print significance codes.
cat("Signif. codes:  0 '***',  0.001 '**',  0.01 '*',  0.05 '.',  0.1 ' ' , 1")

# View selected betadisper results.
betadisper_results_view

# A tibble: 3 × 4
# contrast                                 F_bd p_BH_bd sig_code
# <chr>                                   <dbl>   <dbl> <chr>   
# 1 Group_1_pollutant vs Group_1_antibiotic  7.97  0.0222 "*"     
# 2 Group_2_pollutant vs Group_2_antibiotic  1.74  0.227  ""      
#3 Group_3_pollutant vs Group_3_antibiotic  5.12  0.0810 "."  

#Permutation-based tests of multivariate dispersion showed significant difference in dispersion between antibiotics and pollutants in Priority Group 1 (F =7.97, adjusted p = 0.02*), but not in Group 2 (F = 1.74, adjusted p = 0.23) and Group 3 (F = 5.12, adjusted p = 0.08). This suggests that the PERMANOVA results for Priority Group 1 may be influenced by differences in multivariate dispersion, while the results for Priority Groups 2 and were not primarily driven by significant dispersion differences.

## ---------- 2.11) Run overall PERMANOVA ----------

# Run PERMANOVA using adonis2 to test whether physicochemical-property space
# differs by Group_type.

set.seed(123)

permanova_Group_type_all <- adonis2(
  X ~ Group_type,
  data = perm_dat,
  method = "euclidean",
  permutations = n_perm
)

permanova_Group_type_all

# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 10000

# adonis2(formula = X ~ Group_type, data = perm_dat, permutations = n_perm, method = "euclidean")
# Df SumOfSqs      R2      F    Pr(>F)    
# Model      7   818.05 0.46907 26.631 9.999e-05 ***
# Residual 211   925.95 0.53093                     
# Total    218  1744.00 1.00000                     
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#physicochemical-property space differs by Group_type p = 9.999e-05 ***

## ---------- 2.12) Run pairwise PERMANOVA ----------

# Run pairwise PERMANOVA between all Group_type levels.

set.seed(123)

pairwise_results <- pairwise.adonis2(
  X ~ Group_type,
  data = perm_dat,
  method = "euclidean",
  permutations = n_perm
)


## ---------- 2.13) Extract full pairwise PERMANOVA results ----------

# Extract all rows from each pairwise PERMANOVA result into one data frame.

pairwise_full_df <- pairwise_results %>%
  .[names(.) != "parent_call"] %>%
  imap_dfr(
    ~ .x %>%
      rownames_to_column("term") %>%
      mutate(contrast = .y),
    .id = NULL
  )

# Split contrast names into group1 and group2.
pairwise_full_df <- pairwise_full_df %>%
  separate_wider_delim(
    contrast,
    delim = "_vs_",
    names = c("group1", "group2")
  )

# Adjust p-values using the Benjamini-Hochberg method.
pairwise_full_df <- pairwise_full_df %>%
  mutate(
    p_BH = p.adjust(`Pr(>F)`, method = "BH")
  )

# Reorder columns.
pairwise_full_df <- pairwise_full_df %>%
  select(
    group1,
    group2,
    term,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    p_BH
  )


## ---------- 2.14) Extract pairwise PERMANOVA model rows ----------

# Extract the main model row from each pairwise PERMANOVA result.
# This produces a simpler results table 

pairwise_df <- pairwise_results %>%
  .[names(.) != "parent_call"] %>%
  imap_dfr(
    ~ tibble(
      contrast = .y,
      Df = .x$Df[1],
      SumOfSqs = .x$SumOfSqs[1],
      R2 = .x$R2[1],
      F = .x$F[1],
      p = .x$`Pr(>F)`[1]
    )
  ) %>%
  separate_wider_delim(
    contrast,
    delim = "_vs_",
    names = c("group1", "group2")
  )

# Adjust p-values using the Benjamini-Hochberg method.
pairwise_df <- pairwise_df %>%
  mutate(
    p_BH = p.adjust(p, method = "BH")
  )

# Rename PERMANOVA columns.
pairwise_df <- pairwise_df %>%
  select(
    group1,
    group2,
    Df_adonis = Df,
    SumOfSqs_adonis = SumOfSqs,
    R2_adonis = R2,
    F_adonis = F,
    p_adonis = p,
    p_BH_adonis = p_BH
  )

# Add significance codes based on adjusted p-values.
pairwise_df <- pairwise_df %>%
  mutate(
    sig_code = case_when(
      p_BH_adonis <= 0.001 ~ "***",
      p_BH_adonis <= 0.01  ~ "**",
      p_BH_adonis <= 0.05  ~ "*",
      p_BH_adonis <= 0.1   ~ ".",
      TRUE                 ~ ""
    )
  )

# Make a contrast column 
pairwise_df <- pairwise_df %>%
  mutate(
    contrast = paste(group1, "vs", group2)
  ) %>%
  select(
    contrast,
    everything()
  )

## ---------- 2.15) Create PERMANOVA results table for priority groups ----------

# Create a smaller PERMANOVA results table showing only the antibiotic vs
# pollutant comparisons within Priority Groups 1-3.

contrast_order <- c(
  "Group_1_pollutant vs Group_1_antibiotic",
  "Group_2_pollutant vs Group_2_antibiotic",
  "Group_3_pollutant vs Group_3_antibiotic"
)

pairwise_df_view <- pairwise_df %>%
  filter(
    contrast %in% contrast_order
  ) %>%
  select(
    contrast,
    R2_adonis,
    F_adonis,
    p_BH_adonis,
    sig_code,
    Df_adonis,
    SumOfSqs_adonis
  ) %>%
  arrange(
    factor(contrast, levels = contrast_order)
  )

# View selected PERMANOVA results.

head(pairwise_df_view)

# A tibble: 3 × 7
# contrast                                R2_adonis F_adonis p_BH_adonis sig_code Df_adonis SumOfSqs_adonis
# <chr>                                       <dbl>    <dbl>       <dbl> <chr>        <dbl>           <dbl>
# 1 Group_1_pollutant vs Group_1_antibiotic     0.158     7.50    0.000280 ***              1            23.4
# 2 Group_2_pollutant vs Group_2_antibiotic     0.326     9.20    0.001000 ***              1            82.8
# 3 Group_3_pollutant vs Group_3_antibiotic     0.653    20.7     0.00131  **               1            19.0

#Despite visual overlap within the PCA space, pairwise PERMANOVA showed significant differences between antibiotics and pollutants meeting Priority Group 1 (R2 = 0.16, pseudo-F = 7.50, adjusted p = 0.0003***), Group 2 (R2 = 0.33 , pseudo-F = 9.20 , adjusted p =0.001***), and Group 3 (R2 =0.65 , pseudo-F = 20.7 , adjusted p = 0.0013**) criteria. Although, Priority Group 1 and Group 2 showed weaker separation between antibiotics and pollutants than Group 3, based on their smaller PERMANOVA effect sizes (R2 = 0.16 and 0.33 vs 0.65), aligning with the closer proximity observed in the PCA plot. 

## ---------- 2.16) Save PERMANOVA and betadisper results ----------

# Set working directory to where the PERMANOVA outputs will be saved.
setwd(here("Files", "PCA_PERMANOVA_NO_MW_Results"))


# Save pairwise PERMANOVA, full pairwise PERMANOVA and betadisper results
# into one Excel workbook with separate sheets.

write_xlsx(
  list(
    pairwise_PERMANOVA = pairwise_df,
    full_pairwise_PERMANOVA = pairwise_full_df,
    betadisper_results = betadisper_results
  ),
  "PERMANOVA_betadisper_results_no_mw.xlsx"
)
