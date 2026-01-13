# ==============================================================================
# 1. DEPENDENCIES & ENVIRONMENT SETUP
# ==============================================================================
library(phyloseq)    # Microbiome data management
library(ggplot2)     # Graphics
library(reshape2)    # Data reshaping (wide-to-long)
library(patchwork)   # Plot assembly
library(rcartocolor) # Color-blind safe palettes
library(car)         # Statistical testing (Levene's)
library(RColorBrewer)# Color palettes
library(corrplot)    # Base-R correlation heatmaps
library(SpiecEasi)   # SparCC and microbial network inference

# ==============================================================================
# 2. DATA IMPORT & TAXONOMIC CLEANING
# ==============================================================================
beef <- import_biom("Species.biom")

# Remove prefix notation (e.g., "k__", "p__") from taxonomy strings
beef@tax_table@.Data <- substring(beef@tax_table@.Data, 4) 

# Assign standard rank names
colnames(beef@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# --- Normalization (Scaling to Median Depth) ---
total = median(sample_sums(beef))
standf = function(x, t = total) round(t * (x / sum(x)))
beef = transform_sample_counts(beef, standf)

# Filter for Bacteria only
beef_bac <- subset_taxa(beef, Kingdom == "Bacteria")

# ==============================================================================
# 3. TAXONOMIC ABUNDANCE VISUALIZATION (Phylum & Genus)
# ==============================================================================

# [!] DATA PREPARATION STEP [!]
# THE FOLLOWING SECTIONS REQUIRE .CSV FILES EXPORTED FROM MICROBIOMEANALYST.CA. 
# Ensure your BIOM file has been processed and results are saved in the working directory.

# --- Phylum - Location ---
phy_loc <- read.csv("phylum_abundance_location.csv") 
phy_loc <- melt(phy_loc, id = c("Sample", "Location")) 
colnames(phy_loc)[3:4] <- c("Phylum", "Proportion")

# Removes empty/zero rows dynamically
phy_loc <- subset(phy_loc, Proportion > 0 & !is.na(Proportion))

phy_loc$Sample <- as.factor(phy_loc$Sample) 
phy_loc$Location <- factor(phy_loc$Location, levels = c("Hotbox", "Cooler"))

plot1 <- ggplot(phy_loc, aes(x = Sample, y = Proportion, fill = Phylum)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Location, scales = "free_x", strip.position = "top") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() +
  theme(
    legend.title = element_text(face = "bold"), 
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(color = "black", fill = "white")
  )
A <- plot1 + scale_fill_carto_d(palette = "Safe")

# --- Phylum - Year ---
phy_time <- read.csv("phylum_abundance_time.csv")
phy_time <- melt(phy_time, id = c("Sample", "Year"))
colnames(phy_time)[3:4] <- c("Phylum", "Proportion")

phy_time <- subset(phy_time, Proportion > 0 & !is.na(Proportion))
phy_time$Sample <- as.factor(phy_time$Sample)
phy_time$Year <- factor(phy_time$Year, levels = c("2019", "2021"))

plot2 <- ggplot(phy_time, aes(x = Sample, y = Proportion, fill = Phylum)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Year, scales = "free_x", strip.position = "top") + 
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.title = element_text(face = "bold"), 
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(color = "black", fill = "white")
  ) 
B <- plot2 + scale_fill_carto_d(palette = "Safe")

# --- Genus - Location ---
gen_loc <- read.csv("genus_abundance_location.csv", header = TRUE)
gen_loc <- melt(gen_loc, id = c("Sample", "Location"))
colnames(gen_loc)[3:4] <- c("Genus", "Proportion")

gen_loc <- subset(gen_loc, Proportion > 0 & !is.na(Proportion))
gen_loc$Sample <- as.factor(gen_loc$Sample)
gen_loc$Location <- factor(gen_loc$Location, levels = c("Hotbox", "Cooler"))

plot3 <- ggplot(gen_loc, aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Location, scales = "free_x", strip.position = "top") + 
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "italic"), # ITALICIZES GENUS NAMES
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(color = "black", fill = "white")
  )
C <- plot3 + scale_fill_carto_d(palette = "Safe")

# --- Genus - Year ---
gen_time <- read.csv("genus_abundance_time.csv", header = TRUE)
gen_time <- melt(gen_time, id = c("Sample", "Year"))
colnames(gen_time)[3:4] <- c("Genus", "Proportion")

gen_time <- subset(gen_time, Proportion > 0 & !is.na(Proportion))
gen_time$Sample <- as.factor(gen_time$Sample)
gen_time$Year <- factor(gen_time$Year, levels = c("2019", "2021"))

plot4 <- ggplot(gen_time, aes(x = Sample, y = Proportion, fill = Genus)) +
  geom_bar(stat = "identity") + 
  facet_wrap(~Year, scales = "free_x", strip.position = "top") + 
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(
    legend.title = element_text(face = "bold"), 
    legend.text = element_text(face = "italic"), # ITALICIZES GENUS NAMES
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(color = "black", fill = "white")
  ) 
D <- plot4 + scale_fill_carto_d(palette = "Safe")


Final_Plot <- wrap_plots(A, B, C, D, ncol = 2, nrow = 2)

ggsave("Taxonomic_Abundance_Grid.png", 
       plot = Final_Plot, 
       width = 16, 
       height = 12, 
       units = "in", 
       dpi = 300)

# ==============================================================================
# 4. ALPHA DIVERSITY GRAPHICS & STATISTICS
# ==============================================================================

# --- Alpha Diversity Plots ---

alpha1 <- plot_richness(beef_bac, x = "Location", color = "Location", 
                        measures = c("Observed", "Shannon"), title = "Alpha Diversity - Hotbox vs. Cooler") + 
  geom_boxplot(alpha = 0.3) + 
  scale_x_discrete(limits = c("Hotbox", "Cooler")) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_carto_d(palette = "Safe")

alpha2 <- plot_richness(beef_bac, x = "Year", color = "Year", 
                        measures = c("Observed", "Shannon"), title = "Alpha Diversity - 2019 vs. 2021") + 
  geom_boxplot(alpha = 0.3) + 
  scale_x_discrete(limits = c("2019", "2021")) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_carto_d(palette = "Safe")

# Arrange the two alpha diversity plots vertically
alpha1_2 <- alpha1 / alpha2 

# --- Mann-Whitney U (Wilcoxon Rank-Sum) Test ---

# Calculate diversity metrics
alpha_metrics <- estimate_richness(beef_bac, measures = c("Shannon", "Observed"))

# Prepare metadata and merge with metrics
metadata <- data.frame(sample_data(beef_bac))
alpha_stats <- cbind(alpha_metrics, metadata)
alpha_stats$Location <- as.factor(alpha_stats$Location)
alpha_stats$Year <- as.factor(alpha_stats$Year)

# --- Statistical analysis for 'Shannon' index ---

# Test difference by location
levene_shannon_location  <- leveneTest(Shannon ~ Location, data = alpha_stats)
wilcox_shannon_location <- wilcox.test(Shannon ~ Location, data = alpha_stats)
# Test difference by year
levene_shannon_year  <- leveneTest(Shannon ~ Year, data = alpha_stats)
wilcox_shannon_year <- wilcox.test(Shannon ~ Year, data = alpha_stats)


# --- Statistical analysis for 'Observed' ---

# Test difference by location
levene_observed_location  <- leveneTest(Observed ~ Location, data = alpha_stats)
wilcox_observed_location <- wilcox.test(Observed ~ Location, data = alpha_stats)

# Test difference by year
levene_observed_year  <- leveneTest(Observed ~ Year, data = alpha_stats)
wilcox_observed_year <- wilcox.test(Observed ~ Year, data = alpha_stats)


# --- PRINT RESULTS ---
print(levene_shannon_location)
print(levene_shannon_year)
print(wilcox_shannon_location)
print(wilcox_shannon_year)

print(levene_observed_location)
print(levene_observed_year)
print(wilcox_observed_location)
print(wilcox_observed_year)

# ==============================================================================
# 5. CORRELATION PLOT OF TOP 20 GENERA
# ==============================================================================

# Agglomerate (merge) all OTUs belonging to the same Genus
beef_bac_genus <- tax_glom(beef_bac, taxrank = "Genus")

# Extract the OTU table as a matrix to calculate total abundance
beef_bac_abund_genus <- otu_table(beef_bac_genus)
beef_bac_abund_matrix_genus <- as.matrix(beef_bac_abund_genus)

# Calculate the sum of reads for each genus across all samples
beef_genus_total_abundance <- rowSums(beef_bac_abund_matrix_genus)

# Sort genera by total abundance and extract the names (IDs) of the Top 20
sorted_beef_bac_genera <- sort(beef_genus_total_abundance, decreasing = TRUE)
top20_beef_bac_genera_ids <- names(sorted_beef_bac_genera)[1:20]

# Filter the phyloseq object to keep only these Top 20 taxa
beef_bac_top20_genera <- subset_taxa(beef_bac_genus, taxa_names(beef_bac_genus) %in% top20_beef_bac_genera_ids)


# --- DATA RE-FORMATTING & NAMING ---


# Extract abundance data and transpose (SparCC requires samples in rows, taxa in columns)
df_top20_beef_bac_abund_genus <- as.data.frame(as.matrix(otu_table(beef_bac_top20_genera)))
df_top20_beef_bac_abund_genus <- t(df_top20_beef_bac_abund_genus)

# Create a mapping vector to translate OTU IDs into human-readable Genus names
tax_info <- tax_table(beef_bac_top20_genera)
genus_names <- as.character(tax_info[, "Genus"])
genus_name_map <- setNames(genus_names, rownames(tax_info))

# Clean up Genus names for R compatibility and assign them to column names
column_names_to_use <- genus_name_map[colnames(df_top20_beef_bac_abund_genus)]
column_names_to_use <- make.names(column_names_to_use, unique = TRUE)
colnames(df_top20_beef_bac_abund_genus) <- column_names_to_use

# --- STATISTICAL ANALYSIS: SparCC & BOOTSTRAPPING ---

# Run the SparCC algorithm to compute the correlation matrix
sparcc_results <- sparcc(df_top20_beef_bac_abund_genus)
beef_cor_matrix <- sparcc_results$Cor

# Set seed for reproducibility and run bootstrapping (R = 1000 permutations)
# This step calculates p-values to determine if correlations are significant
set.seed(13)
sparccboot_results <- sparccboot(df_top20_beef_bac_abund_genus, R = 1000)
beef_sig_sparcc <- pval.sparccboot(sparccboot_results)
beef_pval_vector <- beef_sig_sparcc$pvals

# Reconstruct a symmetrical P-value matrix from the SparCC results
p_value_matrix <- matrix(NA, 
                         nrow = nrow(beef_cor_matrix), 
                         ncol = ncol(beef_cor_matrix))

# Fill the upper and lower triangles of the matrix
p_value_matrix[upper.tri(p_value_matrix)] <- beef_pval_vector
p_value_matrix[lower.tri(p_value_matrix)] <- t(p_value_matrix)[lower.tri(p_value_matrix)]
p_value_matrix[is.na(p_value_matrix)] <- 0  # Replace diagonal NAs with 0

# Re-attach Genus names to the resulting matrices
colnames(beef_cor_matrix) <- column_names_to_use
row.names(beef_cor_matrix) <- column_names_to_use
colnames(p_value_matrix) <- column_names_to_use
row.names(p_value_matrix) <- column_names_to_use

# --- VISUALIZATION (CORRPLOT) ---

# Generate the heatmap
beef_cor_plot <- corrplot(beef_cor_matrix,
                          cex.main = 2,  
                          font.main = 2,  
                          type = "upper",      # Only show the top half for clarity
                          order = "hclust",     # Reorder taxa based on similar correlation patterns
                          hclust.method = "ward.D2",
                          addrect = 2,          # Highlight the 2 most distinct clusters
                          p.mat = p_value_matrix, 
                          sig.level = c(0.001, 0.01, 0.05), # Visual significance thresholds
                          insig = "label_sig",  # Use * symbols for significant correlations
                          method = "color",     # Map correlation strength to colors
                          font = 4,             # Bold-Italic font for labels
                          pch.col = "azure",    # Color of the significance markers
                          pch.cex = 1.5,
                          tl.col = "black",     # Text Label color
                          tl.srt = 90,          # Vertical rotation for labels
                          tl.cex = 0.9, 
                          cl.cex = 0.8,         # Color Legend font size
                          cl.length = 11,
                          col = brewer.pal(n = 10, name = "RdBu"), # Red-to-Blue scale
                          diag = FALSE)         # Don't show the 1:1 diagonal

# Add Title and Legend labels using Margin Text
mtext("Taxonomic Correlation of Top 20 Genera Within \n Beef Processing Facility Drain Samples", 
      line = -1, adj = 0.55, cex = 2.5)  

mtext(text = "Correlation Coefficient", 
      at = 11, side = 4, line = -20, cex = 1.5, las = 0, srt = 270)


