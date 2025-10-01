# Updated Bird Project Script 30/09/25

# Load libraries
library(tidyverse)
library(vegan)
library(broom)
library(grid)
library(gridExtra)


# Bird Data ----

# Load data
bird <- read_csv("your_path")
View(bird)

# Reshape to long format
birds_long <- bird %>%
  pivot_longer(cols = starts_with("Site"),
               names_to = "Site",
               values_to = "Count")

# Assign habitat type based on site number
birds_long <- birds_long %>%
  mutate(Habitat = case_when(
    Site %in% c("Site_1", "Site_5", "Site_6") ~ "Plantation",
    Site %in% c("Site_2", "Site_3", "Site_4") ~ "Heathland",
    Site %in% c("Site_7", "Site_8", "Site_9") ~ "Caledonian"
  ))

view(birds_long)

# Summarise total abundance per site
site_summary <- birds_long %>%
  group_by(Site, Habitat) %>%
  summarise(Total_Count = sum(Count), .groups = "drop")

view(site_summary)

# Boxplot of abundance by habitat (DO NOT USE BOXPLOT - not enough sites per habitat)
(ggplot(site_summary, aes(x = Habitat, y = Total_Count, fill = Habitat)) +
    geom_boxplot(size = 1) +
    scale_fill_manual(values = c("green4", "darkorange", "dodgerblue3")) +
    theme_test() +
    theme(axis.text = element_text(size = 12, face = "bold"), 
          axis.title = element_text(size = 12, face = "bold"), legend.position = "none") +
    labs(y = "Bird Species Richness", x = "Habitat"))

# Grouping by total count
site_summary_size <- site_summary %>% 
  group_by(Habitat, Total_Count) %>% 
  summarise(count = n())

# Dot plot - preferred over boxplot
ggplot(site_summary_size, aes(x = Habitat, y = Total_Count, color = Habitat)) +
  geom_point(shape = 15, aes(size = as.factor(count))) +
  scale_color_manual(values = c("green4", "darkorange", "dodgerblue3")) +
  scale_size_manual(values = c(3, 5)) +
  theme_test() +
  theme(
    axis.text   = element_text(size = 12, face = "bold"),
    axis.title  = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")  # make "Count of points" bigger
  ) +
  labs(y = "Bird Species Richness", x = "Habitat", size = "Count of points") +
  guides(color = "none")


# ANOVA (DO NOT USE)
anova_result <- aov(Total_Count ~ Habitat, data = site_summary)
summary(anova_result) 
# Not valid as we don't have enough sites per habitat (minimum 10 required)
# Ignore all anova and tukey's post-hoc test

# Checking for normality
resid_anova <- residuals(anova_result)
shapiro.test(resid_anova)

# Checking for equal variances
bartlett.test(Count~Habitat, data = birds_long)

# Running Tukey's HSD post-hoc test on the anova output and setting the confidence level to 95%
(bird_test <- TukeyHSD(anova_result, conf.level=.95))

# To convert results into a better presented format of the summary table you can use the broom package
(tukey_results <- broom::tidy(bird_test)) # Creating tidy data frame using broom

# Plotting Tukey's test result
plot(bird_test)


# NMDS ----

# NMDS preparation

# Transpose so rows = sites, cols = species
bird_matrix <- t(bird)

# Remove the Species column so only counts remain
bird_counts <- bird[ , -1]

# Transpose so rows = sites, cols = species
bird_matrix <- t(bird_counts)

# Make sure everything is numeric
bird_matrix <- apply(bird_matrix, 2, as.numeric)

# Convert back to matrix
bird_matrix <- as.matrix(bird_matrix)

# Create site-by-species matrix
bird_matrix <- birds_long %>%
  dplyr::select(Site, Species, Count) %>%
  tidyr::pivot_wider(  
    names_from = Species,
    values_from = Count,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Site")

view(bird_matrix)

# Run NMDS
nmds_result <- metaMDS(bird_matrix, distance = "bray", k = 2, trymax = 100)

# Kruskal's stress test
nmds_result$stress
# 0.103152 - acceptable fit

# Keep habitat info separate
habitats <- birds_long %>%
  distinct(Site, Habitat) %>%
  arrange(match(Site, rownames(bird_matrix))) %>%
  pull(Habitat)

# Bray dissimilarity
bray_dist <- vegdist(bird_matrix, method = "bray")

# Run PERMANOVA
permanova <- adonis2(bray_dist ~ habitats)
print(permanova)
# 0.027 - significant

# Check if group dispersions differ
disp <- betadisper(bray_dist, habitats)
anova(disp)
# 0.7561 - not significant which is good as dispersions are not different between groups

# Plot dispersion ellipses (optional)
plot(disp)


# Preparation for plotting NMDS

# Create site scores
site_scores <- scores(nmds_result, display = "sites") %>% as.data.frame() %>% rownames_to_column("Site")
site_scores$Habitat <- habitats
cols_points <- c("Plantation"="dodgerblue3","Heathland"="darkorange","Caledonian"="green4")
cols_fill   <- c("Plantation"="dodgerblue3","Heathland"="darkorange","Caledonian"="green4")

# Ensure factor levels match your colour vectors
site_scores$Habitat <- factor(site_scores$Habitat, levels = names(cols_fill))

# Compute SE (standard error) ellipses
ell <- ordiellipse(
  nmds_result,
  groups      = site_scores$Habitat,
  kind        = "se",
  conf        = 0.95,
  draw        = "none")

# Convert ellipse parameters to coordinates (robust to near-singular covariances)
ellipse_df <- imap_dfr(ell, function(e, g) {
  if (is.null(e$cov) || anyNA(e$cov)) return(NULL)
  theta  <- seq(0, 2*pi, length.out = 200)
  circle <- cbind(cos(theta), sin(theta))
  ev     <- eigen(e$cov, symmetric = TRUE)
  axes   <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), 2, 2)
  coords <- t(apply(circle, 1, function(v) e$center + e$scale * (axes %*% v)))
  tibble(NMDS1 = coords[,1], NMDS2 = coords[,2], Habitat = g)})

# Plotting NMDS
ggplot(site_scores, aes(NMDS1, NMDS2, color = Habitat, fill = Habitat)) +
  geom_point(size = 3) +
  geom_polygon(
    data = ellipse_df,
    aes(group = Habitat),
    alpha = 0.2,
    color = NA   # no extra border to avoid second legend
  ) +
  scale_color_manual(values = cols_points) +
  scale_fill_manual(values = cols_fill, guide = "none") +  # remove fill legend
  coord_equal() +
  theme_test(base_size = 14) +
  theme(legend.position = "bottom", axis.title.x = element_text(size = 13, face = "bold"),  # X-axis title size
        axis.title.y = element_text(size = 13, face = "bold"), legend.text  = element_text(size = 13, face = "bold"),  # change legend labels
        legend.title = element_text(size = 15, face = "bold")) +
  labs(x = "NMDS1", y = "NMDS2")


# SIMPER analysis ----

# Run SIMPER analysis
sim <- simper(bird_matrix, group = habitats, permutations = 999)

# Look at the summary across all group comparisons
summary(sim)

# Convert all pairwise results into a tidy dataframe (stacking all pairwise comparisons)
sim_summary <- map_df(names(sim), function(comp) {
  df <- as.data.frame(summary(sim, ordered = TRUE)[[comp]])
  df <- tibble::rownames_to_column(df, var = "Species")
  df$Comparison <- comp
  df})

# Summarise across ALL comparisons
species_contrib <- sim_summary %>%
  group_by(Species) %>%
  summarise(
    mean_contribution = mean(average, na.rm = TRUE),
    sum_contribution  = sum(average, na.rm = TRUE),
    n_contrasts       = n()
  ) %>%
  arrange(desc(mean_contribution))

# Look at the top 10 species across all habitats
head(species_contrib, 10)

ggplot(species_contrib %>% slice_max(mean_contribution, n = 10),
       aes(x = reorder(Species, mean_contribution), y = mean_contribution)) +
  geom_col(fill = "goldenrod1") +
  coord_flip() +
  labs(x = "Species", y = "Mean Contribution Across All Habitat Comparisons") +
  theme_minimal(base_size = 14)


# Vegetation Data ----

# Load vegetation data
veg <- read_csv("your_path")
# Inspect
glimpse(veg)

# Convert all data to numeric
comm_matrix <- as.data.frame(lapply(veg, as.numeric))

# Assign rownames as site IDs
rownames(comm_matrix) <- paste0("Site_", seq_len(nrow(comm_matrix)))

# Define habitat factor
habitat <- factor(c(rep("Plantation", 3),
                    rep("Heathland", 3),
                    rep("Caledonian", 3)))

# Check
comm_matrix
habitat

# Run PERMANOVA
permanova <- adonis2(comm_matrix ~ habitat, method = "bray")
print(permanova)

# Run NMDS
nmds_result <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)
nmds_result$stress   # check stress
# 0.0156 - good fit

# Extract site scores
site_scores <- scores(nmds_result, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Site")

# Add habitat info
site_scores$Habitat <- factor(habitat, levels = c("Plantation", "Heathland", "Caledonian"))

# Define colors
cols_points <- c("Plantation"="dodgerblue3",
                 "Heathland"="darkorange",
                 "Caledonian"="green4")
cols_fill   <- cols_points

# Compute ellipses
ell <- ordiellipse(
  nmds_result,
  groups = site_scores$Habitat,
  kind   = "se",
  conf   = 0.95,
  draw   = "none")

# Convert ellipses to coordinates
ellipse_df <- imap_dfr(ell, function(e, g) {
  if (is.null(e$cov) || anyNA(e$cov)) return(NULL)
  theta  <- seq(0, 2*pi, length.out = 200)
  circle <- cbind(cos(theta), sin(theta))
  ev     <- eigen(e$cov, symmetric = TRUE)
  axes   <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), 2, 2)
  coords <- t(apply(circle, 1, function(v) e$center + e$scale * (axes %*% v)))
  tibble(NMDS1 = coords[,1], NMDS2 = coords[,2], Habitat = g)})

# Plot NMDS with ellipses
ggplot(site_scores, aes(NMDS1, NMDS2, color = Habitat, fill = Habitat)) +
  geom_point(size = 3) +
  geom_polygon(
    data = ellipse_df,
    aes(group = Habitat),
    alpha = 0.2,
    color = NA
  ) +
  scale_color_manual(values = cols_points) +
  scale_fill_manual(values = cols_fill, guide = "none") +
  coord_equal() +
  theme_test(base_size = 14) +
  theme(legend.position = "bottom") +
  labs(x = "NMDS1", y = "NMDS2")
