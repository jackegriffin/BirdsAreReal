# Our Code

# Bird Project 08/09/25

# Load libraries

library(tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)
library(broom)
library(grid)
library(gridExtra)

# Bird Data ----

# Import data
bird <- read_csv("Data/Birdresults.csv")
View(bird)


# Reshape to long format
birds_long <- bird %>%
  pivot_longer(cols = starts_with("Site"),
               names_to = "Site",
               values_to = "Count")

# Assign habitat type based on site number
birds_long <- birds_long %>%
  mutate(Habitat = case_when(
    Site %in% c("Site_1", "Site_2", "Site_3") ~ "Plantation",
    Site %in% c("Site_4", "Site_5", "Site_6") ~ "Heathland",
    Site %in% c("Site_7", "Site_8", "Site_9") ~ "Caledonian"
  ))

# Summarise total abundance per site
site_summary <- birds_long %>%
  group_by(Site, Habitat) %>%
  summarise(Total_Count = sum(Count), .groups = "drop")

# Boxplot of abundance by habitat
(ggplot(site_summary, aes(x = Habitat, y = Total_Count, fill = Habitat)) +
  geom_boxplot() +
  scale_fill_manual(values = c("green4", "darkorange", "dodgerblue3")) +
  theme_test() +
  theme(axis.text = element_text(size = 12, face = "bold"), 
                    axis.title = element_text(size = 12, face = "bold")) +
  labs(y = "Abundance"))


# ANOVA
anova_result <- aov(Total_Count ~ Habitat, data = site_summary)
summary(anova_result)

# ANOVA results
resid_anova <- residuals(anova_result)
shapiro.test(resid_anova)

bartlett.test(Count~Habitat, data = birds_long)

# Running Tukey's HSD post-hoc test on the anova output and setting the confidence level to 95%
(bird_test <- TukeyHSD(anova_result, conf.level=.95))

# To convert results into a better presented format of the summary table you can use the broom package
(tukey_results <- broom::tidy(bird_test)) # Creating tidy data frame using broom

# Plotting Tukey's test result
plot(bird_test)

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

# Assuming your data is called birds_long (Species, Site, Count, Habitat)
# Convert to wide format: rows = sites, cols = species
bird_matrix <- birds_long %>%
  select(Site, Species, Count) %>%
  pivot_wider(names_from = Species, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Site")

nmds_result <- metaMDS(bird_matrix, distance = "bray", k = 2, trymax = 100)
nmds_result$stress


# Keep habitat info separate
habitats <- birds_long %>%
  distinct(Site, Habitat) %>%
  arrange(match(Site, rownames(bird_matrix))) %>%
  pull(Habitat)


bray_dist <- vegdist(bird_matrix, method = "bray")


# Run PERMANOVA
permanova <- adonis2(bray_dist ~ habitats)
print(permanova)


# Check if group dispersions differ
disp <- betadisper(bray_dist, habitats)
anova(disp)

# Optional: plot dispersion ellipses
plot(disp)


# bird_matrix: rows = sites, columns = species
# habitats: vector of group labels corresponding to rows of bird_matrix

simper_result <- simper(bird_matrix, group = habitats, permutations = 999)

# Overview of SIMPER results
summary(simper_result) # Get SIMPER for a specific comparison, e.g., Plantation vs Heathland
simper_result[["Plantation_Heathland"]]

# Sort by contribution
top_species <- simper_result[["Plantation_Heathland"]] %>%
  as.data.frame() %>%
  arrange(desc(average))  # 'avg' = average contribution

head(top_species, 10)  # top 10 contributing species



# site_scores already built in your session:
 site_scores <- scores(nmds_result, display = "sites") %>% as.data.frame() %>% rownames_to_column("Site")
 site_scores$Habitat <- habitats
 cols_points <- c("Plantation"="dodgerblue3","Heathland"="darkorange","Caledonian"="green4")
 cols_fill   <- c("Plantation"="dodgerblue3","Heathland"="darkorange","Caledonian"="green4")

# Ensure factor levels match your colour vectors
site_scores$Habitat <- factor(site_scores$Habitat, levels = names(cols_fill))

# (Optional) require at least 3 points per group to avoid singular covariances
valid_groups <- names(which(table(site_scores$Habitat) >= 3))

# Compute SE ellipses (no drawing)
ell <- ordiellipse(
  nmds_result,
  groups      = site_scores$Habitat,
  show.groups = if (length(valid_groups)) valid_groups else levels(site_scores$Habitat),
  kind        = "se",
  conf        = 0.95,
  draw        = "none"
)

# Convert ellipse params to coordinates (robust to near-singular covariances)
ellipse_df <- imap_dfr(ell, function(e, g) {
  if (is.null(e$cov) || anyNA(e$cov)) return(NULL)
  theta  <- seq(0, 2*pi, length.out = 200)
  circle <- cbind(cos(theta), sin(theta))
  ev     <- eigen(e$cov, symmetric = TRUE)
  axes   <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), 2, 2)
  coords <- t(apply(circle, 1, function(v) e$center + e$scale * (axes %*% v)))
  tibble(NMDS1 = coords[,1], NMDS2 = coords[,2], Habitat = g)
})


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
  theme(legend.position = "bottom") +
  labs(x = "NMDS1", y = "NMDS2")



  
# Vegetation Data ----

# Choose own path
veg <- read_csv("your_path")
View(veg)


# The first row contains vegetation type -> transpose
vegtype <- as.factor(as.character(unlist(veg[1, ])))
comm <- veg[-1, ]   # remove that row to keep community data only
comm <- apply(comm, 2, as.numeric) # convert to numeric

# PERMANOVA (using Bray-Curtis dissimilarity)
adonis2(comm ~ vegtype, method = "bray")

# NMDS ordination
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)

# Plot NMDS with vegetation types
ordiplot(nmds, type = "n")
points(nmds, display = "sites", col = vegtype, pch = 19)
legend("topright", legend = levels(vegtype), col = 1:length(levels(vegtype)), pch = 19)



# Community matrix: rows = sites, cols = species
comm <- veg[, -1]
comm <- as.data.frame(lapply(comm, as.numeric))
rownames(comm) <- veg$Site

# Define vegetation type factor by site ranges
vegtype <- factor(c(rep("Plantation", 3),
                    rep("Heathland", 3),
                    rep("Caledonian", 3)))

# Check alignment
nrow(comm)        # should be 9
length(vegtype)   # should also be 9

# PERMANOVA
adonis2(comm ~ vegtype, method = "bray")

# NMDS
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)

# Plot with colors by vegtype
ordiplot(nmds, type = "n")
points(nmds, display = "sites", col = vegtype, pch = 19, cex = 1.5)
legend("topright", legend = levels(vegtype),
       col = 1:length(levels(vegtype)), pch = 19)




library(vegan)
library(tibble)
library(dplyr)
library(purrr)
library(ggplot2)

# Run NMDS
nmds_result <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)

# Extract site scores
site_scores <- scores(nmds_result, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Site")

# Define habitats manually to match row order of comm
habitats <- factor(c(rep("Plantation", 3),
                     rep("Heathland", 3),
                     rep("Caledonian", 3)))

# Check alignment
length(habitats)  # should equal nrow(comm)
rownames(comm)    # make sure site order matches expectation


# Add habitat info (factor with correct order)
site_scores$Habitat <- factor(habitats,
                              levels = c("Plantation", "Heathland", "Caledonian"))

# Define colors
cols_points <- c("Plantation"="dodgerblue3",
                 "Heathland"="darkorange",
                 "Caledonian"="green4")

cols_fill   <- cols_points  # same for ellipses

# Compute ellipses
valid_groups <- names(which(table(site_scores$Habitat) >= 3))

ell <- ordiellipse(
  nmds_result,
  groups      = site_scores$Habitat,
  show.groups = if (length(valid_groups)) valid_groups else levels(site_scores$Habitat),
  kind        = "se",
  conf        = 0.95,
  draw        = "none"
)

ellipse_df <- imap_dfr(ell, function(e, g) {
  if (is.null(e$cov) || anyNA(e$cov)) return(NULL)
  theta  <- seq(0, 2*pi, length.out = 200)
  circle <- cbind(cos(theta), sin(theta))
  ev     <- eigen(e$cov, symmetric = TRUE)
  axes   <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), 2, 2)
  coords <- t(apply(circle, 1, function(v) e$center + e$scale * (axes %*% v)))
  tibble(NMDS1 = coords[,1], NMDS2 = coords[,2], Habitat = g)
})

# Plot
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

