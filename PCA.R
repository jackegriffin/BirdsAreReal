#input data
bird_data <- read_csv("Data/bird_data.csv")
View(bird_data)

# make a vegetation matrix
veg_matrix <- bird_data[, c("tree_density","g","h","f","dw","bg","ml")]
veg_dist <- dist(veg_matrix, method = "euclidean")

# Test for differences among habitats
install.packages("vegan")
library(vegan)

# Habitat must be a factor
bird_data$habitat <- factor(bird_data$habitat)

adonis2(veg_dist ~ habitat, data = bird_data)

# Run a PCA 

veg_matrix <- as.data.frame(lapply(veg_matrix, as.numeric))

pca_res <- prcomp(veg_matrix, scale. = TRUE)
summary(pca_res)

# Extract PC1 for each site
bird_data$veg_score <- pca_res$x[,1]

# get site scores (PC1, PC2, etc.)
scores <- as.data.frame(pca_res$x)

# add back habitat and site info
scores$habitat <- bird_data$habitat
scores$site_number <- bird_data$site_number

# look at first few rows
head(scores)

# extract variable loadings (arrows for variables)
loadings <- as.data.frame(pca_res$rotation[,1:2])
loadings
# positive or negative values show the direction of influence
# the bigger the number, the stringer the contribution to that principle component

# make a PCA biplor
# this shows sites and variables together
# sites cluster based on simlarity
# arrows show which variables explain the separation

library(ggplot2)

# Scale loadings so arrows are visible
loadings_scaled <- loadings
loadings_scaled[,1:2] <- loadings_scaled[,1:2] * 3
loadings_scaled$variable <- rownames(loadings)

ggplot(scores, aes(x = PC1, y = PC2, color = habitat)) +
  geom_point(size = 3) +
  geom_segment(data = loadings_scaled,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = loadings_scaled,
            aes(x = PC1, y = PC2, label = variable),
            vjust = 1.5, color = "black") +
  theme_minimal() +
  labs(title = "PCA Biplot: Vegetation Structure")

# test if habitats differ in PC1
# now we ask: are habitats significantly different along the main vegetation gradient (PC1)?

anova_res <- aov(PC1 ~ habitat, data = scores)
summary(anova_res)

# visualize PC1 by habitat
ggplot(scores, aes(x = habitat, y = PC1, fill = habitat)) +
  geom_boxplot() +
  theme_minimal() +
  labs(y = "PC1 Vegetation Score", x = "Habitat")

# which variable is driving the change

