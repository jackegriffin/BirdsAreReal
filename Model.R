#Originally we had planned to build a mixed effect model using our vegetational data as random effects, however our data isn't categorical so this wouldn't work. We also didn't have enough data points to add anymore fixed effects.



#general linear poisson model, of effect of tree density on species richness
glm_trees <- glm(species_richness ~ tree_density, family = poisson, data = bird_datavegscore)
summary(glm_trees)
library(AER)
dispersiontest(glm_trees)

#adding shannon as a factor, shows tree density is a stronger driver of bird species richness than vegetation diversity. However neither are significant and if being critical, two explanatory variables are far too much to fit in due to our sample size. 
glm_trees2 <- glm(species_richness ~ tree_density + shannon, family = poisson, data = bird_datavegscore)
summary(glm_trees2)
dispersiontest(glm_trees2)

#both passed dispersion tests
