# library(BiocInstaller)
# source("http://www.bioconductor.org/biocLite.R")
# useDevel()
# biocLite("microbiome")

library(microbiome)  

data(atlas1006) 
atlas1006
data(dietswap) # Data from http://dx.doi.org/10.1038/ncomms7342
dietswap
dietswap.compositional <- transform(dietswap, "compositional")
g <- global(atlas1006, index = "gini")
sample_data(atlas1006)$diversity <- global(atlas1006, index = "shannon")[,1]

# Compare age and microbiome diversity
plot_regression(diversity ~ age, meta(atlas1006))
p <- prevalence(dietswap, detection = 0, sort = TRUE)

# Taxa with over 50% prevance at .2% relative abundance
dietswap.core <- core(dietswap.compositional, 
                      detection = .2/100, prevalence = 50/100)
dietswap.core

library(ggplot2, quiet = TRUE)
p <- plot_core(transform(dietswap.core, "compositional"), 
               plot.type = "heatmap", 
               colours = gray(seq(0,1,length=5)),
               prevalences = seq(.05, 1, .05), 
               detections = 10^seq(log10(1e-3), log10(.2), length = 10), 
               horizontal = TRUE) +
  xlab("Detection Threshold (Relative Abundance (%))") 
print(p)  
