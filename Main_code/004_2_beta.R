library("lme4")
library("hrbrthemes")
library("gcookbook")
library("tidyverse")
library("devtools")
library("phytools")
library("phyloseq")
library("structSSI")
library("microbiomeViz");
library("gtools")
##################################################
#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
# mainDir <- "/home/csosa"
# chapter <- "chapter_1"
mainDir <- "E:/DADA2"
chapter <- "chapter_1"
levels <- c("kingdom","phylum","class","order","family","genus","species")
#Defining workspace folders
dat_dir <- paste0(mainDir,"/",chapter,"/","data"); if(!file.exists(dat_dir)){dir.create(dat_dir)}
production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)} # Copy here the production data
taxadb_dir <- paste0(mainDir,"/",chapter,"/","db"); if(!file.exists(taxadb_dir)){dir.create(taxadb_dir)}
#production <- read.csv(paste0(production_dir,"/","production_data.csv"),header=T,sep="|")
out_dir <- paste0(mainDir,"/",chapter,"/","outcomes"); if(!file.exists(out_dir)){dir.create(out_dir)}
graph_dir <- paste0(out_dir,"/","graphics"); if(!file.exists(graph_dir)){dir.create(graph_dir)}
csv_dir <- paste0(out_dir,"/","csv"); if(!file.exists(csv_dir)){dir.create(csv_dir)}
seq_dir <- paste0(dat_dir,"/","sequences"); if(!file.exists(seq_dir)){dir.create(seq_dir)}
seq_raw_dir <- paste0(seq_dir,"/","raw"); if(!file.exists(seq_raw_dir)){dir.create(seq_raw_dir)} #Copy here raw reads

seq_proc_dir <- paste0(seq_dir,"/","proc"); if(!file.exists(seq_proc_dir)){dir.create(seq_proc_dir)}
seq_proc_dir_p <- paste0(seq_proc_dir,"/","paired"); if(!file.exists(seq_proc_dir_p)){dir.create(seq_proc_dir_p)}
seq_proc_dir_up <- paste0(seq_proc_dir,"/","unpaired"); if(!file.exists(seq_proc_dir_up)){dir.create(seq_proc_dir_up)}
seq_filt_dir <- paste0(seq_dir,"/","filter"); if(!file.exists(seq_filt_dir)){dir.create(seq_filt_dir)}
qual_pl_dir <- paste0(graph_dir,"/","qual_plot"); if(!file.exists(qual_pl_dir)){dir.create(qual_pl_dir)}
program_dir <- paste0(mainDir,"/","software")
script_dir <-  paste0(mainDir,"/","scripts")#Load de novo cleaned sequences
production_dir <- paste0(dat_dir,"/","production"); if(!file.exists(production_dir)){dir.create(production_dir)}



ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))

ps3 <- readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))



# functions #
st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

tree = phy_tree(ps)
tree <- midpoint.root(tree)
ps@phy_tree <- tree
##### 1) Alpha diversity analyses #####
alpha.1<-as.data.frame(estimate_richness(ps,measures=c("Observed","Shannon","InvSimpson")))
alpha.1.sd<-as(sample_data(ps),'data.frame')
alpha.2<-cbind(alpha.1.sd,alpha.1)
alpha.2$E2<-alpha.2$InvSimpson/alpha.2$Observed
alpha.2$log.reads<-log(alpha.2$UsableReads)
a2<-droplevels((alpha.2))

par(mfrow=c(1, 1))
boxplot(a2$Shannon~a2$Treatment+a2$Phase)


dist.bray <- phyloseq::distance(ps3, method="bray")
dist.uu <- phyloseq::distance(ps3, method="unifrac")
dist.wu <- phyloseq::distance(ps3, method="wunifrac")

# PERMANOVA #
md<-data.frame(sample_data(ps3))
md$UsableReads<-log(md$UsableReads)
adonis(dist.wu ~ Treatment*Phase, data = md)

tree$tip.label


a2_HE <- a2[which(a2$Treatment=="High feed efficiency"),]
alpha_df <- data.frame(
  groups = aggregate(Observed ~ Treatment, mean,data=a2)[,1],
  observation = aggregate(Observed ~ Treatment, mean,data=a2)[,2],
  observation_st = aggregate(  Observed ~ Treatment, st.err,data=a2)[,2],
  log.reads = aggregate(  log.reads ~ Treatment, mean,data=a2)[,2],
  log.reads_st = aggregate(  log.reads ~ Treatment, st.err,data=a2)[,2],
  shannon = aggregate(  Shannon ~ Treatment, mean,data=a2)[,2],
  shannon_st = aggregate(  Shannon ~ Treatment, st.err,data=a2)[,2],
  InvSimpson = aggregate(  InvSimpson ~ Treatment, mean,data=a2)[,2],
  InvSimpson_st = aggregate(  InvSimpson ~ Treatment, st.err,data=a2)[,2],
  E2 = aggregate(  E2 ~ Treatment, mean,data=a2)[,2],
  E2_st = aggregate(  E2 ~ Treatment, st.err,data=a2)[,2],
  FCR = aggregate(  FCR ~ Treatment, mean,data=a2)[,2],
  FCR_st = aggregate(  FCR ~ Treatment, st.err,data=a2)[,2],
  pH = aggregate(  pH ~ Treatment, mean,data=a2)[,2],
  pH_st = aggregate(  pH ~ Treatment, st.err,data=a2)[,2]
  
  
)

alpha_df

#phylosig(tree = tree, x = merge(tree$tip.label,ps@sam_data$FCR,"row.names")$y, method="K",  test=TRUE, nsim=999)

# Transform to compositional abundances
pseq.rel <- microbiome::transform(ps, "compositional")

pseq.core <- core(pseq.rel, detection = 0.1/100, prevalence = 10/100)


# Core with compositionals:
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-3), log10(.2), length = 10)


p <- plot_core(pseq.rel, plot.type = "heatmap",
               prevalences = prevalences, detections = detections) +
  xlab("Detection Threshold (Relative Abundance)")
print(p)  


tab <- diversities(ps, index = "all")
tab <- dominance(ps, index = "all")

tab <- coverage(ps, threshold = 0.5)
tab <- inequality(ps)


##############
taxic <- as.data.frame(ps@tax_table) 
otu.df <- abundances(ps)
# make a dataframe for OTU information.
otu.df <- as.data.frame(otu.df)

# check the rows and columns
# head(otu.df) 

# Add the OTU ids from OTU table into the taxa table at the end.
taxic$OTU <- row.names(otu.df)  

taxmat <- as.matrix(taxic)
new.tax <- tax_table(taxmat)  
tax_table(pseq.rel) <- new.tax 
pseq.fam <- aggregate_taxa(pseq.rel, "Family")
pseq.famrel <- transform(pseq.fam, "compositional")
# 
# p.famrel <- plot_composition(pseq.famrel, sample.sort = NULL, otu.sort = NULL,
#                              x.label = "ibd_subtype", plot.type = "barplot", verbose = FALSE)
ps@otu_table <- t(ps@otu_table)

plot_taxa_prevalence(ps, "Phylum")

ps_HE <- subset_samples(ps, Treatment=="High feed efficiency")
ps_LE <- subset_samples(ps, Treatment=="Low feed efficiency")
par(mfrow=c(1, 2))
plot_taxa_prevalence(ps_HE, "Phylum")
plot_taxa_prevalence(ps_LE, "Phylum")

require(phyloseqGraphTest)

gt <- graph_perm_test(physeq = ps3, sampletype = "Treatment",
                      distance = "jaccard", type = "mst",nperm=100)

plotNet2=plot_test_network(gt) + theme(legend.text = element_text(size = 8),
                                       legend.title = element_text(size = 9))
plotPerm2=plot_permutations(gt)
grid.arrange(ncol = 2,  plotNet2, plotPerm2)



library(PMA)
X <- otu_table(ps3)
#metab <- log(1 + metab, base = 10)
cca_res <- CCA(X,  ps3@sam_data[,c(5:9)], penaltyx = .15, penaltyz = .15)
cca_res

out.dpcoa.log <- ordinate(ps3, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(ps3, out.dpcoa.log, color = "Treatment", label= "SampleID",
                shape = "Phase") +
  labs(col = "Treatment", shape = "Phase")+
  coord_fixed(sqrt(evals[2] / evals[1]))


plot_ordination(ps3, out.dpcoa.log, type = "species", color = "Kingdom",shape="Treatment") +
  coord_fixed(sqrt(evals[2] / evals[1]))


abund <- otu_table(ps)
abund_ranks <- t(apply(abund, 1, rank))

abund_ranks <- abund_ranks - 10
abund_ranks[abund_ranks < 1] <- 1
library(dplyr)
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")


ps_ccpna <- ordinate(ps3, "CCA", formula = ps3 ~ Treatment + Phase)
plot_ordination

tax <- tax_table(ps3)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
tax$otu_id <- row.names(tax)
ss <- sample_data(ps3)
ss$Subject <- as.character(ss$Subject)
library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$Subject <- rownames(sites)
sites <- sites %>%
  left_join(ss)

species <- data.frame(ps_scores$species)
species$otu_id <- colnames(t(otu_table(ps)))
species <- species %>%
  left_join(tax)


evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -10),
                  aes(x = CCA1, y = CCA2, label = Phylum),
                  size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ Treatment) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))

##########################
###########################


library(caret)
dataMatrix <- data.frame(Treatment = sample_data(ps3)$Treatment, otu_table(ps3))
# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(ps3)$Subject), size = 8)
inTrain <- which(sample_data(ps3)$Subject %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]
plsFit <- train(Treatment ~ ., data = training,
                method = "pls", preProc = "center")

plsClasses <- predict(plsFit, newdata = testing)
table(plsClasses, testing$Treatment)

library(randomForest)
rfFit <- train(Treatment ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$Treatment)

####

pls_biplot <- list("loadings" = loadings(plsFit$finalModel),
                   "scores" = scores(plsFit$finalModel))
class(pls_biplot$scores) <- "matrix"

pls_biplot$scores <- data.frame(sample_data(ps3)[inTrain, ],
                                pls_biplot$scores)

tax <- tax_table(ps3)@.Data %>%
  data.frame(stringsAsFactors = FALSE)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
class(pls_biplot$loadings) <- "matrix"
pls_biplot$loadings <- data.frame(tax, pls_biplot$loadings)


ggplot() +
  geom_point(data = pls_biplot$scores,
             aes(x = Comp.1, y = Comp.2), shape = 2) +
  geom_point(data = pls_biplot$loadings,
             aes(x = 25 * Comp.1, y = 25 * Comp.2, col = Order),
             size = 0.3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Axis1", y = "Axis2", col = "Binned Age") +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  facet_grid( ~ Treatment) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
#############


rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(ps3)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = Treatment,shape=Phase),
             size = 6, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")
