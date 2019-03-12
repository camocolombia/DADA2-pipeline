
##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##ASSIGNING OTUs                                 #
##################################################

##################################################
#Calling libraries
#source("https://bioconductor.org/biocLite.R")
#biocLite("microbiomeSeq")
library("ggplot2")
library("gridExtra")
library("reshape")
library("plyr")
library("phyloseq")
library("microbiome")
#library("microbiomeSeq")
library("vegan")
library("ggvegan")
library("dplyr")
library("ggrepel")
library("ggtree")
library("data.table")
library("DESeq2")
library("BSDA")
library("adegenet")
#install.packages("devtools")
# library("devtools")
# library("phytools")
# library("massMap")
# library("structSSI")
# library("microbiomeViz");
# library("stringr")
# library("StructFDR")
library("gtools")
#install_github("JiyuanHu/massMap")

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


##################################################
##################################################
#Defining functions to be use in the FDR approaches
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#permutation function
perm.func <- function (X, Y) {
  return(list(X=X, Y=sample(Y)))
}



#T.TEST
#   test.func <- function (X, Y) {
#     obj <- apply(X, 1, function(x) {
#       ttest.obj <- t.test(x ~ Y)#wilcox.test(x ~ Y,paired=F)
#       c(ttest.obj$p.value, sign(ttest.obj$statistic))
#     })
#     return(list(p.value=obj[1, ], e.sign=obj[2, ]))
#   }


#Wilcoxon 

test.func <- function (X, Y) {  
  Y <- as.numeric(factor(Y))
  obj <- apply(X, 1, function(x) {                
    p.value <- suppressWarnings(wilcox.test(x ~ Y)$p.value)
    e.sign <- sign(mean(x[Y == 2]) - mean(x[Y == 1]))
    c(p.value, e.sign)          
  })
  return(list(p.value=obj[1, ], e.sign=obj[2, ])) 
}

##################################################
##################################################
#Loading Phyloseq object
ps <- readRDS(paste0(csv_dir,"/","Phyloseq_object",".RDS"))

levels(ps@sam_data$Treatment) <- c("Low feed efficiency","High feed efficiency")
sample_data(ps)$UsableReads<-sample_sums(ps)


##################################################
##################################################
#Exploring differences in libraries per group

HL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="liquid"))))
HS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="solid"))))
LL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="liquid"))))
LS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="solid"))))





lib_df <- data.frame(samples = c(names(HL),names(HS),names(LL),names(LS)),
                        sampleType = c(rep("liquid",length(HL)),rep("solid",length(HS)),rep("liquid",length(LL)),rep("solid",length(LS))),
                        Treatment = c(rep("High feed efficiency",length(HL)),rep("High feed efficiency",length(HS)),rep("Low feed efficiency",length(LL)),rep("Low feed efficiency",length(LS))),  
                          libSize = as.numeric(c(colSums(HL),colSums(HS),colSums(LL),colSums(LS)))
                        )
lib_df$libSize_rel <- NA; lib_df$libSize_rel <- lib_df$libSize/sum(lib_df$libSize)

#lib__df_df2 <- lib_df[which(lib_df$libSize>=100000),]

####


#ggsave(paste0(graph_dir,"/","Production_BP_scaled",".pdf"),p1bp,dpi=300,width =90,height=80,units = "cm",scale=1.2,limitsize = FALSE)
#p1bp
##################################################
##################################################
#Prevalence cleaning up approach



#ps = 
ps <- subset_taxa(ps, Phylum!=is.na(Phylum))

prev0 = apply(X = otu_table(ps),
              MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 5)]
prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
# Execute prevalence filter, using `prune_taxa()` function
ps1 = prune_taxa((prev0 > prevalenceThreshold), ps)
ps1
# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps1, Phylum %in% names(keepPhyla))
saveRDS(ps2,paste0(csv_dir,"/","Phyloseq_object_filter",".RDS"))
ps <- ps2
############after prevalence treatment


HL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="liquid"))))
HS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="High feed efficiency" & Phase=="solid"))))
LL <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="liquid"))))
LS <-as.data.frame(t(otu_table(subset_samples(ps, Treatment=="Low feed efficiency" & Phase=="solid"))))

lib_df$libSize_T <- NA;lib_df$libSize_T <- as.numeric(c(colSums(HL),colSums(HS),colSums(LL),colSums(LS)))
lib_df$libSize_T_rel <- NA; lib_df$libSize_T_rel <- lib_df$libSize_T/sum(lib_df$libSize_T)

rm(HL,HS,LL,LS)
lib_df2 <- melt(lib_df)
lib_df2 <- lib_df2[which(lib_df2$variable=="libSize" | lib_df2$variable=="libSize_T"),]
lib_df2$variable <- gsub("libSize","Before filter",lib_df2$variable)
lib_df2$variable <- gsub("Before filter_T","After filter",lib_df2$variable)

lib_df2$variable[which(lib_df2$variable=="libSize")] <- "Before filter"
lib_df2$variable[which(lib_df2$variable=="libSize_T")] <- "After filter"
lib_df2$Type <- NA;lib_df2$Type <- paste0(lib_df2$Treatment,"-",lib_df2$sampleType)

lib_size <- ggplot(data=lib_df2,aes(x=Type,y=value,fill=variable))+
  geom_boxplot(position = "dodge",outlier.colour = NA,outlier.fill=NA,outlier.alpha=1,outlier.size =NA,na.rm = T) +
  stat_boxplot(geom ='errorbar') +
  #stat_summary(fun.y=mean, geom="line", aes(group=1))  +
  #stat_summary(fun.y=mean, geom="point")+
  #guides(fill=FALSE)+
  #ylim(c(0,100))+
  ylab("Scaled values")+ xlab("")+
  ggtitle("")+
    scale_fill_manual("", values = c("blue","lightblue","orange","red"))+
  theme(panel.background = element_rect(fill = "gray95"),
        text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle = 45, hjust = 1),
        axis.text.y  = element_text(size=60,colour="black")+
          scale_y_continuous(limits = c(50000,140000),breaks=seq(50000,140000,5000))
        
          #,
        #legend.title=element_blank(),legend.position="none"
  )+ 
  coord_flip()


lib_df$libSize - lib_df$libSize_T
lib_size


ps3 <- transform_sample_counts(ps,function(x) {x/sum(x)} )
saveRDS(ps3,paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))





# p <- plot_density(ps3, variable='Treatment')
# p


ps_kingdom <- aggregate_taxa(ps, 'Kingdom');ps_kingdom2 <- transform_sample_counts(ps,function(x) {x/sum(x)} ) 
ps_phylum <- aggregate_taxa(ps, 'Phylum');ps_phylum2 <- transform_sample_counts(ps,function(x) {x/sum(x)} ) 
ps_class <- aggregate_taxa(ps, 'Class');ps_class2 <- transform_sample_counts(ps,function(x) {x/sum(x)} ) 
ps_order <- aggregate_taxa(ps, 'Order');ps_order2 <- transform_sample_counts(ps,function(x) {x/sum(x)} ) 
ps_family <- aggregate_taxa(ps, 'Family');ps_family <- transform_sample_counts(ps,function(x) {x/sum(x)} ) 
ps_genus <- aggregate_taxa(ps, 'Genus');ps_genus2 <- transform_sample_counts(ps,function(x) {x/sum(x)} ) 

dist_otu = phyloseq::distance(ps3, method="bray")

ps_list <-list(ps_kingdom2,ps_phylum2,ps_class2,ps_order2,ps_family,ps_genus2)

mantel_levels <- lapply(1:length(ps_list),function(i){
cat(i,"\n")  
dist1 = phyloseq::distance(subset_samples(ps_list[[i]], Treatment=="High feed efficiency" & Phase=="solid"), method="bray")
dist2 = phyloseq::distance(subset_samples(ps_list[[i]], Treatment=="High feed efficiency" & Phase=="liquid"), method="bray")
dist3 = phyloseq::distance(subset_samples(ps_list[[i]], Treatment=="Low feed efficiency" & Phase=="solid"), method="bray")
dist4 = phyloseq::distance(subset_samples(ps_list[[i]], Treatment=="Low feed efficiency" & Phase=="liquid"), method="bray")
dist_H = phyloseq::distance(subset_samples(ps_list[[i]], Treatment=="High feed efficiency"), method="bray")
dist_L = phyloseq::distance(subset_samples(ps_list[[i]], Treatment=="Low feed efficiency"), method="bray")

d1_2 <- vegan::mantel(dist1, dist2, permutations=10000)
d1_3 <- vegan::mantel(dist1, dist3, permutations=10000)
d1_4 <- vegan::mantel(dist1, dist4, permutations=10000)
d2_3 <- vegan::mantel(dist2, dist3, permutations=10000)
d2_4 <- vegan::mantel(dist2, dist4, permutations=10000)
d3_4 <- vegan::mantel(dist3, dist4, permutations=10000)
dH_L <- vegan::mantel(dist_H, dist_L, permutations=10000)

x <- data.frame(level   =  levels[[i]],
                HS_HL_R =  d1_2$statistic,
                HS_HL_P =  d1_2$signif,
                HS_LS_R =  d1_3$statistic,
                HS_LS_P =  d1_3$signif,
                HS_LL_R =  d1_4$statistic,
                HS_LL_P =  d1_4$signif,
                HL_LS_R =  d2_3$statistic,
                HL_LS_P =  d2_3$signif,
                LS_LL_R =  d3_4$statistic,
                LS_LL_P =  d3_4$signif,
                H_L_R   =  dH_L$statistic,
                H_L_P   =  dH_L$signif
)

return(x)
cat("   ","\n");cat("DONE!","\n");cat("   ","\n")
})
mantel_levels <- do.call(rbind,mantel_levels)
mantel_levels2 <- mantel_levels

cols <- c(3,5,7,9,11,13); for(i in 1:length(cols)){
  p_i <- as.numeric(cols[i])
mantel_levels2[,p_i][which(mantel_levels2[,p_i] >0.05)] <- NA};rm(i)

# PERMANOVA #
md<-data.frame(sample_data(ps3))
md$UsableReads<-log(md$UsableReads)
adonis(dist_otu ~ Treatment*Phase*pH+UsableReads+ADG, data = md)



#Sort the Phyla by abundance and pick the top 5

top20P.names = sort(tapply(taxa_sums(ps_phylum), tax_table(ps_phylum)[, "Phylum"], sum), TRUE)[1:10]
#Cut down the physeq.tree data to only the top 10 Phyla
top5P = subset_taxa(ps_phylum, Phylum %in% names(top20P.names))
#Plot
pp <- plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = (~Kingdom)) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

bar_to_merge <- merge(taxa_sums(top5P),tax_table(top5P),"row.names") #(taxa_sums(top5P)/sum(taxa_sums(top5P))*100)

top_rel <- ggplot(data=bar_to_merge,aes(x=reorder(Family, x),y=x,fill="red"))+
  geom_bar(stat="identity",show.legend = F)+
  xlab("")+
  ylab("Relative abundance (%)")+
  # ylim(0,0.1)+
  geom_hline(yintercept=0.05, linetype="dashed", color = "black")+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60),
        axis.text.x  = element_text(size=60,colour="black",angle=90),
        axis.text.y  = element_text(size=60,colour="black"))+
  coord_flip()#+
  #scale_y_continuous(limits = c(0,0.8),breaks=seq(0,0.8,0.05))
top_rel


#mantel_levels[,2][which()]

HL.correlog <- mantel.correlog(dist_H, dist_L,nperm=1000,mult="BH",n.class=50)
summary(HL.correlog)
HL.correlog   
plot(HL.correlog)


require(corrplot)
corrplot(as.matrix(dist1))

comp_list <- list(
  # SH = (df$SH/(sum(df$SH,df$SL))),
  # SL = (df$SL/(sum(df$SH,df$SL))),
  # LH = (df$LH/(sum(df$LH,df$LL))),
  # LL = (df$LL/(sum(df$LH,df$LL)))
  SH = rowSums(otu_table(subset_samples(s, Treatment=="High feed efficiency" & Phase=="solid")))/
    sum(taxa_sums(subset_samples(s, Treatment=="High feed efficiency"))),
  LH = rowSums(otu_table(subset_samples(s, Treatment=="High feed efficiency" & Phase=="liquid")))/
    sum(taxa_sums(subset_samples(s, Treatment=="High feed efficiency"))),
  
  SL = rowSums(otu_table(subset_samples(s, Treatment=="Low feed efficiency" & Phase=="solid")))/
    sum(taxa_sums(subset_samples(s, Treatment=="Low feed efficiency"))),
  LL = rowSums(otu_table(subset_samples(s, Treatment=="Low feed efficiency" & Phase=="liquid")))/
    sum(taxa_sums(subset_samples(s, Treatment=="Low feed efficiency")))
  
  )

comp_list  <- melt(comp_list)
#comp_list$value_rel  <- NA; comp_list$value_rel  <-  comp_list$value/sum(comp_list$value)

ggplot(comp_list, aes(value, fill = L1,colour = L1)) +
  geom_freqpoly(bins=400)#+
#  scale_x_continuous(limits = c(0,0.0001),breaks=seq(0,0.0001,0.00001))

#kruskal.test(L1 ~ value,data=comp_list)
kruskal(comp_list$value, comp_list$L1, console = TRUE)
out<-with(comp_list,kruskal(y=value,trt=L1,mian="comp_list",group=F, p.adj="BH"))
out
plot_bar(ps,x="Treatment", fill="Phylum")

##################################################
##################################################
#####DESEQ2

#Converting to DESeq2object
diagdds = phyloseq_to_deseq2(ps, ~ Treatment)
#Applying geometric means
geoMeans = apply(counts(diagdds), 1, gm_mean)
#Estimating size factors
diagdds = estimateSizeFactors(diagdds)#, geoMeans = geoMeans)
diagdds <- estimateDispersions(diagdds)
#varianceStabilizingTransformation(diagdds, blind = TRUE, fitType = "parametric")
abundds2 <- getVarianceStabilizedData(diagdds)
diagdds2 <- DESeq2::DESeq(diagdds,test="Wald",fitType="local") 
res2 <- as.data.frame(results(diagdds2)); 
res2 <- res2[,c("baseMean", "log2FoldChange", "lfcSE", "padj")]

sigtab = res2[(res2$padj <= 0.05), ]; sigtab <- sigtab[complete.cases(sigtab),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
sigtabgen = subset(sigtab, !is.na(Genus))
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ss_Simp <- simper(otu_table(ps3), ps@sam_data$Treatment, permutations=100)




#StructFDR
#https://cran.r-project.org/web/packages/StructFDR/vignettes/StructFDR.html





ps_family = tax_glom(ps, taxrank = "Family")
#sps_family = tax_glom(ps, taxrank = "Family")
ntaxa(ps_family); ntaxa(ps)
colnames(tax_table(ps_family))

# Phylogenetic agglomeration
# Phylogenetic agglomeration



















#Prevalence
ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = "black")) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7,show.legend = F) +
  scale_y_log10() + scale_x_log10() +
  xlab("Total Abundance") +
  facet_wrap(~Phylum)+ 
  scale_color_manual("Phylum", values = c("black"))




# short_names <- substr(rownames(abundds), 1, 5)%>%
#   make.names(unique = TRUE)
#rownames(abundds) <- short_names
#diagdds = DESeq(diagdds,test="Wald",fitType="local") 


abund_sums <- rbind(data.frame(sum = colSums(abundds2),
                               sample = colnames(abundds2),
                               type = "DESeq2"),
                    data.frame(sum = rowSums(otu_table(ps3)),
                               sample = rownames(otu_table(ps3)),
                               type = "log(1 + x)"))
ggplot(abund_sums) +
  geom_histogram(aes(x = sum), binwidth = 50) +
  facet_grid(type ~ .) +
  xlab("Total abundance within sample")
##
#plot(hfdr_res, height = 5000) # opens in a browser

#row.names(tax2) <-row.names(tax)

tax2 <- tax2[which(tax2$adj.significance=="**" | tax2$adj.significance=="***" ),]

tax2 <- merge(tax2,res2,by.x="seq",by.y="row.names")
# data("GlobalPatterns")
# GP = GlobalPatterns
otu_table = t(otu_table(ps3))

ps3@otu_table <- otu_table

# # # Use the OTUs that make up 99% of the total
# lib_size = colSums(otu_table)
# mat = sapply(1:ncol(otu_table), function(i)
#   otu_table[,i]/lib_size[i] >= 0.01)
# ind = rowSums(mat)>=1
#
ps3 = fix_duplicate_tax(ps3)#

tr = parsePhyloseq(ps3)

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)
sm_fam <- table(ps2@tax_table[,5])[order(table(ps2@tax_table[,5]),decreasing = T)]
p <- tree.backbone(tr, size=1)
anno.data <- data.frame(node= c("f__Lachnospiraceae",
                                "g__Butyrivibrio",
                                "f__Prevotellaceae",
                                "g__Prevotella",
                                "f__Ruminococcaceae",
                                "g__Ruminococcaceae_UCG-014",
                                "f__Succinivibrionaceae",
                                "g__Succinivibrionaceae_UCG-001",
                                "f__Synergistaceae",
                                "g__Fretibacterium",
                                "f__Muribaculaceae",
                                "f_Family_XIII",
                                "g__Family_XIII_UCG-001",
                                "f__Veillonellaceae",
                                "f__Methanobacteriaceae"), 
                          
                          #paste0("f__",
                               #      as.character(unique(tax2[which(tax2$adj.significance=="**" & !is.na(tax2$Family)),]$Family))),
                          color=c("red",
                                  "red",
                                  "blue",
                                  "blue",
                                  "orange",
                                  "orange",
                                  "purple",
                                  "purple",
                                  "yellow",
                                  "yellow",
                                  "green",
                                  "brown",
                                  "brown",
                                  "gray",
                                  "gray90"),#mypalette[1:13],
                        stringsAsFactors = FALSE)
p <- clade.anno(p,anno.data = anno.data)
p
#diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

#diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
#diagdds = DESeq(diagdds, fitType="local")
#https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
#http://joey711.github.io/phyloseq-extensions/DESeq2.html
#Saving results

res = results(diagdds2)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
#head(sigtab)

posigtab = sigtab
#posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]

posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]

sigtabgen = subset(sigtab, !is.na(Genus))
sigtabgen = sigtab


# Phylum order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# faimly order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Family, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Family = factor(as.character(sigtabgen$Family), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
# Species order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Species, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Species = factor(as.character(sigtabgen$Species), levels=names(x))



#exp_plot <- ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
exp_plot <- ggplot(sigtabgen, aes(y=log2FoldChange, x=reorder(Family, -log2FoldChange), color=Phylum, fill=Phylum)) + 
   #geom_point(size=6) + 
  geom_bar(stat="identity")+
  geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
  xlab("")+
    #geom_line()+
  theme(text=element_text(size=20),
        legend.text=element_text(face="italic",size=30,colour="black"),
        axis.text.x  = element_text(face="italic",angle = -90,size=25, hjust = 0, vjust=0.5),
        axis.text.y  = element_text(face="italic",size=25,colour="black"))+
  scale_y_continuous(limits = c(-30,5),breaks=seq(-30,5,1))+
coord_flip()
exp_plot  
ggsave(paste0(graph_dir,"/","log2_change",".pdf"),exp_plot,dpi=300,width =70,height=40,units = "cm",scale=0.8,limitsize = FALSE)


library(FSA)
skrus <- kruskal.test(list(df$SH,df$SL,df$LH,df$LL))


map = as.data.frame(sample_data(ps))

#map = as.data.frame(sample_data(phy))
X = NULL
#X = map[,c('Phase')] # a data frame, each element could either numeric or factor.
Y = map[['Treatment']]#a vector instead of a data frame.
Y <- as.numeric((Y))
Y[which(Y==1)]<-0; Y[which(Y==2)]<-1;
obstime <- NULL
otu.tab = otu_table(ps)
tax.tab = tax_table(ps)
tree = phy_tree(ps)
tree <- midpoint.root(tree)

# Not run:
res.Surv = massMap(X=X,Y=Y, obstime = obstime,delta= NULL,otu.tab=otu.tab,tax.tab=tax.tab,
                  tree=tree,outcome.trait="binary",screening.rank = "Phylum",
                  target.rank="Family",alpha=1,n.perm=1e4)
# res.Surv$res.screening ##The group association test at the screening rank
# res.Surv$res.target ##The microbial association test at the target rank
row.names(df) <- df$Row.names
df3 <- merge(df,res.Surv$res.target,by="row.names",all=T)
df3 <- df3[,-c(1,2)]
df3 <- df3[which(df3$status.HBH==1),]


# exp_plot2 <- ggplot(df3, aes(y=change_Treat, x=reorder(Family, -change_Treat), color=Phylum, fill=Phylum)) + 
#   #geom_point(size=6) + 
#   geom_bar(stat="identity")+
#   geom_hline(yintercept = 0.0, color = "black", size = 1.5) +
#   xlab("")+
#   #geom_line()+
#   theme(text=element_text(size=20),
#         legend.text=element_text(face="italic",size=30,colour="black"),
#         axis.text.x  = element_text(face="italic",angle = -90,size=25, hjust = 0, vjust=0.5),
#         axis.text.y  = element_text(face="italic",size=25,colour="black"))+
#   scale_y_continuous(limits = c(-5,5),breaks=seq(-5,5,0.5))+
#   coord_flip()
# exp_plot2  
#ggsave(paste0(graph_dir,"/","log2_change",".pdf"),exp_plot,dpi=300,width =70,height=40,units = "cm",scale=0.8,limitsize = FALSE)
#https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/doc/MicrobiomeWorkshopII.html#sparse-discriminant-analyses-using-a-phylogenetic-tree








#Plotting tree
ps_tree2 <- plot_tree(ps,nodelabf=nodeplotboot(),# nodelabf=nodeplotboot(80,0,3),
                      color="Treatment", 
                      size = "Abundance",
                      #label.tips="taxa_names",
                      justify = "yes please", 
                      ladderize="left",
                      plot.margin=0.4#,
                      #shape="Phase"
)+
  scale_size_continuous(range = c(1, 3)) +
  scale_color_manual("Treatment", values = c("blue","red"))+
  theme(text=element_text(size=60),
        legend.text=element_text(size=60)) #+
  #coord_polar(theta="y")    


ggsave(paste0(graph_dir,"/","ps_tree",".pdf"),ps_tree2,dpi=300,width =100,height=150,units = "cm",scale=0.8,limitsize = FALSE)


#plot_richness(ps, x="FCR", measures=c("Shannon", "Simpson"), color="Phase") + theme_bw()
#Phylum
pserie_phylum <- ps %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Phylum)     
  #Family
  pserie_family <- ps %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Family)     
  
  tsfamily <-table(pserie_family$Family)
  tsfamily <- tsfamily[order(tsfamily,decreasing = T)]
    #Genus
pserie_genus <- ps %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort d
#Species
pserie_sp <- ps %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Species)  
