

##################################################
##Chrystian C. Sosa 2018 Chapter 1               #
##2018-06-29                                     #
##NUIG - Teagasc Athenry                         #
##Testing beta diversity PERMANOVA               #
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
library("missMDA")
# library("adegenet")
# #install.packages("devtools")
library("devtools")
library("phytools")
# library("massMap")
# library("structSSI")
# library("microbiomeViz");
library("stringr")
# library("StructFDR")
# library("gtools")
library("agricolae")
library("car")
library("ggdendro")
library("rgl")

#install_github("JiyuanHu/massMap")

##################################################
#Defining paths to be analyzed
# mainDir <- "E:/Dropbox/Dropbox/Paper_PhD"
# chapter <- "chapter_1
# mainDir <- "/home/csosa"
# chapter <- "chapter_1"
mainDir <- "E:/DADA2"
chapter <- "chapter_2"
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
ps_object2 <- readRDS(paste0(csv_dir,"/","Phyloseq_taxa_list_trams",".RDS"))

#ps@sam_data$Treatment <- factor(ps@sam_data$Treatment,levels = c("High feed efficiency","Low feed efficiency"))
ps3 <-  readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))
metadata <- as(sample_data(ps3), "data.frame")
#GP.chl = prune_samples(ps@sam_data$sampleID==10707,ps)

taxa_df <- data.frame(seqId=row.names(ps3@tax_table),ps3@tax_table)
taxa_df$seqId <- as.character(taxa_df$seqId)

outliers <- read.csv(paste0(seq_dir,"/","Outlier.csv"),header=T)

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#RANDOM FORSEST
# if(!file.exists(paste0(out_dir,"/","csv/","rF_list.RDS"))){
# library(pROC)
# library(ROCR)
# library(caret)
# set.seed(825)
# 
# dataMatrix <- data.frame(Treatment = sample_data(ps3)$Treatment, 
#                          Phase= sample_data(ps3)$Phase,
#                          #pH = sample_data(ps3)$pH,
#                          otu_table(ps3))
# 
# # take 8 mice at random to be the training set, and the remaining 4 the test set
# #ps3@sam_data$Treatment <- factor(ps3@sam_data$Treatment,levels=c("High feed efficiency","Low feed efficiency"))
# testingMice <- sample(unique(sample_data(ps3)$Subject), 
#                       size = ceiling(sum(sample_sums(ps3))/5))
# inTest <- which(sample_data(ps)$Subject %in% testingMice)
# testing <- dataMatrix[inTest,]
# training <- dataMatrix[-inTest,]
# 
# table(training$Treatment)
# table(testing$Treatment)
# 
# lis<- list(testing=testing,
#            training=training)
# 
# #saveRDS(lis,paste0(out_dir,"/","csv/","FDR_list.RDS"))
# 
# control <- trainControl(method='repeatedcv',
#                         number=10,
#                         repeats=3,
#                         search = 'random')
# 
# # specify that the resampling method is 
# # fit_control <- trainControl(## 10-fold CV
# #   method = "cv",
# #   number = 10)
# 
# # plsFit <- train(Treatment ~ ., data = training,
# #                 method = "pls", prePrquireoc = "center")
# # 
# # plsClasses <- predict(plsFit, newdata = testing)
# # #table(plsClasses, testing$Treatment)
# # varImp(plsFit)
# # confusionMatrix(plsClasses, testing$Treatment)
# # 
# 
# rfFit <- train(Treatment ~ ., data = training,
#                method = "rf", preProc = "center",
#                metric="Accuracy",proximity = TRUE,
#                tuneLength  = 50, 
#                trControl = control
# )
# 
#   rfClasses <- predict(rfFit, newdata = testing)
# confusionMatrix(rfClasses, testing$Treatment)
# varRF <- varImp(rfFit)
# importance <- varRF$importance#
# importance$OTU <- row.names(importance)
# importance <- importance[,c(2,1)]
# importance <- importance[order(varRF$importance$Overall,decreasing = TRUE),]
# #importance <- importance[1:25,]
# #importance <- merge(importance,ps3@tax_table,by.x="OTU",by.y="row.names",all=F)
# colnames(taxa_df)[1] <- "OTU"
# importance <- left_join(importance,taxa_df,by="OTU")
# importance <- importance[order(importance$Overall,decreasing = TRUE),]
# 
# #View(importance)
# summary(rfFit)
# plot(rfFit)
# var_plot <- plot(varRF,top=20)
# var_plot
# #https://topepo.github.io/caret/variable-importance.html
# 
# ggplot(rfFit)
# varImp(rfFit)
# 
# 
# 
# rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
#   data.frame(sample_data(ps3)[-inTest, ])
# ggplot(rf_prox) +
#   geom_point(aes(x = X1, y = X2, col = Treatment),
#              size = 4, alpha = 0.6) +
#   scale_color_manual(values = c("red", "blue")) +
#   guides(col = guide_legend(override.aes = list(size = 3))) +
#   labs(col = "Groups", x = "Axis1", y = "Axis2")
# require(randomForest)
# as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])
# require(randomForest)
# impOtu <- as.vector(otu_table(ps3)[,which.max(importance(rfFit$finalModel))])
# maxImpDF <- data.frame(sample_data(ps3), abund = impOtu)
# ggplot(maxImpDF) +
#   geom_histogram(aes(x = abund)) +
#   facet_grid(Treatment ~ .) +
#   labs(x = "Abundance of discriminative microbe", y = "Number of samples")
# #filterVarImp(training,rfClasses)
# 
# lis2 <-  list(dataMatrix = dataMatrix,
#               testing=testing,
#               training=training,
#               rfFit= rfFit,
#               rfClasses = rfClasses,
#               importance = importance,
#               rf_prox = rf_prox)
# 
# saveRDS(lis2,paste0(out_dir,"/","csv/","rF_list.RDS"))
# } else {
#   lis2 <- readRDS(paste0(out_dir,"/","csv/","rF_list.RDS"))
#     
#   }
# 
# rfClasses <- predict(lis2$rfFit, newdata = lis2$testing)
# require(caret)
# confusionMatrix(rfClasses, lis2$testing$Treatment)

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#Beta diversity


 ####BRAY CURTIS

dist_b <- phyloseq::distance(ps3,"bray")
# dist_j <- phyloseq::distance(ps3, method = "gower", binary = TRUE)

# hc <- hclust(dist_b,method="ward.D")
# 
# 
# labels <-   cbind(as.character(ps@sam_data$Subject),
#                         as.character(ps@sam_data$Breed),
#                   as.character(ps@sam_data$sampleID),
#                   as.character(ps@sam_data$Phase)
#                   #paste0(as.character(ps@sam_data$sampleID),"-",as.character(ps@sam_data$Phase))
#                   )
# 
# row.names(labels) <- as.character(ps@sam_data$Subject)
# labels <- as.data.frame(labels); labels$x <- row.names(labels)
# labels$V1 <- as.character(labels$V1)
# #We will color the labels according to countries(group_info[,1])
# hc_d <- dendro_data(as.dendrogram(hc))
# hc_d$labels$V1 <- as.character(hc_d$labels$x)
# # hc_d$labels$Type <- merge(hc_d$labels,labels,by.y="row.names",by.x="label")[,8]
# # hc_d$labels$Treatment <- merge(hc_d$labels,labels,by.y="row.names",by.x="label")[,6]
# 
# hc_d$labels$Type <-left_join(hc_d$labels,labels,"V1")[,8]
# hc_d$labels$Treatment <-left_join(hc_d$labels,labels,"V1")[,6]
#   
# hc_d$labels$color <- NA
# hc_d$labels$color[which(hc_d$labels$Treatment=="Low feed efficiency")] <- "red"
# hc_d$labels$color[which(hc_d$labels$Treatment=="High feed efficiency")] <- "blue"
# hc_d$labels$color[which(hc_d$labels$Treatment=="High feed efficiency")] <- "green"
# #hc_d$labels$color=cols[hc_d$labels$Type]
# 
# ## Plot clusters
# p1 <- ggplot(data = segment(hc_d)) +
#   geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
#   coord_flip() +
#   scale_x_discrete(labels=label(hc_d)$Type) +
#   geom_point(data=hc_d$label, aes(x = x, y = y, color = color), alpha = 1) +
#   #p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
#   #scale_color_manual(values = cols)+
#   scale_color_manual(name="Groups",labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
#   ylab("Distance (beta diversity: Bray-Curtis)") + theme_bw()+
#   theme(axis.text.y = element_text(color = hc_d$labels$Treatment),
#         axis.title.y = element_blank())
# p1 <- p1 + geom_text(data=label(hc_d),
#                 aes(label=Type, x=x, y=-0.5, colour=hc_d$labels$color))
# p1

out.wuf.log <- ordinate(ps3, method = "PCoA", distance = "bray")
evals <- out.wuf.log$values$Eigenvalues #MDS
po3 <- plot_ordination(ps3, out.wuf.log, color = "Breed") +
  labs(col = "Breed") +
  coord_fixed(sqrt(evals[2] / evals[1]))

po3

scatter3d(x = out.wuf.log$vectors[,1], y = out.wuf.log$vectors[,2], z = out.wuf.log$vectors[,3], 
          groups = factor(paste0(po3$data$Breed),
                          levels=c("Cheviot",
                                   "Connemara",
                                   "Lanark")),
          xlab= paste("Axis 1",round(out.wuf.log$values$Relative_eig[1]*100,2),"%"),
          ylab= paste("Axis 2",round(out.wuf.log$values$Relative_eig[2]*100,2),"%"),
          zlab= paste("Axis 3",round(out.wuf.log$values$Relative_eig[3]*100,2),"%"),
          surface=FALSE, grid = FALSE, ellipsoid = F,
labels = paste0(po3$data$sampleID," - ",po3$data$Phase),
point.col =  c("black", "black", "black"),id.n=length(po3$data$sampleID),
surface.col = c("darkblue", "red", "green"),
axis.col = c("black", "black", "black")
)

legend3d("topright", legend = c("Cheviot",
                                   "Connemara",
                                "Lanark"),
         pch = 16, col = c("darkblue", "red", "green"), cex=2, inset=c(0.01))

###############################################
# 
# dist_b_r <-(as.matrix(dist_b))
# diag(dist_b_r)<- NA
# 
# dist_br_f  <- as.data.frame(cbind(ps@sam_data$Subject,
#                                   ps@sam_data$sampleID,
#                                   ps@sam_data$Treatment,
#                                   ps@sam_data$pH,
#                                   ps@sam_data$FCR,
#                                   colMedians(dist_b_r,na.rm = T),
#                                   colSds(dist_b_r,na.rm = T)
#                                   ))
# dist_br_f[,2] <- as.factor(ps@sam_data$sampleID)
# dist_br_f[,3] <- as.factor(ps@sam_data$Treatment)
# colnames(dist_br_f) <- c("Subject","Sample","Treatment","pH","FCR","Dist","Dist_SD")
# #plot(ps@sam_data$FCR,colMedians(dist_b_r,na.rm = T),col= ps@sam_data$Treatment,pch=)
# 
# 
# 
# p_dist<- ggplot(data = dist_br_f) +
#   # geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
#   # coord_flip() +
#   # scale_x_discrete(labels=label(hc_d)$Type) +
#   geom_point(data=dist_br_f, aes(x = pH , y = Dist, color = Treatment), alpha = 1) +
#   #p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
#   #scale_color_manual(values = cols)+
#   scale_color_manual(name="Groups",labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
#   ylab("Distance (beta diversity: Bray-Curtis)") + theme_bw()+
#   theme(axis.text.y = element_text(color = dist_br_f$Treatment),
#         axis.title.y = element_blank())
# p_dist <- p_dist + geom_text(data=dist_br_f,
#                      aes(label=Sample, x=pH, y=Dist, colour=Treatment))
# p_dist
# #require(animation);require(rgl);require(magick)
# # dir.create("animation")
# # p_images <- list()
# # for (i in 1:90) {
# #   p_images[[i]] <- view3d(userMatrix=rotationMatrix(2*pi * i/90, 1, -1, -1))
# #   p_images[[i]] <- rgl.snapshot(filename=paste("animation/frame-",
# #                               sprintf("%03d", i), ".png", sep=""))
# #   p_images[[i]] <- image_read(paste("animation/frame-", sprintf("%03d", i), ".png", sep=""))
# #   
# # }
# # library(magick)
# # 
# # p_images_c <- do.call(c,p_images)
# # A.img <- image_scale(p_images_c)
# # A.ani <- image_animate(A.img, fps = 20, dispose = "previous")
# # image_write(A.ani, paste0("A_animation.gif"))
# # rm(p_images,p_images_c,A.img,A.ani);gc()
# 
# 
# 


#PERMANOVA BRAY CURTIS

#http://thebiobucket.blogspot.com/2011/04/assumptions-for-permanova-with-adonis.html

#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
library(vegan)
permanova <- adonis(phyloseq::distance(ps3, method = "bray") ~
                      sample_data(ps3)$Breed ,
  #otu_table(ps3) ~ sample_data(ps3)$Treatment,
                     permutations=10000, method = "bray")
beta <- betadisper(vegdist(otu_table(ps3), method = "bray"), metadata$Breed)
perm_beta <- permutest(beta,permutations = 10000)
plot(beta)

# coef <- as.coefficients(permanova)
# top.coef <- coef[rev(order(abs(coef)))[1:20]]
# par(mar = c(3, 14, 2, 1))
# barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
#   
# 
# ssCor <- as.data.frame(t(cor(sample_data(ps3)$FCR,otu_table(ps3))))
# ssCor$OTU <- NA;ssCor$OTU <- row.names(ssCor)
# ssCor <- ssCor[which(ssCor[,1] >= 0.5 | ssCor[,1] <= -0.5),]
# library(vegan)


#https://bioconductor.org/packages/devel/bioc/vignettes/microbiome/inst/doc/vignette.html
# Cross correlate data sets and return a table



##################################################################################
##################################################################################
# #RDA
# 
 bray_not_na <- phyloseq::distance(ps3, method = "bray")
ps4 <- ps
# for(i in 1:ncol(ps4@otu_table)){
#   cat(i,"\n")
#   ps4@otu_table[,i][which(ps4@otu_table[,i]>0)] <- 1
# };rm(i)
bray_jacc <- phyloseq::distance(ps4, method = "jaccard",binary = TRUE)
hclust_Jacc <- hclust(bray_jacc,"ward.D")
plot(hclust(bray_jacc))
# s_df <- data.frame(sample_data(ps3))
# s_df$Subject <- factor(s_df$Subject)



########################################################
###############################################
# 
# dist_j_r <-(as.matrix(bray_jacc))
# diag(dist_j_r)<- NA
# 
# dist_j_f  <- as.data.frame(cbind(ps@sam_data$Subject,
#                                   ps@sam_data$sampleID,
#                                   ps@sam_data$Treatment,
#                                   ps@sam_data$pH,
#                                   ps@sam_data$FCR,
#                                   colMedians(dist_j_r,na.rm = T),
#                                   colSds(dist_j_r,na.rm = T),
#                                  global(ps, index = "shannon"),
#                                  global(ps, index = "simpson")
# ))
# dist_j_f[,2] <- as.factor(ps@sam_data$sampleID)
# dist_j_f[,3] <- as.factor(ps@sam_data$Treatment)
# colnames(dist_j_f) <- c("Subject","Sample","Treatment","pH","FCR","Dist","Dist_SD","Shannon","Simpson")
# #plot(ps@sam_data$FCR,colMedians(dist_b_r,na.rm = T),col= ps@sam_data$Treatment,pch=)
# 
# (dist_j_f$Dist/max(dist_j_f$Dist))/(dist_j_f$Shannon/max(dist_j_f$Shannon))
# 
# p_dist_j <- ggplot(data = dist_j_f) +
#   # geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
#   # coord_flip() +
#   # scale_x_discrete(labels=label(hc_d)$Type) +
#   geom_point(data=dist_j_f, aes(x = Shannon , y = Dist, color = Treatment), alpha = 1) +
#   #p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
#   #scale_color_manual(values = cols)+
#   scale_color_manual(name="Groups",labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
#   ylab("Distance (beta diversity: Jaccard)") + theme_bw()+
#   theme(axis.text.y = element_text(color = dist_j_f$Treatment),
#         axis.title.y = element_blank())
# p_dist_j <- p_dist_j + geom_text(data=dist_j_f,
#                              aes(label=Sample, x=Shannon, y=Dist, colour=Treatment))
# p_dist_j
# 
# 

#######################################################



# 
# #hc_d$labels$label <- as.character(hc_d$labels$label)
# ###CLuster
# require(ggdendro)
# #We will color the labels according to countries(group_info[,1])
# hc_d <- dendro_data(as.dendrogram(hclust_Jacc))
# hc_d$labels$Subject <- hc_d$labels$label
# 
# hc_d$labels$Type <- left_join(hc_d$labels,s_df,by="Subject")[,"Phase"]
# hc_d$labels$Treatment <-  left_join(hc_d$labels,s_df,by="Subject")[,"Treatment"]
# hc_d$labels$SubjectID <-  merge(hc_d$labels,s_df,by.y="Subject",by.x="label")[,"sampleID"]
# 
# 
#   #base::merge(hc_d$labels,s_df,by.y="Subject",by.x="label")[,"Phase"]#
# 
# # 
# # hc_d$labels$Treatment <-  merge(hc_d$labels,s_df,by.y="Subject",by.x="label")[,"Treatment"]
# # hc_d$labels$SubjectID <-  merge(hc_d$labels,s_df,by.y="Subject",by.x="label")[,""]
# 
# hc_d$labels$SubjectID <-  paste0(hc_d$labels$SubjectID," - ",hc_d$labels$Type)
# #  merge(hc_d$labels,labels,by.y="row.names",by.x="label")[,5]
# 
# hc_d$labels
# hc_d$labels$color <- NA
# hc_d$labels$color[which(hc_d$labels$Treatment=="Low feed efficiency")] <- "red"
# hc_d$labels$color[which(hc_d$labels$Treatment=="High feed efficiency")] <- "blue"
# #hc_d$labels$color=cols[hc_d$labels$Type]
# 
# ## Plot clusters
# p1 <- ggplot(data = segment(hc_d)) +
#   geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
#   coord_flip() +
#   scale_x_discrete(labels=label(hc_d)$Type) +
#   geom_point(data=hc_d$label, aes(x = x, y = y, color = color), alpha = 1) +
#   #p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
#   #scale_color_manual(values = cols)+
#   scale_color_manual(name="Groups",labels=c("High feed efficiency","Low feed efficiency"),values = c("blue","red"))+
#   ylab("Distance: Jaccard") + theme_bw()+
#   theme(axis.text.y = element_text(color = hc_d$labels$Treatment),
#         axis.title.y = element_blank())
# p1 <- p1 + geom_text(data=label(hc_d),
#                      aes(label=SubjectID, x=x, y=-0.8, colour=hc_d$labels$color))
# p1


ordu_j = ordinate(ps4, "PCoA",distance = "jaccard", binary=T)
p_jac <- plot_ordination(ps4, ordu_j, color="Breed", shape="Phase")
p_jac <- p_jac + 
  scale_color_manual(values = c("blue","red","green"))+
  geom_point(size=9)+
  geom_text_repel(aes(label=sampleID), size=9,colour="black")

p_jac


permanova_j <- adonis(bray_jacc ~
                      sample_data(ps3)$Breed,
                    #otu_table(ps3) ~ sample_data(ps3)$Treatment,
                    permutations=10000, method = "bray")
beta_j <- betadisper(bray_jacc, metadata$Breed)
perm_beta_j <- permutest(beta_j,permutations = 10000)
plot(beta_j)

ordu_UF = ordinate(ps3, "PCoA", "bray")
plot_UF <- plot_ordination(ps3, ordu_UF, color="Breed", shape="Phase")
plot_UF <- plot_UF + 
  scale_color_manual(values = c("blue","red","green"))+
  geom_point(size=9)+
  geom_text_repel(aes(label=sampleID), size=9,colour="black")
plot_UF

   ################


out.wuf.log_J <- ordinate(ps4, method = "PCoA", distance = "jaccard",binary=T)
evals <- out.wuf.log_J$values$Eigenvalues #MDS
po3_j <- plot_ordination(ps4, out.wuf.log_J, color = "Treatment") +
  labs(col = "Treatment") +
  coord_fixed(sqrt(evals[2] / evals[1]))

po3_j


scatter3d(x = out.wuf.log_J$vectors[,1], y = out.wuf.log_J$vectors[,2], z = out.wuf.log_J$vectors[,3], 
          groups = factor(paste0(po3$data$Breed),
                          levels=c("Cheviot","Connemara","Lanark")),
          xlab= paste("Axis 1",round(out.wuf.log$values$Relative_eig[1]*100,2),"%"),
          ylab= paste("Axis 2",round(out.wuf.log$values$Relative_eig[2]*100,2),"%"),
          zlab= paste("Axis 3",round(out.wuf.log$values$Relative_eig[3]*100,2),"%"),
          surface=FALSE, grid = FALSE, ellipsoid = F,
          labels = paste0(po3$data$sampleID," - ",po3$data$Phase),
          point.col =  c("black", "black", "black"),id.n=length(po3$data$sampleID),
          surface.col = c("darkblue", "red", "green"),
          axis.col = c("black", "black", "black")
)

legend3d("topright", legend = c("Cheviot",
                                "Connemara",
                                "Lanark"),
         pch = 16, col = c("darkblue", "red", "green"), cex=2, inset=c(0.01))




permanova_results <- 
data.frame(bray_df =permanova$aov.tab$Df[1],
           bray_SumsOfSqs = permanova$aov.tab$SumsOfSqs[1],
           bray_MeanSqs = permanova$aov.tab$MeanSqs[1],
           bray_R2=permanova$aov.tab$R2[1], 
           bray_F=permanova$aov.tab$F.Model[1],
           bray_Pr_F=permanova$aov.tab$`Pr(>F)`[1],
           bray_perm_bet_df =perm_beta$tab$Df[1],
           #beta
           bray_perm_bet_SumsOfSqs =perm_beta$tab$`Sum Sq`[1],
           bray_perm_bet_df_MeanSq =perm_beta$tab$`Mean Sq`[1],
           bray_perm_bet_df_F =perm_beta$tab$F[1],
           bray_perm_bet_df_F =perm_beta$tab$N.Perm[1],
           bray_perm_bet_df_Pr_F =perm_beta$tab$`Pr(>F)`[1],
           ##################################################
           jaccard_df =permanova_j$aov.tab$Df[1],
           jaccard_SumsOfSqs = permanova_j$aov.tab$SumsOfSqs[1],
           jaccard_MeanSqs = permanova_j$aov.tab$MeanSqs[1],
           jaccard_R2=permanova_j$aov.tab$R2[1], 
           jaccard_F=permanova_j$aov.tab$F.Model[1],
           jaccard_Pr_F=permanova_j$aov.tab$`Pr(>F)`[1],
           #beta
           jaccard_perm_bet_SumsOfSqs =perm_beta_j$tab$`Sum Sq`[1],
           jaccard_perm_bet_df_MeanSq =perm_beta_j$tab$`Mean Sq`[1],
           jaccard_perm_bet_df_F =perm_beta_j$tab$F[1],
           jaccard_perm_bet_df_F =perm_beta_j$tab$N.Perm[1],
           jaccard_perm_bet_df_Pr_F =perm_beta_j$tab$`Pr(>F)`[1]
           )

write.csv(permanova_results,paste0(out_dir,"/","csv/","permanova_results.csv"),quote = F,row.names = F)



lis_beta<- list(permanova=permanova,
                beta=beta,
                perm_beta=perm_beta)
saveRDS(lis_beta,paste0(out_dir,"/","csv/","permanova.RDS"))
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
#https://oliviarata.wordpress.com/2014/09/06/rda-in-ggplot2/
df_to <- as(ps3@sam_data[,-c(1,2,3,4,5,13,15)], "matrix")
ps3a <- ps3

nb <- missMDA::estim_ncpPCA(df_to,ncp.max=10,scale = T) ## Time consuming, nb = 2
res.comp <- missMDA::imputePCA(df_to,ncp=as.numeric(nb[1]))
ps3a@sam_data <- ps3@sam_data[,-c(1,2,3,4,5,13,15)]
  ps3a@sam_data$Subject <- ps3@sam_data$Subject
  ps3a@sam_data$sampleID <- ps3@sam_data$sampleID
  ps3a@sam_data$ADG_Treatment <- ps3@sam_data$ADG_Treatment
  ps3a@sam_data$Treatment <- ps3@sam_data$Treatment
  ps3a@sam_data$Liquid <- ps3@sam_data$Liquid
  ps3a@sam_data$Breed <-   ps3@sam_data$Breed
  ps3a@sam_data$Phase <-   ps3@sam_data$Phase
  #
  ps3a@sam_data$pH <- res.comp$completeObs[,1]
  ps3a@sam_data$ADG <- res.comp$completeObs[,2]
  ps3a@sam_data$FCR <- res.comp$completeObs[,3]
  ps3a@sam_data$R_weight_full <-res.comp$completeObs[,4]
  ps3a@sam_data$R_weight_empty <- res.comp$completeObs[,5]
  ps3a@sam_data$Live_weight <-   res.comp$completeObs[,6]
  #
  ps3a@sam_data$Dead_weight <- res.comp$completeObs[,7]
  ps3a@sam_data$LW.Gain <- res.comp$completeObs[,8]
  
  
  # ############################################################################################################################
  # # CAP ordinate
  # cap_ord <- ordinate(
  #   physeq = ps3a, 
  #   method = "RDA",
  #   distance = bray_jacc,
  #   formula =  OTU ~ 
  #     
  #     pH + ADG + R_weight_full + R_weight_empty + Live_weight + Dead_weight + LW.Gain  
  #   # pH + 
  #   # FCR + #R_weight_full + R_weight_empty + 
  #   #ADG +Live_weight + Dead_weight + LW.Gain #use this
  #   
  # )
  # scor = vegan::scores(cap_ord, display=c("sp", "cn", "bp"), scaling=2) 
  # 
  # # CAP plot
  # cap_plot <- plot_ordination(
  #   physeq = ps3a, 
  #   ordination = cap_ord, 
  #   color = "Breed", 
  #   axes = c(1,2),
  #   label="Subject"
  # ) + 
  #   aes(shape = Phase) + 
  #   geom_point(aes(colour = Breed), alpha = 0, size = 8) + 
  #   #geom_text(size = 15)
  #   # geom_point(colour = "grey90", size = 1.5) + 
  #   scale_color_manual(values = c("blue","red","green")
  #   )
  # ##################################
  # # cap_plota <- phyloseq::plot_ordination(
  # #   physeq = ps3, 
  # #   ordination = cap_ord, 
  # #   type = "sites",
  # #   color = "Breed", 
  # #   axes = c(1,2),
  # #   shape="Phase",
  # #   label="Subject"
  # # )
  # 
  # cap_plot2 <- phyloseq::plot_ordination(
  #   physeq = ps3a, 
  #   ordination = cap_ord, 
  #   type = "taxa",
  #   color = "Phylum", 
  #   axes = c(1,2)#,
  #   #shape="Breed"#,
  #   #label="Subject"
  # )
  # # cap_plot
  # # cap_plot2
  # #)# +
  # # aes(shape = Phase) + 
  # # geom_point(aes(colour = Breed), alpha = 0, size = 8) + 
  # # geom_point(colour = "grey90", size = 1.5) + 
  # # scale_color_manual(values = c("blue","red","green","orange")
  # 
  # ##################################
  # RsquareAdj(cap_ord)$adj.r.squared
  # # anova(cap_ord)
  # # Now add the environmental variables as arrows
  # arrowmat <- vegan::scores(cap_ord, display = "bp")
  # 
  # # Add labels, make a data.frame
  # arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  # 
  # # Define the arrow aesthetic mapping
  # arrow_map <- aes(xend = RDA1, 
  #                  yend = RDA2, 
  #                  x = 0, 
  #                  y = 0, 
  #                  shape = NULL, 
  #                  color = NULL, 
  #                  label = labels)
  # 
  # label_map <- aes(x = 1.3 * RDA1, 
  #                  y = 1.3 * RDA2, 
  #                  shape = NULL, 
  #                  color = NULL, 
  #                  label = labels)
  # 
  # arrowhead = arrow(length = unit(0.02, "npc"))
  # 
  # # Make a new graphic
  # cap_plot <- cap_plot  +
  #   geom_segment(
  #     mapping = arrow_map, 
  #     size = 1, 
  #     data = arrowdf, 
  #     color = "gray", 
  #     arrow = arrowhead
  #   ) + 
  #   geom_text(
  #     mapping = label_map, 
  #     size = 4,  
  #     data = arrowdf, 
  #     show.legend = FALSE
  #   )
  # 
  # ##################################
  # cap_plot2 <- cap_plot2 + geom_segment(
  #   mapping = arrow_map, 
  #   size = 1, 
  #   data = arrowdf, 
  #   color = "gray", 
  #   arrow = arrowhead
  # ) + 
  #   geom_text(
  #     mapping = label_map, 
  #     size = 4,  
  #     data = arrowdf, 
  #     show.legend = FALSE
  #   )
  # # cap_plot <- cap_plot +
  # #   
  # #   theme(text=element_text(size=69),
  # #                 legend.text=element_text(size=60,colour="black"),
  # #                 axis.text.x  = element_text(size=60, hjust = 0, vjust=0.5),
  # #                  axis.text.y  = element_text(size=60,colour="black"))#+
  # # scale_y_continuous(limits = c(-7,7),breaks=seq(-7,7,0.5))+
  # 
  # cap_plot
  # cap_plot2
  # ggsave(paste0(graph_dir,"/","RDA_PLOT",".pdf"),cap_plot,dpi=300,width =100,height=100,units = "cm",limitsize = FALSE)
  # ggsave(paste0(graph_dir,"/","RDA_PLOT_species",".pdf"),cap_plot2,dpi=300,width =100,height=100,units = "cm",limitsize = FALSE)
  # 
  # 
  # 
  
  
  
  ############################################################################################################################
  ps4 <- ps_object2$ps_phylum2
  ps4@phy_tree <- NULL

  
  row.names(ps4@tax_table) <- ps4@tax_table[,2] 
  colnames(ps4@otu_table) <- ps4@tax_table[,2] 
  ps4@sam_data$Subject <- ps3@sam_data$Subject
  ps4@sam_data$sampleID <- ps3@sam_data$sampleID
  ps4@sam_data$ADG_Treatment <- ps3@sam_data$ADG_Treatment
  ps4@sam_data$Treatment <- ps3@sam_data$Treatment
  ps4@sam_data$Liquid <- ps3@sam_data$Liquid
  ps4@sam_data$Breed <-   ps3@sam_data$Breed
  ps4@sam_data$Phase <-   ps3@sam_data$Phase
  #
  ps4@sam_data$pH <- res.comp$completeObs[,1]
  ps4@sam_data$ADG <- res.comp$completeObs[,2]
  ps4@sam_data$FCR <- res.comp$completeObs[,3]
  ps4@sam_data$R_weight_full <-res.comp$completeObs[,4]
  ps4@sam_data$R_weight_empty <- res.comp$completeObs[,5]
  ps4@sam_data$Live_weight <-   res.comp$completeObs[,6]
  #
  ps4@sam_data$Dead_weight <- res.comp$completeObs[,7]
  ps4@sam_data$LW.Gain <- res.comp$completeObs[,8]
  
  bray_jacc_2 <- phyloseq::distance(ps4, method = "jaccard",binary = TRUE)
  ############################################################################################################################
  # CAP ordinate
  

  cap_ord_ph <- ordinate(
    physeq = ps4, 
    method = "RDA",
    distance = bray_jacc_2,
    formula =  OTU ~ 
      
      pH + ADG + R_weight_full + R_weight_empty + Live_weight + Dead_weight + LW.Gain  
    # pH + 
    # FCR + #R_weight_full + R_weight_empty + 
    #ADG +Live_weight + Dead_weight + LW.Gain #use this
    
  )
  scor = vegan::scores(cap_ord_ph, display=c("sp", "cn", "bp"), scaling=2) 
  
  # CAP plot
  cap_plot_ph <- plot_ordination(
    physeq = ps4, 
    ordination = cap_ord_ph, 
    color = "Breed", 
    axes = c(1,2),
    label="Subject"
  ) + 
    aes(shape = Phase) + 
    geom_point(aes(colour = Breed), alpha = 0, size = 10) + 
    
    # geom_point(colour = "grey90", size = 1.5) + 
    scale_color_manual(values = c("blue","red","green")
    )
  
  cap_plot_ph <- cap_plot_ph +
    geom_point(aes(size=50),show.legend = F)#+
    #geom_text(size=50)
  # cap_plot_ph <- cap_plot_ph + 
  #       geom_point(size=10)
  ##################################
  # cap_plota <- phyloseq::plot_ordination(
  #   physeq = ps3, 
  #   ordination = cap_ord, 
  #   type = "sites",
  #   color = "Breed", 
  #   axes = c(1,2),
  #   shape="Phase",
  #   label="Subject"
  # )
  
  cap_plot_ph2 <- phyloseq::plot_ordination(
    physeq = ps4, 
    ordination = cap_ord_ph, 
    type = "taxa",
    color = "Phylum", 
    axes = c(1,2)#,
    #shape="Breed"#,
    #label="Subject"
  )
  cap_plot_ph2 <- cap_plot_ph2 +
    geom_point(size=10)+
    scale_color_brewer(palette="Set3")
  # cap_plot
  # cap_plot2
  #)# +
  # aes(shape = Phase) + 
  # geom_point(aes(colour = Breed), alpha = 0, size = 8) + 
  # geom_point(colour = "grey90", size = 1.5) + 
  # scale_color_manual(values = c("blue","red","green","orange")
  
  ##################################
  RsquareAdj(cap_ord_ph)$adj.r.squared
  # anova(cap_ord)
  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_ord_ph, display = "bp")
  
  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  
  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = RDA1, 
                   yend = RDA2, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  label_map <- aes(x = 1.3 * RDA1, 
                   y = 1.3 * RDA2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  arrowhead = arrow(length = unit(0.02, "npc"))
  
  # Make a new graphic
  cap_plot_ph <- cap_plot_ph  +
    geom_segment(
      mapping = arrow_map, 
      size = 1, 
      data = arrowdf, 
      color = "gray", 
      arrow = arrowhead
    ) + 
    geom_text(
      mapping = label_map, 
      size = 4,  
      data = arrowdf, 
      show.legend = FALSE
    )
  
  ##################################
  cap_plot_ph2 <- cap_plot_ph2 + geom_segment(
    mapping = arrow_map, 
    size = 1, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
    geom_text(
      mapping = label_map, 
      size = 4,  
      data = arrowdf, 
      show.legend = FALSE
    )
  # cap_plot <- cap_plot +
  #   
  #   theme(text=element_text(size=69),
  #                 legend.text=element_text(size=60,colour="black"),
  #                 axis.text.x  = element_text(size=60, hjust = 0, vjust=0.5),
  #                  axis.text.y  = element_text(size=60,colour="black"))#+
  # scale_y_continuous(limits = c(-7,7),breaks=seq(-7,7,0.5))+
  
  cap_plot_ph
  cap_plot_ph2
  ggsave(paste0(graph_dir,"/","RDA_PLOT_Phylum",".pdf"),cap_plot_ph,dpi=300,width =100,height=100,units = "cm",limitsize = FALSE,scale = -2)
  ggsave(paste0(graph_dir,"/","RDA_PLOT_species_Phylum",".pdf"),cap_plot_ph2,dpi=300,width =100,height=100,units = "cm",limitsize = FALSE)
  
 