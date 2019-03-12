library("microbiomeViz");
library("stringr")
ps3 <- readRDS(paste0(csv_dir,"/","Phyloseq_object_trans",".RDS"))
ps3 = microbiomeViz::fix_duplicate_tax(ps3)#
ps3@otu_table <- t(ps3@otu_table)
tr = parsePhyloseq(ps3)

cool = rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols = c(rev(cool), rev(warm))
mypalette <- colorRampPalette(cols)(255)
sm_fam <- table(ps3@tax_table[,5])[order(table(ps3@tax_table[,5]),decreasing = T)]
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

ggsave(paste0(graph_dir,"/","microbiome_viz",".pdf"),p,dpi=600,width =90,height=80,units = "cm",scale=0.8,limitsize = FALSE)
