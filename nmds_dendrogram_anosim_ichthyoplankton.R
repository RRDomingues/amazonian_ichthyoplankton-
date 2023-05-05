##################################################
##################################################
############ Dendrogram and NMDS #################
##################################################
##################################################

https://jkzorz.github.io/2019/06/06/NMDS.html
https://www.youtube.com/watch?v=OMrtxobDhrM

###### Install packages #####
library (vegan)
library (ggplot2)
library (factoextra)
library (dendextend)
library (igraph)

# load script_estatistica_metabarcoding
load("ictioplancton.RData")

################################ EGGS ##################################


### input file
fish_eggs <- read.csv('csv_nmds_especies_eggs.csv', header = T)
str (fish_eggs)
head(fish_eggs)


### transform data into a matrix w/ presence absence information

# fish.matrix<-as.matrix(fish)
com_eggs = fish_eggs[,3:ncol(fish_eggs)]
m_com_eggs = as.matrix(com_eggs)


###########################
###########################
##### dendrogram_eggs #####
###########################
###########################

clust_eggs<-vegdist (m_com_eggs, method = "jaccard")
hc_eggs <- hclust(clust_eggs, method = "average")

#rename labels of dendrogram
labels_egg = c("3", "4A", "4B", "5A", "6A","6B", "7A", "7B")
labels_egg

hc_eggs$labels_egg <- labels_egg
rm (labels)
#raw dendrogram graph
tiff(file = "dendogram_eggs.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
fviz_dend(hc_eggs, cex = 1.5, main = "", lwd = 1.5, ylab = "")
dev.off()

#customized dendrogram graph
#tiff(file = "color_dendogram_eggs.tiff", width=10, heigh=8, 
#     unit="in",res = 300, compression="lzw")
#fviz_dend(hc_eggs, k = 2, # cortando em 2 grupos
#          cex = 1.5, # tamanho do rótulo
#          ylab = "",
#          main = "",
#          k_colors = c("#211B15", "#211B15"),
#          color_labels_by_k = TRUE, # cores por grupo
#          rect = TRUE, # Adicionar retângulo ao redor dos grupos
#          rect_border = c("#C8C5BE", "#C8C5BE"),
#          rect_fill = TRUE)
#dev.off()

#colored dendrogram
fviz_dend(hc_eggs, cex = 1.5, k = 2, # corte em 2 grupos
          k_colors = "jco", main = "", lwd = 1.5, ylab = "")

#phylogenetic dendrogram
fviz_dend(hc_eggs, k = 2, k_colors = "jco",
          type = "phylogenic", repel = TRUE, main = "", lwd = 1.5, ylab = "")


#####################
#####################
##### NMDS_eggs #####
#####################
#####################

# run NMDS analysis
set.seed(123)
nmds = metaMDS(m_com_eggs, distance = "jaccard")
nmds

# See data distribution 
plot(nmds)
ordiplot (nmds, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds)$sites)

# add site information
data.scores$Sites = fish_eggs$Site
data.scores$Local = fish_eggs$Local
head (data.scores)

#plot NMDS in ggplot2

tiff(file = "nmds_eggs.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
nmds_eggs_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 12, aes( colour = Sites))+ 
  scale_color_viridis_d()+
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 20), 
        legend.text = element_text(size = 20, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 20), 
        axis.title.x = element_text(face = "bold", size = 20, colour = "black"), 
        legend.title = element_text(size = 20, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Sites", y = "NMDS2") 
nmds_eggs_plot
dev.off()

#######################
#######################
##### ANOSIM_eggs #####
#######################
#######################

#ANOSIM #anosim necessita de replicas por localidades
ano = anosim(m_com_eggs, fish_eggs$Local, distance = "jaccard", permutations = 9999)
ano


################################ LARVAE ##################################


### input file
fish_larvae <- read.csv('csv_nmds_especies_larvas.csv', header = T)
str (fish_larvae)
head(fish_larvae)


### transform data into a matrix w/ presence absence information

#fish.matrix<-as.matrix(fish_larvae)
com_larvae = fish_larvae[,4:ncol(fish_larvae)]
m_com_larvae = as.matrix(com_larvae)

######################
######################
##### dendrogram #####
######################
######################

clust_larvae<-vegdist (m_com_larvae, method = "jaccard")
hc_larvae <- hclust(clust_larvae, method = "average")

#rename labels of dendrogram
labels = c("1", "2", "3", "4A", "4B", "5A", "6A","6B", "7A", "7B")
labels

hc_larvae$labels <- labels

#raw dendrogram graph
tiff(file = "dendogram_larvae.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
fviz_dend(hc_larvae, cex = 1.5, main = "", lwd = 1.5, ylab = "")
dev.off()

#customized dendrogram graph
#tiff(file = "color_dendogram_larvae.tiff", width=10, heigh=8, 
#     unit="in",res = 300, compression="lzw")
#fviz_dend(hc_larvae, k = 2, # cortando em 2 grupos
#          cex = 1.5, # tamanho do rótulo
#          ylab = "",
#          main = "",
#          k_colors = c("#211B15", "#211B15"),
#          color_labels_by_k = TRUE, # cores por grupo
#          rect = TRUE, # Adicionar retângulo ao redor dos grupos
#          rect_border = c("#C8C5BE", "#C8C5BE"),
#          rect_fill = TRUE)
#dev.off()

#colored dendrogram
fviz_dend(hc_larvae, cex = 1.5, k = 2, # corte em 2 grupos
          k_colors = "jco", main = "", lwd = 1.5, ylab = "")

#phylogenetic dendrogram
fviz_dend(hc_larvae, k = 2, k_colors = "jco",
          type = "phylogenic", repel = TRUE, main = "", lwd = 1.5, ylab = "")


#####################
#####################
##### NMDS_eggs #####
#####################
#####################

# run NMDS analysis
set.seed(123)
nmds_larvae = metaMDS(m_com_larvae, distance = "jaccard")
nmds_larvae

# See data distribution 
plot(nmds_larvae)
ordiplot (nmds_larvae, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(nmds_larvae)$sites)

# add site information
data.scores$Sites = fish_larvae$Site
data.scores$Local = fish_larvae$Local
data.scores$Type = fish_larvae$Type
head (data.scores)

#plot NMDS in ggplot2

tiff(file = "nmds_larvae.tiff", width=10, heigh=8, 
     unit="in",res = 300, compression="lzw")
nmds_larvae_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 12, aes( shape = Type, colour = Sites))+
  scale_color_viridis_d()+
  theme(axis.text.y = element_text(colour = "black", size = 20, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 20), 
        legend.text = element_text(size = 20, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 20), 
        axis.title.x = element_text(face = "bold", size = 20, colour = "black"), 
        legend.title = element_text(size = 20, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Sites", y = "NMDS2", shape = "Type") 
nmds_larvae_plot
dev.off()

#######################
#######################
##### ANOSIM_eggs #####
#######################
#######################

ano_larvae = anosim(m_com_larvae, fish_larvae$Type, distance = "jaccard", permutations = 9999)
ano_larvae

ano_larvae1 = anosim(m_com_larvae, fish_larvae$Local, distance = "jaccard", permutations = 9999)
ano_larvae1

# save R environment
save.image(file = "ictioplancton.RData")
