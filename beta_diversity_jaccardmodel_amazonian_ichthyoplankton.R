
#https://rpubs.com/Bruno_Vilela/770763

#install packages
install.packages("tidyverse")
install.packages("betapart")
install.packages("wesanderson")
install.packages("reshape2")

# carregar pacotes
library(tidyverse)
library(betapart)
library(wesanderson)
library(reshape2)
library (ggpubr)

# dataset including reservoir 
#carregar o inputfile

data <- read.csv("diversidade_a_y_b.csv" , fill = TRUE, sep = ";")
data

#beta diversidade total
data1 <- ifelse(data > 0, 1, 0)

beta_total <- beta.multi(data1,
                         index.family = "jaccard")
beta_total
#beta.JTU = turnover
#beta.JNE = nestedess
#beta.JAC = overall beta diversity

#beta diversidade par a par
beta_par <- beta.pair(data1,
                      index.family = "jaccard")

beta_par 


#grafico beta diversidade par a par
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

jtu <- get_upper_tri(as.matrix(beta_par$beta.jtu))

jac <- get_upper_tri(as.matrix(beta_par$beta.jac))

jne <- get_upper_tri(as.matrix(beta_par$beta.jne))

melted_cormat <- reshape2::melt(jtu)
melted_cormat2 <- reshape2::melt(jac)
melted_cormat3 <- reshape2::melt(jne)

melted_cormat <- rbind(melted_cormat2, melted_cormat, melted_cormat3)

melted_cormat$metric <- factor(rep(c("Total", "Turnover", "Nestedness"), 
                                   each = nrow(melted_cormat2)), 
                               c("Total", "Turnover", "Nestedness"))

pal <- wes_palette("Zissou1", 100, type = "continuous")


g <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = pal, 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        text = element_text(size = 16)) +
  scale_x_continuous(name = "", breaks = 1:7) + 
  scale_y_continuous(name = "", breaks = 1:7) + 
  coord_fixed() +
  xlab("") +
  ylab("") +
  facet_wrap(metric ~.) +
  ggtitle("")
g


## dados sem reservatórios

data_sem_reservatorios <- read.csv("diversidade_a_y_b_sem_reservatorios.csv", fill = TRUE, sep = ";")
data_sem_reservatorios

#beta diversidade total
data_sem_reservatorios_1 <- ifelse(data_sem_reservatorios > 0, 1, 0)

beta_total_sem_reservatorios <- beta.multi(data_sem_reservatorios_1,
                         index.family = "jaccard")
beta_total_sem_reservatorios
#beta.JTU = turnover
#beta.JNE = nestedess
#beta.JAC = overall beta diversity

#beta diversidade par a par
beta_par_sem_reservatorios <- beta.pair(data_sem_reservatorios_1,
                      index.family = "jaccard")
beta_par_sem_reservatorios

##################################################################

#grafico beta diversidade par a par
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

jtu1 <- get_upper_tri(as.matrix(beta_par_sem_reservatorios$beta.jtu))

jac1 <- get_upper_tri(as.matrix(beta_par_sem_reservatorios$beta.jac))

jne1 <- get_upper_tri(as.matrix(beta_par_sem_reservatorios$beta.jne))

melted_cormat <- reshape2::melt(jtu1)
melted_cormat2 <- reshape2::melt(jac1)
melted_cormat3 <- reshape2::melt(jne1)

melted_cormat <- rbind(melted_cormat2, melted_cormat, melted_cormat3)

melted_cormat$metric <- factor(rep(c("Total", "Turnover", "Nestedness"), 
                                   each = nrow(melted_cormat2)), 
                               c("Total", "Turnover", "Nestedness"))

pal <- wes_palette("Zissou1", 100, type = "continuous")


g1 <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours = pal, 
                       name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0),
        text = element_text(size = 16)) +
  scale_x_continuous(name = "", labels = c("4", "5", "6", "7"), breaks = 1:4) + 
  scale_y_continuous(name = "", labels = c("4", "5", "6", "7"), breaks = 1:4) + 
  coord_fixed() +
  xlab("") +
  ylab("") +
  facet_wrap(metric ~.) +
  ggtitle("")

g1


### two figures together
tiff('_jaccard_pairwise_beta_diversity_ichthyoplanckton_amazon.tiff', 
     units="in", width=10, height=8, res=300, compression = 'lzw')
ggarrange(g, g1, 
          labels = c("A", "B"),
          hjust = -5, vjust = 5,
          ncol = 1, nrow = 2)
dev.off ()
