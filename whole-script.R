# clean the environment 

rm(list=ls())

# load pacakges 

library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(EnvStats)
library(ggpubr)
library(ggrepel)
library(heatmaply)
library(heatmap3)
library(data.table)
library(viridis)
library(RColorBrewer)
library(ggpval)  
library(corrplot)
library("factoextra")
library("FactoMineR")
library(ggplot2)
library(tidyr)
library(psych)
library(heatmap3)
library(gplots)
library(heatmaply)

# Load data

## antibody reactivity data
az_proteins <-  read.csv("az_Proteinswithmetadata.csv")

## protein families data 
protein_families <- read_tsv("no.antigens.participants.tsv", 
                             show_col_types = FALSE)

# make colum names for protein families unique 
colnames(protein_families) <- make.names(colnames(protein_families), 
                                         unique = T)


## Make protein names unique

protein_families$GeneID.and.Domain <- make.names(protein_families$GeneID.and.Domain,
                                                   unique = T)

## convert to data frame 

protein_families <- as.data.frame(protein_families)

########################################################################
# Determining seroprevalance based on protein family 
################################################################################

# make new matrix

data_new <- matrix(data = NA, nrow = 698, ncol = 5)

# assign matrix Column names

colnames(data_new) <- c(" antigen_number", "participants_no","seroprevalence", 
                        "Family", "Proteins")


# Loop for determining number of proteins with value greater than threshold across participants 

for (i in 1:698){
  data_new[i] <- sum(protein_families[i,5:57]>1.5625)
}

# Fill second column of new dataframe  

data_new[,2]=53

# Determine seroprevalence

data_new[,3]=(data_new[,1]/data_new[,2])*100

# include protein families in the new dataframe 

data_new[,4] <- protein_families$Family

data_new[,5] <- protein_families$GeneID.and.Domain

# Convert to dataframe 

data_new <- as.data.frame(data_new)

# Convert variables to become factors or numerials

data_new$Family <- as.factor(data_new$Family)

data_new$seroprevalence <- as.numeric(data_new$seroprevalence)

data_new$Protein_IDs <- protein_families$GeneID.and.Domain


# Make a plot for the seroprevalence 

seroplot <- ggplot(data_new, aes(Family, seroprevalence))+
  geom_boxplot()+
  theme_classic()+
  geom_dotplot(binaxis = "y", stackdir = "center", stackratio = 1.5, 
               dotsize = 1.5,binwidth = 0.5)+
  geom_hline(yintercept = 10, linetype="dashed", colour="red")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title=element_text(colour = "black", face="bold"))+
  xlab(" Protein groups")+
  ylab("Seroprevalence (%)")+
  scale_y_continuous(breaks = seq(0,100,10))

seroplot
ggsave("seroplot.tiff", height = 4, width = 8, dpi = 300)

################################################################################

# Relationship between gravidity and Number of antigens identified 

################################################################################



# Select specific rows and columns 

az_proteins <- az_proteins[5:57,]

rownames(az_proteins) <- NULL
# convert columns to numeric form 

az_proteins[,33:730] <- as.numeric(unlist(az_proteins[,33:730]))

# assign protein names to columns 

colnames(az_proteins)[33:730] <- protein_families$GeneID.and.Domain

# Matrix of participants 

gravidity_matrix=matrix(data=NA, nrow = 53, ncol =2 )

# Assign column names 

colnames(gravidity_matrix) <- c("Gravidity", "Antigen_number")

# Loop 

for(i in 1:53){
  gravidity_matrix[,1] <- az_proteins[1:53,4]
  gravidity_matrix[i,2] <- sum(az_proteins[i,33:730]>1.5625)
}


# Convert to become dataframe 

gravidity_matrix <- as.data.frame(gravidity_matrix)

# Convert gravidity to become a factor

gravidity_matrix$Gravidity <- as.factor(gravidity_matrix$Gravidity)

# Plot a curve for gravidity versus number of antigens 

gravidity.vs.antibody_breadth <- ggplot(gravidity_matrix, aes(Gravidity, Antigen_number))+
  geom_boxplot()+
  geom_jitter(width = 0.05)+
  theme_classic()+
  ylab("Antibody breadth")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))

gravidity.vs.antibody_breadth

ggsave("gravidity.vs.antibody_breadth.tiff", height = 6, width = 8, dpi = 300)

###############################################################################
# Scatterplot for antibody breadth versus age of participants 
###############################################################################

# Make matrix 
age_matrix <- matrix(data = NA, nrow = 53, ncol = 2)

# Assign column names to matrix 

colnames(age_matrix) <- c("Age", "Antibody_breadth")

# enter data

age_matrix[,1] <- az_proteins[,3]
age_matrix[,2] <- gravidity_matrix[,2]


# Convert to data frame
age_matrix <- as.data.frame(age_matrix)


# Plot scatter plot for age versus antibody breadth

antibody_breadth.vs.age <- ggplot(age_matrix, aes(Age, Antibody_breadth))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="black")+
  theme_classic()+
  ylab("Antibody breadth")+
  theme(axis.text = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  stat_regline_equation(label.x =34)+
  stat_cor(label.x = 34, label.y = 580)

antibody_breadth.vs.age

ggsave("antibody_breadth.vs.age.tiff", height = 6, width = 8, dpi = 300)

##############################################################################
# Plot volcano plot to compare antibody responses between primigravida and 
# multigravida 
################################################################################

# Clean up data to only use multigravida 2 and 3 

az_proteins_volcano <- az_proteins%>%
  filter(Gravida %in% c(1,2,3))



# Classify gravidity into multigravid and primigravid 

az_proteins_volcano$Gravida[az_proteins_volcano$Gravida==1]="Primigravida"
az_proteins_volcano$Gravida[az_proteins_volcano$Gravida==2]="Multigravida"
az_proteins_volcano$Gravida[az_proteins_volcano$Gravida==3]="Multigravida"


# make Gravida a factor

az_proteins_volcano <- az_proteins_volcano%>%
  mutate(Gravida=factor(Gravida, 
                        levels = c("Primigravida", "Multigravida")))

# Identify proteins whose seropositivity was less than 10%

seroposite_proteins <- data_new %>%
  filter(seroprevalence>10)

# Maake the dataframe longer

full_data <- az_proteins_volcano%>%
  pivot_longer(cols = PF3D7_0100100_DBLa0.11:PF3D7_1116100)

# select only the seroreactive proteins for downstream analysis

 full_data <- full_data%>%
   filter(name %in% c(seroposite_proteins$Protein_IDs))


# select primagravida data

primigravida <- full_data%>%
  filter(Gravida=="Primigravida")%>%
  dplyr::select(Study.Number, name, value)%>%
  pivot_wider(names_from = Study.Number, 
              values_from = value) %>%
  as.data.frame()

# Ensure data is in numeric form

primigravida[,2:13] <- as.data.frame(sapply(primigravida[,2:13], as.numeric))


# Assign row names to primigravida 

rownames(primigravida) <- primigravida$name


# Remove the first name

primigravida <- primigravida[-1]


# Assign all rows to become primigravida 

colnames(primigravida)[1:12] <- "primigravida"


# Check if all values are greater than background 

primagravida.check <-matrix(data = NA, ncol = 1, nrow = 678)

colnames(primagravida.check)="sum"

for (i in 1:678){
  primagravida.check[i,] <- sum(primigravida[i,]>1.562500)
}



# Select multigravidas

multigravida <- full_data%>%
  filter(Gravida=="Multigravida")%>%
  dplyr::select(Study.Number, name, value)%>%
  pivot_wider(names_from = Study.Number, 
              values_from = value) %>%
  as.data.frame()



multigravida[,2:31] <- as.data.frame(sapply(multigravida[,2:31], as.numeric))

# Assign rownames to multigravida set

rownames(multigravida) <- multigravida$name


# Remove first row 

multigravida <- multigravida[-1]

# Assign all rows to become multigravida 

colnames(multigravida)[1:30] <- "multigravida"

# Check reactivity 


multigravida.check <-matrix(data = NA, ncol = 1, nrow = 678)

colnames(multigravida.check)="sum"

for (i in 1:678){
  multigravida.check[i,] <- sum(multigravida[i,]>1.562500)
}



# Combine two dataset

combined.dataset <- cbind(primigravida,multigravida)


# Create a matrix

raw_values=matrix(data=NA, nrow = 678, ncol = 1)

colnames(raw_values) <- "P-values"

# Determine raw p-values 

for (i in 1:678){
  raw_values[i] <- t.test(combined.dataset[i,1:12], 
                          combined.dataset[i,13:42], 
                          var.equal=T)$p.value
}



# Find log of all values 

combined.dataset <- log(combined.dataset)

# Find the log mean of each variable 
mean.primagravida <- apply(combined.dataset[1:12],1, geoMean)
mean.multigravida <- apply(combined.dataset[13:42],1, geoMean)

# Calculate the log fold change 

fold.change <- mean.primagravida-mean.multigravida 

# make a dataframe containing log fold change values and p-values 
final.dataset <- as.data.frame(cbind(fold.change,raw_values))

# Indicate protein name as a column 
final.dataset$Protein_ID <- rownames(final.dataset)

# Remove rownames
rownames(final.dataset) <- NULL

# Check column names
colnames(final.dataset)[2] <- "raw_pvalues"

# classify p-values 
final.dataset$Significance[final.dataset$raw_pvalues < 0.05] <- "p-values < 0.05"
final.dataset$Significance[final.dataset$raw_pvalues > 0.05] <- "p-values > 0.05"

# volcano plot
vol.plot <- ggplot(final.dataset, aes(x=fold.change, y=-log10(raw_pvalues), 
                                      col=Significance))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red")+
  xlab("log 2 fold change")+
  ylab("- log 10(p-values)")+
  theme_classic()+
  scale_color_manual(values=c('red', "black"))+
  theme(legend.position = "bottom")+
  geom_text_repel(aes(label=Protein_ID), 
                  subset(final.dataset, 
                         raw_pvalues < 0.05), size=3)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title = element_text(face = "bold"))

vol.plot
ggsave("vol.plot.tiff", height = 6, width = 10, dpi = 300)

###############################################################################

# Select protein that were significant for mutigravida 

data_multigravida <- final.dataset%>%
  filter(fold.change<0.0 & raw_pvalues <0.05)

# filter the proteins from the original dataset

multigravida_data <- full_data%>%
  filter(name %in% c(data_multigravida$Protein_ID))

# Convert to required class

multigravida_data$Gravida <- as.factor(multigravida_data$Gravida)
multigravida_data$value <- as.numeric(multigravida_data$value)

# Plot significant proteins in multi-gravida

multigravida_significant <- ggboxplot(multigravida_data, 
                                      x="Gravida", 
                                      y="value", 
                                      color= "Gravida", 
                                      add="jitter", 
                                      palette = "jco")+
  facet_wrap(~name)+
  ylab("Antibody response")+
  labs(fill="Protein")+
  theme(axis.text = element_text(colour = "black"))+
  stat_compare_means(aes(group=Gravida), method = "t.test")+
  theme(axis.title = element_text(colour = "black", face = "bold"))

multigravida_significant

ggsave("multigravida_significant.tiff", height = 8, width = 10, dpi = 300)


###############################################################################

# Extract data for significant proteins for primigravida

data_primigrivida <- final.dataset%>%
  filter(fold.change>0.0 & raw_pvalues <0.05)

# Extract heatmap data for primigravida

heatmap_data <- full_data%>%
  filter(name %in% c(data_primigrivida$Protein_ID))%>%
  arrange(Gravida)%>%
  select(Study.Number, name, value, Gravida)%>%
  filter(Gravida=='Primigravida')%>%
  as.data.frame()


# Develop heatmap using geom tiles
 Heatmap_primigravida <- ggplot(heatmap_data, aes(Study.Number, name, 
                                                  fill=value))+
   geom_tile(color="grey")+
   ylab("Proteins")+
   xlab("")+
   scale_fill_gradient(low = "white", high = "red")+
   theme(axis.text.y = element_text(colour = "black"))+
   theme(axis.title = element_text(colour = "black", face = "bold"))+
   theme(axis.text.x = element_blank())+
   theme(legend.position = "none")+
   ggtitle("Primigravida")+
   theme(plot.title= element_text(hjust = 0.5, size = 12))
   
Heatmap_primigravida


# Extract heatmap data for primigravida

heatmap_data_multigravida <- full_data%>%
  filter(name %in% c(data_primigrivida$Protein_ID))%>%
  arrange(Gravida)%>%
  select(Study.Number, name, value, Gravida)%>%
  filter(Gravida=="Multigravida")%>%
  as.data.frame()


# Develop heatmap using geom tiles
Heatmap_multigravida <- ggplot(heatmap_data_multigravida, aes(Study.Number, name, 
                                                 fill=value))+
  geom_tile(color="grey")+
  ylab("")+
  xlab("")+
  scale_fill_gradient(low = "white", high = "red")+
  theme(axis.text.y = element_text(colour = "black"))+
  theme(axis.title = element_text(colour = "black", face = "bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  ggtitle("Multigravida")+
  theme(plot.title= element_text(hjust = 0.5, size = 12))

Heatmap_multigravida


library(cowplot)

plot_grid(Heatmap_primigravida, Heatmap_multigravida)

##############################################################################
# PCA Analysis
###################################################################################

# select data for PCA analysis 

pca_data <- az_proteins_volcano[,33:730]

# Assign row names for PCA data 

rownames(pca_data) <- az_proteins_volcano$Study.Number


# Change entire dataframe to become numeric

pca_data <- as.data.frame(sapply(pca_data, as.numeric))

# Run PCA analysis

res.pca <- PCA(pca_data, scale.unit = TRUE, ncp = 6)

# plot the scree plot for PCA analysis 

pca.screeplot <- fviz_eig(res.pca, addlabels =T, ylim=c(0,25))+
  theme_classic()+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

pca.screeplot

# categorize pca plot based on gravidity 

dim1.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                geom.ind = "point", 
                label = "ind", 
                palette=c("red", "green"), 
                addEllipses = TRUE, 
                legend.title="Gravidity", 
                select.var = list(contrib = 5), 
                pointshape=19, 
                col.var = "black")+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim1.biplot

dim2.3.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                 geom.ind = "point", 
                                 label = "ind", 
                                 palette=c("red", "green"), 
                                 addEllipses = TRUE, 
                                 legend.title="Gravidity", 
                                 select.var = list(contrib = 5), 
                                 pointshape=19, 
                                 col.var = "black", axes = c(2,3))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim2.3.biplot


dim3.4biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                geom.ind = "point", 
                                label = "ind", 
                                palette=c("red", "green"), 
                                addEllipses = TRUE, 
                                legend.title="Gravidity", 
                                select.var = list(contrib = 5), 
                                pointshape=19, 
                                col.var = "black",
                                axes = c(3,4))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim3.4biplot

dim4.5.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                                 geom.ind = "point", 
                                 label = "ind", 
                                 palette=c("red", "green"), 
                                 addEllipses = TRUE, 
                                 legend.title="Gravidity", 
                                 select.var = list(contrib = 5), 
                                 pointshape=19, 
                                 col.var = "black",
                                 axes = c(4,5))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim4.5.biplot

dim5.6.biplot <- fviz_pca_biplot(res.pca, col.ind= az_proteins_volcano$Gravida, 
                geom.ind = "point", 
                label = "ind", 
                palette=c("red", "green"), 
                addEllipses = TRUE, 
                legend.title="Gravidity", 
                select.var = list(contrib = 5), 
                pointshape=19, 
                col.var = "black",
                axes = c(5,6))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(colour = "black"))

dim5.6.biplot


library(cowplot)
combined.pca <- plot_grid(pca.screeplot, dim1.biplot, 
                          dim2.3.biplot, dim3.4biplot, dim4.5.biplot, 
                          dim5.6.biplot, nrow = 2, labels = c("A", "B", "C", "D", 
                                                              "E", "F"))


combined.pca
ggsave("combined.pca.tiff", height = 12, width = 20)


###################################################################

 #VAR2CSA 
  VAR2CSA <- full_data%>%
   filter(name %in% c("PF3D7_1200600_DBLpam1",
                      "PF3D7_1200600_DBLpam2", 
                      "PF3D7_1200600_CIDRpam",
                      "PF3D7_1200600_DBLpam3",
                      "PF3D7_1200600_DBLepam4",
                      "PF3D7_1200600_DBLepam5", 
                      "PF3D7_1200600_DBLe10"))
 PF3D7_1200600_DBLe10 <- VAR2CSA%>%
   filter(name=="PF3D7_1200600_DBLe10")
 
var2csa_plot <- ggboxplot(VAR2CSA, 
          x="Gravida", 
          y="value", 
          add = "jitter", 
          palette = "jco", 
          color = "Gravida", 
          xlab="Gravida", 
          ylab = "Reactivity")+
  facet_wrap(~name, scales = "free", ncol = 4)+
  stat_compare_means(aes(group=Gravida), method = "t.test", label.y = 30)+
  scale_y_continuous(breaks = seq(0,30,5))+
  theme(axis.text = element_text(size = 10))
  
 ggsave("var2csa_plot.pdf", height = 8, width = 12, dpi = 300)
 
 ###############################################################################
