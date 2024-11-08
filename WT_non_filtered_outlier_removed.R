library(pacman)
p_load(tidyverse, magrittr, ggplot2, tools, flowCore, FlowSOM, dplyr, Rtsne, 'Biobase', 'gplots', pheatmap, Seurat, reticulate, caret, umap, readxl, patchwork, vsn)


names <- c("B72", "CD3", "CD4", "CD8", "CD11c", "CD206", "F480", "FOXP3", "Klrb1c", "Ly.6G", "sample", "condition")
conditions <- read_csv("raw_data/conditions.csv", col_names = T)
raw_data <- c("WT_intensities")

intensities <- as.data.frame(matrix(ncol = length(names), nrow = 0))

#read in the intensities and normalize per sample within range 1 and 100000
for (i in 1:length(excel_sheets("raw_data/COMET_cell_counts_WT.xlsx"))) {
  temp <- read_xlsx("raw_data/COMET_cell_counts_WT.xlsx", sheet = i, col_names = F)
  sampl <- unlist(c(temp[1,1]))
  if (sampl == "18CW25.4") {
    next
  }else {
  cells <- temp[3:dim(temp)[1],1]
  temp <- temp[3:dim(temp)[1],2:dim(temp)[2]]
  temp <- as.matrix(temp) 
  temp <- matrix(as.numeric(temp), ncol = ncol(temp))
  rownames(temp) <- cells$...1
  temp <- na.omit(temp)
  colnames(temp) <- c('B72', 'CD3', 'CD4', 'CD8', 'CD11c', 'CD206','F480', 'FOXP3', 'Klrb1c', 'Ly.6G')
  processed <- preProcess(temp, method=c("range"), rangeBounds = c(1, 100000))
  temp <- predict(processed, temp)
  temp <- as.data.frame(temp)
  temp$sample <- sampl
  temp$condition <- c(conditions[,sampl])
  intensities <- rbind(intensities, temp)
  rm(temp,sampl)}
}

head(intensities)

#convert the table to work-able object

r_data <- intensities[,1:(dim(intensities)[2]-2)]

head(r_data)

# Check data and data column names -- for this script to work, the first row must be the column names
fcsfilename <- paste(raw_data,".fcs", sep = "")

# Create FCS file metadata - column names with descriptions
metadata <- data.frame(name=dimnames(r_data)[[2]],
                       desc=paste(dimnames(r_data)[[2]],'positive cells'),
                       range = apply(apply(r_data,2,range),2,diff),
                       minRange = apply(r_data,2,min),
                       maxRange = apply(r_data,2,max)
)

head(metadata)

# Create flowframe with intensity data
r_data <- as.matrix(r_data)
data.ff <- new("flowFrame",
               exprs=r_data,
               parameters=AnnotatedDataFrame(metadata))
write.FCS(data.ff, paste("analyses/",fcsfilename, sep = ""))

#read in created flowframe
data <- flowCore::exprs(flowCore::read.FCS(paste("analyses/",fcsfilename, sep = ""), transformation = FALSE, truncate_max_range = FALSE))

# create flowFrame object (required input format for FlowSOM)

data_FlowSOM <- flowCore::flowFrame(data)

###################
### RUN FLOWSOM ###
###################

# set seed for reproducibility

set.seed(4321)

# run FlowSOM (initial steps prior to meta-clustering)

transformList <- flowCore::estimateLogicle(data_FlowSOM, channels = colnames(as.data.frame(data)))
out1 <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE, transformList = transformList)
out1 <- FlowSOM::BuildSOM(out1)
out1 <- FlowSOM::BuildMST(out1)


# extract cluster labels (pre meta-clustering) from output object

labels_pre <- out1$map$mapping[, 1]

# specify final number of clusters for meta-clustering (can also be selected automatically, but this often does not perform well)

k <- 9

# run meta-clustering


seed <- 4321
out <- FlowSOM::metaClustering_consensus(out1$map$codes, k = k, seed = seed)


# extract cluster labels from output object

labels <- out[labels_pre]

# summary of cluster sizes and number of clusters

table(labels)
length(table(labels))

# save cluster labels

res <- data.frame(cluster = labels)
r_data <- as.data.frame(r_data)
head(r_data)
r_data$cluster <- res$cluster
r_data$sample <- intensities$sample
r_data$condition <- intensities$condition
r_data$condition <- sapply(r_data$condition, toString)
head(r_data)

write.csv(r_data, file = paste0("analyses/",raw_data,"_cluster.txt"))

################
### HEAT MAP ###
################

#raw_data = read.table(paste0(raw_data,"_cluster.txt"), header = T, sep = "\t")
#Select all markers and cluster columns


p = r_data[r_data$cluster == "1", ]
w = apply(p[,1:10],2,mean)

for (i in 2:9){
  p = r_data[r_data$cluster == i, ]
  x = apply(p[,1:10],2,mean)
  w = rbind(w,x)
}
rownames(w) <- c(1:9)
head(w)

heatmap <- pheatmap(w, scale = "column", color = bluered(100), legend = F, cluster_rows = F, cluster_cols = F, 
                    border_color = NA)

head(r_data)
r_data <- cbind(r_data, read.table(text = as.character(r_data$condition), sep = " "))
colnames(r_data) <- c("B72","CD3","CD4","CD8","CD11c","CD206","F480","FOXP3","Klrb1c","Ly.6G","cluster","sample","condition","sex","treatment")
head(r_data)

percent <- r_data %>%
  select(sample, condition, cluster, treatment) %>%
  group_by(cluster)%>%
  count(treatment) 
percent <- as.data.frame(percent)

head(percent)

percentV <- percent[percent$treatment=="(vehicle)",]
percentV$percent <- percentV$n/sum(percentV$n)*100
head(percentV)

percentT <- percent[percent$treatment=="(AOM/DSS)",]
percentT$percent <- percentT$n/sum(percentT$n)*100
head(percentT)


percent2 <- rbind(percentV, percentT)
percent2$treatment <- factor(percent2$treatment, levels = c("(vehicle)", "(AOM/DSS)"))
percent2$cluster <- factor(percent2$cluster)
head(percent2)

full <- ggplot(data=percent2, aes(x=percent, y= treatment, fill=cluster))+
  geom_bar(stat = "identity")+
  ylab(NULL)+
  scale_fill_manual(name = "cluster",
                    values = c("1" = "#fd7f6f",
                               "2" = "#7eb0d5",
                               "3" = "#b2e061",
                               "4" = "#bd7ebe",
                               "5" = "#ffb55a",
                               "6" = "#ffee65",
                               "7" = "#beb9db",
                               "8" = "#fdcce5",
                               "9" = "#8bd3c7"))+
  theme_classic()+
  scale_y_discrete(limits=rev)
full

half <- ggplot(data=percent2, aes(x=percent, y= treatment, fill=cluster))+
  geom_bar(stat = "identity", show.legend = FALSE)+
  ylab(NULL)+
  scale_fill_manual(name = "Cluster",
                    values = c("1" = "#fd7f6f",
                               "2" = "#7eb0d5",
                               "3" = "#b2e061",
                               "4" = "#bd7ebe",
                               "5" = "#ffb55a",
                               "6" = "#ffee65",
                               "7" = "#beb9db",
                               "8" = "#fdcce5",
                               "9" = "#8bd3c7"))+
  theme_classic()+
  xlim(c(0,3))+
  scale_y_discrete(limits=rev)



full /
(half | plot_spacer()) + plot_layout(guides = 'collect')




zoom <- ggplot(data=percent2, aes(x=cluster, y=percent, fill=treatment))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_classic()+
  ylim(c(0,10))+
  scale_fill_manual(name="Treatment",
                    values = c("(vehicle)"="#aaccff",
                               "(AOM/DSS)"="#ffee65"))
zoom

top <- ggplot(data=percent2, aes(x=cluster, y=percent, fill=treatment))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme_classic()+
  scale_fill_manual(name="Treatment",
                    values = c("(vehicle)"="#aaccff",
                               "(AOM/DSS)"="#ffee65"))


top





#####################
### dim-reduction ###
#####################





###### only immune cells #####


subset <-  as_tibble(data) %>%
  add_column(cluster = r_data$cluster) %>%
  dplyr::filter(cluster != "1")


#################
### RUN t-SNE ###
#################


set.seed(1234)
# prepare data for Rtsne (matrix format required)
data_Rtsne <- subset[,1:dim(subset)[2]-1]
data_Rtsne <- as.matrix(data_Rtsne)

head(data_Rtsne)
dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]


dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm; runtime: 2-3 min)


set.seed(1234)
out_Rtsne <- Rtsne(data_Rtsne, pca = TRUE, verbose = TRUE)




###################
### CREATE PLOT ###
###################

# load cluster labels (if not still loaded)


labels <- subset$cluster

# select points used by Rtsne

labels_plot <- labels[!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne
dim(data_Rtsne)
# prepare Rtsne output data for plot

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot$Cluster_type <- as.factor(labels_plot)

head(data_plot)

as.tibble(data_plot) %>%
  group_by(Cluster_type)%>%
  count(Cluster_type)

# plot 2-dimensional t-SNE projection

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = data_plot$Cluster_type)) + 
  geom_point(aes(colour = as.factor(Cluster_type)),
             size = 0.5) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#fd7f6f",
                                "2" = "#7eb0d5",
                                "3" = "#b2e061",
                                "4" = "#bd7ebe",
                                "5" = "#ffb55a",
                                "6" = "#ffee65",
                                "7" = "#beb9db",
                                "8" = "#fdcce5",
                                "9" = "#8bd3c7")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))



##############
#### UMAP ####
##############

library(umap)

data_umap <- umap(subset[,1:dim(subset)[2]-1])

data_umap


head(data_umap$layout)

plot_umap <- as.data.frame(data_umap$layout)

colnames(plot_umap) <- c("UMAP1", "UMAP2")

plot_umap$cluster <- as.factor(labels)



ggplot(plot_umap, aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(aes(colour = as.factor(cluster)),
             size = 0.5) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#fd7f6f",
                                "2" = "#7eb0d5",
                                "3" = "#b2e061",
                                "4" = "#bd7ebe",
                                "5" = "#ffb55a",
                                "6" = "#ffee65",
                                "7" = "#beb9db",
                                "8" = "#fdcce5",
                                "9" = "#8bd3c7")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))








###### Full data ##### 

#################
### RUN t-SNE ###
#################


set.seed(1234)
# prepare data for Rtsne (matrix format required)
data_Rtsne <- data
data_Rtsne <- as.matrix(data_Rtsne)

head(data_Rtsne)
dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]


dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm; runtime: 2-3 min)


set.seed(1234)
out_Rtsne <- Rtsne(data_Rtsne, pca = TRUE, verbose = TRUE)




###################
### CREATE PLOT ###
###################

# load cluster labels (if not still loaded)

labels <- r_data$cluster

# select points used by Rtsne

labels_plot <- labels[!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne
dim(data_Rtsne)

# prepare Rtsne output data for plot

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot$Cluster_type <- as.factor(labels_plot)

head(data_plot)


# plot 2-dimensional t-SNE projection

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = data_plot$Cluster_type)) + 
  geom_point(aes(colour = as.factor(Cluster_type)),
             size = 0.5) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#fd7f6f",
                                "2" = "#7eb0d5",
                                "3" = "#b2e061",
                                "4" = "#bd7ebe",
                                "5" = "#ffb55a",
                                "6" = "#ffee65",
                                "7" = "#beb9db",
                                "8" = "#fdcce5",
                                "9" = "#8bd3c7")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))






##############
#### UMAP ####
##############

library(umap)

#data_umap <- umap(data)
load("analyses/umap.Rdata")
head(data_umap)

plot_umap <- as.data.frame(data_umap$layout)

colnames(plot_umap) <- c("UMAP1", "UMAP2")
samples <- r_data$sample
groups <- r_data$condition
plot_umap$cluster <- as.factor(labels)
plot_umap$sample <- as.factor(samples)
plot_umap$group <- as.factor(groups)
head(plot_umap)

ggplot(plot_umap, aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(aes(colour = as.factor(cluster)),
             size = 0.5) +
  scale_color_manual(name = "cluster",
                     values = c("1" = "#fd7f6f",
                                "2" = "#7eb0d5",
                                "3" = "#b2e061",
                                "4" = "#bd7ebe",
                                "5" = "#ffb55a",
                                "6" = "#ffee65",
                                "7" = "#beb9db",
                                "8" = "#fdcce5",
                                "9" = "#8bd3c7")) + 
  guides(colour = guide_legend(override.aes = list(size=5))) +
  coord_fixed(ratio = 1) + 
  theme(plot.title = element_text(size = 10, face = "bold"),panel.background = element_rect(fill='white', colour='white'))



