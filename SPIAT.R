library(pacman)
p_load(readxl, SPIAT, tidyverse, ggplot2)

#coordinate data not supplied due to the size 

#SAMPLE PROCESSING
for (i in 1:length(excel_sheets("raw_data/COMET cell counts WT with coordinates.xlsx"))) {
  temp <- read_xlsx("raw_data/COMET cell counts WT with coordinates.xlsx", sheet = i, col_names = F)
  sample <- unlist(c(temp[1,1]))
  print(paste0(sample))
  if (sample == "18CW25.4") {
    next
  }
  cells <- temp[3:dim(temp)[1],1]
  temp <- temp[3:dim(temp)[1],2:dim(temp)[2]]
  temp <- as.matrix(temp)
  colnames(temp) <- c('CD86', 'CD3', 'CD4', 'CD8', 'CD11c', 'CD206','F4/80', 'FOXP3', 'Klrb1c', 'Ly-6G', 'X', 'Y', 'DAPI')
  rownames(temp) <- cells$...1
  intensity_sheet <- t(temp)
  
  phenotypes <- read.table("analyses/WT_intensities_unfiltered_sample_removed_cluster.txt", sep = ",", row.names = 1, header = T)
  phenotypes2 <- phenotypes[phenotypes$sample == sample,]
  
  for (n in 1:length(rownames(phenotypes2))) {
    if (phenotypes2$cluster[n]==1) {
      phenotypes2$phenotype[n] <- c("OTHER")
    }
    if (phenotypes2$cluster[n]==2) {
      phenotypes2$phenotype[n] <- c("FOXP3,CD4")
    }
    if (phenotypes2$cluster[n]==3) {
      phenotypes2$phenotype[n] <- c("CD11c")
    }
    if (phenotypes2$cluster[n]==4) {
      phenotypes2$phenotype[n] <- c("Ly.6G")
    }
    if (phenotypes2$cluster[n]==5) {
      phenotypes2$phenotype[n] <- c("CD3,CD8")
    }
    if (phenotypes2$cluster[n]==6) {
      phenotypes2$phenotype[n] <- c("F4.80,CD4")
    }
    if (phenotypes2$cluster[n]==7) {
      phenotypes2$phenotype[n] <- c("F4.80,CD206")
    }
    if (phenotypes2$cluster[n]==8) {
      phenotypes2$phenotype[n] <- c("CD86")
    }
    if (phenotypes2$cluster[n]==9) {
      phenotypes2$phenotype[n] <- c("CD86,CD3,CD4,CD11c,CD206,F4.80,Klrb1c")
    }
  }
  
  phenotypes2 <- phenotypes2[colnames(intensity_sheet),]
  
  coord_x <- as.numeric(intensity_sheet['X',])
  coord_y <- as.numeric(intensity_sheet['Y',])
  
  intensity_matrix <- intensity_sheet[c(1:10,13),]
  intensity_matrix <- apply(intensity_matrix, 2, as.numeric)
  general_format_image <-  format_image_to_spe(format = "general",
                                               intensity_matrix = intensity_matrix,
                                               phenotypes = phenotypes2$phenotype,
                                               coord_x = coord_x,
                                               coord_y = coord_y)
  
  
  
  formatted_image <- define_celltypes(
    general_format_image, 
    categories = c("OTHER", "FOXP3,CD4", "CD11c", "Ly.6G", "CD3,CD8", "F4.80,CD4", "F4.80,CD206", "CD86", "CD86,CD3,CD4,CD11c,CD206,F4.80,Klrb1c"), 
    category_colname = "Phenotype", 
    names = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9"),
    new_colname = "Cell.Type")
  
  my_colors <- c("#7eb0d5","#b2e061","#bd7ebe","#ffb55a","#ffee65","#beb9db","#fdcce5")
  
  plot <- plot_cell_categories(spe_object = formatted_image, 
                       categories_of_interest = 
                         c("Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8"), 
                       colour_vector = my_colors, 
                       feature_colname = "Cell.Type",
                       cex = 0.5)+ggtitle(paste0(sample))
  
  
  distances <- calculate_pairwise_distances_between_celltypes(
    spe_object = formatted_image, 
    cell_types_of_interest = c("Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8"),
    feature_colname = "Cell.Type")
  
  line <- data.frame(
    comparison = unique(distances$Pair),
    threshold = 19.536
  )
  
  distance_plot <- plot_cell_distances_violin(distances)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle(paste0(sample))
  
  total_interactions <-  distances %>%
    group_by(Pair)%>%
    count(Pair, name = "total_interactions")
  
  cutoff <- distances %>%
    filter(Distance < 19.536)%>%
    group_by(Pair)%>%
    count(Pair, name = "interactions_within_20_um")
    
  distance_counts <- left_join(total_interactions, cutoff, by = "Pair") %>%
    mutate(proportion = interactions_within_20_um/total_interactions)
  
  write.csv(distance_counts, file=paste0("analyses/SPIAT/",sample,".csv"))
  
  save(general_format_image, formatted_image, my_colors, plot, p_cells, distances, line, distance_plot, distance_counts, file=paste0("analyses/SPIAT/",sample,".R"))
  
  print(paste0(sample," sample DONE"))

}



### tissue mapping 

load("analyses/SPIAT/18CW7.3.R")


plot <- plot_cell_categories(spe_object = formatted_image, 
                             categories_of_interest = 
                               c("Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9"), 
                             colour_vector = c("#7eb0d5","#b2e061","#bd7ebe","#ffb55a","#ffee65","#beb9db","#fdcce5", "#8bd3c7"), 
                             feature_colname = "Cell.Type",
                             cex = 0.5)+scale_y_reverse()
plot$data$Cell.Type <-factor(plot$data$Cell.Type, levels = c("OTHER", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8",  "Cluster 9"))


plot$data <- plot$data %>%
  arrange(Cell.Type)


plot+ ggtitle(NULL)+ xlab("Cell X Position") + ylab("Cell Y Position")
plot+ ylim(4150,4650)+xlim(1990,3000)



### Violin plots ###

load("analyses/SPIAT/155 colon.R")

test <- data.frame(matrix(ncol = dim(distances)[2], nrow = 0))
colnames(test) <- colnames(distances)

samples <- c()

files <- list.files(path = "analyses/SPIAT/", pattern = ".R")

for (i in files) {
  load(paste0("analyses/SPIAT/",i))
  test <- rbind(test,distances)
  samples <- c(samples, rep(i, dim(distances)[1]))
}

test$sample <- samples


line <- data.frame(
  comparison = unique(distances$Pair),
  threshold = 19.536
)


plot_cell_distances_violin <- function (cell_to_cell_dist) 
{
  Pair <- Distance <- NULL
  ggplot(cell_to_cell_dist, aes(x = Pair, y = Distance)) + 
    geom_violin() + facet_wrap(~Pair, scales = "free") + 
    theme_bw()
}

plot_cell_distances_violin(test)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ ylim(0,10100)

conditions$treatment <- ifelse(grepl("AOM/DSS", conditions$condition), "AOM/DSS", "Vehicle")
conditions$samples2 <- paste0(conditions$sam,".R")

AOM_DSS <- c(conditions$samples2[conditions$treatment=="AOM/DSS"])
AOM_DSS

test$treatment <- ifelse(test$sample %in% AOM_DSS, "AOM/DSS", "Vehicle")

test_treated <- test[test$treatment == "AOM/DSS",]
test_untreated <- test[test$treatment == "Vehicle",]




plot_cell_distances_violin(test_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Treated tissue")+ ylim(0,10100)
plot_cell_distances_violin(test_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Vehicle tissue")+ ylim(0,10100)

plot_cell_distances_violin(test_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Treated tissue")+ ylim(0,20)
plot_cell_distances_violin(test_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Vehicle tissue")+ ylim(0,20)



test_untreated %>%
  group_by(Pair)%>%
  count(Pair)%>%
  print(n=50)


test_treated %>%
  group_by(Pair)%>%
  count(Pair)%>%
  print(n=50)

test_untreated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  count(Pair)%>%
  print(n=37)


test_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  count(Pair)%>%
  print(n=50)
