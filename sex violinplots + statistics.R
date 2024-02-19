





### Violin plots ###
#females
load("analyses/SPIAT/155 colon.R")

test <- data.frame(matrix(ncol = dim(distances)[2], nrow = 0))
colnames(test) <- colnames(distances)

samples <- c()

files <- list.files(path = "analyses/SPIAT/", pattern = ".R")
files_F <- c("18CW7.3.R", "17CW47.2.R", "155 colon.R", "167 colon.R", "18CW17.1.R", "18CW17.3.R", "18CW13.1.R", "17CW50.1.R")

  
pb <- txtProgressBar(min = 0,     
                       max = length(files), 
                       style = 3,    
                       width = 50,   
                       char = "=") 
x=0
  
  for (i in files_F) {
    x = x+1
    load(paste0("analyses/SPIAT/",i))
    test <- rbind(test,distances)
    samples <- c(samples, rep(i, dim(distances)[1]))

    setTxtProgressBar(pb, x)
  }
close(pb)

test$sample <- samples


line <- data.frame(
  comparison = unique(distances$Pair),
  threshold = 19.536
)


test_F <- test

### males
load("analyses/SPIAT/155 colon.R")

test <- data.frame(matrix(ncol = dim(distances)[2], nrow = 0))
colnames(test) <- colnames(distances)

samples <- c()

files <- list.files(path = "analyses/SPIAT/", pattern = ".R")
files_M <- c("17CW49.11.R", "18CW13.9.R", "17CW49.9.R", "17CW49.5.R", "17CW50.3.R", "18CW27.5.R", "18CW27.6.R")


pb <- txtProgressBar(min = 0,     
                     max = length(files), 
                     style = 3,    
                     width = 50,   
                     char = "=") 
x=0

for (i in files_M) {
  x = x+1
  load(paste0("analyses/SPIAT/",i))
  test <- rbind(test,distances)
  samples <- c(samples, rep(i, dim(distances)[1]))
  
  setTxtProgressBar(pb, x)
}
close(pb)

test$sample <- samples


line <- data.frame(
  comparison = unique(distances$Pair),
  threshold = 19.536
)


test_M <- test




plot_cell_distances_violin <- function (cell_to_cell_dist) 
{
  Pair <- Distance <- NULL
  ggplot(cell_to_cell_dist, aes(x = Pair, y = Distance)) + 
    geom_violin() + facet_wrap(~Pair, scales = "free") + 
    theme_bw()
}



conditions <- read.csv("raw_data/conditions.csv", header = FALSE)
for (i in 1:length(colnames(conditions))) {
  conditions[3,i] <- c(paste0(conditions[1,i],".csv"))
  
}
rownames(conditions) <- c("sam","condition","samples")
conditions <- t(conditions)
conditions <- as.data.frame(conditions)



conditions$treatment <- ifelse(grepl("AOM/DSS", conditions$condition), "AOM/DSS", "Vehicle")
conditions$samples2 <- paste0(conditions$sam,".R")

AOM_DSS <- c(conditions$samples2[conditions$treatment=="AOM/DSS"])
AOM_DSS

test_F$treatment <- ifelse(test_F$sample %in% AOM_DSS, "AOM/DSS", "Vehicle")

test_F_treated <- test_F[test_F$treatment == "AOM/DSS",]
test_F_untreated <- test_F[test_F$treatment == "Vehicle",]

test_M$treatment <- ifelse(test_M$sample %in% AOM_DSS, "AOM/DSS", "Vehicle")

test_M_treated <- test_M[test_M$treatment == "AOM/DSS",]
test_M_untreated <- test_M[test_M$treatment == "Vehicle",]




plot_cell_distances_violin(test_F_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Females treated")+ ylim(0,10100)
plot_cell_distances_violin(test_F_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Females untreated")+ ylim(0,10100)
plot_cell_distances_violin(test_M_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Males treated")+ ylim(0,10100)
plot_cell_distances_violin(test_M_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Males untreated")+ ylim(0,10100)

plot_cell_distances_violin(test_F_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Females treated")+ ylim(0,20)
plot_cell_distances_violin(test_F_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Females untreated")+ ylim(0,20)
plot_cell_distances_violin(test_M_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Males treated")+ ylim(0,20)
plot_cell_distances_violin(test_M_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Males untreated")+ ylim(0,20)


plot_cell_distances_violin(test_F_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Females treated")+ ylim(0,50)
plot_cell_distances_violin(test_F_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Females untreated")+ ylim(0,50)
plot_cell_distances_violin(test_M_treated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Males treated")+ ylim(0,50)
plot_cell_distances_violin(test_M_untreated)+geom_abline(data = line, aes(intercept = threshold, slope = 0, colour='#E41A1C'), show.legend = FALSE)+ggtitle("Males untreated")+ ylim(0,50)


test_F_treated %>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)

write_csv(test_F_untreated %>%
            group_by(Pair)%>%
            count(Pair), file = "figures/SPIAT/sex_differences/F_V_full.csv")


write_csv(test_M_treated %>%
              group_by(Pair)%>%
              count(Pair), file = "figures/SPIAT/sex_differences/M_T_full.csv")

write_csv(test_M_untreated %>%
            group_by(Pair)%>%
            count(Pair), file = "figures/SPIAT/sex_differences/M_V_full.csv")



write_csv(test_F_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  count(Pair), file = "figures/SPIAT/sex_differences/F_T_20.csv")

write_csv(test_F_untreated %>%
            filter(Distance<19.536)%>%
            group_by(Pair)%>%
            count(Pair), file = "figures/SPIAT/sex_differences/F_V_20.csv")

write_csv(test_M_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  count(Pair), file = "figures/SPIAT/sex_differences/M_T_20.csv")
  
write_csv(test_M_untreated %>%
              filter(Distance<19.536)%>%
              group_by(Pair)%>%
              count(Pair), file = "figures/SPIAT/sex_differences/M_V_20.csv")

  
  
write_csv(test_F_treated %>%
  filter(Distance<50)%>%
  group_by(Pair)%>%
  count(Pair), file = "figures/SPIAT/sex_differences/F_T_50.csv")

write_csv(test_F_untreated %>%
            filter(Distance<50)%>%
            group_by(Pair)%>%
            count(Pair), file = "figures/SPIAT/sex_differences/F_V_50.csv")


  write_csv(test_M_treated %>%
  filter(Distance<50)%>%
  group_by(Pair)%>%
  count(Pair), file = "figures/SPIAT/sex_differences/M_T_50.csv")
  
  write_csv(test_M_untreated %>%
              filter(Distance<50)%>%
              group_by(Pair)%>%
              count(Pair), file = "figures/SPIAT/sex_differences/M_V_50.csv")


ggplot(test_untreated %>%
         filter(Distance < 100, Pair == "Cluster 2/Cluster 5"), aes(x=Distance))+
  geom_histogram()



test_F %>%
  group_by(Type1)%>%
  count(treatment)


test_M %>%
  group_by(Type1)%>%
  count(treatment)


source("statistics function.R")

#### SIGNIFICANCE testing

list_comp <- c("Cluster 2/Cluster 2", "Cluster 2/Cluster 3","Cluster 2/Cluster 4", "Cluster 2/Cluster 5", "Cluster 2/Cluster 6", "Cluster 2/Cluster 7", "Cluster 2/Cluster 8",
               "Cluster 3/Cluster 3", "Cluster 3/Cluster 4", "Cluster 3/Cluster 5", "Cluster 3/Cluster 6", "Cluster 3/Cluster 7", "Cluster 3/Cluster 8", 
               "Cluster 4/Cluster 4", "Cluster 4/Cluster 5", "Cluster 4/Cluster 6", "Cluster 4/Cluster 7", "Cluster 4/Cluster 8", 
               "Cluster 5/Cluster 5", "Cluster 5/Cluster 6", "Cluster 5/Cluster 7", "Cluster 5/Cluster 8", 
               "Cluster 6/Cluster 6", "Cluster 6/Cluster 7", "Cluster 6/Cluster 8",
               "Cluster 7/Cluster 7", "Cluster 7/Cluster 8",
               "Cluster 8/Cluster 8")

#T vs V Females
comp_F_treated <- test_F_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_F_treated$cond <- "F_treated"

comp_F_untreated <- test_F_untreated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_F_untreated$cond <- "F_vehicle"


comp_F <- rbind(comp_F_treated,comp_F_untreated)


cond_F <- c("F_treated", "F_vehicle")
sig_F <- statistics(comparison = comp_F, list = list_comp, conditions = cond_F)


#T vs V Males

comp_M_treated <- test_M_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_M_treated$cond <- "M_treated"

comp_M_untreated <- test_M_untreated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_M_untreated$cond <- "M_vehicle"


comp_M <- rbind(comp_M_treated,comp_M_untreated)

cond_M <- c("M_treated", "M_vehicle")

sig_M <- statistics(comparison = comp_M, list = list_comp, conditions = cond_M)



#M vs F vehicle
comp_M_vehicle <- test_M_untreated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_M_vehicle$cond <- "M_vehicle"

comp_F_vehicle <- test_F_untreated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_F_vehicle$cond <- "F_vehicle"


comp_V <- rbind(comp_M_vehicle,comp_F_vehicle)


cond_V <- c("M_vehicle", "F_vehicle")
sig_V <- statistics(comparison = comp_V, list = list_comp, conditions = cond_V)


#M vs F treated

comp_M_treated <- test_M_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_M_treated$cond <- "M_treated"

comp_F_treated <- test_F_treated %>%
  filter(Distance<19.536)%>%
  group_by(Pair)%>%
  group_by(sample)%>%
  count(Pair)
comp_F_treated$cond <- "F_treated"


comp_T <- rbind(comp_M_treated,comp_F_treated)

cond_T <- c("M_treated", "F_treated")

sig_T <- statistics(comparison = comp_T, list = list_comp, conditions = cond_T)


sig_F$p_adj <- p.adjust(sig_F$p_val, method = "BH")
sig_M$p_adj <- p.adjust(sig_M$p_val, method = "BH")
sig_V$p_adj <- p.adjust(sig_V$p_val, method = "BH")
sig_T$p_adj <- p.adjust(sig_T$p_val, method = "BH")


write.csv(sig_F, "analyses/SPIAT/treat_sign_females.csv")
write.csv(sig_M, "analyses/SPIAT/treat_sign_males.csv")
write.csv(sig_V, "analyses/SPIAT/sex_sign_vehicle.csv")
write.csv(sig_T, "analyses/SPIAT/sex_sign_treated.csv")


# treated vs untreated 

Treated <- rbind(comp_F_treated, comp_M_treated)
Treated$cond <- c("treated")

vehicle <- rbind(comp_F_vehicle, comp_M_vehicle)
vehicle$cond <- c("vehicle")



comp <- rbind(Treated, vehicle)


condit <- c("treated", "vehicle")

signif <- statistics(comparison = comp, list = list_comp, conditions = condit)

signif$p_adj <- p.adjust(signif$p_val, method = "BH")
write.csv(signif, "analyses/SPIAT/treat_sign.csv")
