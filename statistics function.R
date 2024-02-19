statistics <- function(comparison, list, conditions){
  conditions <- sort(conditions)
  sign <- matrix(nrow = 0, ncol = 5)
  colnames(sign) <- c("Cluster", paste0("Mean ", conditions[1]), paste0("Mean ", conditions[2]), "p_val", "test")
  sign <- as.data.frame(sign)

for (i in 1:length(list)) {
if(list[i] %in% comparison$Pair){  
  n <- comparison[comparison$Pair==list[i],]
  
  
  
  if(length(unique(n$cond))==1){
    n[length(rownames(n))+1,] <- data.frame(
      sample = "dummy",
      Pair = n$Pair[1],
      n=0,
      cond=conditions[which(conditions!=unique(n$cond))]
    )
  }
  
  m <- n%>%
    group_by(cond)%>%
    count(cond)  
  
  
  
  if(m$n[1]<3 | m$n[2]<3){
    if(m$n[1]<3){
      while(m$n[1]<3){
        
        n[length(rownames(n))+1,] <- data.frame(
          sample = "dummy",
          Pair = n$Pair[1],
          n=0,
          cond=m$cond[1]
        )
        m <- n%>%
          group_by(cond)%>%
          count(cond)
      }
    }
    if(m$n[2]<3){
      while(m$n[2]<3){
        n[length(rownames(n))+1,] <- data.frame(
          sample = "dummy",
          Pair = n$Pair[1],
          n=0,
          cond=m$cond[2]
        )
        m <- n%>%
          group_by(cond)%>%
          count(cond)
      }
    }
  }
  
  t_W <- wilcox.test(n ~ cond, data = n, paired = F)
  t_T <- try(t.test(n~cond,data=n), silent = T)
  
  norm <- shapiro.test(n$n)
  if (norm$p.value <0.05) {
    sign[i,] <- list(list[i],
                     if (is(t_T, "try-error")) mean(n$n[n$cond==conditions[1]]) else unname(t_T$estimate[1]),
                     if (is(t_T, "try-error")) mean(n$n[n$cond==conditions[2]]) else unname(t_T$estimate[2]),
                     t_W$p.value,
                     t_W$method)
  }else{
    sign[i,] <- list(list[i],
                     if (is(t_T, "try-error")) mean(n$n[n$cond==conditions[1]]) else unname(t_T$estimate[1]),
                     if (is(t_T, "try-error")) mean(n$n[n$cond==conditions[2]]) else unname(t_T$estimate[2]),
                     if (is(t_T, "try-error")) "Inf" else t_T$p.value,
                     if (is(t_T, "try-error")) "t-test" else t_T$method)
  }
} else{
  sign[i,] <- list(list[i],
                   "NA",
                   "NA",
                   "NA",
                   "NA")
}
  
  cat("✓", i, "out of", length(list),"\n")

}
return(sign)
}
