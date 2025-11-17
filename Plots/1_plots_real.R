#%%
library(ggplot2)
library(dplyr)
library(tidyr)

dataset <- 'Photoswitch' 
vers <- 'ori' # syn / real
nr_steps = 25
nr_tests = 100
nr_train <- 30
nr_validate <- 100


metrics <- matrix(0,6,8)
coll <- c("black", "cadetblue", "darkred",'darkorange', 'chocolate4', 'green', 'red')
linetyy <- c("solid", "dashed", 'dotted', "dotdash", 'longdash', 'twodash')
pointtyy <- c(4, 16, 17, 1, 2, 5, 9, 6)

seq_str <- c('Random', 'twCRPS1', 'twCRPS2', 'ItwCRPS1', 'ItwCRPS2', 'TMSE','Entropy', 'TIMSE', 'IBV')

nr_steps <- 25
rmse <- matrix(NA,nr_tests,length(seq_str))
crps <- matrix(NA,nr_tests,length(seq_str))
twcrps1 <- matrix(NA,nr_tests,length(seq_str))
twcrps2 <- matrix(NA,nr_tests,length(seq_str))

for(ss in 1:length(seq_str)){
  for(r in 1:nr_tests){
  
  sequential = seq_str[ss]
  
    res <- tryCatch(
      read.csv(paste(sequential, vers, nr_train, nr_validate, r - 1, 'results.txt', sep = '_')),
      error = function(e) {
        cat(paste(sequential, vers, nr_train, nr_validate, r - 1, 'results.txt', sep = '_'), conditionMessage(e), "\n")
        return(matrix(NA, nr_steps,length(seq_str)) )
      }       )
    
    if(sum(is.na(res))==0){
      res <- res[nr_steps,]
      rmse[r,ss] <- res$RMSE
      crps[r,ss] <- res$CRPS
      twcrps1[r,ss] <- res$twCRPS1
      twcrps2[r,ss] <- res$twCRPS2 }
    
  }
  
}

seq_str <- c('Random', 'CRPSv1', 'CRPSv2', 'ICRPSv1', 'ICRPSv2', 'TMSE','Entropy','TIMSE', 'IBV')
coll <- c('0'="grey",'1'="darkred", '2'="cadetblue")
seq_str <- factor(seq_str, levels = seq_str)


colnames(twcrps1) <- seq_str
long_data <- as.data.frame(twcrps1) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                             ifelse(long_data$Method %in% c('CRPSv1', 'CRPSv2', 'ICRPSv1', 'ICRPSv2'),'1','2'))

p1 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = expression(bar(CRPS)[gamma[1]]), x = "Method \n (a)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

median(twcrps2[,5])/median(twcrps2[,1])

colnames(twcrps2) <- seq_str
long_data <- as.data.frame(twcrps2) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSv1', 'CRPSv2', 'ICRPSv1', 'ICRPSv2'),'1','2'))

p4 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  ylim(c(0.05,0.12)) +
  labs(title = expression(bar(CRPS)[gamma[2]]), x = "Method \n (b)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

median(twcrps1[,2])/median(twcrps1[,1])

colnames(crps) <- seq_str
long_data <- as.data.frame(crps) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSv1', 'CRPSv2', 'ICRPSv1', 'ICRPSv2'),'1','2'))

p7 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = expression(bar(CRPS)), x = "Method \n (c)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

median(crps[,5])/median(crps[,1])

colnames(rmse) <- seq_str
long_data <- as.data.frame(rmse) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSv1', 'CRPSv2', 'ICRPSv1', 'ICRPSv2'),'1','2'))

p8 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = 'RMSE', x = "Method \n (h)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

median(twcrps2[,5])/median(twcrps2[,1])
median(twcrps2[,8])/median(twcrps2[,1])
median(twcrps2[,9])/median(twcrps2[,1])

library(ggpubr)
jpeg(file='Real_comp.png', width=17, height=7, units="cm", res=300)
ggarrange(p1,p4,p7, ncol = 3, nrow = 1, heights=c(0.9,1))
dev.off()