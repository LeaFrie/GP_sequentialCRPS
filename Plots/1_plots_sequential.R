#%%
library(ggplot2)
library(dplyr)
library(tidyr)

dataset <- 'Photoswitch' 
vers <- 'syn' # syn / real
nr_steps = 50
nr_tests = 100
nr_train <- 30
nr_validate <- 100

colMedians <- function(x, na.rm = FALSE) {
  apply(x, 2, median, na.rm = na.rm)
}


metrics <- matrix(0,6,8)
coll <- c("black", "cadetblue", "darkred",'darkorange', 'chocolate4', 'green', 'red')
linetyy <- c("solid", "dashed", 'dotted', "dotdash", 'longdash', 'twodash')
pointtyy <- c(4, 16, 17, 1, 2, 5, 9, 6)
pointtyy <- c(4, 16, 17, 15, 18)

seq_str <- c('Random', 'twCRPS1', 'twCRPS2', 'ItwCRPS1', 'ItwCRPS2')

rmse <- matrix(NA,nr_tests,nr_steps)
crps <- matrix(NA,nr_tests,nr_steps)
twcrps1 <- matrix(NA,nr_tests,nr_steps)
twcrps2 <- matrix(NA,nr_tests,nr_steps)
tpr <- matrix(NA,nr_tests,nr_steps)
tnr <- matrix(NA,nr_tests,nr_steps)
rmsep <- matrix(NA,nr_tests,nr_steps)
rmsegamma <- matrix(NA,nr_tests,nr_steps)

rmseM <- matrix(NA,length(seq_str),nr_steps)
crpsM <- matrix(NA,length(seq_str),nr_steps)
twcrps1M <- matrix(NA,length(seq_str),nr_steps)
twcrps2M <- matrix(NA,length(seq_str),nr_steps)
tprM <- matrix(NA,length(seq_str),nr_steps)
tnrM <- matrix(NA,length(seq_str),nr_steps)
rmsepM <- matrix(NA,length(seq_str),nr_steps)
rmsegammaM <- matrix(NA,length(seq_str),nr_steps)

for(ss in 1:length(seq_str)){
  
  sequential = seq_str[ss]
  
  for(r in 1:nr_tests){
    
    res <- tryCatch(
      read.csv(paste(sequential, vers, nr_train, nr_validate, r - 1, 'results.txt', sep = '_')),
      error = function(e) {
        cat(paste(sequential, vers, nr_train, nr_validate, r - 1, 'results.txt', sep = '_'), conditionMessage(e), "\n")
        return(matrix(NA, nr_steps,length(seq_str)) )
      }       )
    
    if(sum(is.na(res))==0){
    res <- res[1:nr_steps,]
    rmse[r,] <- res$RMSE 
    crps[r,] <- res$CRPS
    twcrps1[r,] <- res$twCRPS1
    twcrps2[r,] <- res$twCRPS2
    rmsep[r,] <- res$RMSElevel
    rmsegamma[r,] <- res$RMSElevel1
    tpr[r,] <- res$TP/(res$TP + res$FN)
    tnr[r,] <- res$TP/(res$TP + res$FP) }
    
  }
  
  rmseM[ss,] <- colMeans(rmse, na.rm = TRUE) 
  rmsepM[ss,] <- colMeans(rmsep, na.rm = TRUE)
  rmsegammaM[ss,] <- colMeans(rmsegamma, na.rm = TRUE)
  crpsM[ss,] <- colMeans(crps, na.rm = TRUE)
  twcrps1M[ss,] <- colMeans(twcrps1, na.rm = TRUE)
  twcrps2M[ss,] <- colMeans(twcrps2, na.rm = TRUE)
  tprM[ss,] <- colMeans(tpr, na.rm = TRUE)
  tnrM[ss,] <- colMeans(tnr, na.rm = TRUE)
  
}

seq_str <- c('Random', 'CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2')
seq_str <- factor(seq_str, levels = seq_str)

# crps
df <- as.data.frame(crpsM[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p7 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (g)", y = 'CRPS') +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)

# rmse
df <- as.data.frame(rmseM[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p8 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (h)", y = 'RMSE') +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)


# rmse p
df <- as.data.frame(rmsepM[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str  
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p2 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (b)", y = expression(RMSE[p])) +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)

# rmse gamma
df <- as.data.frame(rmsegammaM[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str  
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p5 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  # geom_line(size = 1) +  # Line for each row
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (e)", y = expression(RMSE[Gamma])) +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)

# twcrps1
df <- as.data.frame(twcrps1M[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str  
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p1 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  # geom_line(size = .7) +  # Line for each row
  geom_point(aes(shape=row_id),size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (a)", y = expression(bar(CRPS)[gamma[1]])) +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)

# twcrps2
df <- as.data.frame(twcrps2M[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str  
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p4 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  # geom_line(size = 1) +  # Line for each row
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy) +
  labs(title = "", x = "Step \n (d)", y = expression(bar(CRPS)[gamma[2]]))

# 
df <- as.data.frame(tprM[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2)
df$row_id <- seq_str  
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p6 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (f)", y = "Sensitivity") +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)

# 
df <- as.data.frame(tnrM[,seq(1,nr_steps,2)])
colnames(df) <- seq(1,nr_steps,2) # 1:nr_steps
df$row_id <- seq_str  
df_long <- pivot_longer(df, cols = -row_id, names_to = "step", values_to = "value")
df_long$step <- as.numeric(df_long$step)

p3 <- ggplot(df_long, aes(x = step, y = value, group = row_id, color = row_id)) +
  geom_point(aes(shape=row_id), size = 1.5) + # Points for each value
  scale_color_discrete(name = "Criterion") +  # Label the legend
  labs(title = "", x = "Step \n (c)", y = "Precision") +
  theme_minimal() +
  scale_colour_manual(values = coll) +
  scale_linetype_manual(values = linetyy) +
  scale_shape_manual(values = pointtyy)

setwd("~/GitHub/GP_sequentialCRPS/Results")
library(ggpubr)
jpeg(file='Photo_syn_seq.png', width=17, height=17, units="cm", res=300)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 3, nrow = 3, heights=c(0.9,1),
          common.legend = TRUE, legend="bottom")
dev.off()




###################################
# comparison

setwd("~/GitHub/GP_sequentialCRPS/Results/synthetic_sequential")

seq_str <- c('Random', 'twCRPS1', 'twCRPS2', 'ItwCRPS1', 'ItwCRPS2', 'TMSE','Entropy' ,'TIMSE', 'IBV')
# seq_str <- c('Random', 'Entropy' ,'Margin')


nr_steps <- 25
rmse <- matrix(NA,nr_tests,length(seq_str))
crps <- matrix(NA,nr_tests,length(seq_str))
twcrps1 <- matrix(NA,nr_tests,length(seq_str))
twcrps2 <- matrix(NA,nr_tests,length(seq_str))
tpr <- matrix(NA,nr_tests,length(seq_str))
tnr <- matrix(NA,nr_tests,length(seq_str))
rmsep <- matrix(NA,nr_tests,length(seq_str))
rmsegamma <- matrix(NA,nr_tests,length(seq_str))

for(ss in 1:length(seq_str)){
  
  sequential = seq_str[ss]
  
  for(r in 1:nr_tests){
    
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
      twcrps2[r,ss] <- res$twCRPS2
      rmsep[r,ss] <- res$RMSElevel
      rmsegamma[r,ss] <- res$RMSElevel1
      tpr[r,ss] <- res$TP/(res$TP + res$FN)
      tnr[r,ss] <- res$TP/(res$TP + res$FP) }
    
  }
  
}

seq_str <- c('Random', 'CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2', 'TMSE','Entropy','TIMSE', 'IBV')
# seq_str <- c('Random', 'Entropy' ,'Margin')
coll <- c('0'="grey",'1'="darkred", '2'="cadetblue")
seq_str <- factor(seq_str, levels = seq_str)

colnames(rmsegamma) <- seq_str
long_data <- as.data.frame(rmsegamma) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))
p5 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = expression(RMSE[Gamma]), x = "Method \n (e)", y = "") +
  theme_minimal() +
  ylim(c(15,75)) + 
  scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

# median(rmsegamma[,6])/median(rmsegamma[,1])
# median(rmsegamma[,4])/median(rmsegamma[,1])


colnames(tpr) <- seq_str
long_data <- as.data.frame(tpr) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))
p6 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = "Sensitivity", x = "Method \n (f)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

# median(tpr[,3])/median(tpr[,1])

colnames(tnr) <- seq_str
long_data <- as.data.frame(tnr) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))

p3 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = "Precision", x = "Method \n (c)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

# median(tnr[,2])/median(tnr[,1])


colnames(rmsep) <- seq_str
long_data <- as.data.frame(rmsep) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))
p2 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = expression(RMSE[p]), x = "Method \n (b)", y = "") +
  theme_minimal() +
  ylim(c(10,68)) +
  scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

# median(rmsep[,2])/median(rmsep[,1])

colnames(twcrps1) <- seq_str
long_data <- as.data.frame(twcrps1) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                             ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))

p1 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  ylim(c(1.5,7.5)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = expression(bar(CRPS)[gamma[1]]), x = "Method \n (a)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"  )

# median(twcrps1[,2])/median(twcrps1[,1])

colnames(twcrps2) <- seq_str
long_data <- as.data.frame(twcrps2) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))

p4 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  ylim(c(0.045,0.105)) + 
  labs(title = expression(bar(CRPS)[gamma[2]]), x = "Method \n (d)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

# median(twcrps2[,5])/median(twcrps2[,1])

colnames(crps) <- seq_str
long_data <- as.data.frame(crps) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))

p7 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title =  expression(bar(CRPS)), x = "Method \n (g)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

# median(crps[,5])/median(crps[,1])

colnames(rmse) <- seq_str
long_data <- as.data.frame(rmse) %>%
  pivot_longer(cols = everything(), names_to = "Method", values_to = "Value")
long_data$Method <- factor(long_data$Method, levels = seq_str)
long_data$Method1 <- ifelse(long_data$Method == 'Random','0',
                            ifelse(long_data$Method %in% c('CRPSγ1', 'CRPSγ2', 'ICRPSγ1', 'ICRPSγ2'),'1','2'))

p8 <- ggplot(long_data, aes(x = Method, y = Value, fill = Method1)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 16) +
  labs(title = 'RMSE', x = "Method \n (h)", y = "") +
  theme_minimal() + scale_fill_manual(values = coll) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    legend.position = "none"
  )

library(ggpubr)
jpeg(file='Syn_comp.png', width=17, height=17, units="cm", res=300)
ggarrange(p1,p2,p3,p4,p5,p6,p7,p8, ncol = 3, nrow = 3, heights=c(0.9,1))
dev.off()