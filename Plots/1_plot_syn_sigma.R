library(ggplot2)
library(dplyr)
library(tidyr)

dataset <- 'Photoswitch' 
vers <- 'syn' # syn / real
nr_steps = 25
nr_tests = 50
nr_train <- 30
nr_validate <- 100

sequential <- 'twCRPS2'
seq_str <- factor(c(0.2,0.4,0.6,0.8,1.0,1.2), levels = c(0.2,0.4,0.6,0.8,1.0, 1.2))

crps <- matrix(NA, nr_tests, length(seq_str))
sens <- matrix(NA, nr_tests, length(seq_str))
prec <- matrix(NA, nr_tests, length(seq_str))
crpsM <- NA
sensM <- NA
precM <- NA

for(ss in 1:length(seq_str)){
  
  for(r in 1:nr_tests){
    
    res <- tryCatch(
      read.csv(paste('sigma', sequential, vers, nr_train, nr_validate, r, ss, 'results.txt', sep='_')),
      error = function(e) {
        cat(paste('sigma', sequential, vers, nr_train, nr_validate, r, ss, 'results.txt', sep='_'), conditionMessage(e), "\n")
        return(matrix(NA, nr_steps,length(seq_str)) )
      }       )
    
    if(sum(is.na(res))==0){
    res <- res[nr_steps,]
    crps[r,ss] <- res$CRPS
    sens[r,ss] <- res$TP/(res$TP + res$FN)
    prec[r,ss] <- res$TP/(res$TP + res$FP)  }
    
    
  }
  
  crpsM <- colMeans(crps, na.rm = TRUE)
  sensM <- colMeans(sens, na.rm = TRUE)
  precM <- colMeans(prec, na.rm = TRUE)
  
}



crps_df1 <- data.frame(Value = crpsM, Index = 66.05*c(0.2,0.4,0.6,0.8,1.0,1.2))
crps_df1$method = 'Pointwise'
tpr_df1 <- data.frame(Value = sensM, Index = 66.05*c(0.2,0.4,0.6,0.8,1.0,1.2))
tpr_df1$method = 'Pointwise'
tnr_df1 <- data.frame(Value = precM, Index = 66.05*c(0.2,0.4,0.6,0.8,1.0,1.2))
tnr_df1$method = 'Pointwise'

sequential <- 'ItwCRPS2'
seq_str <- factor(c(0.2,0.4,0.6,0.8,1.0,1.2), levels = c(0.2,0.4,0.6,0.8,1.0, 1.2))

crps <- matrix(NA, nr_tests, length(seq_str))
sens <- matrix(NA, nr_tests, length(seq_str))
prec <- matrix(NA, nr_tests, length(seq_str))
crpsM <- NA
sensM <- NA
precM <- NA

for(ss in 1:length(seq_str)){
  
  for(r in 1:nr_tests){
    
    res <- tryCatch(
      read.csv(paste('sigma', sequential, vers, nr_train, nr_validate, r, ss, 'results.txt', sep='_')),
      error = function(e) {
        cat(paste('sigma', sequential, vers, nr_train, nr_validate, r, ss, 'results.txt', sep='_'), conditionMessage(e), "\n")
        return(matrix(NA, nr_steps,length(seq_str)) )
      }       )
    
    if(sum(is.na(res))==0){
      res <- res[nr_steps,]
      crps[r,ss] <- res$CRPS
      sens[r,ss] <- res$TP/(res$TP + res$FN)
      prec[r,ss] <- res$TP/(res$TP + res$FP)  }
    
    
  }
  
  crpsM <- colMeans(crps, na.rm = TRUE)
  sensM <- colMeans(sens, na.rm = TRUE)
  precM <- colMeans(prec, na.rm = TRUE)
  
}



crps_df2 <- data.frame(Value = crpsM, Index = 66.05*c(0.2,0.4,0.6,0.8,1.0,1.2))
crps_df2$method = 'SUR'
tpr_df2 <- data.frame(Value = sensM, Index = 66.05*c(0.2,0.4,0.6,0.8,1.0,1.2))
tpr_df2$method = 'SUR'
tnr_df2 <- data.frame(Value = precM, Index = 66.05*c(0.2,0.4,0.6,0.8,1.0,1.2))
tnr_df2$method = 'SUR'

crps_df <- rbind(crps_df1, crps_df2)
tpr_df <- rbind(tpr_df1, tpr_df2)
tnr_df <- rbind(tnr_df1, tnr_df2)

coll <- c("cadetblue", "darkred")


p1 <- ggplot(crps_df, aes(x = Index, y = Value, shape=method, color=method)) +
  geom_point(size = 3) + # Optional: Add points to highlight data
  labs(
    title = '(a)',
    x = expression(sigma[nu]),
    y = "CRPS"
  ) +
  ylim(c(16,20)) +
  scale_colour_manual(values = coll) +
  theme_minimal(base_size = 15) +       # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 12) # Center the title
  )

p2 <- ggplot(tpr_df, aes(x = Index, y = Value, shape=method,color=method)) +
  geom_point(size = 3) + # Optional: Add points to highlight data
  labs(
    title = '(b)',
    x = expression(sigma[gamma]),
    y = 'Sensitivity'
  ) +
  ylim(c(0.6,0.75)) +
  scale_colour_manual(values = coll) +
  theme_minimal(base_size = 15) +       # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 12) # Center the title
  )

p3 <- ggplot(tnr_df, aes(x = Index, y = Value, shape=method, color=method)) +
  geom_point(size = 3) + # Optional: Add points to highlight data
  labs(
    title = '(c)',
    x = expression(sigma[gamma]),
    y = "Precision"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  ylim(c(0.82,0.9)) +
  scale_colour_manual(values = coll) +
  theme_minimal(base_size = 15) +       # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 12) # Center the title
  )


jpeg(file='sigma_nu.png', width=18, height=7, units="cm", res=300)
ggarrange(p1,p2,p3, ncol = 3, nrow = 1, heights=c(0.9,1),
          common.legend = TRUE, legend="bottom")
dev.off()


