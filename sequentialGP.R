library(foreach)
library(doParallel)
library(kergp)
library(mvtnorm)
library(scoringRules)
library(cubature)

dataset <- 'Photoswitch' 
version <- 'syn' # syn / real
nr <- 0 # for the cluster (15 tests at once)
nr_validate <- 100
nr_tests <- 15
nr_train <- 30
nr_seq <- 50
seq_str <- c('Random', 'twCRPS1', 'twCRPS2', 'ItwCRPS1', 'ItwCRPS2', 'TMSE','Entropy', 'IBV', 'TIMSE')
nr_para <- nr_tests * length(seq_str)

if(dataset == 'Photoswitch' & version == 'syn'){
  data <- read.csv("data_syn.csv")
  data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
  colnames(data) <- c(paste0 ("x", seq(ncol(data)-2)), 'response', 'response_clean')
  nvar <- ncol(data)-2
  ndata <- nrow(data)}

if(dataset == 'Photoswitch' & version == 'real'){
  data <- read.csv("data.csv", row.names = 1)
  data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
  colnames(data) <- c(paste0 ("x", seq(ncol(data)-1)), 'response')
  nvar <- ncol(data)-1
  ndata <- nrow(data)}

ker.tanimoto <- function(x1, x2, par) {
  #x1,x2: input vectors
  #par: parameters. Here par = c(sigma^2)
  scalar.prod <-  sum(x1*x2)
  norm2.x1 <- sum(x1^2) #squared norm of x1
  norm2.x2 <- sum(x2^2) #squared norm of x2
  d1 <- scalar.prod / (norm2.x1 + norm2.x2 - scalar.prod)
  K <- par[1]*d1
  attr(K, "gradient") <- c(sigma2 <- d1) #gradient w.r.t. sigma^2
  return(K)
}

#Transform the custom function into a covariance kernel that can be used by the package 
cov.tanimoto <- covMan(kernel = ker.tanimoto, hasGrad=TRUE, d=nvar,
                       label = "Tanimoto kernel", parNames=c("sigma2"), 
                       par = c(sigma2=100), parLower =c(sigma2=0), parUpper=c(sigma2=Inf))
                       
                       cov.exp <- covRadial(k1Fun1 = k1Fun1Exp, d=nvar, cov="homo", iso=1L, hasGrad = TRUE)

fitGP <- function(formula, train, cov, obs_noise, obs_varNoiseIni, obs_varNoiseLower){
  fit <- try(gp(formula=formula, data=train, cov=cov, 
                noise=obs_noise,
                varNoiseLower=obs_varNoiseLower, varNoiseUpper=obs_varNoiseUpper,
                optimFun = "stats::optim",
                optimMethod="L-BFGS-B", checkNames=FALSE))
  if(is.character(fit)){
    fit <- try(gp(formula = formula, 
                  data = train, cov=cov, noise=obs_noise,
                  varNoiseLower=obs_varNoiseLower, varNoiseUpper=obs_varNoiseUpper,
                  optimFun = "nloptr::nloptr",
                  opts=list(xrel_tol=1e-4)))  }
  return(fit)  }

T <- quantile(data$response, probs=0.8)
cov <- cov.tanimoto 
formula <- as.formula("response ~ 1")  
obs_noise <- TRUE
kernel <- 'Tanimoto'
if(version == 'syn'){
  obs_varNoiseLower <- 77.46
  obs_varNoiseUpper <- 77.46 }
if(version == 'real'){
  obs_varNoiseLower <- 1e-10
  obs_varNoiseUpper <- 100 }


#################################################################################
# functions for twcrps...

twcrps1_scoring <- function(y, pred_mean, pred_sd){
  y_tildeT <- (pmax(y,T)-pred_mean)/pred_sd 
  T_tilde <- (T-pred_mean)/pred_sd
  val <- mean(pred_sd*(-T_tilde*pnorm(T_tilde)^2 
                                 + y_tildeT*(2*pnorm(y_tildeT)-1)
                                 + (2*dnorm(y_tildeT) - 2*dnorm(T_tilde)*pnorm(T_tilde))
                                 - 1/sqrt(pi)*(1-pnorm(T_tilde*sqrt(2)))))
  return(val)  }

twcrps2_scoring <- function(y, pred_mean, pred_sd, k1){

  value <- NA
  for(inn in 1:length(pred_mean)){
    
    m <- c(pred_mean[inn] - T, pred_mean[inn] - T)
    C <- matrix(c(pred_sd[inn]^2 + k1^2, k1^2, k1^2, pred_sd[inn]^2 + k1^2),2,2)
    value1 <- pmvnorm(upper=c(0,0),mean=m,sigma=C)
    
    m <- c(pred_mean[inn] - T, y[inn] - T)
    C <- matrix(c(pred_sd[inn]^2 + k1^2, k1^2, k1^2, k1^2),2,2)
    value2 <- pmvnorm(upper=c(0,0),mean=m,sigma=C)
    
    value[inn] <- value1 - 2*value2 + pnorm((T - y[inn])/k1)  }
  
  return(mean(value))  }

twcrps1_pointwise <- function(pred_mean, pred_sd){
  T_tilde <- (T-pred_mean)/pred_sd
  crit_val <- pred_sd*
      (T_tilde*pnorm(T_tilde)^2 - T_tilde*pnorm(T_tilde) - 1/sqrt(pi)*pnorm(T_tilde*sqrt(2)) +
         2*dnorm(T_tilde)*pnorm(T_tilde) - dnorm(T_tilde) + 1/sqrt(pi))
  return(which.max(crit_val))}

twcrps2_pointwise <- function(pred_mean, pred_sd, k1){

    crit1 <- NA
    for(inn in 1:length(pred_mean)){
      m <- c(pred_mean[inn] - T, T - pred_mean[inn])
      C <- matrix(c(pred_sd[inn]^2 + k1^2, -k1^2, -k1^2, pred_sd[inn]^2 + k1^2),2,2)
      crit1[inn] <- pmvnorm(upper=c(0,0),mean=m,sigma=C) 
    }
    return(which.max(crit1))}

twcrp1_SUR <- function(pred_mean, pred_sd, pred_cov, pred_tau2){
  
  crit1 <- NaN
  zz <- rnorm(100000)
  zzz <- rnorm(100000)
  zzzz <- rnorm(100000)
  for(pp in 1:length(pred_mean)){  # which point to add
    
    unc <- NaN
    for(xin in 1:length(pred_mean)){  # uncertainty in molecule xin
      sigma2_n1 <- pmax(pred_sd[xin]^2 - pred_cov[pp,xin]^2/(pred_sd[pp]^2 + pred_tau2), 1e-20)
      alpha <- pred_cov[pp,xin]/sqrt(pred_sd[pp]^2+pred_tau2)
      mn <- pred_mean[xin]
      unc[xin] <- sqrt(sigma2_n1)*mean(pmax(zz - pmax(zzz,(T-mn-alpha*zzzz)/sqrt(sigma2_n1)),0))  }
    
    crit1[pp] <- mean(unc) }
  
  return(which.min(crit1)) }

twcrps2_SUR <- function(pred_mean, pred_sd, pred_cov, pred_tau2, k1){
  crit1 <- NaN
  for(pp in 1:length(pred_mean)){  # which point to add
    
    unc <- NaN
    for(xin in 1:length(pred_mean)){  # uncertainty in molecule xin
      
      m <- c(pred_mean[xin] - T, T - pred_mean[xin])
      sd2_upd <- pred_sd[xin]^2 - pred_cov[pp,xin]^2/(pred_sd[pp]^2 + pred_tau2)
      alpha <- pred_cov[pp,xin]/sqrt(pred_sd[pp]^2+pred_tau2)
      C <- matrix(c(sd2_upd + alpha^2 + k1^2,
                    -alpha^2 - k1^2,
                    -alpha^2 - k1^2,
                    sd2_upd + alpha^2 + k1^2 ),2,2)
      unc[xin] <- pmvnorm(upper=c(0,0),mean=m,sigma=C) }
    
    crit1[pp] <- mean(unc) }
  
  return(which.min(crit1))    }

#################################################################################


# for(ss in 1:length(seq_str)){
registerDoParallel(cores=nr_para)
foreach(ss = 1:nr_para)%dopar%{

  sequential = seq_str[ceiling(ss/nr_tests)]
  r = ss %% nr_tests + nr*15

  set.seed(r) 
  
  if(dataset == 'Photoswitch' & version == 'syn'){
    data <- read.csv("data_syn.csv")
    data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
    colnames(data) <- c(paste0 ("x", seq(ncol(data)-2)), 'response', 'response_clean')}
  
  if(dataset == 'Photoswitch' & version == 'real'){
    data <- read.csv("data.csv", row.names = 1)
    data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
    colnames(data) <- c(paste0 ("x", seq(ncol(data)-1)), 'response')}
  
  # take out a batch for final testing (not considered to add seq)
  shuffle_index <- sample(seq(nrow(data)), replace=FALSE)
  first_i <- shuffle_index[1:nr_validate]
  test_fin <- data[first_i, ]
  data <- data[-first_i, ]
  
  # initial training set
  shuffle_index <- sample(seq(nrow(data)), replace=FALSE)
  first_i <- shuffle_index[1:nr_train]
  train <- data[first_i, ]
  test <- data[-first_i, ]
  
  # set up result frame
  results <- data.frame(row.names = 1:nr_seq)
  results$RMSE <- NA
  results$CRPS <- NA
  results$twCRPS1 <- NA
  results$twCRPS2_1 <- NA 
  results$twCRPS2_2 <- NA 
  results$RMSEgamma <- NA
  results$RMSEp <- NA
  results$TP <- NA
  results$TN <- NA
  results$FP <- NA
  results$FN <- NA
  
  # for twCRPS2
  k1_all <- sd(data$response)
  k1_start <- sd(train$response)
  
  for(it in 1:nr_seq){
  
      # fit the GP model
      fit <- fitGP(formula, train, cov, obs_noise, obs_varNoiseIni, obs_varNoiseLower)
      if(is.character(fit)){ stop("GP fitting failed. Stopping script.") }
      
      #######################################################
      # evaluation with test_fin
      pred_results <- predict(object=fit, newdata = test_fin)
      
      if(version == 'syn'){
        pred_mean <- pred_results$mean
        pred_sd <- pred_results$sd 
        y <- test_fin$response_clean}
      
      if(version == 'real'){
        pred_mean <- pred_results$mean
        pred_sd <- pred_results$sd + sqrt(fit$varNoise)
        y <- test_fin$response}
      
      
      # metrics for both versions
      results$RMSE[it] <- sqrt(mean((y - pred_mean)^2))
      results$CRPS[it] <- mean(crps_norm(y, pred_mean, pred_sd))
      
      # twCRPS1
      results$twCRPS1[it] <- twcrps1_scoring(y, pred_mean, pred_sd)
      
      # twCRPS2
      results$twCRPS2_1[it] <- twcrps2_scoring(y, pred_mean, pred_sd, 0.5*k1_all)
      results$twCRPS2_2[it] <- twcrps2_scoring(y, pred_mean, pred_sd, k1_all)
      
      # metrics for synthetic version
      if(version == 'syn'){
        
        true_class <- 1*(y >= T)
        exc_prob <- pnorm((pred_mean - T)/pred_sd)
        pred_class <- 1*(exc_prob >= 0.5)
        results$TP[it] <- sum(true_class*pred_class)
        results$TN[it] <- sum((1-true_class)*(1-pred_class))
        results$FP[it] <- sum((1-true_class)*(pred_class))
        results$FN[it] <- sum((true_class)*(1-pred_class))  
        
        results$RMSEp[it] <- ifelse(sum(pred_class) == 0, 0,
                                         sqrt(1/sum(pred_class)*sum(pred_class*(y - pred_mean)^2)))
        results$RMSEgamma[it] <- sqrt(1/sum(true_class)*sum(true_class*(y - pred_mean)^2))

      }
      
      ###############################################################
      # sequential data acquisition
      # predict for test set
      pred_results <- predict(object=fit, newdata = test)
      pred_mean <- pred_results$mean
      pred_sd <- pred_results$sd
      pred_tau2 <- fit$varNoise
      exc_prob <- pnorm((pred_results$mean - T)/pred_results$sd)
      
      # if SUR, then we need the covariance matrix
      if(sequential %in% c('IBV', 'TIMSE','ItwCRPS1', 'ItwCRPS2')){
        pred_results <- predict(object=fit, newdata = test, covCompute = TRUE)
        pred_cov <- pred_results$cov}
      
      # random
      if(sequential == 'Random'){new_i <- sample(seq(nrow(test)), replace=FALSE)[1]}
      
      # TMSE
      if(sequential == 'TMSE'){
        crit_val <- pred_sd^2*dnorm(T,pred_mean,pred_sd) 
        new_i <- which.max(crit_val)   }
      
      # entropy
      if(sequential == 'Entropy'){
        crit_val <- -exc_prob*log(exc_prob) - (1-exc_prob)*log(1-exc_prob)
        new_i <- which.max(crit_val)      }
      
      #twCRPS1: pointwise
      if(sequential == 'twCRPS1'){new_i <- twcrps1_pointwise(pred_mean, pred_sd)      }
      
      # twCRPS2: pointwise
      if(sequential == 'twCRPS2'){ new_i <- twcrps2_pointwise(pred_mean, pred_sd, 0.5*k1_all) }
      
      # TIMSE
      if(sequential == 'TIMSE'){
        weight <- dnorm(T,pred_mean,pred_sd)
        timse_vec <- NaN
        for(pp in 1:nrow(test)){ 
          sd2_upd <- pred_sd^2 - 1/(pred_sd[pp]^2+pred_tau2)*c(pred_cov[pp,])*c(pred_cov[pp,])
          timse_vec[pp] <- mean(sd2_upd*weight)}
        new_i <- which.min(timse_vec)}
      
      # IBV
      if(sequential == 'IBV'){
        sur_vec <- NaN
        for(pp in 1:nrow(test)){ 
          sd2_upd <- pmax(pred_sd^2 - 1/(pred_sd[pp]^2+pred_tau2)*c(pred_cov[pp,])*c(pred_cov[pp,]),1e-20)
          a <- (pred_mean[-pp]-T)/sqrt(sd2_upd[-pp])
          c <- pred_sd[-pp]^2/sd2_upd[-pp]
          int <- NaN
          int <- sapply(seq((nrow(test)-1)), function(j){
            pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0), mean=c(a[j],-a[j]), sigma=matrix(c(c[j], 1-c[j], 1-c[j], c[j]), nrow=2, byrow=TRUE))[1] })
          sur_vec[pp] <- mean(int)  }
        new_i <- which.min(sur_vec)  }
      
      # twcrps1: SUR
      if(sequential == 'ItwCRPS1'){
        new_i <- twcrp1_SUR(pred_mean, pred_sd, pred_cov, pred_tau2)   }
      
      # twcrps2: SUR
      if(sequential == 'ItwCRPS2'){
        new_i <- twcrps2_SUR(pred_mean, pred_sd, pred_cov, pred_tau2, 0.5*k1_all) }

      #### ADAPT TRAIN AND TEST 
      train <- rbind(train, test[new_i, ])
      test <- test[-c(new_i), ] 
      
      
    }
    
  write.csv(results, file=paste(dataset, version, sequential, r, 'results.txt', sep='_'), row.names = FALSE)   
    
}





