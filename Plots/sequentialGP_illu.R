# GP modeling and data acquisition ##############################################

# install and load packages
pack = TRUE

if (pack == TRUE) {
  # 1. Create personal library if missing
  user_lib <- "~/Rlibs"
  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
  .libPaths(c(user_lib, .libPaths()))
  # 2. Packages you want
  packages <- c("kergp", "mvtnorm", "tidyverse", "ggplot2",
                "foreach", "doParallel", "scoringRules")
  # 3. Install Bioconductor/CRAN manager if not present
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  # 4. Install each package if missing
  installed <- rownames(installed.packages())
  for (p in packages) {
    if (!(p %in% installed)) {
      tryCatch(
        install.packages(p, dependencies = TRUE, repos = "https://cloud.r-project.org"),
        error = function(e) {
          message(sprintf("Could not install %s: %s", p, e$message))
        }
      )
    }
  }
  # 5. Load all packages
  lapply(packages, require, character.only = TRUE)    }


if(pack == FALSE){
  library(foreach)
  library(doParallel)
  library(mvtnorm)
  library(scoringRules)
  library(cubature)  
  library(kergp)}


### FUNCTIONS ###

# GP model ######################################################################
# tanimoto kernel
ker.tanimoto <- function(x1, x2, par) {
  #x1,x2: input vectors
  #par: parameters. Here par = c(sigma^2)
  scalar.prod <-  sum(x1*x2)
  norm2.x1 <- sum(x1^2) #squared norm of x1
  norm2.x2 <- sum(x2^2) #squared norm of x2
  d1 <- scalar.prod / (norm2.x1 + norm2.x2 - scalar.prod)
  K <- par[1]*d1
  attr(K, "gradient") <- c(sigma2 <- d1) #gradient w.r.t. sigma^2
  return(K)    }

# fit GP model
fitGP <- function(formula, train, cov, obs_noise, obs_varNoiseIni, obs_varNoiseLower){
  fit <- try(gp(formula=formula, data=train, cov=cov, 
                noise=obs_noise,
                varNoiseLower=obs_varNoiseLower, varNoiseUpper=obs_varNoiseUpper,
                optimFun = "stats::optim",
                optimMethod="L-BFGS-B", checkNames=FALSE))
  if(is.character(fit)){  # in case fitting does not converge
    fit <- try(gp(formula = formula, 
                  data = train, cov=cov, noise=obs_noise,
                  varNoiseLower=obs_varNoiseLower, varNoiseUpper=obs_varNoiseUpper,
                  optimFun = "nloptr::nloptr",
                  opts=list(xrel_tol=1e-4)))  }
  return(fit)  }

# twcrps evaluation #############################################################

# twCRPS1
twCRPS1_e <- function(response, pred_mean, pred_sd, T){
  y_tildeT <- (pmax(response,T)-pred_mean)/pred_sd 
  T_tilde <- (T-pred_mean)/(pred_sd )
  res <- mean(pred_sd*
                                (-T_tilde*pnorm(T_tilde)^2 
                                 + y_tildeT*(2*pnorm(y_tildeT)-1)
                                 + (2*dnorm(y_tildeT) - 2*dnorm(T_tilde)*pnorm(T_tilde))
                                 - 1/sqrt(pi)*(1-pnorm(T_tilde*sqrt(2)))))
  return(res)  }

# twCRPS2
twCRPS2_e <- function(response, pred_mean, pred_sd, T, k1=33){
  value <- NA
    for(inn in 1:length(pred_mean)){
      
      m <- c(pred_mean[inn] - T, pred_mean[inn] - T)
      C <- matrix(c(pred_sd[inn]^2 + k1^2, k1^2, k1^2, pred_sd[inn]^2 + k1^2),2,2)
      value1 <- pmvnorm(upper=c(0,0),mean=m,sigma=C)
      
      m <- c(pred_mean[inn] - T, response[inn] - T)
      C <- matrix(c(pred_sd[inn]^2 + k1^2, k1^2, k1^2, k1^2),2,2)
      value2 <- pmvnorm(upper=c(0,0),mean=m,sigma=C)
      
      value[inn] <- value1 - 2* value2 + pnorm((T - response[inn])/k1)   }
    return(value) }

# twcrps acquisition ############################################################

twCRPS1_a <- function(pred_mean, pred_sd, T){
  T_tilde <- (T-pred_mean)/pred_sd
  crit_val <- pred_sd*
    (T_tilde*pnorm(T_tilde)^2 - T_tilde*pnorm(T_tilde) - 1/sqrt(pi)*pnorm(T_tilde*sqrt(2)) +
       2*dnorm(T_tilde)*pnorm(T_tilde) - dnorm(T_tilde) + 1/sqrt(pi))
  return(crit_val) }

twCRPS2_a <- function(pred_mean, pred_sd, T, k1=33){
  crit1 <- NA
  for(inn in 1:length(pred_mean)){
    m <- c(pred_mean[inn] - T, T - pred_mean[inn])
    C <- matrix(c(pred_sd[inn]^2 + k1^2, -k1^2, -k1^2, pred_sd[inn]^2 + k1^2),2,2)
    crit1[inn] <- pmvnorm(upper=c(0,0),mean=m,sigma=C)   }
  return(crit1) }

twCRPS1_aI <- function(pred_mean, pred_sd, pred_cov, pred_tau2, T){
  crit1 <- NaN
  zz <- rnorm(10000)
  zzz <- rnorm(10000)
  zzzz <- rnorm(10000)
  for(pp in 1:length(pred_mean)){  # which point to add
    unc <- NaN
    for(xin in 1:length(pred_mean)){  # uncertainty in molecule xin
      
      sigma2_n1 <- pred_sd[xin]^2 - pred_cov[pp,xin]^2/(pred_sd[pp]^2 + pred_tau2)
      alpha <- pred_cov[pp,xin]/sqrt(pred_sd[pp]^2+pred_tau2)
      mn <- pred_mean[xin]
      
      unc[xin] <- sqrt(sigma2_n1)*mean(pmax(zz - pmax(zzz,(T-mn-alpha*zzzz)/sqrt(sigma2_n1)),0))  }
    
    crit1[pp] <- mean(unc) }
    return(crit1)  }

twCRPS2_aI <- function(pred_mean, pred_sd, pred_cov, pred_tau2, T, k1=33){
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
  return(crit1) }


# function evolving from initial design ###########################################

GPevolve <- function(sequential, test_ind, data, nr_seq, k1 = 33){
  
  set.seed(test_ind) # make sure to use same seed for the repetitions across the strategies
  
  # take out validation set for final testing 
  shuffle_index <- sample(seq(nrow(data)), replace=FALSE)
  first_i <- shuffle_index[1:nr_validate]
  test_fin <- data[first_i, ]
  
  dataa <- data
  data <- data[-first_i, ]
  dataa$validate <- '0'
  dataa$validate[first_i] <- '1'
  
  # initial training data
  shuffle_index <- sample(seq(nrow(data)), replace=FALSE)
  first_i <- shuffle_index[1:nr_train]
  train <- data[first_i, ]
  test <- data[-first_i, ]
  
  dataa$init <- NA
  dataa$init[which(dataa$response %in% train$response)] <- '1'
  dataa$pred_init <- NA
  
  # set up result frame
  results <- data.frame(row.names = 1:nr_seq)
  results$RMSE <- NA
  results$CRPS <- NA
  results$twCRPS1 <- NA
  results$twCRPS2 <- NA 
  if(vers == 'syn'){
    results$RMSElevel <- NA
    results$RMSElevel1 <- NA
    results$TP <- NA
    results$TN <- NA
    results$FP <- NA
    results$FN <- NA    }
  
  dataa$added <- NA
  for(it in 1:nr_seq){
    
    # fit the GP model
    fit <- fitGP(formula, train, cov, obs_noise, obs_varNoiseIni, obs_varNoiseLower)
    if(is.character(fit)){ stop("GP fitting failed. Stopping script.") }
    
    # evaluation validation set #
    pred_results <- predict(object=fit, newdata = test_fin)
    pred_mean <- pred_results$mean
    pred_sd <- pred_results$sd 
    
    if(it==1){
      for(qq in 1:length(pred_mean)){
        dataa$pred_init[which(dataa$response_clean == test_fin$response_clean[qq])] <- pred_mean[qq]  }
    }
    
    # evaluation original (noisy)
    if(vers == 'ori'){
      pred_sd <- pred_sd + sqrt(fit$varNoise)
      results$RMSE[it] <- sqrt(mean((test_fin$response - pred_results$mean)^2))
      results$CRPS[it] <- mean(crps_norm(test_fin$response, pred_results$mean, pred_sd))
      results$twCRPS1[it] <- mean(twCRPS1_e(test_fin$response, pred_mean, pred_sd, T))
      results$twCRPS2[it] <- mean(twCRPS2_e(test_fin$response, pred_mean, pred_sd, T, k1))    }
    
    # evaluation synthetic
    if(vers == 'syn'){
      results$RMSE[it] <- sqrt(mean((test_fin$response_clean - pred_results$mean)^2))
      results$CRPS[it] <- mean(crps_norm(test_fin$response_clean, pred_results$mean, pred_sd))
      results$twCRPS1[it] <- mean(twCRPS1_e(test_fin$response_clean, pred_mean, pred_sd, T))
      results$twCRPS2[it] <- mean(twCRPS2_e(test_fin$response_clean, pred_mean, pred_sd, T, k1=33))
      
      true_class <- 1*(test_fin$response_clean >= T)
      exc_prob <- pnorm((pred_mean - T)/pred_sd)
      pred_class <- 1*(exc_prob >= 0.5)
      results$TP[it] <- sum(true_class*pred_class)
      results$TN[it] <- sum((1-true_class)*(1-pred_class))
      results$FP[it] <- sum((1-true_class)*(pred_class))
      results$FN[it] <- sum((true_class)*(1-pred_class))
      
      results$RMSElevel[it] <- ifelse(sum(pred_class) == 0, 0,
                                      sqrt(1/sum(pred_class)*sum(pred_class*(test_fin$response_clean - pred_mean)^2)))
      results$RMSElevel1[it] <- sqrt(1/sum(true_class)*sum(true_class*(test_fin$response_clean - pred_mean)^2))   }
    
    # sequential data acquisition
    # predict for candidate set
    pred_results <- predict(object=fit, newdata = test)
    pred_mean <- pred_results$mean
    pred_sd <- pred_results$sd
    pred_tau2 <- fit$varNoise
    exc_prob <- pnorm((pred_mean - T)/pred_sd)
    if(sequential %in% c('IBV', 'TIMSE', 'ItwCRPS1', 'ItwCRPS2')){
      pred_results <- predict(object=fit, newdata = test, covCompute = TRUE)
      pred_cov <- pred_results$cov    }
    
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
      
      # margin-based
    if(sequential == 'Margin'){
      crit_val <- abs(2*exc_prob - 1)
      new_i <- which.min(crit_val)      }
    
    # twcrps1 pointwise
    if(sequential == 'twCRPS1'){
      crit_val <- twCRPS1_a(pred_mean, pred_sd, T)
      new_i <- which.max(crit_val)       }
    
    # twcrps2 pointwise
    if(sequential == 'twCRPS2'){
      crit_val <- twCRPS2_a(pred_mean, pred_sd, T, k1)
      new_i <- which.max(crit_val)       }
    
    # twcrps1 SUR
    if(sequential == 'ItwCRPS1'){
      crit_val <- twCRPS1_aI(pred_mean, pred_sd, pred_cov, pred_tau2, T)
      new_i <- which.min(crit_val)       }
    
    # twcrps2 SUR
    if(sequential == 'ItwCRPS2'){
      crit_val <- twCRPS2_aI(pred_mean, pred_sd, pred_cov, pred_tau2, T, k1)
      new_i <- which.min(crit_val)       }
    
    # IBV
    if(sequential == 'IBV'){
      sur_vec <- NaN
      for(pp in 1:nrow(test)){ 
        sd2_upd <- pmax(pred_sd^2 - 1/(pred_sd[pp]^2+pred_tau2)*c(pred_cov[pp,])*c(pred_cov[pp,]),1e-10)
        a <- (pred_mean[-pp]-T)/sqrt(sd2_upd[-pp])
        c <- pred_sd[-pp]^2/sd2_upd[-pp]
        int <- NaN
        int <- sapply(seq((nrow(test)-1)), function(j){
          pmvnorm(lower=c(-Inf,-Inf), upper=c(0,0), mean=c(a[j],-a[j]), sigma=matrix(c(c[j], 1-c[j], 1-c[j], c[j]), nrow=2, byrow=TRUE))[1] })
        sur_vec[pp] <- mean(int)  }
      new_i <- which.min(sur_vec)  }
    
    # TIMSE
    if(sequential == 'TIMSE'){
      weight <- dnorm(T,pred_mean,pred_sd)
      timse_vec <- NaN
      for(pp in 1:nrow(test)){ 
        sd2_upd <- pred_sd^2 - 1/(pred_sd[pp]^2+pred_tau2)*c(pred_cov[pp,])*c(pred_cov[pp,])
        timse_vec[pp] <- mean(sd2_upd*weight)}
      new_i <- which.min(timse_vec)  }
    
    dataa$added[which(dataa$response == test[new_i, ]$response)] <- '1'
    
    train <- rbind(train, test[new_i, ])
    test <- test[-c(new_i), ] 
  }
  
  # fit the GP model
  fit <- fitGP(formula, train, cov, obs_noise, obs_varNoiseIni, obs_varNoiseLower)
  if(is.character(fit)){ stop("GP fitting failed. Stopping script.") }
  
  # evaluation validation set #
  pred_results <- predict(object=fit, newdata = test_fin)
  pred_mean <- pred_results$mean
  pred_sd <- pred_results$sd 
  
  dataa$pred <- NA
  for(qq in 1:length(pred_mean)){
    dataa$pred[which(dataa$response_clean == test_fin$response_clean[qq])] <- pred_mean[qq]  }

  return(dataa)
}


# run the acquisition ###########################################################

vers <- 'syn'  # syn' or 'ori' dataset

# load data
if(vers == 'syn'){data <- read.csv("Data/data_syn.csv")
  data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
  colnames(data) <- c(paste0 ("x", seq(ncol(data)-2)), 'response', 'response_clean')
  nvar <- ncol(data)-2 }
if(vers == 'ori'){data <- read.csv("Data/data.csv")
  data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
  colnames(data) <- c(paste0 ("x", seq(ncol(data)-1)), 'response')
  nvar <- ncol(data)-1 }
  ndata <- nrow(data)

# specifications
nr_validate <- 100 # validation molecules 
nr_tests <- 1 # repetitions
pppr <- 0.8 # for the threshold
T <- quantile(data$response, probs=pppr) # threshold excursion set
nr_train <- 30 # initial training molecules
nr_seq <- 25 # added molecules

# ordinary kriging with noise
formula <- as.formula("response~1") 
# obs noise fixed (syn) or estimated (ori)
if(vers == 'syn'){
  obs_noise <- TRUE
  obs_varNoiseLower <- 77.46
  obs_varNoiseUpper <- 77.46 
  obs_varNoiseIni <- 77.46 }
if(vers == 'ori'){
  obs_noise <- TRUE
  obs_varNoiseLower <- 1e-10
  obs_varNoiseUpper <- 100 
  obs_varNoiseIni <- 20}
  
# Tanimoto kernel
cov <- covMan(kernel = ker.tanimoto, hasGrad=TRUE, d=nvar,
                       label = "Tanimoto kernel", parNames=c("sigma2"), 
                       par = c(sigma2=100), parLower =c(sigma2=0), parUpper=c(sigma2=Inf))

sequential  = 'twCRPS1' 
dataa1 <- GPevolve(sequential, test_ind = 45, data, nr_seq, k1 = 33)
dataa1 <- dataa1 %>% arrange(response_clean)
sequential  = 'twCRPS2' 
dataa2 <- GPevolve(sequential, test_ind = 45, data, nr_seq, k1 = 33)
dataa2 <- dataa2 %>% arrange(response_clean)
ii <- min(which(dataa$response_clean>=T))

p1 <- ggplot(dataa1, aes(x = seq(1, nrow(dataa1), 1), y = response_clean)) +
  geom_point(color = 'lightgray', size=0.8) +
  geom_hline(yintercept = T) +
  geom_vline(xintercept = ii, linetype='dashed') +
  geom_point(aes(y = ifelse(init == '1', response, NA)), 
             color = 'darkred', shape = 17, size = 2, alpha=0.8) +   
  geom_point(aes(y = ifelse(added == '1', response, NA)), 
             color = 'darkorange', shape = 15, size = 2, , alpha=.8) + 
  geom_point(aes(y = pred), col = 'cadetblue', size=2.5, shape='X', alpha=1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = expression(CRPS[nu[1]])) + xlab('Molecule index \n (b)') + ylab('Response')

p2 <- ggplot(dataa2, aes(x = seq(1, nrow(dataa2), 1), y = response_clean)) +
  geom_point(color = 'lightgray', size=0.8) +
  geom_hline(yintercept = T) +
  geom_vline(xintercept = ii, linetype='dashed') +
  geom_point(aes(y = ifelse(init == '1', response, NA)), 
             color = 'darkred', shape = 17, size = 2, alpha=0.8) +  
  geom_point(aes(y = ifelse(added == '1', response, NA)), 
             color = 'darkorange', shape = 15, size = 2, , alpha=.8) + 
  geom_point(aes(y = pred), col = 'cadetblue', size=2.5, shape='X', alpha=1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = expression(CRPS[nu[2]])) + xlab('Molecule index \n (c)') + ylab('Response')

p3 <- ggplot(dataa2, aes(x = seq(1, nrow(dataa2), 1), y = response_clean)) +
  geom_point(color = 'lightgray', , size=0.8) +
  geom_hline(yintercept = T) +
  geom_vline(xintercept = ii, linetype='dashed') +
  geom_point(aes(y = ifelse(init == '1', response, NA)), 
             color = 'darkred', shape = 17, size = 2, alpha=0.8) +   
  
  geom_point(aes(y = pred_init), col = 'cadetblue', size=2.5, shape='X', alpha=1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = 'Initial state') + xlab('Molecule index \n (a)') + ylab('Response')

setwd("~/GitHub/GP_sequentialCRPS")
jpeg(file=paste0("twcrps_illu_photo.jpeg"), width=17, height=7, units="cm", res=300)
ggarrange(p3,p1,p2, ncol = 3, nrow = 1, heights=c(0.9,1))
dev.off()
 