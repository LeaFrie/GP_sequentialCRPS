# generate synthetic version of the dataset ####################################
set.seed(15)

# install and load packages
packages <- c("kergp", "mvtnorm", "tidyverse","ggplot2")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)}}
lapply(packages, library, character.only = TRUE)

# load original dataset
# https://github.com/Ryan-Rhys/The-Photoswitch-Dataset
data <- read.csv("data.csv", row.names=1)
data[, 1:ncol(data)] <- lapply(data[, 1:ncol(data)], as.numeric)
colnames(data) <- c(paste0 ("x", seq(ncol(data)-1)), 'response')
nvar <- ncol(data)-1
ndata <- nrow(data)

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
cov.tanimoto <- covMan(kernel = ker.tanimoto, hasGrad=TRUE, d=nvar,
                       label = "Tanimoto kernel", parNames=c("sigma2"), 
                       par = c(sigma2=100), parLower =c(sigma2=0), parUpper=c(sigma2=Inf))

# ordinary kriging with noise
formula <- as.formula("response~1") 
obs_noise <- TRUE
obs_varNoiseIni <- 20
obs_varNoiseLower <- 1e-10

# fit a GP to the full dataset
fit <- gp(formula=formula, data=data, cov=cov.tanimoto, 
          noise=obs_noise,
          varNoiseLower=obs_varNoiseLower, varNoiseIni=obs_varNoiseIni,
          optimFun = "stats::optim",
          optimMethod="L-BFGS-B", checkNames=FALSE)

# use predictive mean as ground truth, add noise
pred <- predict(object=fit, newdata = data)
data_syn <- data
data_syn$response_clean <- pred$mean
data_syn$response <- data_syn$response_clean + sqrt(fit$varNoise)*rnorm(nrow(data_syn))

# save synthetic dataset
write.csv(data_syn, file="data_syn.csv", row.names = FALSE)   
