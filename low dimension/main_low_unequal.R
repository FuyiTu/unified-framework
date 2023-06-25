#! /usr/bin/env Rscript
setwd("low dimension")
Rcpp::sourceCpp("vestimator.cpp")
source("estimators.R")
source("Model1.R")
source("Model2.R")
source("Model3.R")
source("Model4.R")
source("tables.R")
source("parallel.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(202301)
mc.reset.stream()

n = 1000
lambda = 0.75
weight = 1
omega = c(0.1,0.1,0.4,0.4)
sigma1 = 3
sigma0 = 1
pi = 2/3
Iternum = 2000

# Model 1
Model = "Model1"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re1 = mclapply(1:Iternum,sim_low,mc.cores = 26)
tv1 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re1,tv1,file = "0624_re123.RData")

# Model 2
Model = "Model2"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re2 = mclapply(1:Iternum,sim_low,mc.cores = 26)
tv2 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re2,tv2,file = "0624_re223.RData")

#Model 3
Model = "Model3"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re3 = mclapply(1:Iternum,sim_low,mc.cores = 26)
tv3 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re3,tv3,file = "0624_re323.RData")

# Model 4
Model = "Model4"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re4 = mclapply(1:Iternum,sim_low,mc.cores = 26)
tv4 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re4,tv4,file = "0624_re423.RData")
