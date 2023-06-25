#! /usr/bin/env Rscript
setwd("high dimension")
Rcpp::sourceCpp("vestimator.cpp")
source("estimators.R")
source("Model1.R")
source("Model2.R")
source("Model3.R")
source("Model4.R")
source("Model5.R")
source("Model6.R")
source("Model7.R")
source("Model8.R")
source("tables.R")
source("DML.R")
source("parallel.R")
RNGkind("L'Ecuyer-CMRG")
set.seed(202209)
mc.reset.stream()

n = 1000
lambda = 0.75
weight = 1
omega = c(0.1,0.1,0.4,0.4)
sigma1 = 3
sigma0 = 1
pi = 1/2
Iternum = 2000
p = 200
#name_methods = c("rf","nn","rpart","lasso","enet","gbm","ensemble","best")

# Model 5
Model = "Model5"
alpha0 = 1
alpha1 = 4
betavec0 = c(75,35,125,80)
betavec1 = c(100,80,60,40)
re5 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv5 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re5,tv5,file = "0909_re5.RData")

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
re6 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv6 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re6,tv6,file = "0909_re6.RData")

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
re7 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv7 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re7,tv7,file = "0909_re7.RData")

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
re8 = mclapply(1:Iternum,sim_DML,mc.cores = 26)
tv8 = trueval(Model,alpha0,alpha1,betavec0,betavec1)
save(re8,tv8,file = "0909_re8.RData")
