#! /usr/bin/env Rscript
setwd("~/semi/202209/0907")
library(microbenchmark)
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
t5<-microbenchmark(sim_DML(1),times = 1)

# Model 6
Model = "Model6"
alpha0 = -3
alpha1 = 0
betavec0 = c(10,24,15,20)
betavec1 = c(20,27,10)
t6<-microbenchmark(sim_DML(1),times = 1)

#Model 7
Model = "Model7"
alpha0 = 5
alpha1 = 2
betavec0 = c(42,83)
betavec1 = c(30,75)
t7<-microbenchmark(sim_DML(1),times = 1)

# Model 8
Model = "Model8"
alpha0 = 5
alpha1 = 5
betavec0 = c(20,30,50)
betavec1 = c(20,30,65)
t8<-microbenchmark(sim_DML(1),times = 1)

save(t5,t6,t7,t8,file = "time_DML.RData")