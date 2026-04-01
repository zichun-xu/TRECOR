source("~/Research_projects/tree_cov_reg/config/config_real_data.R")

# Ensure renv is available
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}


# Load the project-local environment
setwd(path)

# Run the algorithm
library(matrixcalc)
library(dplyr)
library(mvtnorm)
library(statmod)
library(stats)
library(BayesLogit)
library(matrixNormal)
library(GeneralizedHyperbolic)
source("R/helper.R")
source("R/LTNM_gen.R")
source("R/Gibbs_R_LTNM.R")

load(paste0(data_path, "processed_obj"))
n = nrow(Y)
p = ncol(Y)
d = ncol(X)

touse = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
for(rank in ranks){
  print(paste0("Running: rank-",rank," LTNM covariance regression with GL prior"))
  # result = GibbsRCovMeanRegGL(Phi = Phi_est, X = X, R = rank, niter = niter)
  
  # Initialization
  if(chain == 1){
    set.seed(202408181)
    Theta = 0.5*cov(Phi_est)+0.5*diag(p) # The initialization is changed to handle high-dimensional p
    B = array(rnorm(p*d*rank), dim = c(p,d,rank))
    B0 = matrix(rnorm(p*d),p,d)
    set.seed(NULL)
  } else if(chain == 2){
    set.seed(202408182)
    Theta = diag(diag(cov(Phi_est)))
    B = array(rnorm(p*d*rank), dim = c(p,d,rank))
    B0 = matrix(rnorm(p*d),p,d)
    set.seed(NULL)
  }else if(chain == 3){
    set.seed(202408183)
    A = matrix(rnorm(p*p),p,p)
    Theta = diag(sqrt(diag(cov(Phi_est))))%*%cov2cor(A%*%t(A))%*%diag(sqrt(diag(cov(Phi_est))))
    B = array(rnorm(p*d*rank), dim = c(p,d,rank))
    B0 = matrix(rnorm(p*d),p,d)
    set.seed(NULL)
  }
  
  result = GibbsLTNMRCovMeanRegGL(Y = Y, N = N, X = X, R = rank, niter = niter,
                                  Theta = Theta, B = B, B0 = B0)
  lppd <- pwaic2 <- 0
  for( i in 1:n ){
    pred_density = numeric(niter-burnin)
    for( iter in (burnin+1):niter ){
      pred_density[iter-burnin] = DensityLTNM(Theta = result[[1]][[iter]], 
                                              B = result[[2]][[iter]], 
                                              Phii = result[[4]][[iter]][i,],
                                              Xi = X[i,], Yi = Y[i,], Ni = N[i,],
                                              B0 = result[[3]][[iter]])
    }
    lppd = lppd + log(mean(pred_density))
    pwaic2 = pwaic2 + var(log(pred_density)[pred_density != 0])
  }
  WAIC = -1*lppd + pwaic2
  
  date_of_run = Sys.time()
  note = paste0("Real data: rank", rank, " LTNM covariance matrix regression with the graphical LASSO prior")
  
  save(result, WAIC, 
       date_of_run, note, 
       file = paste0(result_path, "LTNM_R_",rank,"_cov_reg_chain_", chain, "_iter_", touse, sep=""))
}
