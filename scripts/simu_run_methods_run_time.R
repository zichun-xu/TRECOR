# config.R for simulations
path = "~/Research_projects/tree_cov_reg/"
data_path = "~/Research_projects/tree_cov_reg/data/synthetic/"
result_path = "~/Research_projects/tree_cov_reg/results/time/"
niter = 8000
burnin = niter/2

n = 500
p = 25
d = 3
R = 3
rank = 1
sparse = 0.5
Theta_structure  = "tree"

prefix = paste("n",n, "p",p, "sparse", sparse, Theta_structure, sep="_")

# Ensure renv is available
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}


# Load the project-local environment
setwd(path)

# Run the algorithm
library(igraph)
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
source("R/Gibbs_R_cov_mean_reg.R")
source("R/Gibbs_cov.R")
load(paste0(data_path, "params_", prefix))
touse = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(20241015+touse*10)

graph_true = as.numeric(Theta != 0) %>% matrix(nrow = p)
graph_true_vec = graph_true[upper.tri(graph_true)]
R = dim(B)[3]

Sigma <- array(NA, dim = c(p,p,n))
for(i in 1:n){
  Sigma[,,i] = Theta
  for(r in 1:R){
    Sigma[,,i] = Sigma[,,i] + B[,,r]%*%X[i,]%*%t(X[i,])%*%t(B[,,r])
  }
}

### Generate data
LTNM_data = LTNMGen(X = X, B = B, Theta = Theta, phy_tree = phy_tree,
                    seq_depth = seq_depth, mu = mu)
Y = LTNM_data$Y
N = LTNM_data$N
Phi_true = LTNM_data$Phi
Phi_est = matrix(NA, nrow=n, ncol=p)
for (i in 1:n){
  for (j in 1:p){
    if (N[i,j]==0){
      Phi_est[i,j]=0
    }
    else{
      tmp=Y[i,j]/N[i,j]
      if (tmp==0){
        Phi_est[i,j]=-15
      }
      else{
        if (tmp==1){
          Phi_est[i,j]=15
        }
        else{
          Phi_est[i,j]=stats::qlogis(tmp)
        }
      }
    }
  }
}

### Run
set.seed(202408181)
Theta_init = 0.5*cov(Phi_est)+0.5*diag(p) # The initialization is changed to handle high-dimensional p
B_init = array(rnorm(p*d*rank), dim = c(p,d,rank))
B0_init = matrix(rnorm(p*d),p,d)
print(paste0("Running: rank-",rank," LTNM covariance regression with GL prior"))

start_time <- Sys.time()
result = GibbsLTNMRCovMeanRegGL(Y = Y, N = N, X = X, R = rank, niter = niter,
                                Theta = Theta_init, B = B_init, B0 = B0_init)
end_time <- Sys.time()
time = end_time - start_time
save(time, 
     file = paste0(result_path, prefix, "_LTNM_R_",rank,"_cov_reg_simulation_time_", touse))
