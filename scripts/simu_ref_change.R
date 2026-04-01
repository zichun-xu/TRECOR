### Simulation LTNM high-dimensional sparse covariance regression model
### This is a small experiment to validate whether changing the reference group changes the result
rm(list = ls())
getwd()
library(matrixcalc)
library(dplyr)
library(igraph)
library(mvtnorm)
library(statmod)
library(stats)
library(BayesLogit)
library(matrixNormal)
library(GeneralizedHyperbolic)
source("/home/dxu/Research_projects/tree_cov_reg/R/LTNM_gen.R")
source("/home/dxu/Research_projects/tree_cov_reg/R/Gibbs_R_LTNM.R")
source("/home/dxu/Research_projects/tree_cov_reg/R/helper.R")
touse = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(20241015+touse)

### Global parameters
niter = 10000
burnin = niter/2
result_folder = "/home/dxu/Research_projects/tree_cov_reg/results/output/simu/ref/"

Theta = readRDS(paste0(result_folder,"Theta"))
p = dim(Theta)[1]
graph_true = as.numeric(Theta != 0) %>% matrix(nrow = p)
graph_true_vec = graph_true[upper.tri(graph_true)]
X = readRDS(paste0(result_folder,"X"))
B = readRDS(paste0(result_folder,"B"))
R = dim(B)[3]
seq_depth = readRDS(paste0(result_folder, "seq_depth"))
phy_tree = readRDS(paste0(result_folder, "phy_tree"))
n = dim(X)[1]
d = dim(X)[2]


mu = numeric(p)
tryCatch(
  {mu = readRDS(paste0(result_folder,"mu"))},
  error = function(e){
    message("Set default to 0")
  }
)

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
        Phi_est[i,j]=-5
      }
      else{
        if (tmp==1){
          Phi_est[i,j]=5
        }
        else{
          Phi_est[i,j]=stats::qlogis(tmp)
        }
      }
    }
  }
}

### Changed reference design matrix
X_changed = matrix(0,n,d)
X_changed[,d] = numeric(n) + 1
X_changed[,1] = as.numeric(X[,1] == 0 & X[,2] == 0)
X_changed[,2] = X[,2]
colnames(X_changed) = c("A", "C", "B_reference")

### Initialization
Theta_intial = 0.5*cov(Phi_est)+0.5*diag(p)
B_intial = array(0, dim = c(p,d,R))
B0_intial = matrix(0,p,d)

### Save results
result = GibbsLTNMRCovMeanRegGL(Y = Y, N = N, X = X, R = R, niter = niter,
                                Theta = Theta_intial, B = B_intial, B0 = B0_intial)
result_changed = GibbsLTNMRCovMeanRegGL(Y = Y, N = N, X = X_changed, R = R, niter = niter,
                                        Theta = Theta_intial, B = B_intial, B0 = B0_intial)
save(result, result_changed,
        file = paste0(result_folder, "output/","LTNM_R_",R,"_cov_reg_simulation_results_", touse, sep=""))

print("Done!")