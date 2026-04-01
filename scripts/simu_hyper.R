source("~/Research_projects/tree_cov_reg/config/config_simu_hyper.R")

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
load(paste0(data_path, "params_", paste("n",n, "p",p, "sparse", sparse, Theta_structure, sep="_")))
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

### Save results

# LTNM covariance regression
parameters = data.frame(ac = c(10,10,5), 
                        ad = c(10,10,5), 
                        bc = c(1,0.5,1), 
                        bd = c(1,0.5,1))
ranks = sample(ranks, size = length(ranks), replace = FALSE)
for(par_set in 1:nrow(parameters)){
  prefix = paste("ac",parameters$ac[par_set], "ad",parameters$ad[par_set], "bc",parameters$bc[par_set], "bd",parameters$bd[par_set], sep="_")
  for(rank in ranks){
    set.seed(202408181)
    Theta_init = 0.5*cov(Phi_est)+0.5*diag(p) # The initialization is changed to handle high-dimensional p
    B_init = array(rnorm(p*d*rank), dim = c(p,d,rank))
    B0_init = matrix(rnorm(p*d),p,d)
    
    
    print(paste0("Running: rank-",rank," LTNM covariance regression with GL prior"))
    result = GibbsLTNMRCovMeanRegGL(Y = Y, N = N, X = X, R = rank, niter = niter,
                                    Theta = Theta_init, B = B_init, B0 = B0_init,
                                    a_c = parameters$ac[par_set], a_d = parameters$ad[par_set], b_c = parameters$bc[par_set], b_d = parameters$bd[par_set])
    
    date_of_run = Sys.time()
    session_info = devtools::session_info()
    
    # evaluation metrics: posterior means, posterior MSEs, WAIC, and ROC
    Theta_pmean = matrix(0,p,p)
    B_pmean = array(0, dim=c(p,d,rank))
    B0_pmean = matrix(0,p,d)
    TFP_pmean = matrix(0,3,1000)
    TPR <- FPR <- PPV <- matrix(NA, nrow = burnin, ncol = 1000)
    Sigma_pmean = array(0, dim = c(p,p,n))
    lppd <- pwaic2 <- 0
    mse_Theta <- mse_Sigma <- numeric(niter-burnin) + NA
    
    # calculate these metrics
    for(iter in (burnin+1):niter){
      
      # ROC
      graph_est = cov2cor(result[[1]][[iter]])[upper.tri(Theta)]
      for(j in 1:1000){
        a = (j-1)/1000
        graph_binary = (abs(graph_est) >= a)%>%as.numeric()
        TPR[iter-burnin,j] = sum(graph_true_vec == 1 & graph_binary == 1)/sum(graph_true_vec == 1)
        FPR[iter-burnin,j] = sum(graph_true_vec == 0 & graph_binary == 1)/sum(graph_true_vec == 0)
        PPV[iter-burnin,j] = sum(graph_true_vec == 1 & graph_binary == 1)/sum(graph_binary == 1)
      }
      
      # posterior means
      Theta_pmean = Theta_pmean + result[[1]][[iter]]/burnin
      B_pmean = B_pmean + result[[2]][[iter]]/burnin
      B0_pmean = B0_pmean + result[[3]][[iter]]/burnin
      
      Sigma_est = array(NA, c(p,p,n))
      for(j in 1:n){
        Sigma_est[,,j] = result[[1]][[iter]]
        for(r in 1:rank){
          Sigma_est[,,j] = Sigma_est[,,j]+result[[2]][[iter]][,,r]%*%X[j,]%*%t(X[j,])%*%t(result[[2]][[iter]][,,r])
        }
      }
      Sigma_pmean = Sigma_pmean + Sigma_est/burnin
      
      # posterior MSEs
      mse_Theta[iter-burnin] = norm(Theta -  result[[1]][[iter]])^2
      mse_Sigma[iter-burnin] = sum((Sigma - Sigma_est)^2)/n
    }
    TFP_pmean[1,] = apply(TPR,2,mean, na.rm = T)
    TFP_pmean[2,] = apply(FPR,2,mean, na.rm = T)
    TFP_pmean[3,] = apply(PPV,2,mean, na.rm = T)
    
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
    
    note = paste0("Simulation: rank", rank, " LTNM covariance matrix regression with the graphical LASSO prior")
    save(TFP_pmean, Theta_pmean, B_pmean, B0_pmean, Sigma_pmean,
         mse_Sigma, mse_Theta, 
         WAIC, 
         date_of_run, session_info, note, 
         file = paste0(result_path, prefix, "_LTNM_R_",rank,"_cov_reg_simulation_results_", touse))
  }
}

print("Done!")