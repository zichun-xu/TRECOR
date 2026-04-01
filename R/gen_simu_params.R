# Extract simulation parameters from real data

library(phyloseq)
library(microbiome)
library(igraph)
library(mvtnorm)
library(matrixcalc)
source("~/Research_projects/tree_cov_reg/R/LTNM_gen.R")
source("~/Research_projects/tree_cov_reg/R/helper.R")

load("~/Research_projects/tree_cov_reg/data/raw/Yatsunenko_QIIME_processed.rda")

### Simulation parameters
n = 150
p = 100
d = 3
R = 3
sparse = 0.5 # remove zero probability > sparse, larger means more sparse
Theta_structure  = "tree"
prefix = paste("n",n, "p",p, "sparse", sparse, Theta_structure, sep="_")

### Real data QC
summarize_phyloseq(subset110_genus)

meta = subset110_genus@sam_data
taxonomy = subset110_genus@tax_table
otu = subset110_genus@otu_table %>% t()

# Check proportion of zeros
apply(otu, 2, function(x){mean(x==0)}) %>% summary()

# Check the total number of reads
colSums(otu) %>% summary()

remove_ind = which(apply(otu, 2, function(x){mean(x==0)}) > sparse | colSums(otu) < 100)
data_processed = otu[,-1*remove_ind]

# Check the mean abudance
data_processed_comp = data_processed/rowSums(data_processed)
hist(colMeans(data_processed_comp), breaks = 50)

# The genus 338437 is extremely abundant, should we keep it?
data_processed = data_processed[,-142]

### Randomly sample 100 genus to estimate parameters in simulation
set.seed(202608) # use seed 202504 to replicate for sparsity = 0.3,0.5
used_ind = sample(1:ncol(data_processed), size = p+1, replace = ifelse(ncol(data_processed) <= p+1, TRUE, FALSE))
data_sampled = data_processed[,used_ind]
seq_depth = rowSums(data_sampled)
summary(seq_depth)
hist(seq_depth)
data_sampled_comp = data_sampled/seq_depth
hist(colMeans(data_sampled_comp), breaks = 20)

### Generate parameters
set.seed(202504)
seq_depth = sample(seq_depth, size = n, replace = TRUE)
names(seq_depth) = NULL

# Add pseudo-number by 1
# data_sampled_imputed = data_sampled
# data_sampled_imputed[data_sampled == 0] = 1

set.seed(202504)
phy_tree = make_tree(2*p+1, children = 2, mode = "out")
Y <- N <- matrix(NA, nrow = nrow(data_sampled), ncol = p)
for( i in 1:nrow(data_sampled) ){
  result = CountLeaf2Int(adj = as_adjacency_matrix(phy_tree), Y = data_sampled[i,])
  Y[i,] = result$Y
  N[i,] = result$N
}
# mu = log(colMeans(Y/N, na.rm=T)/(1-colMeans(Y/N, na.rm=T)))
Phi_est = log(Y/(N-Y))
Phi_est[N==0] = NA
Phi_est[Phi_est == Inf] = 15
Phi_est[Phi_est == -Inf] = -15
mu = colMeans(Phi_est, na.rm=T)
hist(mu, breaks = 20)
summary(mu)

### Generate other parameters
set.seed(20241112)
X = rmvnorm(n = n, sigma = diag(d)*0.9 + matrix(0.1, d, d) )
X[,1] = numeric(n)+1

# B different ranks are non-overlapping
set.seed(20241112)
B_mean <- B_sd <- sqrt(1/2/d)
B1 = rnorm(n=as.integer(p*0.2)*d, mean = B_mean, sd = B_sd) %>% matrix(nrow = as.integer(p*0.2), ncol = d)
B1 = rbind(B1, matrix(0, nrow = as.integer(p*0.8), ncol = d))
B2 = rnorm(n=as.integer(p*0.4)*d, mean = B_mean, sd = B_sd) %>% matrix(nrow = as.integer(p*0.4), ncol = d)
B2 = rbind(matrix(0, nrow = as.integer(p*0.2), ncol = d), B2)
B2 = rbind(B2, matrix(0, nrow = as.integer(p*0.4), ncol = d))
B3 = rnorm(n=as.integer(p*0.4)*d, mean = B_mean, sd = B_sd) %>% matrix(nrow = as.integer(p*0.4), ncol = d)
B3 = rbind(matrix(0, nrow = as.integer(p*0.6), ncol = d), B3)
B = array(NA, dim = c(p,d,R))
B[,,1] = B1
B[,,2] = B2
B[,,3] = B3

# Theta
if(Theta_structure == "off"){
  Theta = diag(p)
  for(i in 2:p){
    j = i-1
    Theta[i,j]<-Theta[j,i]<-0.5
  }  
} else if(Theta_structure == "scale"){
  library(igraph)
  library(matrixcalc)
  set.seed(42)
  g <- sample_pa(n = p, power = 1, m = 1, directed = FALSE)
  plot(g)
  A <- as.matrix(as_adjacency_matrix(g, sparse = FALSE))
  sum(A == 1)/2
  corr_mat <- diag(1, nrow(A))
  set.seed(42)
  non_zero_indices <- which(A == 1, arr.ind = TRUE)
  for (i in 1:nrow(non_zero_indices)) {
    row <- non_zero_indices[i, 1]
    col <- non_zero_indices[i, 2]
    if (row < col) {  # Only assign once due to symmetry
      random_value <- runif(1, min = 0.3, max = 0.9)  # Large random correlations
      random_value <- random_value*sample(c(-1,1), size = 1)
      corr_mat[row, col] <- random_value
      corr_mat[col, row] <- random_value
    }
  }
  is.positive.definite(corr_mat+diag(1,p))
  corr_mat = cov2cor(corr_mat+diag(1,p))
  upper = corr_mat[upper.tri(corr_mat)]
  summary(abs(upper[upper!=0]))
  Theta = corr_mat
} else if(Theta_structure == "tree"){
  dist_mat = igraph::distances(phy_tree)[1:p,1:p]
  prob_graph = exp(-dist_mat)
  A = matrix(0L, p, p)
  set.seed(42)
  for (i in 1:(p - 1)){
    for (j in (i + 1):p){
      A[i,j] = rbinom(1, size = 1, prob = prob_graph[i,j])
      A[j,i] = A[i,j]
    }
  }
  corr_mat <- diag(1, nrow(A))
  set.seed(42)
  non_zero_indices <- which(A == 1, arr.ind = TRUE)
  for (i in 1:nrow(non_zero_indices)){
    row <- non_zero_indices[i, 1]
    col <- non_zero_indices[i, 2]
    if (row < col) {  # Only assign once due to symmetry
      random_value <- runif(1, min = 0.3, max = 0.9)  # Large random correlations
      random_value <- random_value*sample(c(-1,1), size = 1)
      corr_mat[row, col] <- random_value
      corr_mat[col, row] <- random_value
    }
  }
  is.positive.definite(corr_mat+diag(1,p))
  corr_mat = cov2cor(corr_mat+diag(1,p))
  upper = corr_mat[upper.tri(corr_mat)]
  summary(abs(upper[upper!=0]))
  Theta = corr_mat
}


Sigma <- array(NA, dim = c(p,p,n))
for(i in 1:n){
  Sigma[,,i] = Theta
  for(r in 1:R){
    Sigma[,,i] = Sigma[,,i] + B[,,r]%*%X[i,]%*%t(X[i,])%*%t(B[,,r])
  }
  Sigma[,,i] = cov2cor(Sigma[,,i])
}

### Check the off-diagonal entries of the average correlation matrices
average_Sigma = apply(Sigma, c(1,2), mean)
average_Sigma[upper.tri(average_Sigma, diag = F)] %>% hist(breaks = 50)
SST_true = 0
for(i in 1:n){
  SST_true = SST_true+sum((Sigma[,,i]-average_Sigma)^2)
}

### True R-squared
R_squared_true = numeric(R)
for(rank in 1:R){
  Sigma_rank = array(NA, dim = c(p,p,n))
  for(i in 1:n){
    Sigma_rank[,,i] = Theta
    for(r in 1:rank){
      Sigma_rank[,,i] = Sigma_rank[,,i] + B[,,r]%*%X[i,]%*%t(X[i,])%*%t(B[,,r])
    }
    Sigma_rank[,,i] = cov2cor(Sigma_rank[,,i])
  }
  R_squared_true[rank] = 1-sum((Sigma-Sigma_rank)^2)/SST_true
}
print(R_squared_true)

### Check generated data
data_syn = LTNMGen(X = X, B = B, Theta = Theta, phy_tree = phy_tree,
                    seq_depth = seq_depth, mu = mu)
data_syn_internal = data_syn$count
print(mean(data_syn_internal == 0))
Y_syn = data_syn$Y
N_syn = data_syn$N
Phi_est_syn = log(Y_syn/(N_syn-Y_syn))
Phi_est_syn[N_syn==0] = NA
Phi_est_syn[Phi_est_syn == Inf] = 15
Phi_est_syn[Phi_est_syn == -Inf] = -15
mu_syn = colMeans(Phi_est_syn, na.rm=T)
plot(mu, mu_syn)
abline(0, 1, col="red")

save(mu, phy_tree, seq_depth, B, Theta, X, 
     file = paste0("~/Research_projects/tree_cov_reg/data/synthetic/params_", prefix))
