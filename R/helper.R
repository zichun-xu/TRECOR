# Helper functions
# Requires the following packages:
# library(mvtnorm)
ProbInt2Leaf = function(adj, probs){
  #' Transform the probabilities corresponding to the internal nodes of a binary 
  #' tree to the probabilities corresponding to the leaf nodes
  #' @param adj is the adjacency matrix of directed (out) binary tree
  #' @param probs is a vector of probabilities corresponding to the internal nodes
  
  p = length(probs) # number of internal nodes
  result = numeric(p+1) + 1
  for( k in (p+1):(2*p+1) ){
    j = k
    while( sum((adj[,j] == 1)) > 0 ){
      j_new = which(adj[,j] == 1)
      if( j == 2*j_new ){
        result[k-p] = result[k-p]*probs[j_new]
      }
      else{
        result[k-p] = result[k-p]*(1-probs[j_new])
      }
      j = j_new
    }
  }
  return(result)
}

CountLeaf2Int = function(adj, Y){
  #' Transform the counts corresponding to the leaf nodes of a binary 
  #' tree to the counts corresponding to the internal nodes
  #' @param adj is the adjacency matrix of directed (out) binary tree
  #' @param Y is a vector of integers corresponding to the observed counts of each leaf nodes
  
  p = length(Y)-1 # number of internal nodes
  result_N <- result_Y <- numeric(p+p+1)
  result_N[(p+1):(2*p+1)] <- result_Y[(p+1):(2*p+1)] <- Y
  for( k in p:1 ){
    children = which(adj[k,] == 1)
    result_N[k] = sum(result_N[children])
    result_Y[k] = result_N[children[1]]
  }
  return(list(N = result_N[1:p], 
              Y = result_Y[1:p]))
}


Deviance = function(Theta, B, Phi, X){
  n = nrow(Phi)
  p = ncol(Phi)
  d = ncol(X)
  R = dim(B)[3]
  dev = 0
  for(i in 1:n){
    sigma = Theta
    for(r in 1:R){
      sigma = sigma + B[,,r]%*%X[i,]%*%t(B[,,r]%*%X[i,])
    }
    dev = dev + -2*dmvnorm(log = T, x = Phi[i,], sigma = sigma)
  }
  return(dev)
}

Density = function(Theta, B, Xi, Yi, B0 = matrix(0,dim(B)[1],dim(B)[2])){
  R = dim(B)[3]
  sigma = Theta
  for(r in 1:R){
    sigma = sigma + B[,,r]%*%Xi%*%t(B[,,r]%*%Xi)
  }
  return(dmvnorm(log = F, x = Yi, mean = B0%*%Xi, sigma = sigma))
}

DensityLTNM = function(Theta, B, Phii, Xi, Yi, Ni, B0 = matrix(0,dim(B)[1],dim(B)[2])){
  R = dim(B)[3]
  p = dim(B)[1]
  sigma = Theta
  for(r in 1:R){
    sigma = sigma + B[,,r]%*%Xi%*%t(B[,,r]%*%Xi)
  }
  den = dmvnorm(log = F, x = Phii, mean = B0%*%Xi, sigma = sigma)
  for(j in 1:p){
    den = den*dbinom(x = Yi[j], size = Ni[j], prob = exp(Phii[j])/(exp(Phii[j])+1))
  }
  return(den)
}

B2B = function(B){
  #' Constrain the parameter to be identifiable
  p = dim(B)[1]
  d = dim(B)[2]
  R = dim(B)[3]
  B_copy = array(NA, dim=c(p,R,d))
  for(j in 1:d){
    for(r in 1:R){
      B_copy[,r,j] = B[,j,r]
    }
  }
  B_c1 = B_copy[,,1]
  B_c1_svd = svd(B_c1)
  Q = B_c1_svd$v
  for(j in 1:d){
    B_copy[,,j] = B_copy[,,j]%*%Q
  }
  for(j in 1:d){
    for(r in 1:R){
      B[,j,r] = B_copy[,r,j]
    }
  }
  return(B)
}

GraphSelect = function(cor_list, kappa, delta){
  n = dim(cor_list)[1]
  p = dim(cor_list)[2]
  result = matrix(NA, n, p*(p-1)/2)
  for(i in 1:n){
    result[i,] = as.numeric(abs(cor_list[i,,][upper.tri(cor_list[i,,])]) > kappa)
  }
  local_FDR = 1 - apply(result, 2, mean)
  thresholds = cummean(sort(local_FDR))
  threshold = max(thresholds[thresholds <= delta])
  graph_vec = as.numeric(local_FDR < threshold)
  graph = matrix(0,p,p)
  graph[upper.tri(graph)] = graph_vec
  graph = graph + t(graph)
  return(graph)
}

FindLeaf = function(adj, k){
  # Given a binary tree (represented by its adjacency matrix), find all leaves of a node k (as in names of adj)
  leaves = c()
  childern = c(k)
  while(length(childern) > 0){
    child = childern[1]
    grandchildern = names(which(adj[child,] == 1))
    if(length(grandchildern) == 0){
      leaves = c(leaves, child)
    } else{
      childern = c(childern, grandchildern)
    }
    childern = childern[-1]
  }
  return(leaves)
}
