LTNMGen = function(X, B, Theta, phy_tree, seq_depth, mu = numeric(nrow(Theta))){
  #' Generate based on the LTNM model
  
  n = nrow(X)
  p = nrow(Theta)
  d = dim(B)[2]
  R = dim(B)[3]
  
  Phi = matrix(NA,n,p)
  for(i in 1:n){
    sigma = Theta
    for(r in 1:R){
      sigma = sigma + B[,,r]%*%X[i,]%*%t(X[i,])%*%t(B[,,r])
    }
    Phi[i,] = rmvnorm(1, mean = mu, sigma = sigma)
  }
  # Generate LTNM data based on the tree
  cp_data = t(apply(exp(Phi)/(exp(Phi) + 1), 1, ProbInt2Leaf, adj = as_adjacency_matrix(phy_tree)))
  
  count_data = matrix(NA, nrow = n, ncol = p+1)
  for( i in 1:n ){
    count_data[i, ] = rmultinom(n = 1, size = seq_depth[i], prob = cp_data[i,])
  }
  
  # Transform data from leaf nodes to internal nodes
  Y <- N <- matrix(NA, nrow = n, ncol = p)
  for( i in 1:n ){
    result = CountLeaf2Int(adj = as_adjacency_matrix(phy_tree), Y = count_data[i,])
    Y[i,] = result$Y
    N[i,] = result$N
  }
  
  return(list(Y = Y, 
              X = X, 
              N = N, 
              Phi = Phi, 
              count = count_data, 
              tree = phy_tree))
}
