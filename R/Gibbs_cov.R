### High-dimensional sparse covariance estimation
### Two kinds of prior: spike-and-slab(SS), graphical LASSO(GL)

# Requires the following packages:
# library(mvtnorm)
# library(matrixNormal)
# library(GeneralizedHyperbolic)
# library(statsmod)

GibbsCovGL = function(S, n, niter, 
                      Theta, 
                      r = 1, s = 0.01){
  p = ncol(S)
  Tau <- matrix(0, p, p)
  Theta_list <- list()
  lambda= rgamma(1,r,s)
  for (l in 2:p){
    for (k in 1:(l-1)){
      mu_prime=sqrt(lambda^2/Theta[k,l]^2)
      ukl=statmod::rinvgauss(1, mu_prime, lambda^2)
      Tau[k,l]=1/ukl
      Tau[l,k]=1/ukl
    }
  }
  for(iter in 1:niter){
    if(iter%%1000==1){
      print(paste0("gibbs iter: ", iter, sep=""))
    }
    
    # Update the covariance matrix by row/column
    for(j in 1:p){
      Theta11 = Theta[-j,-j]
      Theta11_inv = solve(Theta11)
      theta12 = as.numeric(Theta[-j,j])
      theta22 = Theta[j,j]
      u = as.numeric(theta22-t(theta12)%*%Theta11_inv%*%theta12)
      S11 = S[-j,-j]
      s12 = as.numeric(S[-j,j])
      s22 = S[j,j]
      Dtau=diag(as.numeric(Tau[-j,j]))
      A = (Theta11_inv%*%S11%*%Theta11_inv)/u + lambda*Theta11_inv
      w = Theta11_inv%*%s12/u
      sigma = solve(A+chol2inv(chol(Dtau)))
      theta12 = as.numeric(rmvnorm(n = 1, mean = sigma%*%w, sigma = sigma))
      Theta[j,-j] <- Theta[-j,j] <- theta12
      chi = t(theta12)%*%Theta11_inv%*%S11%*%Theta11_inv%*%theta12-2*t(s12)%*%Theta11_inv%*%theta12 + s22
      u = rgig(n = 1, lambda = 1-n/2, psi = lambda, 
               chi = chi)
      theta22 = u+t(theta12)%*%Theta11_inv%*%theta12
      Theta[j,j] = theta22
    }
    for (l in 2:p){
      for (k in 1:(l-1)){
        mu_prime=sqrt(lambda^2/Theta[k,l]^2)
        ukl=rinvgauss(1, mu_prime, lambda^2)
        Tau[k,l]=1/ukl
        Tau[l,k]=1/ukl
      }
    }
    lambda=rgamma(1,r+p*(p+1)/2,s+sum(abs(Theta))/2)
    Theta_list[[iter]] = Theta
  }
  return(Theta_list)
}