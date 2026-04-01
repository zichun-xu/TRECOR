### Compositional high-dimensional sparse covariance regression
### Higher rank model
### Required specification of initilizations

# Requires the following packages:
# library(mvtnorm)
# library(matrixNormal)
# library(GeneralizedHyperbolic)
# library(statsmod)
# library(BayesLogit)

GibbsLTNMRCovMeanRegGL = function(Y, N, X, R, niter, 
                                  Theta, B, B0,
                                  r_GL=1, s_GL=0.01, 
                                  a_c = 5, a_d = 5, b_c = 0.5, b_d = 0.5){
  n = nrow(X)
  d = ncol(X)
  p = ncol(Y)
  ka = Y-N/2
  
  # Initialize Phi as the log ratio
  Phi = matrix(NA, nrow=n, ncol=p)
  for (i in 1:n){
    for (j in 1:p){
      if (N[i,j]==0){
        Phi[i,j]=0
      }
      else{
        tmp=Y[i,j]/N[i,j]
        if (tmp==0){
          Phi[i,j]=-5
        }
        else{
          if (tmp==1){
            Phi[i,j]=5
          }
          else{
            Phi[i,j]=stats::qlogis(tmp)
          }
        }
      }
    }
  }
  
  # Initialize Pólya-Gamma
  Al=matrix(0, nrow=n, ncol=p)
  for (i in 1:n){
    for (j in 1:p){
      if (N[i,j]!=0){
        Al[i,j]=BayesLogit::rpg(1,N[i,j]+0.0001,0)
        while (is.na(Al[i,j])){
          Al[i,j]=BayesLogit::rpg(1,N[i,j]+0.0001,0)
          warning('NA in PG variable')
        }
      }
    }
  }
  
  C <- D <- C_inv <- D_inv <- diag(d)
  lambda= rgamma(1,r_GL,s_GL)
  Tau <- matrix(0, p, p)
  for (l in 2:p){
    for (k in 1:(l-1)){
      mu_prime=sqrt(lambda^2/Theta[k,l]^2)
      ukl=statmod::rinvgauss(1, mu_prime, lambda^2)
      Tau[k,l]=1/ukl
      Tau[l,k]=1/ukl
    }
  }
  ga = rnorm(n*R)%>%matrix(nrow = n, ncol = R)
  
  Theta_list <- B_list <- B0_list <- Phi_list <- ga_list <- list()
  
  for(iter in 1:niter){
    if(iter%%1000==1){
      print(paste0("gibbs iter: ", iter, sep=""))
    }
    res = Phi-X%*%t(B0)
    S = B0%*%D_inv%*%t(B0)
    for(r in 1:R){
      S = S + B[,,r]%*%C_inv%*%t(B[,,r])
      res = res - diag(ga[,r])%*%X%*%t(B[,,r])
    }
    for(i in 1:n){
      S = S + res[i,]%*%t(res[i,])
    }
    rm(res)
    
    ### Update Theta
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
      sigma = sigma/2+t(sigma)/2
      theta12 = as.numeric(rmvnorm(n = 1, mean = sigma%*%w, sigma = sigma))
      Theta[j,-j] <- Theta[-j,j] <- theta12
      chi = t(theta12)%*%Theta11_inv%*%S11%*%Theta11_inv%*%theta12-2*t(s12)%*%Theta11_inv%*%theta12 + s22
      u = rgig(n = 1, lambda = 1-(n+d*R+d)/2, psi = lambda, 
               chi = chi)
      theta22 = u+t(theta12)%*%Theta11_inv%*%theta12
      Theta[j,j] = theta22
    }
    Theta_inv = solve(Theta)
    
    ### Update Tau and lambda
    for (l in 2:p){
      for (k in 1:(l-1)){
        mu_prime=sqrt(lambda^2/Theta[k,l]^2)
        ukl=rinvgauss(1, mu_prime, lambda^2)
        Tau[k,l]=1/ukl
        Tau[l,k]=1/ukl
      }
    }
    lambda=rgamma(1,r_GL+p*(p+1)/2,s_GL+sum(abs(Theta))/2)
    
    ### Update ga
    for(r in 1:R){
      res = Phi-X%*%t(B0)
      for(r_p in c(1:R)[-r]){
        res  = res- diag(ga[,r_p])%*%X%*%t(B[,,r_p])
      }
      for(i in 1:n){
        a = 1/(1+t(X[i,])%*%t(B[,,r])%*%Theta_inv%*%B[,,r]%*%X[i,])
        b = t(res[i,])%*%Theta_inv%*%B[,,r]%*%X[i,]*a
        ga[i,r] = rnorm(n = 1, mean = b, sd = sqrt(a))
      } 
    }
    rm(res)
    
    ### Update B
    for(r in 1:R){
      res = Phi-X%*%t(B0)
      for(r_p in c(1:R)[-r]){
        res  = res- diag(ga[,r_p])%*%X%*%t(B[,,r_p])
      }
      Ga = diag(ga[,r])
      V = solve(t(X)%*%t(Ga)%*%Ga%*%X+C_inv)
      Bn = t(res)%*%Ga%*%X%*%V
      B[,,r] = rmatnorm(M = Bn, U = Theta, V = V) 
    }
    rm(res)
    
    ### Update B0
    res = Phi
    for(r in 1:R){
      res  = res- diag(ga[,r])%*%X%*%t(B[,,r])
    }
    V0 = solve(t(X)%*%X+D_inv)
    B0 = rmatnorm(M = t(res)%*%X%*%V0, U = Theta, V = V0)
    rm(res)
    
    ### Update C
    for(j in 1:d){
      a = 0
      for(r in 1:R){
        a = a + t(B[,j,r])%*%Theta_inv%*%B[,j,r]
      }
      C[j,j]= 1/rgamma(n=1, shape=(R*p)/2+a_c, rate=a/2+b_c)
      C_inv[j,j] = 1/C[j,j]
    }
    
    ### Update D
    for(j in 1:d){
      a = t(B0[,j])%*%Theta_inv%*%B0[,j]
      D[j,j]= 1/rgamma(n=1, shape=p/2+a_d, rate=a/2+b_d)
      D_inv[j,j] = 1/D[j,j]
    }
    
    ### Update Phi
    ### Check this!!!!!!!
    for(i in 1:n){
      D = diag(Al[i,])
      Sigma_inv = chol2inv(chol(Theta_inv + D))
      mu = B0%*%X[i,]
      for(r in 1:R){
        mu = mu + ga[i,r]*B[,,r]%*%X[i,]
      }
      Phi[i,] = rmvnorm(n = 1, mean = Sigma_inv%*%(ka[i,]+Theta_inv%*%mu), sigma = Sigma_inv)
    }
    
    ### Update Polya-Gamma
    for(i in 1:n){
      for(j in 1:p){
        if (N[i,j]!=0){
          if(N[i,j] < 200){
            Al[i,j] = rpg.devroye(1, N[i,j], Phi[i,j]) 
          }
          else{
            Al[i,j] = rpg(1, N[i,j], Phi[i,j]) 
          }
        }
      }
    }
    
    Theta_list[[iter]] = Theta
    B_list[[iter]] = B2B(B)
    B0_list[[iter]] = B0
    Phi_list[[iter]] = Phi
    ga_list[[iter]] = ga
  }
  return(list(Theta_list, B_list, B0_list, Phi_list, ga_list))
}
