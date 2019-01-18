FA_tridiag<- function (dataset, curr_rank, Phi_init = c(1), Max_iter = 1000, tol = 10^-4, lb = 10^-3, Threshold_p = 10^-4) 
{
  l = length(as.vector(Phi_init))
  nsample = nrow(dataset)
  dim = ncol(dataset)
  S = cov(dataset)
  if(l == 1){
    Phi_init = diag(0.2 + runif(dim, min = lb, max = 1))
  }
  f = rep(1, Max_iter)
  Phi = Phi_init
  k = 2
  check = 1
  # m = colMeans(dataset)
  # stddataset = (1/sqrt(nsample - 1)) * (dataset - matrix(rep(m, nsample), nrow = nsample, ncol = dim, byrow = TRUE))
  # S = cov(stddataset)
  #S_half = sqrtm(S)
  
  s_eigen = eigen(S, symmetric = TRUE)
  D = diag(sqrt(s_eigen$values))
  S_half = s_eigen$vectors%*%D%*%t(s_eigen$vectors)
  
  diff_norm = rep(1, Max_iter)
  
  while(check){
    
    # This is the main while loop
    
    Old_Phi = Phi # old assignmet 
    
    # calculating S_Phi_half for mazorization
    
    s1 = S_half%*%Phi%*%S_half
    
    # The S_star matrix 
    
    s1 = (s1 + t(s1))/2
    vd = eigen(s1,curr_rank); d = vd$values
    v = vd$vectors
    
    # calculating subgradient
    
    diff_d = pmax(1 - 1/d, 0)
    A1 = v%*%diag(diff_d)%*%t(v)
    A2 = (S_half%*%A1%*%S_half)
    Subgrad = triu(A2) - triu(A2,2) + tril(A2,-1) - tril(A2,-2);
    
    # updating optimal values of Phi
    
    
    
    Phi = solve_Phi_PGD_mod(Phi, S, Subgrad, 1000, 10^-6)
    
    
    # convergence criterion
    
    check =((norm((Phi-Old_Phi))/norm((dim*(Old_Phi))) > Threshold_p) & (k < Max_iter));
    k = k + 1
    
    diff_norm[k] = norm((Phi-Old_Phi))/norm((dim*(Old_Phi)));
    
  }
  
  # Obtaining Psi by inverting Phi 
  
  Psi = inv(Phi)
  
  # Lambda = calculate_Lambda(Psi, S, rnk)
  out = list(Psi = Psi, norm_difference = diff_norm[c(2:(k - 1))])
  return(out)
}
