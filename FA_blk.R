FA_blk<- function(dataset, curr_rank, blk_number = 0, Phi_init = c(1), Max_iter = 1000, lb = 10^-1, tol = 10^-4)
{
  # packages that we need: 
  # pracma -- for linspace. 
  # rARPACK for eigs
  # psych for tr()
  
  l = length(as.vector(Phi_init))
  nsample = nrow(dataset)
  dim = ncol(dataset)
  S = cov(dataset)
  S = (S + t(S))/2 
  
  # initialization
  if(l == 1){
    Phi_init = diag(0.2 + runif(dim, min = lb, max = 1))
  }
  
  # assigning equally sized blocks from blk_number
  
  if(blk_number){
    
    blk_start = seq(from = 1, to = dim, by = dim/blk_number)
    blk_end =  seq(from = dim/blk_number, to = dim, by = dim/blk_number)
    
  }
  
  # f = rep(1, Max_iter)
  Phi = Phi_init
  k = 2
  check = 1
  
  # m = colMeans(dataset)
  # stddataset = (1/sqrt(nsample - 1)) * (dataset - matrix(rep(m, nsample), nrow = nsample, ncol = dim, byrow = TRUE))
  # S = cov(stddataset)
  # S_half = sqrtm(S)
  
  s_eigen = eigen(S, symmetric = TRUE)
  D = diag(sqrt(s_eigen$values))         ## add non-negative eigenvalue condition
  S_half = s_eigen$vectors%*%D%*%t(s_eigen$vectors)
  
  diff_norm = rep(1, Max_iter)
  
  while(check){
    
    # This is the main while loop
    
    Old_Phi = Phi # old assignmet 
    
    # calculating S_Phi_half for mazorization
    
    s1 = S_half%*%Phi%*%S_half
    
    # The S_star matrix 
    
    s1 = (s1 + t(s1))/2
    vd = eigs(s1, curr_rank); d = vd$values
    v = vd$vectors
    
    # calculating subgradient
    
    diff_d = pmax(1 - 1/d, 0)
    A1 = v%*%diag(diff_d)%*%t(v)
    Subgrad = (S_half%*%A1%*%S_half)
    
    # updating optimal vaue of Phi
    
    for(ind in 1:blk_number){
      
      Phi[blk_start[ind]:blk_end[ind],blk_start[ind]:blk_end[ind]] = 
        solve(S[blk_start[ind]:blk_end[ind],blk_start[ind]:blk_end[ind]] - Subgrad[blk_start[ind]:blk_end[ind],blk_start[ind]:blk_end[ind]])
      
    } # end of for loop 
    
    # convergence criterion
    
    check =((norm((Phi-Old_Phi))/norm((dim*(Old_Phi))) > tol) & (k < Max_iter));
    k = k + 1
    diff_norm[k] = norm((Phi-Old_Phi))/norm((dim*(Old_Phi)));
    
    
  } # end of while loop 
  
  # Obtaining Psi by inverting Phi 
  
  Psi = solve(Phi)
  # Lambda = calculate_Lambda(Psi, S, rnk)
  out = list(Psi = Psi, Psi_inverse = Phi, norm_difference = diff_norm[c(2:(k - 1))])
  return(out)   
  
} # end of function 