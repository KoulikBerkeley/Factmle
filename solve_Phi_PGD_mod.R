solve_Phi_PGD_mod <-  function(Phi_init,S,Subgrad,maxiter,TOL)
{
  #min -logdet(Phi) + < SS, Phi> s.t. Phi = tri-diag
  
  objvals = rep(1000,maxiter);
  
  # finding indices of zero entries of Subgrad
  temp = which(Subgrad == 0, arr.ind = TRUE)
  idszero = matrix(rep(FALSE, nrow(Subgrad)*ncol(Subgrad)), nrow = nrow(Subgrad))
  idszero[temp] = TRUE
  
  
  
  phik = Phi_init; 
  
  #phik(idszero) = 0;  stepsize = 0.01;
  # a good initializatin is needed
  
  SS = S - Subgrad;
  
  for(iter in 2:maxiter){
    
    # gradient calculation and stepsize choose  
    grad = -solve(as.matrix(phik)) + SS;
    
    phik_svd = svd(phik, nu =0, nv = 0) 
    vals = max(1/phik_svd$d);
    
    stepsize = 0.5/(vals^2);
    
    # gradient descent update
    
    phiknew = phik - stepsize*grad; 
    phiknew[idszero] = 0;
    phiknew = (phiknew + t(phiknew) )/2;
    
    
    # tuning stepsize to ensure that phi_new id psd
    
    while (min(eigen(phiknew)$values) < -0.0001){
      
      print('feasibility-viol')
      stepsize = stepsize/2; 
      phiknew = phik - stepsize*grad;  phiknew[idszero] = 0; phiknew = (phiknew + t(phiknew))/2;
      
      if (stepsize <= 10^-8){
        break;
      }
      
    }
    
    phik = phiknew;
    
    # storing objective value
    svd_phiknew = svd(phiknew, nu =0, nv = 0)
    objvals[iter] = -sum(log( svd_phiknew$d )) + tr(as.matrix(t(SS)%*%phiknew)); 
    
    # tuning stepsize to ensure function value decrease                                                
    while ( (iter>3)&&(objvals[iter] > objvals[iter-1]) ){
      
      print('feasibility-viol')
      stepsize = stepsize/2; 
      
      phiknew = phik - stepsize*grad;  phiknew[idszero] = 0; phiknew = (phiknew + t(phiknew))/2;
      
      svd_phiknew2 = svd(phiknew, nu =0, nv = 0)
      objvals[iter] = -sum(log( svd_phiknew2$d )) + tr(as.matrix(t(SS)%*%phiknew));
      
      
      if (stepsize <= 10^-8){
        break
      }
      
    }
    
    phik = phiknew;
    if ((iter > 10)&( abs(objvals[iter] - objvals[iter-1]) <= TOL*abs(objvals[iter-1])) ){
      
      break;
      
    }
    
    
    
  }
  
  
  #objvals = objvals[c(1:iter)];
  Phi_out = phiknew;
  return(Phi_out)                                          
  
  
}