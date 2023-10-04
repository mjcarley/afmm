function Q=gafunc(n, chi)

  ## DLMF https://dlmf.nist.gov/14.15#iii
  
  eta = acosh(chi) ;
  nu = n-1/2 ;
  mu = 0 ;

  al = mu/(nu+1/2) ;
  bt = exp(mu)*((nu-mu+1/2)/(n+mu+1/2))^((nu/2)+1/4) ;
  bt *= ((nu+1/2)^2-mu^2)^(-mu/2) ;
  
  K = besselk(mu, (nu+1/2)*eta) ;
  Q = nu^mu/gamma(nu+mu+1)*sqrt(eta./sinh(eta)).*K ;
