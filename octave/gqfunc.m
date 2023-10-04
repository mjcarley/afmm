function Q=gqfunc(n, chi)

  ## Gradshteyn and Ryzhik 8.852
  
  eta = acosh(chi) ;

  Q = gamma(n+1/2)*sqrt(pi)/gamma(n+1)*exp(-(n+1/2)*eta).*...
      hypergf(1/2, n+1/2, n+1, exp(-2*eta)) ;
