function Q=qleghg(nu, x)

  a = (nu+2)/2 ;
  b = (nu+1)/2 ;
  c = nu+3/2 ;
  z = 1./x.^2 ;

  [a b c]
  
  if 0
  Q = 1 ;
  tm = 1 ;
    for q=0:8
    tm .*= (a+q)*(b+q)/(c+q)/(q+1)*z ;
    Q += tm ;
  endfor
  endif

  ##Q = hypergf(a, b, c, 1./x.^2) ;
  ##Q .*= gamma(nu+1)*gamma(0.5)/2^(nu+1)/gamma(nu+3/2)*x.^(-nu-1) ;
  ##z = 1./x.^2 ;

  al = acosh(x) ;

  Q = hypergf(1/2, 1/2, nu+3/2, 1./(1-exp(2*al))) ;
  Q .*= sqrt(pi)*gamma(nu+1)/gamma(nu+3/2).*...
	exp(-(nu+1)*al)./sqrt(1-exp(-2*al)) ;
  
