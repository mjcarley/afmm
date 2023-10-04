function F=hypergf(a, b, c, z)

  F = 1 ;
  tm = 1 ;
  tol = 1e-14 ;
  
  for q=0:120
    tm .*= (a+q)*(b+q)/(c+q)/(q+1)*z ;
    F += tm ;
    if ( abs(tm) < tol ) break ; endif
  endfor
