function P=legp(nu, m, x)

  P = zeros(size(x)) ;
  ii = find(x < 2) ;
  ii = 1:length(x) ;
  P(ii) = hypergf(nu+m+1, m-nu, m+1, 0.5 - x(ii)/2) ;
  P(ii) .*= gamma(nu+m+1)/2^m/gamma(nu-m+1)*(x(ii).^2-1).^(m/2)/gamma(m+1) ;

  return ;
  ii = find(x >= 2) ;
  F1 = hypergf(m/2-nu/2, m/2-nu/2+1/2, 1/2-nu, 1./x(ii).^2) ;
  F2 = hypergf(nu/2+m/2+1, nu/2+mu/2+1/2, nu+3/2, 1./x(ii).^2) ;

  

