function [G,bt]=gfuncp(m, n, r, r1, z)

  
  nt = 2048 ;
  G = zeros(size(z)) ;

  t = (0:nt-1)*2*pi/nt ;
  C = cos(t) ;
  Cn = exp(j*n*t) ;

  rp = sqrt((r+r1).^2 + z.^2) ;
  rm = sqrt((r-r1).^2 + z.^2) ;
  bt = (r.^2 + r1.^2 + z.^2)./rp./rm ;

  for i=1:length(z)
    R = sqrt(bt(i) + sqrt(bt(i)^2-1)*C) ;
    ##S = sqrt(r(i)^2 + r1(i)^2 + 2*r(i)*r1(i)*C + z(i)^2) ;
    f = Cn./R.^m ;
    G(i) = sum(f)*diff(t(1:2))/4/pi ;
  endfor

  G ./= (-1)^n*sqrt(rp.*rm) ;
