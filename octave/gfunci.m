function [G,dGz,dGr,chi]=gfunci(N, m, n, r, r1, z)

  ##N = 8 ;
  G = zeros(length(z), N+1) ;

  chi = (r.^2 + r1.^2 + z.^2)/2./r1./r ;
  [max(chi) min(chi)]
  mu = sqrt(2./(1+chi)) ;

  ##[K,E] = ellipke(mu.^2) ;
  ## initialize recursion
  ##G(:,1) = mu.*K./sqrt(r.*r1)/2/pi ;
  ##G(:,2) = (chi.*mu.*K - (chi + 1).*mu.*E)./sqrt(r.*r1)/2/pi ;

  if 0
    [K, E] = ellipke(mu.^2) ;
    G(:, 1) = mu.*K./sqrt(r.*r1)/2/pi ;
    G(:, 2) = (chi.*mu.*K - (chi + 1).*mu.*E)./sqrt(r.*r1)/2/pi ;
    for n=2:N
      ##G(:,n+1) = 4*(n-1)/(2*n-1)*chi.*G(:,n) - (2*n-3)./(2*n-1).*G(:,n-1) ;
      G(:, n+1) = gqfunc(n, chi)/2/pi./sqrt(r.*r1) ;
    endfor
  endif

  Gt = zeros(length(z), 2*N+1) ;
  if 1
    Gt(:,end) = 1e-9 ; 
    Gt(:,end-1) = 1e-9 ; 
    for n=2*N:-1:2
      Gt(:,n+1-2) = (4*(n-1)*chi.*Gt(:,n+1-1) - (2*n-1)*Gt(:,n+1))/(2*n-3) ;
    endfor
    [K,E] = ellipke(mu.^2) ;
    sc = mu.*K./sqrt(r.*r1)/2/pi./Gt(:,1) ;
    for n=1:N+1
      Gt(:,n) .*= sc ;
    endfor
    G += Gt(:, 1:N+1) ;
  endif
  
  dGz  = zeros(length(z), N+1) ;
  dGr = zeros(length(z), N+1) ;

  rho1 = sqrt((r-r1).^2 + z.^2) ;
  rho2 = sqrt((r+r1).^2 + z.^2) ;

  ##A1 = 2*r1.*(r.^2 - r1.^2 - z.^2)./(rho1.^2.*rho2.^2) ;
  ##A2 = (r.^4 - (r1.^2 + z.^2).^2)./(rho1.^2.*rho2.^2)./r ;
  ##F1 = 4*r.*r1.*z./(rho1.^2.*rho2.^2) ;
  ##F2 = 2*z.*(r.^2 + r1.^2 + z.^2)./(rho1.^2.*rho2.^2) ;
  n = 0 ;
  ##dGr(:,n+1) = -G(:,n+1)/2./r + ...
  ##(n+1/2)*(A1.*G(:,n+2) - A2.*G(:,n+1)) ;
  ##dGz(:,n+1) = (n+1/2)*(F1.*G(:,n+2) - F2.*G(:,n+1)) ;
  dGr(:,n+1) = -G(:,n+1)/2./r + ...
	       (n-1/2)/2./r.*(r.^2-r1.^2-z.^2).*...
	       ((G(:,n+1) + G(:,n+1+1))./rho2.^2 + ...
		(G(:,n+1) - G(:,n+1+1))./rho1.^2) ;
  
  dGz(:, n+1) = (n-1/2)*((G(:, n+1) - G(:, n+1+1))./rho1.^2 + ...
			 (G(:, n+1) + G(:, n+1+1))./rho2.^2).*z ;

  ##F1 = (rho2.^2 + rho1.^2)./rho1.^2./rho2.^2 ;
  ##F2 = (rho2.^2 - rho1.^2)./rho1.^2./rho2.^2 ;

  for n=1:N
    ##g = (2*n*chi.*G(:,n+1) - (n-1/2)*G(:,n-1+1))/(n+1/2) ;
    ##dGz(:,n+1) = (n+1/2)*(F1.*g - F2.*G(:,n+1)) ;
    ##dGz(:, n+1) = (n-1/2)*(F1.*G(:, n+1) - F2.*G(:,n-1+1)).*z ;
    dGz(:, n+1) = (n-1/2)*((G(:, n+1) - G(:, n-1+1))./rho1.^2 + ...
			   (G(:, n+1) + G(:, n-1+1))./rho2.^2).*z ;
    ##dGr(:,n+1) = -G(:,n+1)/2./r + ...
    ##		 (n+1/2)*(A1.*g - A2.*G(:,n+1)) ;
    dGr(:,n+1) = -G(:,n+1)/2./r + ...
		 (n-1/2)/2./r.*(r.^2-r1.^2-z.^2).*...
		 ((G(:,n+1) + G(:,n-1+1))./rho2.^2 + ...
		  (G(:,n+1) - G(:,n-1+1))./rho1.^2) ;
  endfor
