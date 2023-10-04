r = 1.3 ; r1 = 2.5 ; z = linspace(-2, 2, 65)' ;
n = 0 ; m = 1 ;

N = 64 ;

r2 = (r+r1)^2 + z.^2 ;
rho = sqrt(r2) ;
lm = 4*r*r1./r2 ;

rp = sqrt((r+r1)^2 + z.^2) ;
rm = sqrt((r-r1)^2 + z.^2) ;

b = 0.5*(rp.^2 + rm.^2)./rp./rm ;

Gr = zeros(length(z), N+1) ;

if 1
Gr(:, end) = 1e-6 ;

for n=N-2:-1:0
  Gr(:,n+1) = ((m-2*n-4)*Gr(:,n+2+1) + ...
	       4*(n+1)*b./sqrt(b.^2-1).*Gr(:,n+1+1))/(m+2*n) ;
endfor
endif

stop

[K,E] = ellipke(lm) ;
Gr(:, 1) = K./rho/pi ;
Gr(:, 2) = (2*(K-E)./lm - K)/pi./rho ;

for n=0:N
  Gr(:,n+2+1) = (-4*(n+1)*b./sqrt(b.^2-1).*Gr(:,n+1+1) + ...
		 (m+2*n)*Gr(:,n+1))/(m-2*n-4) ;  
endfor
