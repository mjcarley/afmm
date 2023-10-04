r = 1.3 ; r1 = 2.5 ; z = linspace(-2, 2, 65)' ;
n = 32 ; m = 1 ;

nt = 4096 ;
t = (0:nt-1)/nt*2*pi ;

G = zeros(length(z), 3) ;
Gc = G ;
Gp = G ;

r2 = (r+r1)^2 + z.^2 ;
rho = sqrt(r2) ;
lm = 4*r*r1./r2 ;

rp = sqrt((r+r1)^2 + z.^2) ;
rm = sqrt((r-r1)^2 + z.^2) ;

b = 0.5*(rp.^2 + rm.^2)./rp./rm ;

for i=1:length(z)
  R = sqrt(r^2 + r1^2 - 2*r*r1*cos(t) + z(i)^2) ;
  f = cos(n*t)./R.^m/4/pi ;
  G(i,1) = sum(f)*t(2) ;
  if 0
    f = cos((n+1)*t)./R.^m/4/pi ;
    G(i,2) = sum(f)*t(2) ;
    f = cos((n+2)*t)./R.^m/4/pi ;
    G(i,3) = sum(f)*t(2) ;
  else
    f = cos(n*t)./R.^(m+2)/4/pi ;
    G(i,2) = sum(f)*t(2) ;
    f = cos(n*t)./R.^(m+4)/4/pi ;
    G(i,3) = sum(f)*t(2) ;
  endif
endfor

## check on recursion
d = (m-2*n-4)*G(:,3) + 4*(n+1)*b./sqrt(b.^2-1).*G(:,2) - (m+2*n)*G(:,1) ;

##stop

d = (m/2-n+1)*(rp.*rm).^2*gamma(m/2+2)/gamma(m/2+2-n).*G(:,3) - ...
    (m+1)*b.*(rp.*rm)*gamma(m/2+1)/gamma(m/2+1-n).*G(:,2) + ...
    (m/2+n)*gamma(m/2)/gamma(m/2-n)*G(:,1) ;

d = (m-2*n+2)*(rp.*rm).^2*(m/2+1)*m/2/gamma(m/2+2-n).*G(:,3) - ...
    2*(m+1)*b.*(rp.*rm)*m/2/gamma(m/2+1-n).*G(:,2) + ...
    (m+2*n)/gamma(m/2-n)*G(:,1) ;

d = m*(m+2)*(rp.*rm).^2.*G(:,3) - 2*(m+1)*m*b.*(rp.*rm).*G(:,2) + ...
    (m+2*n)*(m-2*n)*G(:,1) ;

[K,E] = ellipke(lm) ;
G0 = K./rho/pi ;
G1 = (2*(K-E)./lm - K)/pi./rho ;

##stop

t = (0:nt-1)/nt*pi ;
for i=1:length(z)
  R = sqrt(1 - lm(i)*cos(t).^2) ;
  Gc(i) = sum(cos(2*n*t)./R.^m)*t(2) ;
endfor

Gc ./= rho.^m*2*pi ;

##P = hypergf(m/2-1+1, -(m/2-1), 1-n, 0.5 - b/2) ;
##P .*= ((b+1)./(b-1)).^(n/2) ;

nu = m/2-1 ;
## Gradshteyn and Ryzhik 8.751.1
P = hypergf(m/2-1+n+1, n-(m/2-1), n+1, 0.5 - b/2) ;
P .*= gamma(m/2-1+n+1)/2^n/gamma(m/2-1-n+1)*(b.^2-1).^(n/2)/gamma(n+1) ;

P = legp(m/2-1, n, b) ;

Gp = P/2./sqrt(rp.*rm).^m*gamma(m/2-n)/gamma(m/2) ;

