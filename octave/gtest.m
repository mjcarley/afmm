np = 65 ;
n = 0 ; m = 1 ;
N = 256 ;

r = 1.8 ;
r1 = 1.6 ;
##r = r1/sqrt(2)*1.1 ;
ee = 3/sqrt(8)-1 ;
zmin = sqrt(2*ee*r*r1 - (r-r1)^2) ;
##3/sqrt(8)*2*r*r1 - r^2 - r1^2) ;
r = r1 ;
zmin = 1e-2 ;
z = linspace(1.0*zmin, 1, np)' ;
r  *= ones(np, 1) ;
r1 *= ones(np, 1) ;
ee = 1e-6 ;

G = [] ;
dGz = [] ;
dG1 = [] ;
Gp = [] ;
for n=0:N
  G = [G gfunc(m, n, r, r1, z)] ;
  Gp = [Gp gfuncp(m, n, r, r1, z)] ;
  dG = [dGz (gfunc(m, n, r, r1, z+ee/2) - gfunc(m, n, r, r1, z-ee/2))/ee] ;
  dG1 = [dG1 (gfunc(m, n, r, r1+ee/2, z) - gfunc(m, n, r, r1-ee/2, z))/ee] ;
endfor

[Gi,dGi,dGri,chi] = gfunci(N, m, n, r, r1, z) ;
stop

[Gt, G2z, Grt] = gfunci(m, n, r, r1, z+ee/2) ;
[Gt, Gtz, Grt] = gfunci(m, n, r, r1, z-ee/2) ;

dG2z = (G2z - Gtz)/ee ;

[Gt, Gtz, G2r] = gfunci(m, n, r, r1+ee/2, z) ;
[Gt, Gtz, Grt] = gfunci(m, n, r, r1-ee/2, z) ;

dG2r = (G2r - Grt)/ee ;

