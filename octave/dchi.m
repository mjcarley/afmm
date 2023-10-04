## check on polynomial evaluation of Taylor series terms in expansion
## of Q

np = 8 ;

r  = rand(np, 1) ;
r1 = rand(np, 1) ;
z  = rand(np, 1) ;
dr  = rand(np, 1) ;
dr1 = rand(np, 1) ;
dz  = rand(np, 1) ;
##dr1 = 0 ;
##dz = 0 ;

q = 12 ;
f = 2*r.*r1.*((r+dr).^2 + (r1+dr1).^2 + (z+dz).^2) - ...
    2*(r+dr).*(r1+dr1).*(r.^2 + r1.^2 + z.^2) ;

f = f.^q ;

g = multinom(q, r, r1, z, dr, dr1, dz) ;

if 0 
## coefficients of polynomial form
p = [02 2 1 0   1 0 0 ; ...
     -2 3 0 0   0 1 0 ; ...
     02 1 1 0   2 0 0 ; ...
     -2 2 0 0   1 1 0 ; ...
     02 1 2 0   0 1 0 ; ...
     -2 0 3 0   1 0 0 ; ...
     02 1 1 0   0 2 0 ; ...
     -2 0 2 0   1 1 0 ; ...
     04 1 1 1   0 0 1 ; ...
     02 1 1 0   0 0 2 ; ...
     -2 0 1 2   1 0 0 ; ...
     -2 1 0 2   0 1 0 ; ...
     -2 0 0 2   1 1 0] ;

q = polymul(p, p) ;
##q = polymul(q, q) ;
#q = polymul(q, q) ;
#q = polymul(q, p) ;

f = f.^2 ;

g = przval(q, r, r1, z, dr, dr1, dz) ;
endif
