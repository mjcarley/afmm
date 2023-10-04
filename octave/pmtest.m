np = 16 ;

r  = rand(np, 1) ;
r1 = rand(np, 1) ;
z  = rand(np, 1) ;
dr  = rand(np, 1)*4 ;
dr1 = rand(np, 1)*3 ;
dz  = rand(np, 1)*2 ;

a = [(0:4)' 1+(0:4)' (4:-1:0)' (0:4)'+2 1+(0:4)'+1 (4:-1:0)'+2] ;
a = [5*rand(size(a,1),1) a] ;

b = a + 1 ;

fa = przval(a, r, r1, z, dr, dr1, dz) ;
fb = przval(b, r, r1, z, dr, dr1, dz) ;

f = fa.*fb ;

c = polymul(a, b) ;
g = przval(c, r, r1, z, dr, dr1, dz) ;
