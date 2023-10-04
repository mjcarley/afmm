rmax =  2 ;
zmin = -1 ;
zmax =  1 ;

nb = 8 ;
zbs = zmin + (zmax - zmin)*linspace(0, 1, nb+1) ;
rbs = rmax*linspace(0, 1, nb+1).^3 ;
zbt = zmin + (zmax - zmin)*linspace(0, 1, 2*nb+1) ;
rbt = rmax*linspace(0, 1, 2*nb+1) ;

ir = 6 ; iz = 5 ;
[rs, zs] = sbox(rbs, zbs, ir, iz) ;
##jr = 2*(ir-1)-2 ; jz = 2*(iz-1)-3 ;
##[rt, zt] = sbox(rbt, zbt, jr, jz) ;
jr = ir + 2 ; jz = iz + 2 ;
[rt, zt] = sbox(rbs, zbs, jr, jz) ;

t = 1/2 ; s = 1/2 ;
r  = (1-s)*rt(1) + s*rt(2) ;
z  = mean(zt(1:4)) ;
r1  = (1-t)*rs(1) + t*rs(2) ;
z1 = mean(zs(1:4)) ;

chi0 = (r^2 + r1^2 + (z-z1)^2)/2/r/r1 ;

chi = ij = [] ;
for i=1:4
  for j=1:4
    c = (rt(i)^2 + rs(j)^2 + (zt(i)-zs(j))^2)/2/rt(i)/rs(j) ;
    chi = [chi; c] ;
    ij = [ij; i j] ;
  endfor
endfor

rc = chi0 - 1 ;
