r0 = 0.8 ; r10 = 0.2 ; z0 = 0.0 ; N = 256 ;

##z = z0*ones(65, 1) + linspace(-0.1, 0.1, 65)' ;
##z = z0*ones(65, 1) + 0.025 ;
##r = r0 + linspace(-0.1, 0.1, 65)' ;
##r = r0 + 0.005 ;
##ones(size(z)) ;

[r,z] = meshgrid(r0+linspace(-0.1,0.1,65), z0+linspace(-0.1,0.1,65)) ;

r1 = r10*ones(size(z)) ;

[G,dGz,dGr,chi] = gfunci(N, 0, 0, r(:), r1(:), z(:)) ;

if 0
  ee = 1e-6 ;
  [gt,gz,Grp] = gfunci(N, 0, 0, r0+ee/2, r10, z0) ;
  [gt,gz,Grm] = gfunci(N, 0, 0, r0-ee/2, r10, z0) ;
  
  Gr = (Grp-Grm)/ee ;

  [dG, nz, dF1, dF2, dr] = dgfunc(N, 48, 16, 48, r0, r10, z0) ;

  stop
endif

[dG, nz, dF1, dF2] = dgfunc(N, 24, 8, 24, r0, r10, z0) ;

##stop

if 0
  rp2 = (r + r1).^2 + z.^2 ;
  rm2 = (r - r1).^2 + z.^2 ;

  A1 = 2*r1.*(r.^2 - r1.^2 - z.^2)./rp2./rm2 ;
  A2 = (r.^4 - (r1.^2 + z.^2).^2)./r./rp2./rm2 ;
  
  a1 = a2 = 0 ;
  for i=0:48
    a1 += dF1(i+1)*(r-r0).^i/gamma(i+1) ;
    a2 += dF2(i+1)*(r-r0).^i/gamma(i+1) ;
  endfor
  
  stop
endif

g = [] ;
for n=0:N
  g = [g gfeval(dG(n*nz+1:(n+1)*nz, :), 16, 0, 16, r0, r, r10, r1, z0, z)] ;
endfor
