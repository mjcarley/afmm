function G=gfeval(dG, I, J, K, r0, r, r10, r1, z0, z)

  G = 0 ;

  dr = 1 ;
  ##for j=0:J
  for i=0:I
    dz = 1 ;
    for k=0:K
      ##G += dG(k+1, j+1)*(r1 - r10).^j/gamma(j+1).*(z - z0).^k/gamma(k+1) ;
      G += dG(k+1, i+1).*dr.*dz ;
      dz .*= (z - z0) ; # /(k+1) ;
    ##(r1 - r10).^j/gamma(j+1).*(z - z0).^k/gamma(k+1) ;
    endfor
    dr .*= (r - r0) ; #/(i+1) ;
  endfor
