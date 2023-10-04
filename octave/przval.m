function f=przval(p, r, r1, z, dr, dr1, dz)

  f = 0 ;
  for n=1:size(p, 1)
    f += p(n,1)*r.^p(n,2).*r1.^p(n,3).*z.^p(n,4).*...
	 dr.^p(n,5).*dr1.^p(n,6).*dz.^p(n,7) ;
  endfor
