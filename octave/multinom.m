function f=multinom(q, r, r1, z, dr, dr1, dz)

  A = r.*r1.*dr.^2  - r1.*(r1.^2 - r.^2 + z.^2).*dr ;
  B = r.*r1.*dr1.^2 - r.* (r.^2 - r1.^2 + z.^2).*dr1 ;
  C = r.*r1.*dz.^2 + 2*r.*r1.*z.*dz ;
  D = -(r.^2 + r1.^2 + z.^2).*dr.*dr1 ;

  f = 0 ;
  for n1=0:q
    for n2=0:q-n1
      for n3=0:(q-n1-n2)
	n4 = q - n1 - n2 - n3 ;
	##[n1 n2 n3 n4 n1+n2+n3+n4]
	f += nchoosek4(q, n1, n2, n3, n4)*A.^n1.*B.^n2.*C.^n3.*D.^n4 ;	
      endfor
    endfor
  endfor

  f *= 2^q ;
