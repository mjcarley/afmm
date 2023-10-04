Q = 5 ;

idx = [] ;

id = 0 ;
for q = 0:Q
  for i = 0:q
    for j = 0:q-i
      k = q - i - j ;
      idx = [idx; i j k id] ;
      id ++ ;
    endfor
  endfor
endfor

ii = [] ;
for i=1:size(idx, 1)
  ijk = idx(i,1:3) ;
  q = sum(ijk) ;
  off = q*(q+1).*(q+2)/6 ;
  ##ii = [ii; off + ijk(1)/2*(2*q+3) - ijk(1)^2/2 + ijk(2)] ;
  ii = [ii; off + ijk(1)/2*(2*q+3-ijk(1)) + ijk(2)] ;  
    
endfor
