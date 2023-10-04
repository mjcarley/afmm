function p=polymul(a, b)

  ## multiply two polynomials in r, r1, z, provided in the format
  ##
  ## [c I J K i j k]
  ##
  ## where p(r,r_{1},z) = \sum c r^I r1^J z^K (dr^i dr1^j dz^k)

  c = [] ;
  for n=1:size(a, 1)
    for m=1:size(b, 1)
      c = [c; a(n,1)*b(m,1) a(n, 2:7) + b(m, 2:7)] ;
    endfor
  endfor

  [cc,ii] = sort(sum(c(:,5:7),2)) ;
  c = c(ii,:) ;

  p = c(1,:) ;
  for i=2:size(c, 1)
    if ( max(abs(p(end,2:end)-c(i,2:end))) == 0 )
      p(end,1) += c(i,1) ;
    else
      p = [p; c(i,:)] ;
    endif
  endfor
