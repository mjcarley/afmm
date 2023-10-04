function [dG,nz,drmr,drpr]=dgfunc(N, I, J, K, r0, r10, z0)

  ## N: maximum azimuthal order of modal Green's functions
  ## I: maximum derivative in r
  ## J: maximum derivative in r1
  ## K: maximum derivative in z

  ## size of data block for each azimuthal mode
  nz = (K+1) ;

  ## not dealing with r1 derivatives yet
  dG = zeros((N+1)*nz, 2*(I+1)) ;

  ## fill in basic modal quantities
  [G,dGz,dGr]=gfunci(N+1, 0, 0, r0, r10, z0) ;

  ##size(dG)
  ##size(dGz)
  ##size(dGz(1:N+1))
  dG(1:nz:end,1) = G(1:N+1)' ;
  dG(2:nz:end,1) = [dGz(1:N+1).'] ;
  dG(1:nz:end,2) = [dGr(1:N+1).'] ;

  rm2 = (r0 - r10)^2 + z0^2 ;
  rp2 = (r0 + r10)^2 + z0^2 ;  

  ## derivatives w.r.t. z
  drpz = drmz = zeros(2*K+1, 1) ;
  drpz(1) = z0./rp2 ; drmz(1) = z0./rm2 ;
  drpz(2) = (1 - 2*z0*drpz(1))/rp2 ;
  drmz(2) = (1 - 2*z0*drmz(1))/rm2 ;
  for k=2:2*K
    drpz(k+1) = -k*(2*z0*drpz(k-1+1) + (k-1)*drpz(k-2+1))/rp2 ;
    drmz(k+1) = -k*(2*z0*drmz(k-1+1) + (k-1)*drmz(k-2+1))/rm2 ;
  endfor

  ##return

  ## derivatives w.r.t. r and z
  drpr = drmr = zeros(2*K+1, 1) ;
  drpr(1) = (r0^2 - r10^2 - z0^2)/rp2 ;
  drmr(1) = (r0^2 - r10^2 - z0^2)/rm2 ;
  drpr(2) = drmr(2) = -2*z0 ;
  drpr(3) = drmr(3) = -2 ;
  dpp = [rp2 2*z0 2 0 0] ;
  dpm = [rm2 2*z0 2 0 0] ;

  drpr(2) = -2*z0*(1 + drpr(1))/rp2 ;
  drmr(2) = -2*z0*(1 + drmr(1))/rm2 ;
  drpr(3) = -(2 + 4*z0*drpr(2) + 2*drpr(1))/rp2 ;
  drmr(3) = -(2 + 4*z0*drmr(2) + 2*drmr(1))/rm2 ;
  for k=3:2*K
    drpr(k+1) = -k*(2*z0*drpr(k-1+1) + (k-1)*drpr(k-2+1))/rp2 ;
    drmr(k+1) = -(k*2*z0*drmr(k-1+1) + k*(k-1)*drmr(k-2+1))/rm2 ;
  endfor

  
  if 0
    for k=1:2*K
      for q=1:min(k,2)
	drpr(k+1) -= nchoosek(k, q)*dpp(q+1)*drpr(k-q+1) ;
	drmr(k+1) -= nchoosek(k, q)*dpm(q+1)*drmr(k-q+1) ;
      endfor
      drpr(k+1) /= rp2 ;
      drmr(k+1) /= rm2 ;
    endfor
  endif
  
  for k=1:K-1
    di = 1 ;
    for n=0:N
      cft = 1 ;
      offn  = n*nz + 1 ;
      offm1 = (n+di)*nz + 1 ;
      for q=0:k
	dG(offn+k+1,1) += ...
	(n-1/2)*cft*(drpz(k-q+1)*(dG(offn+q,1) + dG(offm1+q,1)) + ...
		     drmz(k-q+1)*(dG(offn+q,1) - dG(offm1+q,1))) ;
	cft *= (k-q)/(q+1) ;
      endfor
      di = -1 ;
    endfor
  endfor

  for k=1:K
    di = 1 ;
    cft = 1 ;
    for n=0:N
      cft = 1 ;
      offn  = n*nz + 1 ;
      offm1 = (n+di)*nz + 1 ;
      dG(offn+k,2) = -dG(offn+k,1)/2/r0 ;
      for q=0:k
	dG(offn+k,2) += ...
	(n-1/2)*cft*(drpr(q+1)*(dG(offn+k-q,1) + dG(offm1+k-q,1)) + ...
		     drmr(q+1)*(dG(offn+k-q,1) - dG(offm1+k-q,1)))/2/r0 ;
	cft *= (k-q)/(q+1) ;
      endfor
      di = -1 ;
    endfor
  endfor
  
  for k=0:K
    dG(k+1:nz:end,1:2) /= gamma(k+1) ;
  endfor
  ##return
  
  for i=0:I-2
    for k=0:K-2
      for n=0:N
	off = n*nz + 1;
	dG(off+k, i+2+1) = -(k+2)*(k+1)/(i+2)/(i+1)*dG(off+k+2, i+1) ;
	tt = 0 ;
	c = 1/r0 ;
	for u=0:i
	  dt = n^2*(u+1)/(i+1)/(i+2)/r0*dG(off+k, i-u+1) ;
	  dt -= (i-u+1)/(i+2)/(i+1)*dG(off+k, i+1-u+1) ;
	  tt += c*dt ;
	  c /= -r0 ;
	endfor
	dG(off+k, i+2+1) += tt ;
      endfor
    endfor
  endfor
