## check axial derivatives of 1/(\chi\pm1)
z0 = 0.5 ; r = 0.9 ; r10 = 0.93 ;

z  = z0 + linspace(-0.2, 0.2, 65) ;
r1 = r10 + linspace(-0.2, 0.2, 65) ;
[z, r1] = meshgrid(z, r1) ;

chi = (r.^2 + r1.^2 + z.^2)/2./r1./r ;

f = 1./(chi.^2 -1) ;

K = 32 ; J = 32 ;
chi0 = (r.^2 + r10.^2 + z0.^2)/2./r10./r ;
df = zeros(K+1, J+1) ;
df(1) = 1./(chi0.^2 - 1) ;

## tables of elementary derivatives
dr = zeros(K+1,J+1) ;
## derivatives w.r.t. z
dr(1,1) = ((r-r10)^2 + z0^2)*((r+r10)^2 + z0^2) ;
dr(2,1) = 4*z0^3 + 2*z0*((r-r10)^2+(r+r10)^2) ;
dr(3,1) = 12*z0^2 + 2*((r-r10)^2+(r+r10)^2) ;
dr(4,1) = 24*z0 ; 
dr(5,1) = 24 ;

## derivatives w.r.t. r1
dr(1,2) = -2*(r-r10)*((r+r10)^2 + z0^2) + 2*(r+r10)*((r-r10)^2+z0^2) ;
dr(1,3) = 4*z0^2 + 2*(r+r10)^2 + 2*(r-r10)^2 - 8*(r^2-r10^2) ;
dr(1,4) = 24*r10 ;
dr(1,5) = 24 ;

## cross derivatives w.r.t. z and r1
dr(2,2) = 8*z0*r10 ;
dr(2,3) = 8*z0 ;
dr(3,2) = 8*r10 ;
dr(3,3) = 8 ;

dc = zeros(K+1,J+1) ;
dc(1,1) = 4*r^2*r10^2 ;
dc(1,2) = 8*r^2*r10 ;
dc(1,3) = 8*r^2 ;

for k=0:K
  df(k+1, 1) = dc(k+1, 1) ;
  for q=1:min(k, 4)
    df(k+1, 1) -= nchoosek(k, q)*dr(q+1, 1)*df(k-q+1, 1) ;
  endfor
  df(k+1, 1) /= dr(1,1) ;
endfor

for j=0:J
  df(1, j+1) = dc(1,j+1) ;
  for u=1:min(j, 4)
    df(1, j+1) -= nchoosek(j, u)*dr(1, u+1)*df(1, j-u+1) ;
  endfor
  df(1, j+1) /= dr(1,1) ;
endfor

for k=1:K
  for j=1:J
    df(k+1,j+1) = dc(k+1, j+1) ;
    q = 0 ;
    for u=0:min(4,j)
      df(k+1,j+1) -= nchoosek(k,q)*nchoosek(j, u)*...
		     dr(q+1, u+1)*df(k-q+1, j-u+1) ;
    endfor
    for q=1:min(4,k)
      for u=0:min(4,j)
	df(k+1,j+1) -= nchoosek(k,q)*nchoosek(j, u)*...
		       dr(q+1, u+1)*df(k-q+1, j-u+1) ;
      endfor
    endfor
    df(k+1,j+1) /= dr(1,1) ;
  endfor
endfor

g = 0 ;
for k=0:K
  for j=0:J
    g += df(k+1,j+1)*(z-z0).^k.*(r1-r10).^j/gamma(k+1)/gamma(j+1) ;
  endfor
endfor
