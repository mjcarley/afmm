rmax = 2 ;

d = 3 ;
ee = 0.1 ;

cc = 1 + ee ;

cc = 3/2/sqrt(2) ;

nb = 2^d ;
dr = rmax/nb ;

i0 = 31 ;
i0 = 21 ;
i1 = i0+3 ; j = 0 ;

##i1 = ((1+ee) - sqrt(ee)*sqrt(2+ee))*i0 ;
##i1 = (cc + sqrt(cc^2-1)*[-1 1])*i0 ;
ii1 = cc*[i0 i0+1] + sqrt(cc^2-1)*[-i0 i0+1] ;

j = sqrt(2*cc - 2)*i0 ; i1 = i0 ;

r1 = i0*dr ; r = i1*dr ; z = j*dr ;

chi = (r.^2 + r1.^2 + z^2)/2./r./r1 ;
dc = chi - 1 ;

dd = ((i0-i1).^2+j^2)/2/i0./i1 ;



