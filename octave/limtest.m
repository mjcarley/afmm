r = 0.3 ; r1 = 0.125 ;

chi0 = (r^2 + r1^2)/2/r/r1 ;

dr1 = linspace(-0.125, 0.125, 65)' ;

sgn = -1 ;

##rdr = sqrt((r^2 - 3*r*r1 + r1^2)*(r^2 + r*r1 + r1^2)) ;
##rdr = sqrt((r^2 - 2*r*r1 + r1^2)*(r^2 + r1^2)) ;
rdr = abs(r-r1)*sqrt(r^2 + r1^2) ;
##rdr = r^2 - r*r1 + r1^2 + sgn*rdr ;
rdr = r^2 - r*r1 + r1^2 + [rdr -rdr];
rdr = (r1 + [dr1 dr1])/r/r1.*rdr ;

dr = rdr - r ;

##dr = dr - sgn*1e-3 ;

chi = (r + dr).^2 + (r1 + dr1).^2 ;
chi ./= 2*(r+dr).*(r1+dr1) ;

[rr1,rr] = meshgrid(dr1, linspace(min(min(dr)), max(max(dr)), 65)) ;

dchi = (r + rr).^2 + (r1 + rr1).^2 ;
dchi ./= 2*(r+rr).*(r1+rr1) ;

dchi -= chi0 ;
dchi -= chi0-1 ;
