## scratch space for checking formulae
r = 0.4 ; r1 = 0.9 ; z = linspace(-1, 1, 65) ; n = 3 ;

rp2 = (r + r1).^2 + z.^2 ;
rm2 = (r - r1).^2 + z.^2 ;

A = (r.^4 - (r1.^2 + z.^2).^2)./rm2./rp2./r ;

f = -1/2./r - (n+1/2)*A ;

g = r.^2.*(r.^2 - r1.^2 + z.^2) + n*(r.^4 - (r1.^2 + z.^2).^2) ;
g = -g./r./rp2./rm2 ;


