z0 = 0.5 ; r = 0.9 ; r10 = 0.93 ;
z  = z0 + linspace(-0.2, 0.2, 65) ;
r1 = r10 + linspace(-0.2, 0.2, 65) ;

##[z, r1] = meshgrid(z, r1) ;

ee = 1e-6 ;

chi = (r.^2 + r1.^2 + z0.^2)/2./r1./r ;
dc  = (r.^2 + (r1+ee/2).^2 + z0.^2)/2./(r1+ee/2)./r ;
dc -= (r.^2 + (r1-ee/2).^2 + z0.^2)/2./(r1-ee/2)./r ;
dc /= ee ;

df = (r1.^2 - r^2 - z0^2)/2/r./r1.^2 ;
