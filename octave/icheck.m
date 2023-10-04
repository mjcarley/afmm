depth = 4 ;

load ilist.dat
xi = ilist(:,1) ; yi = ilist(:,2) ;

[xg,yg] = meshgrid(0:(2^depth)) ;
[xp,yp] = meshgrid(2*(0:(2^(depth-1)))) ;

