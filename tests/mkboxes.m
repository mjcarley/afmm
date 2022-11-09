function [rb,zb] = mkboxes(rmin, rmax, zmin, zmax, d)

  if ( nargin < 2 )
    d = rmin ;
    rmin = zmin = 0 ;
    rmax = zmax = 1 ;
  endif

  nb = 2^d+1 ;
  
  [rb,zb] = meshgrid(linspace(rmin, rmax, nb), linspace(zmin, zmax, nb)) ;
