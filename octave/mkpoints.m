function [r,z]=mkpoints(file, n, rmin, rmax, zmin, zmax)

  fid = fopen(file, "w") ;

  r = rmin + (rmax - rmin)*rand(n,1) ;
  z = zmin + (zmax - zmin)*rand(n,1) ;

  fprintf(fid, "%d\n", n) ;

  dat = [r(:) z(:)]' ;
  fprintf(fid, "%e %e\n", dat) ;
  
  fclose(fid) ;
