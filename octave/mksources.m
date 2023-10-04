function [r,z]=mksources(file, N, ns, n, rmin, rmax, zmin, zmax)

  fid = fopen(file, "w") ;

  r = rmin + (rmax - rmin)*rand(n,1) ;
  z = zmin + (zmax - zmin)*rand(n,1) ;

  fprintf(fid, "%d %d %d\n", n, N, ns) ;

  C = rand(n, 2*(N+1)*ns) ;
  C(:, 2:(2*N+2):end) = 0 ;

  ##C *= 0 ;
  ##C(:, 1:(2*N+2):end) = 1 ;
  
  for i=1:n
    fprintf(fid, "%e %e", r(i), z(i)) ;
    fprintf(fid, " %e", C(i,:)) ;
    fprintf(fid, "\n") ;
  endfor
  
  fclose(fid) ;
