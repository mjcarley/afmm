function mkdpinput(tcase, N, ns, depth, i, j, di, dj)

  ## generate source and field files to test source-to-local
  ## translations on downward pass

  rmin = 0  ; rmax = 1 ;
  zmin = -1 ; zmax = 1 ;
  
  nb = 2^depth ;
  
  dr = (rmax - rmin)/nb ;
  dz = (zmax - zmin)/nb ;

  ## source box centre
  rc = rmin + (i+0.5)*dr ;
  zc = zmin + (j+0.5)*dz ;

  ## field box centre
  rf = rmin + (i+di+0.5)*dr ;
  zf = zmin + (j+dj+0.5)*dz ;
  
  fid = fopen([tcase "-source.dat"], "w") ;
  fprintf(fid, "%d %d %d\n", 1, N, ns) ;

  r = rc - 0.5*dr + dr*rand(1,1) ;
  z = zc - 0.5*dz + dz*rand(1,1) ;

  fprintf(fid, "%f %f", r, z) ;
  for s=1:ns
    dat = 2*(rand(1, 2*(N+1))-0.5) ;
    ##dat(2) = 0 ;
    fprintf(fid, " %f", dat) ;
    fprintf(fid, "\n") ;
  endfor
  
  fclose(fid) ;

  fid = fopen([tcase "-field.dat"], "w") ;
  fprintf(fid, "%d %d %d\n", 1, 0, 0) ;

  r = rf - 0.5*dr + dr*rand(1,1) ;
  z = zf - 0.5*dz + dz*rand(1,1) ;

  fprintf(fid, "%f %f\n", r, z) ;

  fclose(fid) ;
