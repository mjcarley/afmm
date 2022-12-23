dcty = "Data" ;

nsrc = 65536 ;
dfile = [dcty "/direct-" int2str(nsrc) ".dat"] ;
direct = load(dfile) ;

for d=4:6

  err = [] ;
  T = [] ;

  for o=6:2:16
    tfile = [dcty "/tree-" int2str(d) "-" int2str(o) "-" int2str(nsrc) ".dat"] ;
    
    tree = load(tfile) ;
    
    ee = abs(tree(:,3:end)-direct(:,3:end))./max(abs(direct(:,3:end))) ;
  
    err = [err ; d o nsrc max(ee)] ;
    
    tfile = [dcty "/time-" int2str(d) "-" int2str(o) ".dat"] ;
    tdat = load(tfile) ;
    ii = find(tdat(:,1) == nsrc) ;
    T = [T; tdat(ii,end)] ;
  endfor

  t = (min(T)-5):1:(max(T)+5) ;
  p = polyfit(T, log10(err(:,4)), 2) ;

  fid = fopen(["time-error-" int2str(d) "-" int2str(nsrc) ".dat"], "w") ;
  
  dat = [T log10(err(:,4))]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ;

  f = polyval(p, t) ;
  dat = [t(:) f(:)]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ;
  
  fclose(fid) ;
endfor
