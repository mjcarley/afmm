nsrc = 1024 ;
dcty = "Data" ;

nsrc = 65536 ;
##nsrc = 16384 ;
for d=5:7
##d = 7 ;
##n = 0 ;

  err = [] ;
  for i=1:length(nsrc)
    ##i = 7 ; ##length(nsrc) ;
    dfile = [dcty "/direct-" int2str(nsrc(i)) ".dat"] ;
    direct = load(dfile) ;
    for o=6:2:16 ##6:2:16
      tfile = [dcty "/tree-" int2str(d) "-" int2str(o) "-" ...
		    int2str(nsrc(i)) ".dat"] ;
      
      tree = load(tfile) ;
      
      ee = abs(tree(:,3:end)-direct(:,3:end))./max(abs(direct(:,3:end))) ;
      
      err = [err ; d o nsrc(i) max(ee)] ;
    endfor  
  endfor
  
  for n=[0 8]
    fid = fopen(["error-order-" int2str(nsrc) "-" int2str(d) "-" int2str(n) ...
				".dat"], "w") ;
    
    dat = [err(:,2), log10(err(:,4+2*n))]' ;
    fprintf(fid, "%f %f\n", dat) ;
    fprintf(fid, "\n") ;
    
    p = polyfit(err(:,2), log10(err(:,4+2*n)), 1) ;
    
    dat = [err([1 end],2), polyval(p, err([1 end],2))]' ;
    fprintf(fid, "%f %f %% exponent = %f\n", dat(:,1), p(1)) ;
    fprintf(fid, "%f %f\n", dat(:,2:end)) ;
    fprintf(fid, "\n") ;
    
    fclose(fid) ;
    
    A = exp(p(2)) ; a = p(1) ;
  endfor

endfor
