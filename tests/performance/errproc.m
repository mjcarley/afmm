if 1
  nsrc = [1024 2048 4096 8192 16384 32768 65536] ;
nsrc = [1024 2048 4096 8192 16384 24576 32768 40960 49152 57344 65536] ;
nsrc = [1024 2048 4096 8192 16384 32768 65536] ;## 131072] ;
##nsrc = 16384 ;
##nsrc = 8192 ;
dcty = "Data" ;
##dcty = "." ;

err = [] ;

d = 6 ;
o = 6 ;

nsrc = 65536 ;
for i=1:length(nsrc)
##i = 7 ; ##length(nsrc) ;
dfile = [dcty "/direct-" int2str(nsrc(i)) ".dat"] ;
direct = load(dfile) ;
for o=6:2:16 ##6:2:16
  tfile = [dcty "/tree-" int2str(d) "-" int2str(o) "-" int2str(nsrc(i)) ".dat"] ;

  tree = load(tfile) ;

  ee = abs(tree(:,3:end)-direct(:,3:end))./max(abs(direct(:,3:end))) ;

  err = [err ; d o nsrc(i) max(ee)] ;
endfor  
endfor
endif

fid = fopen(["error-" int2str(nsrc) "-" int2str(d) ".dat"], "w") ;

for i=1:length(6:2:16)
  dat = [(0:16)' log10(err(i,[4:2:end])')]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ;
endfor

fclose(fid) ;
