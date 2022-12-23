##o = 16 ;
d = 4 ;

tt = [] ;

for d=4:7
for o=16
  dat = load(["Data/time-" int2str(d) "-" int2str(o) ".dat"]) ;
  tt = [tt dat(:,5)] ;
endfor
endfor

tr = min(min(tt)) ;
tr = 1 ;

tt /= tr ;

np = dat(:,1) ;
td = dat(:,4)/tr ;

fid = fopen(["time-" int2str(o) ".dat"], "w") ;

dat = [log2(np) log10(td)]' ;
fprintf(fid, "%f %f\n", dat) ;
fprintf(fid, "\n") ;

for i=1:4
  dat = [log2(np) log10(tt(:,i))]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ;
endfor

fclose(fid) ;
