o = 6 ;

pf = [] ;
pp = [] ;
td = [] ;
tp = [] ;
tr = [] ;
tpr = [] ;
ppb = [] ;

for d=4:7
  dat = load(["Data/breakdown-" int2str(d) "-" int2str(o) ".dat"]) ;
  td = [td dat(:,4)] ;
  tp = [tp dat(:,5)-dat(:,4)] ;
  tr = [tr dat(:,4)./(dat(:,1).^2/4^d)] ;
  tpr = [tpr tr(:,end)./dat(:,1)] ;
  ppb = [ppb dat(:,1)/4^d] ;
  ##p = polyfit(dat(:,1).^2/4^d, log(td(:,end)), 1) ;
  ii = find(ppb(:,end)>=1) ;
  p = polyfit(log2(ppb(ii,end)), log10(td(ii,end)), 1) ;
  pf = [pf; p] ;
  p = polyfit(dat(ii,1), tp(ii,end), 1) ;
  pp = [pp; p] ;
endfor
##stop
np = dat(:,1) ;

for i=1:3
  fid = fopen(["bdown-direct-" int2str(i+3) ".dat"], "w") ;
  dat = [log2(ppb(:,i)) log10(td(:,i))]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ; 

  x = dat(1,1):dat(1,end) ;
  y = polyval(pf(i,:), x) ;
  dat = [x(:) y(:)]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ; 
	  
  fclose(fid) ;

  fid = fopen(["bdown-pass-" int2str(i+3) ".dat"], "w") ;
  dat = [np/1e2 tp(:,i)]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ; 

  x = linspace(0, 800,8) ;
  y = polyval(pp(i,:), x*1e2) ;
  dat = [x(:) y(:)]' ;
  fprintf(fid, "%f %f\n", dat) ;
  fprintf(fid, "\n") ; 
	  
  fclose(fid) ;
endfor
