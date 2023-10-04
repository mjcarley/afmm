di = [0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3;
      2, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]' ;

N = 6 ; ns = 2 ; level = 5 ;
stub = "ilist-test" ;

for i=1:12
  file = [stub "-p" int2str(di(i,1)) "-p" int2str(di(i,2))] ;
  mkdpinput(file, N, ns, level, 4, 3, +di(i,1), +di(i,2)) ;

  file = [stub "-p" int2str(di(i,1)) "-m" int2str(di(i,2))] ;
  mkdpinput(file, N, ns, level, 4, 3, +di(i,1), -di(i,2)) ;

  file = [stub "-m" int2str(di(i,1)) "-p" int2str(di(i,2))] ;
  mkdpinput(file, N, ns, level, 4, 3, -di(i,1), +di(i,2)) ;

  file = [stub "-m" int2str(di(i,1)) "-m" int2str(di(i,2))] ;
  mkdpinput(file, N, ns, level, 4, 3, -di(i,1), -di(i,2)) ;
endfor
