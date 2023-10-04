function [r,z]=sbox(rb, zb, i, j)

  r = [rb(i) rb(i+1) rb(i+1) rb(i  ) rb(i)] ;
  z = [zb(j) zb(j  ) zb(j+1) zb(j+1) zb(j)] ;

