/* This file is part of AFMM, a Wide-Band Fast Multipole Method code
 *
 * Copyright (C) 2022 Michael Carley
 *
 * AFMM is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  AFMM is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AFMM.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <fftw3.h>

#include <blaswrap.h>

#include <afmm.h>

#include "afmm-private.h"

gint AFMM_FUNCTION_NAME(afmm_laplace_gfunc)(gint N,
					    AFMM_REAL r,
					    AFMM_REAL r1,
					    AFMM_REAL z,
					    AFMM_REAL *G, gint str)

{
  /*
   * generate modal Green's functions and derivatives
   *
   * Cohl and Tohline, THE ASTROPHYSICAL JOURNAL, 527:86--101, 1999
   * December 10
   *
   */
  AFMM_REAL chi, sc, Gn, Gnm1, mu, K, E ;
  gint n, dN ;

  if ( r == 0 ) {
    for ( n = 1 ; n <= N ; n ++ ) G[n*str] = 0.0 ;
    G[0] = 0.5/SQRT(r1*r1+z*z) ;
    return 0 ;
  }
  
  dN = 32 ;
  chi = 0.5*(r*r + r1*r1 + z*z)/r1/r ;
  /* if ( chi < 1.01 ) { */
  /*   g_error("%s: possible divergence r=%lg; r1=%lg; z=%lg; chi=%lg", */
  /* 	    __FUNCTION__, r, r1, z, chi) ; */
  /* } */

  /*descending recursion to seed the recursion for the array*/
  Gn = 1.0 ; Gnm1 = 1.0 ;
  for ( n = N + dN ; n >= N+1 ; n -- ) {
    sc = Gnm1 ;
    Gnm1 = (4.0*(n-1)*chi*Gnm1 - (2.0*n-1)*Gn)/(2.0*n-3) ;
    Gn = sc ;
  }
  G[N*str] = 1.0 ; G[(N-1)*str] = Gnm1/Gn ;
  /* G[ N   *str] *= 1e-6 ; G[(N-1)*str] *= 1e-6 ; */
  for ( n = N ; n >= 2 ; n -- ) {
    G[(n-2)*str] = (4.0*(n-1)*chi*G[(n-1)*str] - (2.0*n-1)*G[n*str])/(2.0*n-3) ;
  }

  /*value for G[0] from elliptic integral*/
  mu = SQRT(2.0/(1.0+chi)) ;
  AFMM_FUNCTION_NAME(afmm_elliptic_KE)(mu, &K, &E) ;
  sc = mu*K*0.5*M_1_PI/SQRT(r*r1)/G[0] ;

  /*scale Green's functions*/
  for ( n = 0 ; n <= N ; n ++ ) G[n*str] *= sc ;

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_gfunc_comp)(gint N,
						 AFMM_REAL r,
						 AFMM_REAL r1,
						 AFMM_REAL z,
						 AFMM_REAL *G)
/*
 * numerical evaluation of modal Green's functions
 *
 */

{
  gint nt = 8192, i, n ;
  AFMM_REAL R, t ;

  for ( n = 0 ; n <= N ; n ++ ) G[n] = 0.0 ;

  for ( i = 0 ; i < nt ; i ++ ) {
    t = 2.0*M_PI*i/nt ;
    R = SQRT(r*r + r1*r1 - 2*r*r1*COS(t) + z*z) ;
    for ( n = 0 ; n <= N ; n ++ ) {
      G[n] += COS(n*t)/R ;
    }
  }
  
  for ( n = 0 ; n <= N ; n ++ ) G[n] *= 0.5/nt ;

  return 0 ;
}

/* static void laplace_recursion_r(AFMM_REAL *dG, gint nd, gint N, */
/* 				AFMM_REAL r, gint i, gint j, gint k) */

/* { */
/*   gint idxo, idx1, idx2, n, u ; */
/*   AFMM_REAL rp ; */
  
/*   g_assert(i >= 2) ; */
  
/*   /\*output index*\/ */
/*   idxo = afmm_derivative_offset(i+j+k) + afmm_derivative_index_ijk(i,j,k) ; */

/*   for ( n = 0 ; n <= N ; n ++ ) { */
/*     idx1 = afmm_derivative_offset(i+j+k) + */
/*       afmm_derivative_index_ijk(i-2,j,k+2) ; */
/*     dG[n*nd + idxo] = -dG[n*nd+idx1]*(k+2)*(k+1)/(i-1)/i ; */
/*     rp = 1.0/r ; */
/*     for ( u = 0 ; u <= i-2 ; u ++ ) { */
/*       idx1 = afmm_derivative_offset(i+j+k-u-2) + */
/* 	afmm_derivative_index_ijk(i-u-2,j,k) ; */
/*       idx2 = afmm_derivative_offset(i+j+k-u-1) + */
/* 	afmm_derivative_index_ijk(i-u-1,j,k) ; */
/*       /\* dG[n*nd + idxo] += (n*n*dG[n*nd+idx1]*rp/r*(u+1) -  *\/ */
/*       /\* 			  dG[n*nd+idx2]*rp*(i-u-1))/(i-1)/i ; *\/ */
/*       dG[n*nd + idxo] += (n*n*dG[n*nd+idx1]/r*(u+1) -  */
/* 			  dG[n*nd+idx2]*(i-u-1))*rp/(i-1)/i ; */
/*       rp /= -r ; */
/*     } */
/*   } */

/*   return ; */
/* } */

static void laplace_recursion_r(AFMM_REAL *dG, gint nd, gint N,
				AFMM_REAL r, gint i, gint j, gint k)

{
  gint idxo, idx1, idx2, n, u ;
  AFMM_REAL rp, sc ;
  
  g_assert(i >= 2) ;
  
  /*output index*/
  idxo = afmm_derivative_offset(i+j+k) + afmm_derivative_index_ijk(i,j,k) ;

  idx1 = afmm_derivative_offset(i+j+k) +
    afmm_derivative_index_ijk(i-2,j,k+2) ;
  n = N+1 ; sc = -(AFMM_REAL)(k+2)*(k+1)/(i-1)/i ;
#ifdef AFMM_SINGLE_PRECISION
  blaswrap_scopy(n, &(dG[idx1]), nd, &(dG[idxo]), nd) ;
  blaswrap_sscal(n, sc, &(dG[idxo]), nd) ;
#else  /*AFMM_SINGLE_PRECISION*/
  blaswrap_dcopy(n, &(dG[idx1]), nd, &(dG[idxo]), nd) ;
  blaswrap_dscal(n, sc, &(dG[idxo]), nd) ;
#endif  /*AFMM_SINGLE_PRECISION*/

  rp = 1.0/r ;
  for ( u = 0 ; u <= i-2 ; u ++ ) {
    idx1 = afmm_derivative_offset(i+j+k-u-2) +
      afmm_derivative_index_ijk(i-u-2,j,k) ;
    idx2 = afmm_derivative_offset(i+j+k-u-1) +
      afmm_derivative_index_ijk(i-u-1,j,k) ;
    /*this works but is much slower than the loop*/
/* #ifdef AFMM_SINGLE_PRECISION */
/*     g_assert_not_reached() ; */
/* #else /\*AFMM_SINGLE_PRECISION*\/ */
/*     gint i0 = 0, i1 = 1 ; */
/*     AFMM_REAL d1 = 1 ; */
/*     sc = rp/r*(u+1)/(i-1)/i ; */
/*     blaswrap_dgbmv(FALSE, n, n, i0, i0, sc, */
/* 		   (gdouble *)AFMM_N_SQUARED, i1, */
/* 		   &(dG[idx1]), nd, d1, &(dG[idxo]), nd) ; */
/*     sc2 = -rp*(i-u-1)/(i-1)/i ; */
/*     blaswrap_daxpy(n, sc2, &(dG[idx2]), nd, &(dG[idxo]), nd) ; */
    
/* #endif /\*AFMM_SINGLE_PRECISION*\/ */
    for ( n = 0 ; n <= N ; n ++ ) {
      dG[n*nd + idxo] += (n*n*dG[n*nd+idx1]/r*(u+1) -
			  dG[n*nd+idx2]*(i-u-1))*rp/(i-1)/i ;
    }
    rp /= -r ;
  }

  return ;
}

static void laplace_recursion_r1(AFMM_REAL *dG, gint nd, gint N,
				 AFMM_REAL r1,
				 gint i, gint j, gint k)
  
{
  gint idxo, idx1, idx2, n, u ;
  AFMM_REAL r1p ;
  
  g_assert(j >= 2) ;

  /* fprintf(stderr, "b: %d %d %d\n", i, j, k) ; */
  
  /*output index*/
  idxo = afmm_derivative_offset(i+j+k) + afmm_derivative_index_ijk(i,j,k) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    idx1 = afmm_derivative_offset(i+j+k) +
      afmm_derivative_index_ijk(i,j-2,k+2) ;
    dG[n*nd + idxo] = -dG[n*nd+idx1]*(k+2)*(k+1)/(j-1)/j ;
    r1p = 1.0/r1 ;
    for ( u = 0 ; u <= j-2 ; u ++ ) {
      idx1 = afmm_derivative_offset(i+j+k-u-2) +
	afmm_derivative_index_ijk(i,j-u-2,k) ;
      idx2 = afmm_derivative_offset(i+j+k-u-1) +
	afmm_derivative_index_ijk(i,j-u-1,k) ;
      dG[n*nd + idxo] += (n*n*dG[n*nd+idx1]*r1p/r1*(u+1) -
			  dG[n*nd+idx2]*r1p*(j-u-1))/j/(j-1) ;
      r1p /= -r1 ;
    }
  }

  return ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(gint N, gint L,
							AFMM_REAL r,
							AFMM_REAL r1,
							AFMM_REAL x,
							AFMM_REAL *dG)

/*
 * size of dG: (N+1)*afmm_derivative_offset(L+1)
 */
  
{
  gint nd, idxx, idxr, idxr1, n, off, l, q, di, i, j, k ;
  gint idxq ;
  AFMM_REAL rp2, rm2, cft ;
  AFMM_REAL drpx[257], drmx[257] ;
  AFMM_REAL drpr[257], drmr[257] ;
  AFMM_REAL drpr1[257], drmr1[257] ;
  AFMM_REAL drpr1x[257], drmr1x[257] ;
  AFMM_REAL drpxx[257], drmxx[257] ;

  g_assert(L <= 256) ;

  /*number of derivatives per modal order*/
  nd = afmm_derivative_offset(L+1) ;

  /*fill in the modal Green's functions (zero order derivatives)*/
  AFMM_FUNCTION_NAME(afmm_laplace_gfunc)(N, r, r1, x, dG, nd) ;

  /*
   * basic derivatives: 
   * 
   * dG_{n}/dx
   * dG_{n}/dr
   * dG_{n}/dr_{1}
   */

  rm2 = (r-r1)*(r-r1) + x*x ;
  rp2 = (r+r1)*(r+r1) + x*x ; 
  
  n = 0 ;
  off = afmm_derivative_offset(1) ;
  idxr  = off + afmm_derivative_index_ijk(1,0,0) ;
  idxr1 = off + afmm_derivative_index_ijk(0,1,0) ;
  idxx  = off + afmm_derivative_index_ijk(0,0,1) ;
  dG[n*nd+idxr] = -dG[n*nd]/2.0/r +
    (n-0.5)/2.0/r*(r*r - r1*r1 - x*x)*
    ((dG[n*nd] + dG[(n+1)*nd])/rp2 + (dG[n*nd] - dG[(n+1)*nd])/rm2) ;
  dG[n*nd+idxr1] = -dG[n*nd]/2.0/r1 +
    (n-0.5)/2.0/r1*(r1*r1 - r*r - x*x)*
    ((dG[n*nd] + dG[(n+1)*nd])/rp2 + (dG[n*nd] - dG[(n+1)*nd])/rm2) ;
  dG[n*nd+idxx] = (n-0.5)*((dG[n*nd] + dG[(n+1)*nd])/rp2 +
			   (dG[n*nd] - dG[(n+1)*nd])/rm2)*x ;

  for ( n = 1 ; n <= N ; n ++ ) {
    dG[n*nd+idxr] = -dG[n*nd]/2.0/r +
      (n-0.5)/2.0/r*(r*r - r1*r1 - x*x)*
      ((dG[n*nd] + dG[(n-1)*nd])/rp2 + (dG[n*nd] - dG[(n-1)*nd])/rm2) ;
    dG[n*nd+idxr1] = -dG[n*nd]/2.0/r1 +
      (n-0.5)/2.0/r1*(r1*r1 - r*r - x*x)*
      ((dG[n*nd] + dG[(n-1)*nd])/rp2 + (dG[n*nd] - dG[(n-1)*nd])/rm2) ;
    dG[n*nd+idxx] = (n-0.5)*((dG[n*nd] + dG[(n-1)*nd])/rp2 +
			     (dG[n*nd] - dG[(n-1)*nd])/rm2)*x ;
  }

  /*basic higher order derivatives*/
  /*d^k/dx^k (x/\rho_{\pm}^{2})*/
  drpx[0]  = x/rp2 ; 
  drpx[1] = (1.0 - 2.0*x*drpx[0])/rp2 ;
  drmx[0]  = x/rm2 ;
  drmx[1] = (1.0 - 2.0*x*drmx[0])/rm2 ;
  drpx[2] = -2.0*(2.0*x*drpx[1] + drpx[0])/rp2 ;
  drmx[2] = -2.0*(2.0*x*drmx[1] + drmx[0])/rm2 ;

  /*d^k/dr^k [(r^{2}-r_{1}^{2}-x^{2})/\rho_{\pm}^{2}]*/
  drpr[0] = (r*r - r1*r1 - x*x)/rp2 ;
  drpr[1] = -2.0*x*(1.0 + drpr[0])/rp2 ;
  drpr[2] = -(2.0 + 4.0*x*drpr[1] + 2.0*drpr[0])/rp2 ;
  drmr[0] = (r*r - r1*r1 - x*x)/rm2 ;
  drmr[1] = -2.0*x*(1.0 + drmr[0])/rm2 ;
  drmr[2] = -(2.0 + 4.0*x*drmr[1] + 2.0*drmr[0])/rm2 ;

  /*d^k/dr^k [(r_{1}^{2}-r^{2}-x^{2})/\rho_{\pm}^{2}]*/
  drpr1[0] = (r1*r1 - r*r - x*x)/rp2 ;
  drpr1[1] = -2.0*x*(1.0 + drpr1[0])/rp2 ;
  drpr1[2] = -(2.0 + 4.0*x*drpr1[1] + 2.0*drpr1[0])/rp2 ;
  drmr1[0] = (r1*r1 - r*r - x*x)/rm2 ;
  drmr1[1] = -2.0*x*(1.0 + drmr1[0])/rm2 ;
  drmr1[2] = -(2.0 + 4.0*x*drmr1[1] + 2.0*drmr1[0])/rm2 ;

  /*d^k/dx^k d/dr_{1}[(r_{1}^{2}-r^{2}-x^{2})/\rho_{\pm}^{2}]*/
  drpr1x[0] = -2.0*r*((r+r1)*(r+r1)-x*x)/rp2/rp2 ;
  drpr1x[1] = 4*r*x - 4*x*rp2*drpr1x[0] ;
  drpr1x[1] /= rp2*rp2 ;
  drpr1x[2] = 4*r - 2*4*x*rp2*drpr1x[1] - (4*rp2+8*x*x)*drpr1x[0] ;
  drpr1x[2] /= rp2*rp2 ;
  drpr1x[3] = -3*4*x*rp2*drpr1x[2] - 3*(4*rp2+8*x*x)*drpr1x[1] -
    24*x*drpr1x[0] ;
  drpr1x[3] /= rp2*rp2 ;
  drpr1x[4] = -4*4*x*rp2*drpr1x[3] - 6*(4*rp2+8*x*x)*drpr1x[2] -
    4*24*x*drpr1x[1] - 24*drpr1x[0] ;
  drpr1x[4] /= rp2*rp2 ;
  drmr1x[0] = +2.0*r*((r-r1)*(r-r1)-x*x)/rm2/rm2 ;
  drmr1x[1] = -4*r*x - 4*x*rm2*drmr1x[0] ;
  drmr1x[1] /= rm2*rm2 ;
  drmr1x[2] = -4*r - 2*4*x*rm2*drmr1x[1] - (4*rm2+8*x*x)*drmr1x[0] ;
  drmr1x[2] /= rm2*rm2 ;
  drmr1x[3] = -3*4*x*rm2*drmr1x[2] - 3*(4*rm2+8*x*x)*drmr1x[1] -
    24*x*drmr1x[0] ;
  drmr1x[3] /= rm2*rm2 ;
  drmr1x[4] = -4*4*x*rm2*drmr1x[3] - 6*(4*rm2+8*x*x)*drmr1x[2] -
    4*24*x*drmr1x[1] - 24*drmr1x[0] ;
  drmr1x[4] /= rm2*rm2 ;
  
  /*d^k/dx^k [(r_{1}^{2}-r^{2}-x^{2})/\rho_{\pm}^{2}]*/
  drpxx[0] = (r*r - r1*r1 - x*x)/rp2 ;
  drpxx[1] = (-2.0*x - 2.0*x*drpxx[0])/rp2 ;
  drpxx[2] = (-2 - 2*2*x*drpxx[1] - 1*2*drpxx[0])/rp2 ;
  drmxx[0] = (r*r - r1*r1 - x*x)/rm2 ;
  drmxx[1] = (-2.0*x - 2.0*x*drmxx[0])/rm2 ;
  drmxx[2] = (-2 - 2*2*x*drmxx[1] - 1*2*drmxx[0])/rm2 ;

  for ( l = 3 ; l <= L ; l ++ ) {
    drpx [l] = -l*(2.0*x*drpx [l-1] + (l-1)*drpx [l-2])/rp2 ;
    drmx [l] = -l*(2.0*x*drmx [l-1] + (l-1)*drmx [l-2])/rm2 ;
    drpr [l] = -l*(2.0*x*drpr [l-1] + (l-1)*drpr [l-2])/rp2 ;
    drmr [l] = -l*(2.0*x*drmr [l-1] + (l-1)*drmr [l-2])/rm2 ;
    drpr1[l] = -l*(2.0*x*drpr1[l-1] + (l-1)*drpr1[l-2])/rp2 ;
    drmr1[l] = -l*(2.0*x*drmr1[l-1] + (l-1)*drmr1[l-2])/rm2 ;
    drpxx[l] = -l*(2.0*x*drpxx[l-1] + (l-1)*drpxx[l-2])/rp2 ;
    drmxx[l] = -l*(2.0*x*drmxx[l-1] + (l-1)*drmxx[l-2])/rm2 ;
  }

  for ( l = 5 ; l <= L ; l ++ ) {
    drpr1x[l] = -(l*4*x*rp2*drpr1x[l-1] + l*(l-1)/2*(4*rp2+8*x*x)*drpr1x[l-2] +
		  l*(l-1)*(l-2)/6*24*x*drpr1x[l-3] +
		  l*(l-1)*(l-2)*(l-3)/24*24*drpr1x[l-4]) ;
    drpr1x[l] /= rp2*rp2 ;
    drmr1x[l] = -(l*4*x*rm2*drmr1x[l-1] + l*(l-1)/2*(4*rm2+8*x*x)*drmr1x[l-2] +
		  l*(l-1)*(l-2)/6*24*x*drmr1x[l-3] +
		  l*(l-1)*(l-2)*(l-3)/24*24*drmr1x[l-4]) ;
    drmr1x[l] /= rm2*rm2 ;    
  }
  
  /*
   * d^{l}G_{n}/dx^{l}/l!
   */  
  for ( l = 1 ; l <= L-1 ; l ++ ) {
    idxx  = afmm_derivative_offset(l+1) + afmm_derivative_index_ijk(0,0,l+1) ;
    di = 1 ;
    for ( n = 0 ; n <= N ; n ++ ) {
      cft = 1.0/(l+1) ;
      dG[n*nd + idxx] = 0.0 ;
      for ( q = 0  ; q <= l ; q ++ ) {
	idxq  = afmm_derivative_offset(l-q) +
	  afmm_derivative_index_ijk(0,0,l-q) ;
	dG[n*nd + idxx] +=
	  (n-0.5)*cft*(drpx[q]*(dG[n*nd+idxq] + dG[(n+di)*nd+idxq]) +
		       drmx[q]*(dG[n*nd+idxq] - dG[(n+di)*nd+idxq])) ;
	cft /= q+1 ;
      }
      di = -1 ;
    }
  }

  /*
   * d^{l+1}G_{n}/dr     dx^{l}/l!
   * d^{l+1}G_{n}/dr_{1} dx^{l}/l!
   */
  for ( l = 1 ; l <= L-1 ; l ++ ) {
    idxx   = afmm_derivative_offset(l  ) + afmm_derivative_index_ijk(0,0,l) ;
    idxr   = afmm_derivative_offset(l+1) + afmm_derivative_index_ijk(1,0,l) ;
    idxr1  = afmm_derivative_offset(l+1) + afmm_derivative_index_ijk(0,1,l) ;
    di = 1 ;
    for ( n = 0 ; n <= N ; n ++ ) {
      cft = 1.0 ;
      dG[n*nd + idxr ] = -dG[n*nd + idxx]/2.0/r  ;
      dG[n*nd + idxr1] = -dG[n*nd + idxx]/2.0/r1 ;
      for ( q = 0  ; q <= l ; q ++ ) {
	idxq  = afmm_derivative_offset(l-q) +
	  afmm_derivative_index_ijk(0,0,l-q) ;
	dG[n*nd + idxr] +=
	  (n-0.5)*cft*(drpr[q]*(dG[n*nd+idxq] + dG[(n+di)*nd+idxq]) +
		       drmr[q]*(dG[n*nd+idxq] - dG[(n+di)*nd+idxq]))/2.0/r ;
	dG[n*nd + idxr1] +=
	  (n-0.5)*cft*(drpr1[q]*(dG[n*nd+idxq] + dG[(n+di)*nd+idxq]) +
		       drmr1[q]*(dG[n*nd+idxq] - dG[(n+di)*nd+idxq]))/2.0/r1 ;
	cft /= q+1 ;
      }
      di = -1 ;
    }
  }

  /*
   * d^{2+l}G_{n}/dr dr_{1} dx^{l}/l!
   */
  for ( l = 0 ; l <= L-2 ; l ++ ) {
    idxr = afmm_derivative_offset(2+l) + afmm_derivative_index_ijk(1,1,l) ;
    di = 1 ;
    for ( n = 0 ; n <= N ; n ++ ) {
      cft = 1.0 ;
      idxq = afmm_derivative_offset(l+1) + afmm_derivative_index_ijk(0,1,l) ;
      dG[n*nd+idxr] = -dG[n*nd+idxq]/2.0/r ;
      for ( q = 0 ; q <= l ; q ++ ) {
	idxq = afmm_derivative_offset(l-q) +
	  afmm_derivative_index_ijk(0,0,l-q) ;
	dG[n*nd+idxr] += cft*(n-0.5)*
	  ((dG[n*nd+idxq] + dG[(n+di)*nd+idxq])*drpr1x[q] + 
	   (dG[n*nd+idxq] - dG[(n+di)*nd+idxq])*drmr1x[q])/2.0/r ;
	idxq = afmm_derivative_offset(l-q+1) +
	  afmm_derivative_index_ijk(0,1,l-q) ;
	dG[n*nd+idxr] += cft*(n-0.5)*
	  ((dG[n*nd+idxq] + dG[(n+di)*nd+idxq])*drpxx[q] +
	   (dG[n*nd+idxq] - dG[(n+di)*nd+idxq])*drmxx[q])/2.0/r ;
	cft /= q+1 ;
      }
      di = -1 ;
    }
  }

  /*
   * Laplace equation recursion for general (r,r_1,z) derivatives 
   * 
   * all derivatives are scaled on the appropriate factorials
   * 
   * see Strickland and Amos, 1990, A fast solver for systems of
   * axisymmetric ring vortices, SAND90-1925, Sandia Labs for the
   * sequence of evaluations (the recursion itself is not quite the
   * same, since it is based on the Laplace equation rather than the
   * axisymmetric streamfunction)
   */

  for ( l = 2 ; l <= L ; l ++ ) {
    for ( k = l-2 ; k >= 0 ; k -- ) {
      for ( j = 0 ; j <= l-k-2 ; j ++ ) {
	i = l - k - j ;
	laplace_recursion_r(dG, nd, N, r, i, j, k) ;
      }
      /*allow for (1, 1, l-2) case here*/
      for ( j = l-k-1 ; j <= l-k ; j ++ ) {
	i = l - k - j ;
	if ( !(i == 1 && j == 1) ) {
	  laplace_recursion_r1(dG, nd, N, r1, i, j, k) ;
	}
      }
    }
  }

  return nd ;
}

AFMM_REAL AFMM_FUNCTION_NAME(afmm_gfunc_series_eval)(gint L,
						     AFMM_REAL *dG,
						     AFMM_REAL dr,
						     AFMM_REAL dr1,
						     AFMM_REAL dx)

{
  AFMM_REAL G, rp[256], r1p[256], xp[256] ;
  gint i, j, k, l, off ;

  G = 0.0 ;
  rp[0] = r1p[0] = xp[0] = 1.0 ;
  for ( l = 1 ; l <= L ; l ++ ) {
    rp [l] = rp [l-1]*dr  ; 
    r1p[l] = r1p[l-1]*dr1 ; 
    xp [l] = xp [l-1]*dx  ; 
  }
  
  for ( l = 0 ; l <= L ; l ++ ) {
    off = afmm_derivative_offset(l) ;
    for ( i = 0 ; i <= l ; i ++ ) {
      for ( j = 0 ; j <= l-i ; j ++ ) {
	k = l - i - j ;
	G += dG[off+afmm_derivative_index_ijk(i,j,k)]*rp[i]*r1p[j]*xp[k] ;
      }
    }
  }
  
  return G ;
}

gint AFMM_FUNCTION_NAME(afmm_gfunc_coefficients_write)(FILE *f, gint L,
						       AFMM_REAL *dG)

{
  gint i, j, k, l, off ;
  
  for ( l = 0 ; l <= L ; l ++ ) {
    off = afmm_derivative_offset(l) ;
    for ( i = 0 ; i <= l ; i ++ ) {
      for ( j = 0 ; j <= l-i ; j ++ ) {
	k = l - i - j ;
	fprintf(f, "%d %d %d %lg\n", i, j, k, 
		dG[off+afmm_derivative_index_ijk(i,j,k)]) ;
      }
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_field_direct)(AFMM_REAL *rzsrc,
						   AFMM_REAL *src, gint sdist,
						   gint ns,
						   gint nsrc,
						   gint N,
						   AFMM_REAL *rzfld,
						   AFMM_REAL *fld, gint fdist,
						   gint nfld,
						   gboolean zero,
						   AFMM_REAL *work)

/*
 * indexing: (real, imaginary part of) amplitude of mode n of source
 * component j at point i is
 * 
 * src[i*ns*sdist + j*sdist + 2*n+(0,1)]
 *
 * and similar for field
 * 
 * dist must be greater than or equal to 2*N+2 (and ideally 2*N+4 for
 * fftw compatability)
 */
  
{
  gint i, j, k, n ;
  AFMM_REAL *G ;

  /* fprintf(stderr, "%s\n", __FUNCTION__) ; */
  
  if ( zero ) {
    memset(fld, 0, ns*fdist*nfld*sizeof(AFMM_REAL)) ;
  }

  g_assert(sdist >= 2*N+2) ;
  g_assert(fdist >= 2*N+2) ;
  
  G = work ;
  for ( j = 0 ; j < nfld ; j ++ ) {
    for ( i = 0 ; i < nsrc ; i ++ ) {
      AFMM_FUNCTION_NAME(afmm_laplace_gfunc)(N,
 					     rzfld[2*j+0],  rzsrc[2*i+0],
					     rzfld[2*j+1] - rzsrc[2*i+1],
					     G, 1) ;
      for ( k = 0 ; k < ns ; k ++ ) {
	for ( n = 0 ; n <= N ; n ++ ) {
	  fld[j*ns*fdist + k*fdist + 2*n+0] +=
	    G[n]*src[i*ns*sdist + k*sdist + 2*n+0] ;
	  fld[j*ns*fdist + k*fdist + 2*n+1] +=
	    G[n]*src[i*ns*sdist + k*sdist + 2*n+1] ;
	}
      }
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_modes_to_field)(AFMM_REAL *modes,
						     gint ns, gint nm,
						     gint N,
#ifdef AFMM_SINGLE_PRECISION
						     fftwf_plan plan
#else /*AFMM_SINGLE_PRECISION*/
						     fftw_plan plan
#endif /*AFMM_SINGLE_PRECISION*/
						     )

{
#ifdef AFMM_SINGLE_PRECISION
  g_assert_not_reached() ; /*untested code (needs checking for fftwf_plan)*/
  fftwf_execute_dft_c2r(plan, (fftwf_complex *)modes, modes) ;
#else  /*AFMM_SINGLE_PRECISION*/
  fftw_execute_dft_c2r(plan, (fftw_complex *)modes, modes) ;
#endif

  return 0 ;
}

#ifdef AFMM_SINGLE_PRECISION
fftwf_plan
#else /*AFMM_SINGLE_PRECISION*/
fftw_plan
#endif /*AFMM_SINGLE_PRECISION*/
AFMM_FUNCTION_NAME(afmm_laplace_modes_to_field_plan)(AFMM_REAL *modes,
							       gint ns, gint np,
							       gint dist,
							       gint N)

{
  gint rank, n, howmany, istride, idist, ostride, odist ;
  guint flags ;

  rank = 1 ;
  n = 2*N + 2 ;
  /*total number of transforms is number of sources times number of points*/
  howmany = np*ns ; 
  idist = dist/2 ; istride = 1 ;
  odist =   dist ; ostride = 1 ;
  flags = FFTW_MEASURE ;

#ifdef AFMM_SINGLE_PRECISION
  g_assert_not_reached() ; /*untested code (needs checking for fftwf_plan)*/
  fftwf_plan plan = fftwf_plan_many_dft_c2r(rank, &n, howmany,
					    (fftwf_complex *)modes,
					    NULL, istride, idist,
					    modes, NULL, ostride, odist,
					    flags) ;

  return plan ;
#else /*AFMM_SINGLE_PRECISION*/
  fftw_plan plan = fftw_plan_many_dft_c2r(rank, &n, howmany,
					  (fftw_complex *)modes,
					  NULL, istride, idist,
					  modes, NULL, ostride, odist,
					  flags) ;
  return plan ;
#endif /*AFMM_SINGLE_PRECISION*/
  
}

gint AFMM_FUNCTION_NAME(afmm_laplace_field_eval)(gint N, gint L,
						 AFMM_REAL *dG, gint nd,
						 AFMM_REAL *S, gint sdist,
						 gint ns,
						 AFMM_REAL *f, gint fdist,
						 gboolean forward,
						 gboolean outward,
						 gboolean zero)

/*
 * data packing:
 * 
 * dG: see greens_func_derivatives
 *
 * S: see moments.c
 *
 * f: f[s*fdist], fdist >= 2*N+2 (but 2*N+4 is better for fftw)
 */
  
  
{
  gint i, j, k, l, n, off, idx, sgn, ds, di, dj ;
  AFMM_REAL sc ;

  g_assert(fdist > 2*N+1) ;
  
  ds = 1 ;
  /*sign of dz1 terms*/
  if ( forward ) ds = -1 ;
  /*swapping r and r1 derivatives*/
  if ( outward ) { di = 0 ; dj = 1 ; } else { di = 1 ; dj = 0 ; } 

  if ( zero ) {
    for ( n = 0 ; n <= N ; n ++ ) {
      for ( i = 0 ; i < ns ; i ++ ) {
	f[i*fdist+2*n+0] = f[i*fdist+2*n+1] = 0.0 ;
      }
    }
  }

  for ( n = 0 ; n <= N ; n ++ ) {
    for ( l = 0 ; l <= L ; l ++ ) {
      off = afmm_derivative_offset(l) ;
      sgn = 1 ;
      for ( k = 0 ; k <= l ; k ++ ) {
	j = l - k ;
	idx = afmm_derivative_index_ij(j,k) ;
	sc = sgn*dG[n*nd+off+afmm_derivative_index_ijk(j*di,j*dj,k)] ;
	for ( i = 0 ; i < ns ; i ++ ) {
	  f[i*fdist+2*n+0] += sc*S[i*(N+1)*sdist + n*sdist + 2*idx+0] ;
	  f[i*fdist+2*n+1] += sc*S[i*(N+1)*sdist + n*sdist + 2*idx+1] ;
	}
	sgn *= ds ;
      }
    }
  }

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_source_to_local)(gint N, gint L,
						      AFMM_REAL *dG, gint nd,
						      AFMM_REAL *S, gint sdist,
						      gint LS,
						      gint ns,
						      AFMM_REAL *P, gint pdist,
						      gint LP,
						      gboolean forward,
						      gboolean outward)

/*
 * expansion to order LS on source side; LP on potential side;
 * derivatives for shift operation available up to order L
 *
 */
  
{
  gint i, j, k, ls, lp, m, n, s, off, idxs, idxp, ssgn, ds, fsgn, df, di, dj ;
  AFMM_REAL sc ;
  
  g_assert(LS + LP <= L) ;
  /*sign of dz1 terms*/
  ds = ( forward ? -1 : 1 ) ;
  df = -ds ;
  /*swapping r and r1 derivatives*/
  if ( outward ) { di = 0 ; dj = 1 ; } else { di = 1 ; dj = 0 ; } 
  n = N+1 ;
  
  /*loop on source coefficient indices*/
  for ( ls = 0 ; ls <= LS ; ls ++ ) {
    ssgn = 1 ;
    for ( k = 0 ; k <= ls ; k ++ ) {
      j = ls - k ;
      idxs = afmm_derivative_index_ij(j,k) ;
      /*loop on potential coefficient indices*/
      for ( lp = 0 ; lp <= LP ; lp ++ ) {
	fsgn = 1 ;
	for ( m = 0 ; m <= lp ; m ++ ) {
	  i = lp - m ;
	  idxp = afmm_derivative_index_ij(i,m) ;
	  off = afmm_derivative_offset(ls+lp) +
	    afmm_derivative_index_ijk(i*dj+j*di,j*dj+i*di,k+m) ;
	  /* g_assert(k+m < 34) ; */
	  sc = fsgn*ssgn*afmm_binomial(k+m,k) ;
	  
	  for ( s = 0 ; s < ns ; s ++ ) {
	    for ( n = 0 ; n <= N ; n ++ ) {
	      P[s*(N+1)*pdist + n*pdist + 2*idxp+0] +=
		sc*S[s*(N+1)*sdist + n*sdist + 2*idxs+0]*dG[off + n*nd] ;
	      P[s*(N+1)*pdist + n*pdist + 2*idxp+1] +=
		sc*S[s*(N+1)*sdist + n*sdist + 2*idxs+1]*dG[off + n*nd] ;
	    }
	  }
	  fsgn *= df ;
	}
      }
      ssgn *= ds ;
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrix)(gint N, gint L,
						 AFMM_REAL *dG, gint nd,
						 gint LS,
						 gint LP,
						 gboolean forward,
						 gboolean outward,
						 AFMM_REAL *S2L)

/*
 * expansion to order LS on source side; LP on potential side;
 * derivatives for shift operation available up to order L
 *
 */
  
{
  gint i, j, k, ls, lp, m, n, off, idxs, idxp, ssgn, ds, fsgn, df, di, dj ;
  gint ns, np, str ;
  AFMM_REAL sc ;
  
  g_assert(LS + LP <= L) ;
  /*sign of dz1 terms*/
  ds = ( forward ? -1 : 1 ) ;
  df = -ds ;
  /*swapping r and r1 derivatives*/
  if ( outward ) { di = 0 ; dj = 1 ; } else { di = 1 ; dj = 0 ; } 

  /*matrix size*/
  ns = afmm_derivative_offset_2(LS+1) ;
  np = afmm_derivative_offset_2(LP+1) ;
  n = N+1 ; str = ns*np ;
  
  /*loop on source coefficient indices*/
  for ( ls = 0 ; ls <= LS ; ls ++ ) {
    ssgn = 1 ;
    for ( k = 0 ; k <= ls ; k ++ ) {
      j = ls - k ;
      idxs = afmm_derivative_index_ij(j,k) ;
      /*loop on potential coefficient indices*/
      for ( lp = 0 ; lp <= LP ; lp ++ ) {
	fsgn = 1 ;
	for ( m = 0 ; m <= lp ; m ++ ) {
	  i = lp - m ;
	  idxp = afmm_derivative_index_ij(i,m) ;
	  off = afmm_derivative_offset(ls+lp) +
	    afmm_derivative_index_ijk(i*dj+j*di,j*dj+i*di,k+m) ;
	  /* g_assert(k+m < 34) ; */
	  sc = fsgn*ssgn*afmm_binomial(k+m,k) ;

#ifdef AFMM_SINGLE_PRECISION
	  blaswrap_saxpy(n, sc, &(dG[off]), nd, &(S2L[idxp*ns + idxs]), str) ;
#else  /*AFMM_SINGLE_PRECISION*/
	  blaswrap_daxpy(n, sc, &(dG[off]), nd, &(S2L[idxp*ns + idxs]), str) ;
#endif  /*AFMM_SINGLE_PRECISION*/
	  /* for ( n = 0 ; n <= N ; n ++ ) { */
	  /*   S2L[n*ns*np + idxp*ns + idxs] = sc*dG[off+n*nd] ; */
	  /* } */
	  fsgn *= df ;
	}
      }
      ssgn *= ds ;
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(gint N,
						AFMM_REAL *S2L,
						AFMM_REAL *S, gint sdist,
						gint LS,
						gint ns,
						AFMM_REAL *P, gint pdist,
						gint LP)

/*
 * note that P is *incremented* by the shifted coefficients
 */

{
  gint n, s, s2lr, s2lc, i2 = 2 ;
  AFMM_REAL *s2l, d1 = 1 ;

  s2lr = afmm_derivative_offset_2(LP+1) ;
  s2lc = afmm_derivative_offset_2(LS+1) ;
  for ( n = 0 ; n <= N ; n ++ ) {
    s2l = &(S2L[n*s2lr*s2lc]) ;
    for ( s = 0 ; s < ns ; s ++ ) {
#ifdef AFMM_SINGLE_PRECISION
      blaswrap_sgemm(FALSE, FALSE, s2lr, i2, s2lc, d1, s2l, s2lc,
		     &(S[s*(N+1)*sdist + n*sdist]), i2, d1,
		     &(P[s*(N+1)*pdist + n*pdist]), i2) ;
#else /*AFMM_SINGLE_PRECISION*/
      blaswrap_dgemm(FALSE, FALSE, s2lr, i2, s2lc, d1, s2l, s2lc,
		     &(S[s*(N+1)*sdist + n*sdist]), i2, d1,
		     &(P[s*(N+1)*pdist + n*pdist]), i2) ;
#endif /*AFMM_SINGLE_PRECISION*/
    }
  }
  
  return 0 ;
}
