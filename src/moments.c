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

#include <blaswrap.h>

#include <afmm.h>

#include "afmm-private.h"

/*
 * data packing for moments:
 * 
 * moment (i,j) indexed by afmm_derivative_index_ij(i,j)
 * 
 * for (i,j) moment of component s of modal order n:
 *
 * n*ns*sdist + s*sdist + afmm_derivative_index_ij(i,j)
 *
 * sdist = 2*afmm_coefficient_number(L)
 *
 * L = maximum order of moments (i+j <= L)
 *
 * total size of moment array:
 *
 * (N+1)*ns*sdist, 0 <= n <= N, 0 < s < ns
 */

gint AFMM_FUNCTION_NAME(afmm_source_moments)(AFMM_REAL dr1, AFMM_REAL dz1,
					     AFMM_REAL *C,
					     gint N, gint L,
					     AFMM_REAL *S, gint sdist)

{
  gint j, k, l, idx, n, i1 = 1, str ;
  AFMM_REAL r1p[64], z1p[64], sc[] = {0, 0} ;

  g_assert(L < 64) ;
  r1p[0] = z1p[0] = 1.0 ;
  str = sdist/2 ;
  for ( l = 0 ; l <= L ; l ++ ) {
    r1p[l+1] = dr1*r1p[l] ;
    z1p[l+1] = dz1*z1p[l] ;
    for ( j = 0 ; j <= l ; j ++ ) {
      k = l - j ;
      idx = afmm_derivative_index_ij(j,k) ;
      g_assert(2*idx < sdist) ;
      n = N+1 ;
      sc[0] = r1p[j]*z1p[k] ;
#ifdef AFMM_SINGLE_PRECISION
      blaswrap_caxpy(n, sc, &(C[0]), i1, &(S[2*idx+0]), str) ;
      /* blaswrap_saxpy(n, sc, &(C[0]), i2, &(S[2*idx+0]), sdist) ; */
      /* blaswrap_saxpy(n, sc, &(C[1]), i2, &(S[2*idx+1]), sdist) ; */
#else /*AFMM_SINGLE_PRECISION*/
      blaswrap_zaxpy(n, sc, &(C[0]), i1, &(S[2*idx+0]), str) ;
      /* blaswrap_daxpy(n, sc, &(C[0]), i2, &(S[2*idx+0]), sdist) ; */
      /* blaswrap_daxpy(n, sc, &(C[1]), i2, &(S[2*idx+1]), sdist) ; */
#endif /*AFMM_SINGLE_PRECISION*/
      
      /* for ( n = 0 ; n <= N ; n ++ ) { */
      /* 	S[n*sdist + 2*idx+0] += C[2*n + 0]*r1p[j]*z1p[k] ; */
      /* 	S[n*sdist + 2*idx+1] += C[2*n + 1]*r1p[j]*z1p[k] ; */
      /* } */
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_expansion_eval)(AFMM_REAL dr, AFMM_REAL dz,
					     gint N, gint L,
					     AFMM_REAL *P, gint pdist,
					     gint ns,
					     AFMM_REAL *f, gint fdist)

{
  gint i, j, l, idx, n, s, i1 = 1, pstr ;
  AFMM_REAL rp[64], zp[64], sc[2] = {0, 0} ;

  g_assert(L < 64) ;
  rp[0] = zp[0] = 1.0 ;
  pstr = pdist/2 ;
  for ( l = 0 ; l <= L ; l ++ ) {
    rp[l+1] = dr*rp[l] ;
    zp[l+1] = dz*zp[l] ;
    for ( i = 0 ; i <= l ; i ++ ) {
      j = l - i ;
      idx = afmm_derivative_index_ij(i,j) ;
      n = N + 1 ;
      sc[0] = rp[i]*zp[j] ;
      for ( s = 0 ; s < ns ; s ++ ) {
#ifdef AFMM_SINGLE_PRECISION
	blaswrap_caxpy(n, sc, &(P[s*(N+1)*pdist + 2*idx]), pstr,
		       &(f[s*fdist]), i1) ;
	/* blaswrap_saxpy(n, sc, &(P[s*(N+1)*pdist + 2*idx + 0]), pdist, */
	/* 	       &(f[s*fdist+0]), i2) ; */
	/* blaswrap_saxpy(n, sc, &(P[s*(N+1)*pdist + 2*idx + 1]), pdist, */
	/* 	       &(f[s*fdist+1]), i2) ; */
#else /*AFMM_SINGLE_PRECISION*/
	blaswrap_zaxpy(n, sc, &(P[s*(N+1)*pdist + 2*idx]), pstr,
		       &(f[s*fdist]), i1) ;
	/* blaswrap_daxpy(n, sc, &(P[s*(N+1)*pdist + 2*idx + 0]), pdist, */
	/* 	       &(f[s*fdist+0]), i2) ; */
	/* blaswrap_daxpy(n, sc, &(P[s*(N+1)*pdist + 2*idx + 1]), pdist, */
	/* 	       &(f[s*fdist+1]), i2) ; */
#endif /*AFMM_SINGLE_PRECISION*/
	/* for ( n = 0 ; n <= N ; n ++ ) { */
	/*   f[s*fdist + 2*n+0] += */
	/*     P[s*(N+1)*pdist + n*pdist + 2*idx+0]*rp[i]*zp[j] ; */
	/*   f[s*fdist + 2*n+1] += */
	/*     P[s*(N+1)*pdist + n*pdist + 2*idx+1]*rp[i]*zp[j] ; */
	/* } */
      }
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_moments_shift)(gint N,
					    gint Mi,
					    AFMM_REAL *Si, gint idist,
					    gint ns,
					    AFMM_REAL dr, AFMM_REAL dz,
					    gint Mo,
					    AFMM_REAL *So, gint odist)

/*
 * shift moments Si up to order Mi, generated about point (r,z) to
 * increment moments about point (r+dr,z+dz) and store in So up to
 * order Mo
 */
  
{
  gint i, j, idxi, idxo, n, s, q, u, mo ;
  AFMM_REAL drp, dzp, sc ;
  
  for ( mo = 0 ; mo <= Mo ; mo ++ ) {
    for ( i = 0 ; i <= mo ; i ++ ) {
      j = mo - i ;
      idxo = afmm_derivative_index_ij(i,j) ;
      drp = 1.0 ;
      for ( q = 0 ; q <= i ; q ++ ) {
	/* for ( u = 0 ; u <= j ; u ++ ) { */
	dzp = 1.0 ;
	for ( u = 0 ; (u <= j) && (i-q+j-u <= Mi) ; u ++ ) {
	  /* if ( (mi = i-q + j-u) <= Mi ) { */
	  idxi = afmm_derivative_index_ij(i-q,j-u) ;
	  sc = afmm_binomial(i,q)*afmm_binomial(j,u)*drp*dzp ;
	  n = N+1 ;
	  for ( s = 0 ; s < ns ; s ++ ) {
#ifdef AFMM_SINGLE_PRECISION
	    blaswrap_saxpy(n, sc,
			   &(Si[s*(N+1)*idist + 2*idxi+0]), idist,
			   &(So[s*(N+1)*odist + 2*idxo+0]), odist) ;
	    blaswrap_saxpy(n, sc,
			   &(Si[s*(N+1)*idist + 2*idxi+1]), idist,
			   &(So[s*(N+1)*odist + 2*idxo+1]), odist) ;
#else /*AFMM_SINGLE_PRECISION*/
	    blaswrap_daxpy(n, sc,
			   &(Si[s*(N+1)*idist + 2*idxi+0]), idist,
			   &(So[s*(N+1)*odist + 2*idxo+0]), odist) ;
	    blaswrap_daxpy(n, sc,
			   &(Si[s*(N+1)*idist + 2*idxi+1]), idist,
			   &(So[s*(N+1)*odist + 2*idxo+1]), odist) ;
#endif /*AFMM_SINGLE_PRECISION*/
	    /* for ( n = 0 ; n <= N ; n ++ ) { */
	    /*   So[s*(N+1)*odist + n*odist + 2*idxo + 0] += */
	    /* 	Si[s*(N+1)*idist + n*idist + 2*idxi + 0]* */
	    /* 	afmm_binomial(i,q)*afmm_binomial(j,u)*drp*dzp ; */
	    /*   So[s*(N+1)*odist + n*odist + 2*idxo + 1] += */
	    /* 	Si[s*(N+1)*idist + n*idist + 2*idxi + 1]* */
	    /* 	afmm_binomial(i,q)*afmm_binomial(j,u)*drp*dzp ; */
	    /* } */
	  }
	  /* } */
	  dzp *= -dz ;
	}
	drp *= -dr ;
      }
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_expansion_shift)(gint N,
					      gint Mi,
					      AFMM_REAL *Pi, gint idist,
					      gint ns,
					      AFMM_REAL dr, AFMM_REAL dz,
					      gint Mo,
					      AFMM_REAL *Po, gint odist)

/*
 * shift expansion Pi up to order Mi, generated about point (r,z) to
 * increment coefficients of expansion about point (r+dr,z+dz) and
 * store in Po up to order Mo
 */
  
{
  gint i, j, idxi, idxo, n, s, q, u, mo ;
  AFMM_REAL drp, dzp, sc ;
  
  for ( mo = 0 ; mo <= Mo ; mo ++ ) {
    for ( i = 0 ; i <= mo ; i ++ ) {
      j = mo - i ;
      idxo = afmm_derivative_index_ij(i,j) ;
      drp = 1.0 ;
      for ( q = 0 ; q <= Mi-i ; q ++ ) {
	dzp = 1.0 ;
	for ( u = 0 ; (u <= Mi-j) && (i+q+j+u <= Mi) ; u ++ ) {
	  g_assert(i+q+j+u <= Mi) ;
	  idxi = afmm_derivative_index_ij(i+q,j+u) ;
	  sc = afmm_binomial(i+q,i)*afmm_binomial(j+u,j)*drp*dzp ;
	  n = N+1 ;
	  for ( s = 0 ; s < ns ; s ++ ) {
#ifdef AFMM_SINGLE_PRECISION
	    blaswrap_saxpy(n, sc,
			   &(Pi[s*(N+1)*idist + 2*idxi+0]), idist,
			   &(Po[s*(N+1)*odist + 2*idxo+0]), odist) ;
	    blaswrap_saxpy(n, sc,
			   &(Pi[s*(N+1)*idist + 2*idxi+1]), idist,
			   &(Po[s*(N+1)*odist + 2*idxo+1]), odist) ;
#else /*AFMM_SINGLE_PRECISION*/
	    blaswrap_daxpy(n, sc,
			   &(Pi[s*(N+1)*idist + 2*idxi+0]), idist,
			   &(Po[s*(N+1)*odist + 2*idxo+0]), odist) ;
	    blaswrap_daxpy(n, sc,
			   &(Pi[s*(N+1)*idist + 2*idxi+1]), idist,
			   &(Po[s*(N+1)*odist + 2*idxo+1]), odist) ;
#endif /*AFMM_SINGLE_PRECISION*/
	  }
	  dzp *= -dz ;
	}
	drp *= -dr ;
      }
    }
  }
  
  return 0 ;
}
