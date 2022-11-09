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

#include <afmm.h>

#include "afmm-private.h"

static void carlson_R(AFMM_REAL x, AFMM_REAL y, AFMM_REAL *RF, AFMM_REAL *RG,
		      AFMM_REAL tol)

  /* 
   * Carlson's method for complete elliptic integrals
   * 
   * https://arxiv.org/pdf/math/9409227.pdf 
   */
  
{
  gint m ;
  AFMM_REAL c, tmp, two_pm ;
  
  tol = 2.7*SQRT(tol) ;
  x = SQRT(x) ; y = SQRT(y) ;
  
  two_pm = 0.25 ;
  c = 0.25*(x+y)*(x+y) ;
  for ( m = 0 ; m < 16 ; m ++ ) {
    tmp = x ;
    x = 0.5*(x+y) ;
    y = SQRT(y*tmp) ;
    two_pm *= 2.0 ;
    c -= two_pm*(x-y)*(x-y) ;
    if ( ABS(x-y) < tol*ABS(x) ) break ;
  }

  *RF = M_PI/(x+y) ;
  *RG = 0.5*c*(*RF) ;
  
  return ;
}

gint AFMM_FUNCTION_NAME(afmm_elliptic_KE)(AFMM_REAL k,
					  AFMM_REAL *K, AFMM_REAL *E)

{
#ifdef AFMM_SINGLE_PRECISION
  AFMM_REAL tol = 1e-7 ;
#else /*AFMM_SINGLE_PRECISION*/
  AFMM_REAL tol = 1e-14 ;
#endif /*AFMM_SINGLE_PRECISION*/

  carlson_R(1.0 - k*k, 1.0, K, E, tol) ;
  *E *= 2.0 ;
  
  return 0 ;
}

AFMM_REAL AFMM_FUNCTION_NAME(afmm_hypergeometric_2F1)(AFMM_REAL a,
						      AFMM_REAL b,
						      AFMM_REAL c,
						      AFMM_REAL z)

{
#ifdef AFMM_SINGLE_PRECISION
  AFMM_REAL tol = 1e-7, f, cfft ;
#else /*AFMM_SINGLE_PRECISION*/
  AFMM_REAL tol = 1e-14, f, cfft ;
#endif /*AFMM_SINGLE_PRECISION*/
  gint q ;

  f = 1 ; cfft = 1.0 ;
  for ( q = 0 ; q < 40 ; q ++ ) {
    cfft *= (a+q)*(b+q)/(c+q)/(q+1)*z ;
    f += cfft ;
    if ( ABS(cfft) < tol ) break ;
  }
  
  return f ;
}

gint AFMM_FUNCTION_NAME(afmm_legendre_Q)(gint N, gint M,
					 AFMM_REAL chi,
					 AFMM_REAL *Q, gint ldq)
/*
 * associated Legendre polynomials Q_{n-1/2}(\chi) and derivatives
 */
  
{
  AFMM_REAL K, E, mu, eta, cfft, p0, p1, p2, q0, q1, r0, Qt ;
  gint n, m ;

  g_assert(ldq > N) ;
  
  mu = SQRT(2.0/(1.0+chi)) ;
  AFMM_FUNCTION_NAME(afmm_elliptic_KE)(mu, &K, &E) ;
  Q[0] = mu*K ; Qt = chi*mu*K - (1.0 + chi)*mu*E ;
  if ( N > 0 ) {
    Q[1] = Qt ;
    Qt = (2.0*chi*Q[1] - 0.5*Q[0])/1.5 ;    
  }

  if ( N > 1 ) {
    /*
     * recursion is unstable for \chi>1 so use hypergeometric function
     * Gradshteyn and Ryzhik 8.852
     */
    eta = ACOSH(chi) ;
    E = EXP(-eta) ; cfft = EXP(-1.5*eta)*M_PI*0.5 ;
    for ( n = 2 ; n <= N ; n ++ ) {      
      cfft *= 0.5*(2*n-1)/n*E ;
      Q[n] = cfft*AFMM_FUNCTION_NAME(afmm_hypergeometric_2F1)(0.5, n+0.5, n+1,
							      E*E) ;
    }
  }

  if ( M == 0 ) return 0 ;
  /*extra Q to be used in calculating the derivative of Q_{N}*/
  if ( N >= 2 ) {
    n = N+1 ;
    cfft *= 0.5*(2*n-1)/n*E ;
    Qt = cfft*AFMM_FUNCTION_NAME(afmm_hypergeometric_2F1)(0.5, n+0.5, n+1,
							  E*E) ;
  }
  
  n = N ;
  Q[ldq*1+n] = (n+0.5)/(chi*chi - 1.0)*(Qt - chi*Q[n]) ;
  for ( n = 0 ; n < N ; n ++ ) {
    Q[ldq*1+n] = (n+0.5)/(chi*chi - 1.0)*(Q[n+1] - chi*Q[n]) ;
  }
  if ( M == 1 ) return 0 ;
  
  /*Glaser, Liu, Rokhlin, 2007 recursion for higher derivatives*/
  p0 = 1.0 - chi*chi ; p1 = -2.0*chi ; p2 = -2.0 ;
  q0 = -2.0*chi ; q1 = -2.0 ;

  for ( n = 0 ; n <= N ; n ++ ) {
    r0 = n*n - 0.25 ;
    for ( m = 0 ; m <= M-2 ; m ++ ) {
      Q[ldq*(m+2)+n] =
	-(p1*m+q0)*Q[ldq*(m+1)+n] - (0.5*m*(m-1)*p2 + m*q1 + r0)*Q[ldq*m+n] ;
      Q[ldq*(m+2)+n] /= p0 ;      
    }
  }

  return 0 ;
}

/* gint AFMM_FUNCTION_NAME(afmm_legendre)(gint N, AFMM_REAL chi, AFMM_REAL *Q, */
/* 				       gint str) */

/* { */
/*   AFMM_REAL Qn, Qnp1 ; */
/*   gint dN ; */

/*   dN = 32 ; */
/*   /\*descending recursion from N+dN to generate starting values for array*\/ */
  
  
/*   return 0 ; */
/* } */
