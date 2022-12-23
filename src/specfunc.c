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

#ifdef AFMM_SINGLE_PRECISION
#define ELLIP_TOL 1e-7
#define EPS 1.0e-7
#define TINY 1.0e-32
#define NPRE 300
#else /*AFMM_SINGLE_PRECISION*/
#define ELLIP_TOL 1e-14
#define EPS 1.0e-16
#define TINY 1.0e-280
#define NPRE 300
#endif /*AFMM_SINGLE_PRECISION*/

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
 * Legendre functions Q_{n-1/2}(\chi) and derivatives
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

gint AFMM_FUNCTION_NAME(afmm_legendre_nm12)(AFMM_REAL x, gint N,
					    AFMM_REAL *P, gint pstr,
					    AFMM_REAL *Q, gint qstr)

/*
 * Legendre functions P_{\nu}, Q_{\nu}, \nu=n-1/2, n=0, ..., N
 *
 * Based on paper and code of 
 * 
 * Gil and Segura, Evaluation of Legendre functions of argument
 * greater than one, Computer Physics Communications
 * 105(2--3):273--283, 1997.
 * 
 * https://doi.org/10.1016/S0010-4655(97)00076-3
 * 
 */
  
{
  AFMM_REAL alpha, delta, a1, a2, K, E, cant, B, A, fc, C0, D0, PNp1, QNp1 ;
  gint n, m ;
  
  g_assert(x > 1.0) ;

  alpha = LOG(x + SQRT(x*x - 1.0)) ;

  a1 = SQRT((x-1.0)/(x+1.0)) ;
  a2 = SQRT(2.0*SQRT(x*x - 1.0)/(x + SQRT(x*x-1.0))) ;

  AFMM_FUNCTION_NAME(afmm_elliptic_KE)(a1, &K, &E) ;
  P[0*pstr] = 2.0/M_PI*SQRT(2.0/(x+1.0))*K ;
  AFMM_FUNCTION_NAME(afmm_elliptic_KE)(a2, &K, &E) ;
  P[1*pstr] = 2.0/M_PI*SQRT(x+SQRT(x*x-1.0))*E ;

  cant = 0.5*(LOG(2.0*M_PI) + LOG(SINH(alpha))) + LOG(10.0)*NPRE ;

  while ( alpha*N - 0.5*LOG(N - 0.5) > cant ) {
    N -= 3 ;
  }

  /*forward recurrence for P_{n-1/2}*/
  for ( n = 1 ; n < N ; n ++ ) {
    P[(n+1)*pstr] = (2.0*n*x*P[n*pstr] - (n-0.5)*P[(n-1)*pstr])/(n+0.5) ;
  }
  /*spare value of P to avoid array overrun*/
  n = N ;
  PNp1 = (2.0*n*x*P[n*pstr] - (n-0.5)*P[(n-1)*pstr])/(n+0.5) ;

  n = N ;

  m = 0 ;
  B = 2.0*(n+1.0)*x/(n+0.5) ;
  A = 1.0 ;
  fc = TINY ;
  C0 = fc ;
  D0 = 0.0 ;

  do {
    D0 = B + A*D0 ;
    if ( D0 == 0.0 ) D0 = TINY ;
    C0 = B+A/C0 ;
    if ( C0 == 0.0 ) C0 = TINY ;
    D0 = 1.0/D0 ;
    delta = C0*D0 ;
    fc *= delta ;
    m ++ ;
    A = -(1.0 + 1.0/(n + m - 0.5)) ;
    B = 2.0*(n+m+1.0)*x/(n + m + 0.5) ;
  } while (ABS(delta - 1.0) > EPS) ;

  Q[N*qstr] = 1.0/(PNp1 - fc*P[N*pstr])/(N+0.5) ;
  /*spare value of Q*/
  QNp1 = Q[N*qstr]*fc ;

  n = N ;
  Q[(n-1)*qstr] = (2.0*n*x*Q[n*qstr] - (n+0.5)*QNp1)/(n - 0.5) ;
  for ( n = N-1 ; n > 0 ; n -- ) {
    Q[(n-1)*qstr] = (2.0*n*x*Q[n*qstr] - (n+0.5)*Q[(n+1)*qstr])/(n - 0.5) ;
  }
  
  return 0 ;
}

AFMM_REAL AFMM_FUNCTION_NAME(afmm_hyperg_2F1)(AFMM_REAL a, AFMM_REAL b,
					      AFMM_REAL c, AFMM_REAL x)

{
  AFMM_REAL F, tol, cft ;
  gint q ;
  
#ifdef AFMM_SINGLE_PRECISION
  tol = 1e-7 ;
#else /*AFMM_SINGLE_PRECISION*/
  tol = 1e-15 ;
#endif /*AFMM_SINGLE_PRECISION*/

  cft = 1.0/GAMMA(c) ;

  F = cft ;

  for ( q = 0 ; (q < 16) && ( cft > tol) ; q ++ ) {
    cft *= (a+q)/(c+q)*(b+q)/(1.0+q)*x ;
    /* (a+q)*(b+q)/(c+q)/(q+1) ;     */
    F += cft ;
  }

  g_assert(cft < tol) ;
  
  return F ;
}

gint AFMM_FUNCTION_NAME(afmm_legendre_Q_rec)(gint N, AFMM_REAL chi,
					     AFMM_REAL *Q, gint str)

/*
 * Method of
 * 
 * An explicit kernel-split panel-based Nystrom scheme
 * for integral equations on axially symmetric surfaces
 * Johan Helsing, Anders Karlsson, J Comp Phys, 2014
 * 
 * http://dx.doi.org/10.1016/j.jcp.2014.04.053
 * 
 */
  
{
  AFMM_REAL csw, K, E, mu, sc, Qn, Qnm1 ;
  gint n, M ;
  
  csw = 1.008 ;
  g_assert(chi > 1.0) ;
  
  /*forward recursion for \chi<=1.008*/
  if ( chi < csw ) {
    mu = SQRT(2.0/(1.0+chi)) ;
    AFMM_FUNCTION_NAME(afmm_elliptic_KE)(mu, &K, &E) ;

    Q[0*str] = mu*K ;
    if ( N == 0 ) return 0 ;
    Q[1*str] = chi*mu*K - (1.0 + chi)*mu*E ;
    if ( N == 1 ) return 0 ;

    for ( n = 2 ; n <= N ; n ++ ) {
      Q[n*str] = (4.0*n-4)/(2.0*n-1)*chi*Q[(n-1)*str] -
	(2.0*n-3)/(2.0*n-1)*Q[(n-2)*str] ;
    }
    

    sc = chi*mu*K - (1.0+chi)*mu*E ;
    if ( ABS(sc - Q[1*str]) > 1e-9 ) {
      g_error("%s: recursion failure on Q_1=%lg, should be %lg, "
	      "difference = %lg",
	    __FUNCTION__, sc, Q[1*str], ABS(sc-Q[1*str])) ;
    }
    
    return 0 ;
  }
  
  /*backward recursion for \chi>1.008*/
  M = 80 ;

  Qn = 1.0 ; Qnm1 = 1.0 ;
  Qn *= 1e-9 ; Qnm1 *= 1e-9 ; 
  for ( n = N + M ; n >= N+1 ; n -- ) {
    sc = Qnm1 ;
    Qnm1 = (4.0*(n-1)*chi*Qnm1 - (2.0*n-1)*Qn)/(2.0*n-3) ;
    Qn = sc ;
  }
  Q[N*str] = 1.0 ; Q[(N-1)*str] = Qnm1/Qn ;
  /* G[ N   *str] *= 1e-6 ; G[(N-1)*str] *= 1e-6 ; */
  for ( n = N ; n >= 2 ; n -- ) {
    Q[(n-2)*str] = (4.0*(n-1)*chi*Q[(n-1)*str] - (2.0*n-1)*Q[n*str])/(2.0*n-3) ;
  }
  
  mu = SQRT(2.0/(1.0+chi)) ;
  AFMM_FUNCTION_NAME(afmm_elliptic_KE)(mu, &K, &E) ;
  sc = mu*K/Q[0] ;
  for ( n = 0 ; n <= N ; n ++ ) {
    Q[n*str] *= sc ;
  }

  sc = chi*mu*K - (1.0+chi)*mu*E ;
  if ( ABS(sc - Q[1*str]) > 1e-9 ) {
    g_error("%s: recursion failure on Q_1=%lg, should be %lg, "
	    "difference = %lg",
	    __FUNCTION__, sc, Q[1*str], ABS(sc-Q[1*str])) ;
  }
  
  return 0 ;
}
