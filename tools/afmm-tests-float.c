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

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <glib.h>

#include <blaswrap.h>

#include <afmm.h>

#include "afmm-private.h"

gchar *_tests[] = {"elliptic",             /* 0*/
		   "legendre",             /* 1*/
		   "greens",               /* 2*/
		   "series",               /* 3*/
		   "dump",                 /* 4*/
		   "derivative",           /* 5*/
		   "fft",                  /* 6*/
		   "expansion",            /* 7*/
		   "translation",          /* 8*/
		   "tree",                 /* 9*/
		   "shift",                /*10*/
		   "backward",             /*11*/
		   "inward",               /*12*/
		   "morton",               /*13*/
		   "interactions",         /*14*/
		   "s2l",                  /*15*/
		   NULL} ;

static gint parse_test(gchar *test)

{
  gint i ;

  for ( i = 0 ; _tests[i] != NULL ; i ++ ) {
    if ( strcmp(_tests[i], test) == 0 ) return i ;
  }
  
  return -1 ;
}

static gfloat *read_points(FILE *f, gint *np)

{
  gfloat *p ;
  gint i ;
  
  fscanf(f, "%d", np) ;

  p = (gfloat *)g_malloc0((*np)*2*sizeof(gfloat)) ;

  for ( i = 0 ; i < 2*(*np) ; i ++ ) {
    fscanf(f, "%g", &(p[i])) ;
  }
  
  return p ;
}

static void random_coefficients(gfloat *C, gint N)

{
  gint i ;

  C[0] = 2.0*(g_random_double() - 0.5) ;
  C[1] = 0.0 ;
  
  for ( i = 1 ; i <= N ; i ++ ) {
    C[2*i+0] = 2.0*(g_random_double() - 0.5) ;
    C[2*i+1] = 2.0*(g_random_double() - 0.5) ;
  }
  
  return ;
}

static void elliptic_test(void)

{
  gfloat x, K, E ;

  fprintf(stderr, "elliptic integral check\n") ;
  fprintf(stderr, "=======================\n") ;
  
  for ( x = 0 ; x < 1.0 ; x += 1.0/512 ) {
    afmm_elliptic_KE_f(x, &K, &E) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e\n", x, K, E) ;
  }
  
  return ;
}

static void legendre_test(gfloat chi, gint N, gint M, gfloat dx)

{
  gfloat Q0[16384], Q[512], err[64], Qv, x ;
  gint ldq, n, i, j ;

  fprintf(stderr, "associated Legendre function test\n") ;
  fprintf(stderr, "=================================\n") ;
  fprintf(stderr, "N = %d; M = %d\n", N, M) ;
  ldq = N+8 ;

  /*Legendre functions and their derivatives at reference chi*/
  afmm_legendre_Q_f(N, M, chi, Q0, ldq) ;

  for ( n = 0 ; n <= N ; n ++ ) err[n] = 0.0 ;
  
  /* for ( x = - ; x <= 1.0/8 ; x += 1.0/64 ) { */
  for ( j = 0 ; j <= 64 ; j ++ ) {
    x = -dx + 2.0*dx*j/64 ;
    /*Legendre functions at evaluation point*/
    afmm_legendre_Q_f(N, M, chi+x, Q, ldq) ;
    /*Taylor series*/
    fprintf(stdout, "%lg ", x) ;
    for ( n = 0 ; n <= N ; n ++ ) {
      Qv = 0.0 ;
      for ( i = 0 ; i <= M ; i ++ ) {
	Qv += Q0[ldq*i+n]*pow(x, i)/afmm_factorial(i) ;
      }
      fprintf(stdout, "%1.16e ", Qv) ;
      err[n] = MAX(err[n], fabs(Qv-Q[n])) ;
    }
    fprintf(stdout, "\n") ;
  }

  fprintf(stderr, "maximum error:") ;
  for ( n = 0 ; n <= N ; n ++ ) fprintf(stderr, " %lg", err[n]) ;
  fprintf(stderr, "\n") ;
  
  return ;
}

static void greens_func_test(gfloat r, gfloat r1, gfloat z,
			     gint N)

{
  gfloat G[1024], Gc[1024], err ;
  gint n, str ;
  
  fprintf(stderr, "modal Green's function test\n") ;
  fprintf(stderr, "===========================\n") ;
  fprintf(stderr, "(r, r1, z) = %lg, %lg, %lg\n", r, r1, z) ;
  fprintf(stderr, "N = %d\n", N) ;

  str = 3 ;
  afmm_laplace_gfunc_f(N, r, r1, z, Gc, str) ;
  afmm_laplace_gfunc_comp_f(N, r, r1, z, G) ;

  err = 0 ;
  for ( n = 0 ; n <= N ; n ++ ) {
    fprintf(stdout, "%d %lg %lg (%lg)\n", n, G[n], Gc[n*str],
	    ABS(G[n]-Gc[n*str])) ;
    err = MAX(err, fabs(G[n]-Gc[n*str])) ;
  }

  fprintf(stderr, "maximum error: %lg\n", err) ;
  
  return ;
}

static void series_test(gfloat r, gfloat r1, gfloat z,
			gint N, gint L)

{
  gfloat *dG, G[256], Gs, dr, dr1, dz, chi, chimin, chimax, chi1 ;
  gint n, nd, sr, sr1, sz ;
  
  fprintf(stderr, "Green's function series expansion test\n") ;
  fprintf(stderr, "======================================\n") ;
  fprintf(stderr, "(r, r1, z) = %lg, %lg, %lg\n", r, r1, z) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;

  nd = afmm_derivative_offset(L+1) ;
  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z, dG) ;

  dr = dr1 = dz = 0.0 ;
  /* dz = 0.1 ; */
  /* dz =  */
  dz = 0.2 ;
  dr1  = 0.2 ;
  dr = 0.2 ;

  chi = 0.5*(r*r + r1*r1 + z*z)/r1/r ;
  chimin = G_MAXDOUBLE ;
  chimax = 0 ;
  for ( sr = -1 ; sr <= 1 ; sr += 2 ) {
    for ( sr1 = -1 ; sr1 <= 1 ; sr1 += 2 ) {
      for ( sz = -1 ; sz <= 1 ; sz += 1 ) {
	chi1 = 0.5*((r+sr*dr)*(r+sr*dr) + (r1+sr1*dr1)*(r1+sr1*dr1) +
		    (z+sz*dz)*(z+sz*dz))/(r1+sr1*dr1)/(r+sr*dr) ;
	fprintf(stderr, "%d %d %d %lg\n", sr, sr1, sz, chi1) ;
	chimax = MAX(chimax, chi1) ;
	chimin = MIN(chimin, chi1) ;
      }
    }
  }

  fprintf(stderr, "chi: %lg--%lg (%lg)\n", chimin, chimax, chi) ;
  
  afmm_laplace_gfunc_f(N, r+dr, r1+dr1, z+dz, G, 1) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    Gs = afmm_gfunc_series_eval_f(L, &(dG[n*nd]), dr, dr1, dz) ;

    fprintf(stderr, "%d %lg %lg (%lg)\n", n, G[n], Gs, ABS(G[n]-Gs)) ;
  }
  
  return ;
}

static void dump_test(FILE *f, gfloat r, gfloat r1, gfloat z,
		      gint N, gint L, gint n)

{
  gfloat *dG ;
  gint nd ;
  
  fprintf(stderr, "Green's function series coefficients\n") ;
  fprintf(stderr, "====================================\n") ;
  fprintf(stderr, "(r, r1, z) = %lg, %lg, %lg\n", r, r1, z) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "n = %d\n", n) ;
  fprintf(stderr, "L = %d\n", L) ;

  nd = afmm_derivative_offset(L+1) ;

  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;

  
  afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z, dG) ;

  afmm_gfunc_coefficients_write_f(f, L, &(dG[n*nd])) ;
  
  return ;
}

static void derivative_test(gfloat r, gfloat r1, gfloat x,
			    gint N, gint L, gint i, gint j, gint k, gint d)

{
  gfloat *dGp, *dGm, ee, *dGc, dr, dr1, dx, dG, dg, sc ;
  gint nd, idx0, idx, n ;

  g_assert(d >= 1 && d <= 3) ;
  
  fprintf(stderr, "Green's function derivative test\n") ;
  fprintf(stderr, "================================\n") ;
  fprintf(stderr, "(r, r1, x) = %lg, %lg, %lg\n", r, r1, x) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;
  fprintf(stderr, "(i,j,k) = (%d,%d,%d)\n", i, j, k) ;
  fprintf(stderr, "derivative = %d\n", d) ;
  
  nd = afmm_derivative_offset(L+1) ;

  dGp = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  dGm = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  dGc = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;

  ee = 1e-6 ;
  idx0 = afmm_derivative_offset(i+j+k) + afmm_derivative_index_ijk(i,j,k) ;

  dr = dr1 = dx = 0.0 ;
  if ( d == 1 ) {
    /*r derivative*/
    sc = i+1 ;
    dr = ee/2 ;
    idx = afmm_derivative_offset(i+1+j+k) + afmm_derivative_index_ijk(i+1,j,k) ;
  }
  if ( d == 2 ) {
    /*r1 derivative*/
    sc = j + 1 ;
    dr1 = ee/2 ;
    idx = afmm_derivative_offset(i+j+1+k) + afmm_derivative_index_ijk(i,j+1,k) ;
  }
  if ( d == 3 ) {
    /*x derivative*/
    sc = k + 1 ;
    dx = ee/2 ;
    idx = afmm_derivative_offset(i+j+k+1) + afmm_derivative_index_ijk(i,j,k+1) ;
  }

  /* fprintf(stderr, "%lg\n", sc) ; */
  
  afmm_laplace_gfunc_derivatives_f(N, L, r+dr, r1+dr1, x+dx, dGp) ;
  afmm_laplace_gfunc_derivatives_f(N, L, r-dr, r1-dr1, x-dx, dGm) ;
  afmm_laplace_gfunc_derivatives_f(N, L, r, r1, x, dGc) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    dG = (dGp[n*nd+idx0] - dGm[n*nd+idx0])/ee/sc ;
    dg = dGc[n*nd+idx] ;
    fprintf(stderr, "%d %lg %lg (%lg)\n", n, dG, dg, ABS(dG-dg)) ;
  }
  
  return ;
}

static void fft_test(gint ns, gint N)

{
  fftwf_plan plan ;
  gfloat *modes, *C, th, f ;
  gint np, n, i, j, k, dist ;

  np = 2 ;
  dist = 2*N + 4 ;
  /*number of points in physical field*/
  n = 2*N + 2 ;

  fprintf(stderr, "FFT test\n") ;
  fprintf(stderr, "========\n") ;
  fprintf(stderr, "N  = %d\n", N) ;
  fprintf(stderr, "np = %d\n", np) ;
  fprintf(stderr, "ns = %d\n", ns) ;
  fprintf(stderr, "n  = %d\n", n) ;
  fprintf(stderr, "dist = %d\n", dist) ;

  /*randomized modal coefficients*/
  C = (gfloat *)fftw_malloc(np*ns*dist*sizeof(gfloat)) ;
  modes = (gfloat *)fftw_malloc(np*ns*dist*sizeof(gfloat)) ;

  plan = afmm_laplace_modes_to_field_plan_f(C, ns, np, dist, N) ;

  for ( i = 0 ; i < np*ns ; i ++ ) {
    random_coefficients(&(C[i*dist]), N) ;
  }

  /*copy C into modes and generate physical field*/
  memcpy(modes, C, np*ns*dist*sizeof(gfloat)) ;

  afmm_laplace_modes_to_field_f(modes, ns, np, N, plan) ;

  /*evaluate the field from the coefficients*/
  for ( i = 0 ; i < np*ns ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ ) {
      f = C[i*dist+0] ;
      th = (2.0*M_PI*j)/n ;
      for ( k = 1 ; k <= N ; k ++ ) {
	f += 2.0*(C[i*dist+2*k+0]*COS(k*th) -
		  C[i*dist+2*k+1]*SIN(k*th)) ;
      }
      fprintf(stdout, "%1.16e %1.16e %1.16e %1.16e\n",
	      th, f, modes[i*dist+j], ABS(f-modes[i*dist+j])) ;
    }
  }  
  
  return ;
}

static void expansion_test(gfloat r, gfloat z,
			   gfloat r1, gfloat z1,
			   gfloat *rz1, gint nrz1,
			   gint N, gint L,
			   gboolean forward, gboolean outward)

{
  gfloat *C, *S, *dG, rz[2], work[1024], *f, *g ;
  gfloat Z, Z1, R, R1 ;
  gint dist, i, j, n, p, sdist, ns, nd ;
  
  ns = 3 ;
  
  fprintf(stderr, "box source expansion test\n") ;
  fprintf(stderr, "=========================\n") ;
  fprintf(stderr, "(r , z ) = (%g, %g)\n", r, z) ;
  fprintf(stderr, "(r1, z1) = (%g, %g)\n", r1, z1) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;
  fprintf(stderr, "npts = %d\n", nrz1) ;
  fprintf(stderr, "%s shift\n", ( forward ? "forward" : "backward")) ;
  fprintf(stderr, "%s shift\n", ( outward ? "outward" : "inward")) ;
  
  /*local source coefficients*/
  dist = 2*N + 4 ;
  C = (gfloat *)fftw_malloc(nrz1*ns*dist*sizeof(gfloat)) ;
  for ( i = 0 ; i < nrz1*ns ; i ++ ) {
    random_coefficients(&(C[i*dist]), N) ;
  }
  f = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;
  g = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;

  /*local source moments (number of points, number of source terms,
   * number of modes, number of derivatives, times two for complex)
   */
  /*number of points per source point per source component*/
  sdist = 2*afmm_derivative_offset_2(L+1) ;
  S = (gfloat *)g_malloc0(ns*(N+1)*sdist*sizeof(gfloat)) ;
  /*
   * indexing within moment array:
   * sdist = 2*afmm_derivative_offset_2(L+1) entries, 
   * indexed by 2*afmm_derivative_index_ij(j,k) + (0,1)
   *
   * (N+1) blocks of size sdist, n=0, ..., N
   *
   * ns blocks of size (N+1)*sdist, one for each source component
   *
   * (j,k) moment of mode n of source component i is indexed by
   *
   * i*(N+1)*sdist + n*sdist + 2*afmm_derivative_index_ij(j,k)+(0,1)
   */
  
  /*generate local moments*/
  for ( p = 0 ; p < nrz1 ; p ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      afmm_source_moments_f(rz1[2*p+0], rz1[2*p+1],
				&(C[p*ns*dist + i*dist]), N, L,
				&(S[i*(N+1)*sdist]), sdist) ;
    }
  }

  /*modal Green's function derivatives*/
  nd = afmm_derivative_offset(L+1) ;
  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z-z1, dG) ;  

  if ( forward ) {
    Z = z ; Z1 = z1 ;
  } else {
    Z = z1 ; Z1 = z ;
  }

  if ( outward ) {
    R = r ; R1 = r1 ;
  } else {
    R = r1 ; R1 = r ;
  }

  rz[0] = R ; rz[1] = Z ;
  for ( j = 0 ; j < nrz1 ; j ++ ) {
    rz1[2*j+0] += R1 ; rz1[2*j+1] += Z1 ;
  }

  afmm_laplace_field_direct_f(rz1, C, dist, ns, nrz1, N,
				  rz, f, dist, 1, TRUE, work) ;

  afmm_laplace_field_eval_f(N, L, dG, nd, S, sdist, ns, g, dist,
				forward, outward, TRUE) ;
  
  for ( n = 0 ; n <= N ; n ++ ) {
    fprintf(stderr, "%d", n) ;
    for ( i = 0 ; i < ns ; i ++ ) {
      fprintf(stderr, " %g (%g) %g (%g)",
	      f[i*dist+2*n+0],
	      ABS(g[i*dist+2*n+0]-f[i*dist+2*n+0]),
	      f[i*dist+2*n+1],
	      ABS(g[i*dist+2*n+1]-f[i*dist+2*n+1])) ;
    }
    fprintf(stderr, "\n") ;
  }
  
  return ;
}

static void translation_test(gfloat r, gfloat z,
			     gfloat r1, gfloat z1,
			     gfloat *rz1, gint nrz1,
			     gint N, gint L, gint ns,
			     gboolean forward, gboolean outward)

{
  gfloat *C, *P, *S, *dG, rz[32], work[1024], *f, *g ;
  gfloat Z, Z1, R, R1 ;
  gint dist, i, j, n, p, pdist, sdist, nd, LS, LP ;

  LS = L/2 - 1 ; LP = L - LS ;
  rz[0] = 0.1 ; rz[1] = -0.11 ;

  fprintf(stderr, "box-to-box translation test\n") ;
  fprintf(stderr, "===========================\n") ;
  fprintf(stderr, "(r , z ) = (%g, %g)\n", r, z) ;
  fprintf(stderr, "(dr, dz) = (%g, %g)\n", rz[0], rz[1]) ;
  fprintf(stderr, "(r1, z1) = (%g, %g)\n", r1, z1) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;
  fprintf(stderr, "LS = %d\n", LS) ;
  fprintf(stderr, "LP = %d\n", LP) ;  
  fprintf(stderr, "npts = %d\n", nrz1) ;
  fprintf(stderr, "ns = %d\n", ns) ;
  fprintf(stderr, "%s shift\n", ( forward ? "forward" : "backward")) ;
  fprintf(stderr, "%s shift\n", ( outward ? "outward" : "inward")) ;

  /*local source coefficients*/
  dist = 2*N + 4 ;
  C = (gfloat *)fftw_malloc(nrz1*ns*dist*sizeof(gfloat)) ;
  for ( i = 0 ; i < nrz1*ns ; i ++ ) {
    random_coefficients(&(C[i*dist]), N) ;
  }
  f = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;
  g = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;

  /*
   * local source moments (number of points, number of source terms,
   * number of modes, number of derivatives, times two for complex)
   */

  /*number of points per source point per source component*/
  sdist = 2*afmm_derivative_offset_2(LS+1) ;
  S = (gfloat *)g_malloc0(ns*(N+1)*sdist*sizeof(gfloat)) ;
  pdist = 2*afmm_derivative_offset_2(LP+1) ;
  P = (gfloat *)g_malloc0(ns*(N+1)*pdist*sizeof(gfloat)) ;
  
  /*
   * indexing within moment array:
   * sdist = 2*afmm_derivative_offset_2(L+1) entries, 
   * indexed by 2*afmm_derivative_index_ij(j,k) + (0,1)
   *
   * (N+1) blocks of size sdist, n=0, ..., N
   *
   * ns blocks of size (N+1)*sdist, one for each source component
   *
   * (j,k) moment of mode n of source component i is indexed by
   *
   * i*(N+1)*sdist + n*sdist + 2*afmm_derivative_index_ij(j,k)+(0,1)
   */
  
  /*generate local moments*/
  for ( p = 0 ; p < nrz1 ; p ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      afmm_source_moments_f(rz1[2*p+0], rz1[2*p+1],
				&(C[p*ns*dist + i*dist]), N, LS,
				&(S[i*(N+1)*sdist]), sdist) ;
    }
  }

  fprintf(stderr, "  local source moments generated about (r_1, z_1)\n") ;

  /*modal Green's function derivatives*/
  nd = afmm_derivative_offset(L+1) ;
  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z-z1, dG) ;  

  /*translate local source expansion to find local field*/
  afmm_laplace_source_to_local_f(N, L, dG, nd,
				     S, sdist, LS, ns,
				     P, pdist, LP,
				     forward, outward) ;
  fprintf(stderr, "  local field expansion generated about (r,z)\n") ;

  memset(g, 0, ns*dist*sizeof(gfloat)) ;
  afmm_expansion_eval_f(rz[0], rz[1], N, LP, P, pdist, ns, g, dist) ;
  fprintf(stderr, "  local field expansion evaluated at (dr,dz)\n") ;
  
  if ( forward ) {
    Z = z ; Z1 = z1 ;
  } else {
    Z = z1 ; Z1 = z ;
  }

  if ( outward ) {
    R = r ; R1 = r1 ;
  } else {
    R = r1 ; R1 = r ;
  }

  rz[0] += R ; rz[1] += Z ;
  for ( j = 0 ; j < nrz1 ; j ++ ) {
    rz1[2*j+0] += R1 ; rz1[2*j+1] += Z1 ;
  }

  memset(f, 0, ns*dist*sizeof(gfloat)) ;
  afmm_laplace_field_direct_f(rz1, C, dist, ns, nrz1, N,
				  rz, f, dist, 1, TRUE, work) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    fprintf(stderr, "%d", n) ;
    for ( i = 0 ; i < ns ; i ++ ) {
      fprintf(stderr, " %g (%g) %g (%g)",
	      f[i*dist+2*n+0],
	      ABS(g[i*dist+2*n+0]-f[i*dist+2*n+0]),
	      f[i*dist+2*n+1],
	      ABS(g[i*dist+2*n+1]-f[i*dist+2*n+1])) ;
    }
    fprintf(stderr, "\n") ;
  }
  
  return ;
}

static void tree_test(gfloat *rz1, gint nrz1)

{
  gfloat rmin, rmax, zmin, zmax, *p ;
  afmm_box_t *boxes ;
  gint i, j, k ;
  guint depth, nb ;
  afmm_tree_t *tree ;
  
  rmin = 0 ; rmax = 1.0 ;
  zmin = -1 ; zmax = 1 ;
  depth = 5 ;
  
  fprintf(stderr, "tree generation test\n") ;
  fprintf(stderr, "====================\n") ;
  fprintf(stderr, "r = (%g, %g)\n", rmin, rmax) ;
  fprintf(stderr, "z = (%g, %g)\n", zmin, zmax) ;
  fprintf(stderr, "%d points\n", nrz1) ;
  fprintf(stderr, "depth: %u\n", depth) ;
  
  tree = afmm_tree_new_f(rmin, rmax, zmin, zmax, nrz1) ;

  afmm_tree_add_points_f(tree, rz1, 2*sizeof(gfloat), nrz1, FALSE) ;

  for ( i = 1 ; i <= depth ; i ++ ) {
    afmm_tree_refine_f(tree) ;
  }

  nb = 1 << (2*depth) ;
  boxes = tree->boxes[depth] ;
  
  for ( i = 0 ; i < nb ; i ++ ) {
    for ( j = 0 ; j < boxes[i].n ; j ++ ) {
      k = boxes[i].i + j ;
      p = (gfloat *)afmm_tree_point_index(tree, tree->ip[k]) ;      
      fprintf(stdout,
	      "%d %d %g %g\n", i, tree->ip[k], p[0], p[1]) ;
    }
  }
  
  return ;
}

static void shift_test(gfloat r, gfloat z,
		       gfloat r1, gfloat z1,
		       gfloat *rz1, gint nrz1,
		       gint N, gint L, gint ns)

{
  gfloat *C, *P, *Ps, *S, *Ss, *T, *dG, rz[32], work[1024], *f, *g, *gs ;
  gfloat dr, dz ;
  gint dist, i, j, n, p, pdist, psdist, sdist, ssdist, nd, LS, LP, LSs ;
  
  dr = -0.13 ; dz = 0.11 ;
  LS = L/2 - 1 ; LP = L - LS ;
  /*order of shifted moments*/
  LSs = LS + 3 ;
  L += 4 ;
  
  fprintf(stderr, "moment shift test\n") ;
  fprintf(stderr, "=================\n") ;
  fprintf(stderr, "(r , z ) = (%g, %g)\n", r, z) ;
  fprintf(stderr, "(r1, z1) = (%g, %g)\n", r1, z1) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;
  fprintf(stderr, "LS = %d\n", LS) ;
  fprintf(stderr, "LP = %d\n", LP) ;  
  fprintf(stderr, "npts = %d\n", nrz1) ;

  /*local source coefficients*/
  dist = 2*N + 4 ;
  C = (gfloat *)fftw_malloc(nrz1*ns*dist*sizeof(gfloat)) ;
  for ( i = 0 ; i < nrz1*ns ; i ++ ) {
    random_coefficients(&(C[i*dist]), N) ;
  }
  f  = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;
  g  = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;
  gs = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;

  /*local source moments (number of points, number of source terms,
   * number of modes, number of derivatives, times two for complex)
   */
  /*number of points per source point per source component*/
  sdist  = 2*afmm_derivative_offset_2(LS +1) ;
  ssdist = 2*afmm_derivative_offset_2(LSs+1) ;
  S  = (gfloat *)g_malloc0(ns*(N+1)* sdist*sizeof(gfloat)) ;
  Ss = (gfloat *)g_malloc0(ns*(N+1)*ssdist*sizeof(gfloat)) ;
  T = (gfloat *)g_malloc0(ns*(N+1)*ssdist*sizeof(gfloat)) ;
  pdist = 2*afmm_derivative_offset_2(LP+1) ;
  P = (gfloat *)g_malloc0(ns*(N+1)*pdist*sizeof(gfloat)) ;
  psdist = 2*afmm_derivative_offset_2(LP+1) ;
  Ps = (gfloat *)g_malloc0(ns*(N+1)*psdist*sizeof(gfloat)) ;
  
  /*
   * indexing within moment array:
   * sdist = 2*afmm_derivative_offset_2(L+1) entries, 
   * indexed by 2*afmm_derivative_index_ij(j,k) + (0,1)
   *
   * (N+1) blocks of size sdist, n=0, ..., N
   *
   * ns blocks of size (N+1)*sdist, one for each source component
   *
   * (j,k) moment of mode n of source component i is indexed by
   *
   * i*(N+1)*sdist + n*sdist + 2*afmm_derivative_index_ij(j,k)+(0,1)
   */
  
  /*generate local moments*/
  for ( p = 0 ; p < nrz1 ; p ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      afmm_source_moments_f(rz1[2*p+0], rz1[2*p+1],
				&(C[p*ns*dist + i*dist]), N, LS,
				&(S[i*(N+1)*sdist]), sdist) ;
    }
  }

  /*local moments about shifted centre*/
  for ( p = 0 ; p < nrz1 ; p ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      afmm_source_moments_f(rz1[2*p+0]-dr, rz1[2*p+1]-dz,
				&(C[p*ns*dist + i*dist]), N, LSs,
				&(T[i*(N+1)*ssdist]), ssdist) ;
    }
  }

  /*shifted moments*/
  afmm_moments_shift_f(N, LS, S, sdist, ns, dr, dz, LSs, Ss, ssdist) ;

  /*difference in moments*/
  /* gint idx, s ; */
  /* for ( n = 0 ; n <= N ; n ++ ) { */
  /*   for ( p = 0 ; p <= LSs ; p ++ ) { */
  /*     for ( i = 0 ; i <= p ; i ++ ) { */
  /* 	j = p - i ; */
  /* 	idx = afmm_derivative_index_ij(i,j) ; */
  /* 	s = 0 ; */
  /* 	fprintf(stderr, "%d %d %d %g (%g)\n", */
  /* 		n, i, j, T[s*(N+1)*ssdist+n*ssdist+2*idx+0], */
  /* 		ABS(T[s*(N+1)*ssdist+n*ssdist+2*idx+0]- */
  /* 		    Ss[s*(N+1)*ssdist+n*ssdist+2*idx+0])) ; */
  /*     } */
  /*   } */
  /* } */
  
  /* return ; */
  
  /*modal Green's function derivatives*/
  nd = afmm_derivative_offset(L+1) ;
  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z-z1, dG) ;  

  /*translate local source expansion to find local field*/
  rz[0] = -0.1 ; rz[1] = -0.15 ;
  afmm_laplace_source_to_local_f(N, L, dG, nd,
				     S, sdist, LS, ns,
				     P, pdist, LP, TRUE, TRUE) ;

  memset(g, 0, ns*dist*sizeof(gfloat)) ;
  afmm_expansion_eval_f(rz[0], rz[1], N, LP, P, pdist, ns, g, dist) ;

  /*same again with shifted moments*/
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1+dr, z-z1-dz, dG) ;
  /* afmm_laplace_source_to_local_f(N, L, dG, nd, */
  /* 				     Ss, ssdist, LSs, ns, */
  /* 				     Ps, psdist, LP, TRUE, TRUE) ; */
  /* memset(gs, 0, ns*dist*sizeof(gfloat)) ; */
  /* afmm_expansion_eval_f(rz[0], rz[1], N, LP, Ps, psdist, ns, gs, dist) ; */

  /*using shifted field coefficients*/
  afmm_expansion_shift_f(N, LP, P, pdist, ns, dr, dz, LP, Ps, psdist) ;
  
  memset(gs, 0, ns*dist*sizeof(gfloat)) ;
  afmm_expansion_eval_f(rz[0]+dr, rz[1]+dz, N, LP, Ps, psdist, ns, gs, dist) ;
  
  rz[0] += r ; rz[1] += z ;
  for ( j = 0 ; j < nrz1 ; j ++ ) {
    rz1[2*j+0] += r1 ; rz1[2*j+1] += z1 ;
  }

  memset(f, 0, ns*dist*sizeof(gfloat)) ;
  afmm_laplace_field_direct_f(rz1, C, dist, ns, nrz1, N,
				  rz, f, dist, 1, TRUE, work) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    fprintf(stderr, "%d", n) ;
    for ( i = 0 ; i < ns ; i ++ ) {
      fprintf(stderr, " %g (%g) %g (%g)",
	      f[i*dist+2*n+0],
	      ABS(gs[i*dist+2*n+0]-f[i*dist+2*n+0]),
	      f[i*dist+2*n+1],
	      ABS(gs[i*dist+2*n+1]-f[i*dist+2*n+1])) ;
    }
    fprintf(stderr, "\n") ;
  }
  
  return ;
}

static void backward_test(gfloat r, gfloat z,
			  gfloat r1, gfloat z1,
			  gfloat *rz1, gint nrz1,
			  gint N, gint L)

{
  gfloat *C, *P, *S, *dG, rz[32], work[1024], *f, *g ;
  gint dist, i, j, n, p, pdist, sdist, ns, nd, LS, LP ;

  ns = 2 ;

  LS = L/2 - 1 ; LP = L - LS ;
  rz[0] = 0.1 ; rz[1] = -0.11 ;

  fprintf(stderr, "backward translation test\n") ;
  fprintf(stderr, "=========================\n") ;
  fprintf(stderr, "(r , z ) = (%g, %g)\n", r, z) ;
  fprintf(stderr, "(dr, dz) = (%g, %g)\n", rz[0], rz[1]) ;
  fprintf(stderr, "(r1, z1) = (%g, %g)\n", r1, z1) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;
  fprintf(stderr, "LS = %d\n", LS) ;
  fprintf(stderr, "LP = %d\n", LP) ;  
  fprintf(stderr, "npts = %d\n", nrz1) ;
  fprintf(stderr, "ns = %d\n", ns) ;

  /*local source coefficients*/
  dist = 2*N + 4 ;
  C = (gfloat *)fftw_malloc(nrz1*ns*dist*sizeof(gfloat)) ;
  for ( i = 0 ; i < nrz1*ns ; i ++ ) {
    random_coefficients(&(C[i*dist]), N) ;
  }
  f = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;
  g = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;

  /*
   * local source moments (number of points, number of source terms,
   * number of modes, number of derivatives, times two for complex)
   */

  /*number of points per source point per source component*/
  sdist = 2*afmm_derivative_offset_2(LS+1) ;
  S = (gfloat *)g_malloc0(ns*(N+1)*sdist*sizeof(gfloat)) ;
  pdist = 2*afmm_derivative_offset_2(LP+1) ;
  P = (gfloat *)g_malloc0(ns*(N+1)*pdist*sizeof(gfloat)) ;
  
  /*
   * indexing within moment array:
   * sdist = 2*afmm_derivative_offset_2(L+1) entries, 
   * indexed by 2*afmm_derivative_index_ij(j,k) + (0,1)
   *
   * (N+1) blocks of size sdist, n=0, ..., N
   *
   * ns blocks of size (N+1)*sdist, one for each source component
   *
   * (j,k) moment of mode n of source component i is indexed by
   *
   * i*(N+1)*sdist + n*sdist + 2*afmm_derivative_index_ij(j,k)+(0,1)
   */
  
  /*generate local moments*/
  for ( p = 0 ; p < nrz1 ; p ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      afmm_source_moments_f(rz1[2*p+0], rz1[2*p+1],
				&(C[p*ns*dist + i*dist]), N, LS,
				&(S[i*(N+1)*sdist]), sdist) ;
    }
  }

  fprintf(stderr, "  local source moments generated about (r_1, z_1)\n") ;
  
  /*modal Green's function derivatives*/
  nd = afmm_derivative_offset(L+1) ;
  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z-z1, dG) ;  

  /*translate local source expansion to find local field*/
  afmm_laplace_source_to_local_f(N, L, dG, nd,
				     S, sdist, LS, ns,
				     P, pdist, LP, TRUE, TRUE) ;
  fprintf(stderr, "  local field expansion generated about (r,z)\n") ;

  memset(g, 0, ns*dist*sizeof(gfloat)) ;
  afmm_expansion_eval_f(rz[0], rz[1], N, LP, P, pdist, ns, g, dist) ;

  fprintf(stderr, "  local field expansion evaluated at (dr,dz)\n") ;
  
  rz[0] += r ; rz[1] += z ;
  for ( j = 0 ; j < nrz1 ; j ++ ) {
    rz1[2*j+0] += r1 ; rz1[2*j+1] += z1 ;
  }

  memset(f, 0, ns*dist*sizeof(gfloat)) ;
  afmm_laplace_field_direct_f(rz1, C, dist, ns, nrz1, N,
				  rz, f, dist, 1, TRUE, work) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    fprintf(stderr, "%d", n) ;
    for ( i = 0 ; i < ns ; i ++ ) {
      fprintf(stderr, " %g (%g) %g (%g)",
	      f[i*dist+2*n+0],
	      ABS(g[i*dist+2*n+0]-f[i*dist+2*n+0]),
	      f[i*dist+2*n+1],
	      ABS(g[i*dist+2*n+1]-f[i*dist+2*n+1])) ;
    }
    fprintf(stderr, "\n") ;
  }
  
  return ;
}

static void morton_test(gint depth)

{
  guint32 x, y, xd, yd ;
  guint64 z ;
  
  fprintf(stderr, "Morton index test\n") ;
  fprintf(stderr, "=================\n") ;
  fprintf(stderr, "depth = %d\n", depth) ;

  for ( x = 0 ; x < (1 << depth) ; x ++ ) {
    for ( y = 0 ; y < (1 << depth) ; y ++ ) {
      z = afmm_index_encode(x, y) ;
      afmm_index_decode(z, &xd, &yd) ;
      if ( (xd != x) || (yd != y) ) {
	fprintf(stderr,
		"error encoding (%u,%u), code %lu returns (%u,%u)\n",
		x, y, z, xd, yd) ;
	exit(1) ;
      }
    }
  }

  fprintf(stderr, "x = (0, %u)\n", (1 << depth) - 1) ;
  fprintf(stderr, "y = (0, %u)\n", (1 << depth) - 1) ;
  
  return ;
}

static void interaction_test(void)

{
  guint64 b, ilist[64] ;
  guint32 x, y ;
  guint depth ;
  gint ni, i ;

  depth = 4 ;
  x = 2 ; y = 5 ;
  
  fprintf(stderr, "box interaction list\n") ;
  fprintf(stderr, "====================\n") ;
  fprintf(stderr, "(x,y) = (%u, %u)\n", x, y) ;
  fprintf(stderr, "depth = %u\n", depth) ;
  
  b = afmm_index_encode(x, y) ;

  afmm_interaction_list(b, depth, ilist, &ni) ;

  fprintf(stderr, "%d interactions\n", ni) ;

  fprintf(stdout, "%u %u\n", x, y) ;
  for ( i = 0 ; i < ni ; i ++ ) {
    afmm_index_decode(ilist[i], &x, &y) ;
    fprintf(stdout, "%u %u\n", x, y) ;
  }
  
  return ;
}

static void s2l_test(gfloat r, gfloat z,
		     gfloat r1, gfloat z1,
		     gfloat *rz1, gint nrz1,
		     gint N, gint L, gint ns,
		     gboolean forward, gboolean outward)

{
  gfloat *C, *P, *S, *Pm, *dG, *S2L, rz[32], work[1024], *f, *g ;
  gfloat Z, Z1, R, R1 ;
  gint dist, i, j, n, p, pdist, sdist, nd, LS, LP ;
    
  LS = L/2 - 1 ; LP = L - LS ;
  rz[0] = 0.1 ; rz[1] = -0.11 ;

  fprintf(stderr, "source to local matrix test\n") ;
  fprintf(stderr, "===========================\n") ;
  fprintf(stderr, "(r , z ) = (%g, %g)\n", r, z) ;
  fprintf(stderr, "(dr, dz) = (%g, %g)\n", rz[0], rz[1]) ;
  fprintf(stderr, "(r1, z1) = (%g, %g)\n", r1, z1) ;
  fprintf(stderr, "N = %d\n", N) ;
  fprintf(stderr, "L = %d\n", L) ;
  fprintf(stderr, "LS = %d\n", LS) ;
  fprintf(stderr, "LP = %d\n", LP) ;  
  fprintf(stderr, "npts = %d\n", nrz1) ;
  fprintf(stderr, "ns = %d\n", ns) ;
  fprintf(stderr, "%s shift\n", ( forward ? "forward" : "backward")) ;
  fprintf(stderr, "%s shift\n", ( outward ? "outward" : "inward")) ;

  /*local source coefficients*/
  dist = 2*N + 4 ;
  C = (gfloat *)fftw_malloc(nrz1*ns*dist*sizeof(gfloat)) ;
  for ( i = 0 ; i < nrz1*ns ; i ++ ) {
    random_coefficients(&(C[i*dist]), N) ;
  }
  f = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;
  g = (gfloat *)fftw_malloc(ns*dist*sizeof(gfloat)) ;

  /*
   * local source moments (number of points, number of source terms,
   * number of modes, number of derivatives, times two for complex)
   */

  /*number of points per source point per source component*/
  sdist = 2*afmm_derivative_offset_2(LS+1) ;
  S = (gfloat *)g_malloc0(ns*(N+1)*sdist*sizeof(gfloat)) ;
  pdist = 2*afmm_derivative_offset_2(LP+1) ;
  P = (gfloat *)g_malloc0(ns*(N+1)*pdist*sizeof(gfloat)) ;
  Pm = (gfloat *)g_malloc0(ns*(N+1)*pdist*sizeof(gfloat)) ;
  
  /*
   * indexing within moment array:
   * sdist = 2*afmm_derivative_offset_2(L+1) entries, 
   * indexed by 2*afmm_derivative_index_ij(j,k) + (0,1)
   *
   * (N+1) blocks of size sdist, n=0, ..., N
   *
   * ns blocks of size (N+1)*sdist, one for each source component
   *
   * (j,k) moment of mode n of source component i is indexed by
   *
   * i*(N+1)*sdist + n*sdist + 2*afmm_derivative_index_ij(j,k)+(0,1)
   */
  
  /*generate local moments*/
  for ( p = 0 ; p < nrz1 ; p ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      afmm_source_moments_f(rz1[2*p+0], rz1[2*p+1],
				&(C[p*ns*dist + i*dist]), N, LS,
				&(S[i*(N+1)*sdist]), sdist) ;
    }
  }

  fprintf(stderr, "  local source moments generated about (r_1, z_1)\n") ;

  /*modal Green's function derivatives*/
  nd = afmm_derivative_offset(L+1) ;
  dG = (gfloat *)g_malloc0((N+1)*nd*sizeof(gfloat)) ;
  nd = afmm_laplace_gfunc_derivatives_f(N, L, r, r1, z-z1, dG) ;  

  S2L = (gfloat *)g_malloc0((N+1)*
			       afmm_derivative_offset(LS+1)*
			       afmm_derivative_offset(LP+1)*
			       sizeof(gfloat)) ;
  
  /*translate local source expansion to find local field*/
  afmm_laplace_source_to_local_f(N, L, dG, nd,
				     S, sdist, LS, ns,
				     P, pdist, LP,
				     forward, outward) ;
  fprintf(stderr, "  local field expansion generated about (r,z)\n") ;

  /*generate the S2L matrix*/
  afmm_laplace_s2l_matrix_f(N, L, dG, nd, LS, LP,
				forward, outward, S2L) ;
  /*size of S2L*/
  /* s2lr = afmm_derivative_offset_2(LP+1) ; */
  /* s2lc = afmm_derivative_offset_2(LS+1) ; */
  /* for ( n = 0 ; n <= N ; n ++ ) { */
  /*   s2l = &(S2L[n*s2lr*s2lc]) ; */
  /*   for ( s = 0 ; s < ns ; s ++ ) { */
  /*     blaswrap_dgemm(FALSE, FALSE, s2lr, i2, s2lc, d1, s2l, s2lc, */
  /* 		     &(S[s*(N+1)*sdist + n*sdist]), i2, d0, */
  /* 		     &(Pm[s*(N+1)*pdist + n*pdist]), i2) ; */
  /*   } */
  /* } */

  afmm_laplace_shift_s2l_f(N, S2L, S, sdist, LS, ns, Pm, pdist, LP) ;
  
  memset(g, 0, ns*dist*sizeof(gfloat)) ;
  afmm_expansion_eval_f(rz[0], rz[1], N, LP, Pm, pdist, ns, g, dist) ;
  fprintf(stderr, "  local field expansion evaluated at (dr,dz)\n") ;

  /* n = 7 ; */
  /* for ( s = 0 ; s < ns ; s ++ ) { */
  /*   fprintf(stderr, "%lg %lg\n", Pm[(n*ns+s)*pdist+0], P[(n*ns+s)*pdist+0]) ; */
  /*   fprintf(stderr, "%lg %lg\n", Pm[(n*ns+s)*pdist+1], P[(n*ns+s)*pdist+1]) ; */
  /* } */
  
  /* return ; */
  
  if ( forward ) {
    Z = z ; Z1 = z1 ;
  } else {
    Z = z1 ; Z1 = z ;
  }

  if ( outward ) {
    R = r ; R1 = r1 ;
  } else {
    R = r1 ; R1 = r ;
  }

  rz[0] += R ; rz[1] += Z ;
  for ( j = 0 ; j < nrz1 ; j ++ ) {
    rz1[2*j+0] += R1 ; rz1[2*j+1] += Z1 ;
  }

  memset(f, 0, ns*dist*sizeof(gfloat)) ;
  afmm_laplace_field_direct_f(rz1, C, dist, ns, nrz1, N,
				  rz, f, dist, 1, TRUE, work) ;

  for ( n = 0 ; n <= N ; n ++ ) {
    fprintf(stderr, "%d", n) ;
    for ( i = 0 ; i < ns ; i ++ ) {
      fprintf(stderr, " %g (%g) %g (%g)",
	      f[i*dist+2*n+0],
	      ABS(g[i*dist+2*n+0]-f[i*dist+2*n+0]),
	      f[i*dist+2*n+1],
	      ABS(g[i*dist+2*n+1]-f[i*dist+2*n+1])) ;
    }
    fprintf(stderr, "\n") ;
  }
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gfloat r, r1, z, chi, dx, *rz1 ;
  gint M, N, test, n, i, j, k, d, ns, nrz1 ;
  gchar ch, *progname ;
  gboolean forward, outward ;
  
  r = 0.5 ; r1 = 0.9 ; z = 2.5 ; chi = 3.0 ;
  N = 16 ; M = 8 ; dx = 0.25 ;
  n = 0 ; ns = 1 ;
  i = j = k = 0 ; d = 1 ;
  rz1 = NULL ;
  forward = TRUE ; outward = TRUE ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  test = -1 ;
  while ( (ch = getopt(argc, argv,
		       "Bc:D:d:Ii:j:k:M:N:n:R:r:s:T:z:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'B': forward = FALSE ; break ;
    case 'c': chi = atof(optarg) ; break ;
    case 'D': dx  = atof(optarg) ; break ;
    case 'd': d = atoi(optarg) ; break ;
    case 'I': outward = FALSE ; break ;
    case 'i': i = atoi(optarg) ; break ;
    case 'j': j = atoi(optarg) ; break ;
    case 'k': k = atoi(optarg) ; break ;
    case 'M': M = atoi(optarg) ; break ;      
    case 'N': N = atoi(optarg) ; break ;      
    case 'n': n = atoi(optarg) ; break ;
    case 'R': r = atof(optarg) ; break ;
    case 'r': r1 = atof(optarg) ; break ;
    case 's': ns = atoi(optarg) ; break ;
    case 'T': test = parse_test(optarg) ; break ;
    case 'z': z = atof(optarg) ; break ;
    }
  }

  if ( test == -1 ) {
    fprintf(stderr, "%s: undefined test\n", progname) ;

    return 1 ;
  }

  if ( test == 0 ) {
    elliptic_test() ;

    return 0 ;
  }    

  if ( test == 1 ) {
    legendre_test(chi, N, M, dx) ;

    return 0 ;
  }

  if ( test == 2 ) {
    greens_func_test(r, r1, z, N) ;

    return 0 ;
  }
  
  if ( test == 3 ) {
    series_test(r, r1, z, N, M) ;

    return 0 ;
  }

  if ( test == 4 ) {
    dump_test(stdout, r, r1, z, N, M, n) ;

    return 0 ;
  }

  if ( test == 5 ) {
    derivative_test(r, r1, z, N, M, i, j, k, d) ;

    return 0 ;
  }

  if ( test == 6 ) {
    fft_test(ns, N) ;

    return 0 ;
  }

  if ( test == 7 ) {
    rz1 = read_points(stdin, &nrz1) ;

    expansion_test(r, z, r1, z-dx, rz1, nrz1, N, M, forward, outward) ;

    return 0 ;
  }

  if ( test == 8 ) {
    rz1 = read_points(stdin, &nrz1) ;

    translation_test(r, z, r1, z-dx, rz1, nrz1, N, M, ns, forward, outward) ;

    return 0 ;
  }

  if ( test == 9 ) {
    rz1 = read_points(stdin, &nrz1) ;

    tree_test(rz1, nrz1) ;

    return 0 ;
  }

  if ( test == 10 ) {
    rz1 = read_points(stdin, &nrz1) ;

    shift_test(r, z, r1, z-dx, rz1, nrz1, N, M, ns) ;
  }

  if ( test == 11 ) {
    rz1 = read_points(stdin, &nrz1) ;

    backward_test(r, z, r1, z-dx, rz1, nrz1, N, M) ;

    return 0 ;
  }

  if ( test == 12 ) {
    /*should be inward test*/
    rz1 = read_points(stdin, &nrz1) ;

    backward_test(r, z, r1, z-dx, rz1, nrz1, N, M) ;

    return 0 ;
  }
  
  if ( test == 13 ) {
    morton_test(d) ;
    
    return 0 ;
  }

  if ( test == 14 ) {
    interaction_test() ;
    
    return 0 ;
  }

  if ( test == 15 ) {
    rz1 = read_points(stdin, &nrz1) ;

    s2l_test(r, z, r1, z-dx, rz1, nrz1, N, M, ns, forward, outward) ;

    return 0 ;
  }
  
  return 0 ;
}

