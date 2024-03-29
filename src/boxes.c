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

#ifndef AFMM_SINGLE_PRECISION
const guint afmm_ilist_di[] = {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3} ;
const guint afmm_ilist_dj[] = {2, 3, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3} ;
#endif /*AFMM_SINGLE_PRECISION*/

afmm_tree_t *AFMM_FUNCTION_NAME(afmm_tree_new)(AFMM_REAL rmin, AFMM_REAL rmax,
					       AFMM_REAL zmin, AFMM_REAL zmax,
					       guint maxpoints, guint maxfield)

{
  afmm_tree_t *t ;  
  gint i ;

  if ( rmin < 0.0 ) {
    g_error("%s: rmin (%lg) may not be negative", __FUNCTION__, rmin) ;
  }
  if ( rmax <= rmin ) {
    g_error("%s: rmax (%lg) must be greater than rmin (%lg)",
	    __FUNCTION__, rmax, rmin) ;
  }
  if ( zmax <= zmin ) {
    g_error("%s: zmax (%lg) must be greater than zmin (%lg)",
	    __FUNCTION__, zmax, zmin) ;
  }
  
  t = (afmm_tree_t *)g_malloc0(sizeof(afmm_tree_t)) ;

  afmm_tree_point_number_max(t) = maxpoints ;
  afmm_tree_field_number_max(t) = maxfield ;
  afmm_tree_point_number(t) = 0 ;
  afmm_tree_field_number(t) = 0 ;
  afmm_tree_mode_number(t)  = 0 ;
  afmm_tree_source_size(t)  = 0 ;
  t->ip = (guint32 *)g_malloc0(maxpoints*sizeof(guint32)) ;
  t->points = NULL ;
  if ( maxfield != 0 ) 
    t->ifld = (guint32 *)g_malloc0(maxfield*sizeof(guint32)) ;
  t->field  = NULL ;

  afmm_tree_box_separation(t) = 3.0/2.0/sqrt(2.0) - 1.0 ;
  
  for ( i = 0 ; i <= AFMM_TREE_MAX_DEPTH ; i ++ ) t->boxes[i] = NULL ;
  t->boxes[0] = (afmm_box_t *)g_malloc(1*sizeof(afmm_box_t)) ;
  t->boxes[0][0].i = t->boxes[0][0].ip = 0 ;
  t->boxes[0][0].n = t->boxes[0][0].np = 0 ;
  
  afmm_tree_r_min(t) = rmin ; 
  afmm_tree_r_max(t) = rmax ; 
  afmm_tree_z_min(t) = zmin ; 
  afmm_tree_z_max(t) = zmax ; 

  t->size = sizeof(AFMM_REAL) ;

  afmm_tree_source_data(t) = NULL ;

  afmm_tree_sources_sorted(t) = FALSE ;

  t->G = NULL ;
  t->L = 0 ;
  
  return t ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_init_s2l)(afmm_tree_t *t, guint d, guint L)

{
  gint nd, i, j, k, b ;
  guint nb, N ;
  AFMM_REAL *G, *Gb ;
  
  g_assert(t->size == sizeof(AFMM_REAL)) ;

  if ( d == 0 ) d = t->depth ;
  /*maximum number of derivatives*/
  if ( L == 0 ) {
    L = 0 ;
    for ( i = 0 ; i <= d ; i ++ ) {
      L = MAX(L, t->order_s[i] + t->order_f[i]) ;
    }
  }

  N = t->N ;
  nd = afmm_derivative_offset(L+2)*(N+2) ;
  t->L = L ; t->ngd = nd ;

  /*number of boxes per side*/
  nb = 1 << d ;

  /*allocate memory for Green's functions*/
  t->G = g_malloc0(12*nb*nd*sizeof(AFMM_REAL)) ;
  G = (AFMM_REAL *)(t->G) ;
  
  for ( b = 0 ; b < nb ; b ++ ) {
    Gb = &(G[12*b*nd]) ;
    j = 0 ; k = 2 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 0 ; k = 3 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 1 ; k = 2 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 1 ; k = 3 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 2 ; k = 0 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 2 ; k = 1 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 2 ; k = 2 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 2 ; k = 3 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 3 ; k = 0 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 3 ; k = 1 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 3 ; k = 2 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
    j = 3 ; k = 3 ; i = afmm_tree_s2l_index(j, k) ; g_assert(i != -1) ;
    /* fprintf(stderr, "%d %d %d\n", j, k, i) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
						       b+j+0.5, b+0.5, k,
						       &(Gb[i*nd])) ;
  }
  
  return 0 ;
}

guint64 AFMM_FUNCTION_NAME(afmm_point_index_2d)(AFMM_REAL *rz,
						AFMM_REAL rmin, AFMM_REAL rmax,
						AFMM_REAL zmin, AFMM_REAL zmax)

{
  guint32 ri, zi ;
  guint64 i ;

  ri = (guint32)((rz[0] - rmin)/(rmax-rmin)*AFMM_INDEX_SHIFT) ;
  zi = (guint32)((rz[1] - zmin)/(zmax-zmin)*AFMM_INDEX_SHIFT) ;

  i = afmm_index_encode(ri, zi) ;
  
  return i ;
}

static gint compare_field_indexed(gconstpointer a, gconstpointer b,
				  gpointer data)

{
  guint i, j ;
  afmm_tree_t *t = data ;
  guint64 mi, mj ;
  AFMM_REAL *rzi, *rzj ;

  i = *((guint *)a) ; j = *((guint *)b) ;

  rzi = afmm_tree_field_index(t, i) ; 
  rzj = afmm_tree_field_index(t, j) ;
  /*Morton codes*/
  mi = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rzi,
					       afmm_tree_r_min(t),
					       afmm_tree_r_max(t),
					       afmm_tree_z_min(t),
					       afmm_tree_z_max(t)) ;
  mj = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rzj,
					       afmm_tree_r_min(t),
					       afmm_tree_r_max(t),
					       afmm_tree_z_min(t),
					       afmm_tree_z_max(t)) ;

  if ( mi < mj ) return -1 ;
  if ( mi > mj ) return  1 ;

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_add_field)(afmm_tree_t *t,
					     gpointer pts, gsize pstr,
					     guint npts, gboolean sorted)

/*
 * pstr: stride in bytes; components of point are packed in pairs;
 */
  
{
  gint i ;
  AFMM_REAL *rz ;
  
  if ( t->size != sizeof(AFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of tree data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, t->size, sizeof(AFMM_REAL)) ;

  if ( npts > afmm_tree_field_number_max(t) ) 
    g_error("%s: too many field points (%u) for tree (%u)",
	    __FUNCTION__, npts, afmm_tree_field_number_max(t)) ;

  afmm_tree_field_number(t) = npts ;
  t->field = (gchar *)pts ; t->fstr = pstr ;

  for ( i = 0 ; i < npts ; i ++ ) {
    rz = afmm_tree_field_index(t, i) ;
    if ( rz[0] < afmm_tree_r_min(t) ||
	 rz[0] > afmm_tree_r_max(t) ||
	 rz[1] < afmm_tree_z_min(t) ||
	 rz[1] > afmm_tree_z_max(t) ) {
      g_error("%s: point %d (%g,%g) does not lie in bounding box of tree\n"
	      "(r = %g--%g, z = %g--%g)", __FUNCTION__, i,
	      rz[0], rz[1],
	      afmm_tree_r_min(t), afmm_tree_r_max(t),
	      afmm_tree_z_min(t), afmm_tree_z_max(t)) ;
    }
    t->ifld[i] = i ;
  }
  
  g_qsort_with_data(t->ifld, npts, sizeof(guint), compare_field_indexed, 
		    (gpointer)t) ;

  t->boxes[0][0].ip = 0 ; 
  t->boxes[0][0].np = npts ; 

  return 0 ;
}

static gint compare_morton_indexed(gconstpointer a, gconstpointer b,
				   gpointer data)

{
  guint i, j ;
  afmm_tree_t *t = data ;
  guint64 mi, mj ;
  AFMM_REAL *rzi, *rzj ;

  i = *((guint *)a) ; j = *((guint *)b) ;

  rzi = afmm_tree_point_index(t, i) ; 
  rzj = afmm_tree_point_index(t, j) ;
  /*Morton codes*/
  mi = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rzi,
					       afmm_tree_r_min(t),
					       afmm_tree_r_max(t),
					       afmm_tree_z_min(t),
					       afmm_tree_z_max(t)) ;
  mj = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rzj,
					       afmm_tree_r_min(t),
					       afmm_tree_r_max(t),
					       afmm_tree_z_min(t),
					       afmm_tree_z_max(t)) ;

  if ( mi < mj ) return -1 ;
  if ( mi > mj ) return  1 ;

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_add_sources)(afmm_tree_t *t,
					      gpointer pts, gsize pstr,
					      guint npts, gboolean sorted)

/*
 * pstr: stride in bytes; components of point are packed in pairs;
 */
  
{
  gint i ;
  AFMM_REAL *rz ;
  
  if ( t->size != sizeof(AFMM_REAL) )
    g_error("%s: mixed precision not implemented\n"
	    "  (size of tree data type (%lu) not equal to "
	    "size of requested target type (%lu))",
	    __FUNCTION__, t->size, sizeof(AFMM_REAL)) ;

  if ( npts > afmm_tree_point_number_max(t) ) 
    g_error("%s: too many points (%u) for tree (%u)",
	    __FUNCTION__, npts, afmm_tree_point_number_max(t)) ;

  afmm_tree_point_number(t) = npts ;
  t->points = (gchar *)pts ; t->pstr = pstr ;

  afmm_tree_sources_sorted(t) = sorted ;
  
  for ( i = 0 ; i < npts ; i ++ ) {
    rz = afmm_tree_point_index(t, i) ;
    if ( rz[0] < afmm_tree_r_min(t) ||
	 rz[0] > afmm_tree_r_max(t) ||
	 rz[1] < afmm_tree_z_min(t) ||
	 rz[1] > afmm_tree_z_max(t) ) {
      g_error("%s: point %d (%g,%g) does not lie in bounding box of tree\n"
	      "(r = %g--%g, z = %g--%g)", __FUNCTION__, i,
	      rz[0], rz[1],
	      afmm_tree_r_min(t), afmm_tree_r_max(t),
	      afmm_tree_z_min(t), afmm_tree_z_max(t)) ;
    }
    t->ip[i] = i ;
  }
  
  g_qsort_with_data(t->ip, npts, sizeof(guint), compare_morton_indexed, 
		    (gpointer)t) ;

  t->boxes[0][0].i = 0 ; 
  t->boxes[0][0].n = npts ; 

  if ( !sorted ) return 0 ;

  for ( i = 0 ; i < npts ; i ++ ) {
    if ( t->ip[i] != i )
      g_error("%s: points not sorted (called with sorted == TRUE)",
	      __FUNCTION__) ;
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_refine)(afmm_tree_t *t)

{
  guint level = afmm_tree_depth(t) ;
  afmm_box_t *parents, *children ;
  guint np, j ;
  guint64 idx, child, xi, box ;
  AFMM_REAL *x ;

  /* g_assert(t->problem != 0) ; */
  afmm_tree_add_level(t) ;

  /*number of parent boxes to refine*/
  np = 1 << 2*(level) ;

  parents  = t->boxes[level  ] ;
  children = t->boxes[level+1] ;

  /*this could probably be done with binary searches*/
  for ( idx = 0 ; idx < np ; idx ++ ) {
    /*initialize the first child box*/
    child = afmm_box_first_child(idx) ;
    children[child].i = parents[idx].i ;
    children[child].n = 0 ;
    /*start at first parent index*/
    j = parents[idx].i ;
    while ( parents[idx].n != 0 ) {
      /*check if current point is in box*/
      x = afmm_tree_point_index(t, t->ip[j]) ;
      xi = AFMM_FUNCTION_NAME(afmm_point_index_2d)(x,
						   afmm_tree_r_min(t),
						   afmm_tree_r_max(t),
						   afmm_tree_z_min(t),
						   afmm_tree_z_max(t)) ;
      box = afmm_point_locate_box(xi, level+1) ;
      /* g_assert(box >= child && box < child+8) ; */
      if ( box == child ) {
	parents[idx].n -- ;
	parents[idx].i ++ ;
	children[child].n ++ ;
	j ++ ;
      } else {
	child ++ ;
	children[child].n = 0 ; 
	children[child].i = parents[idx].i ;
      }
    }
    g_assert(parents[idx].n == 0) ;
    child = afmm_box_first_child(idx) ;
    parents[idx].i = children[child].i ;
    parents[idx].n =
      children[child+0].n + children[child+1].n +
      children[child+2].n + children[child+3].n ;    
  }

  if ( afmm_tree_field_number(t) == 0 ) return 0 ;

  for ( idx = 0 ; idx < np ; idx ++ ) {
    /*initialize the first child box*/
    child = afmm_box_first_child(idx) ;
    children[child].ip = parents[idx].ip ;
    children[child].np = 0 ;
    /*start at first parent index*/
    j = parents[idx].ip ;
    while ( parents[idx].np != 0 ) {
      /*check if current point is in box*/
      x = afmm_tree_field_index(t, t->ifld[j]) ;
      xi = AFMM_FUNCTION_NAME(afmm_point_index_2d)(x,
						   afmm_tree_r_min(t),
						   afmm_tree_r_max(t),
						   afmm_tree_z_min(t),
						   afmm_tree_z_max(t)) ;
      box = afmm_point_locate_box(xi, level+1) ;
      /* g_assert(box >= child && box < child+8) ; */
      if ( box == child ) {
	parents[idx].np -- ;
	parents[idx].ip ++ ;
	children[child].np ++ ;
	j ++ ;
      } else {
	child ++ ;
	children[child].np = 0 ; 
	children[child].ip = parents[idx].ip ;
      }
    }
    g_assert(parents[idx].np == 0) ;
    child = afmm_box_first_child(idx) ;
    parents[idx].ip = children[child].ip ;
    parents[idx].np =
      children[child+0].np + children[child+1].np +
      children[child+2].np + children[child+3].np ;    
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_coefficient_init)(afmm_tree_t *t,
						    guint l,
						    guint nf, guint ns)

{
  guint nb, nc, i ;
  gint nq, N, dist ;
  afmm_box_t *boxes ;
  AFMM_REAL *c ;

  nq = afmm_tree_source_size(t) ;
  N = afmm_tree_mode_number(t) ;

  if ( nq < 1 )
    g_error("%s: tree has invalid number of components in source terms (%d)",
	    __FUNCTION__, nq) ;
  
  g_assert(l <= afmm_tree_depth(t)) ;

  /*number of boxes at level l*/
  nb = 1 << (2*l) ;

  t->Cs[l] = t->Cf[l] = NULL ;
  t->order_s[l] = ns ; t->order_f[l] = nf ;

  /*coefficients in source expansions*/
  if ( ns != 0 ) {
    nc = afmm_coefficient_number(ns) ;
    dist = nc*2*(N+1) ;
    /* 
     * (number of boxes on level)x(number of coefficients for expansion)x
     * (number of source components)x(elements for N+1 modes)
     */
    t->Cs[l] = fftw_malloc(nb*nq*dist*sizeof(AFMM_REAL)) ;
    memset(t->Cs[l], 0, nb*nq*dist*sizeof(AFMM_REAL)) ;
    c = (AFMM_REAL *)(t->Cs[l]) ;
    /*
     * set box pointers to start of their coefficients
     */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i ++ ) {
      boxes[i].Cs = &(c[i*nq*dist]) ;
    }
  }

  /*coefficients in field expansions*/
  if ( nf != 0 ) {
    nc = afmm_coefficient_number(nf) ;
    dist = nc*2*(N+1) ;
    /* 
     * (number of boxes on level)x(number of coefficients for expansion)x
     * (number of source components)x(elements for N+1 modes)
     */
    t->Cf[l] = fftw_malloc(nb*nq*dist*sizeof(AFMM_REAL)) ;
    memset(t->Cf[l], 0, nb*nq*dist*sizeof(AFMM_REAL)) ;
    /* t->Cf[l] = fftw_malloc(nb*nc*nq*2*(N+1)*sizeof(AFMM_REAL)) ; */
    /* memset(t->Cf[l], 0, nb*nc*nq*2*(N+1)*sizeof(AFMM_REAL)) ; */
    c = (AFMM_REAL *)(t->Cf[l]) ;
    /*
     * set box pointers to start of their coefficients
     */
    boxes = t->boxes[l] ;
    for ( i = 0 ; i < nb ; i ++ ) {
      boxes[i].Cf = &(c[i*nq*dist]) ;
    }
  }

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_leaf_expansions)(afmm_tree_t *t,
						   AFMM_REAL *Cs, gint cdist,
						   gint cstr,
						   gboolean zero_expansions)
{
  guint32 nb, nc, i, j, s, ns, d, idx, sdist ;
  guint64 im ;
  afmm_box_t *boxes ;
  AFMM_REAL *rz, r, rb, z, zb, *S, *C ;
  gint nq = afmm_tree_source_size(t) ;
  gint N = afmm_tree_mode_number(t) ;

  /*depth of leaves*/
  d = afmm_tree_depth(t) ;
  /*order of singular expansions*/
  ns = t->order_s[d] ;
  /*number of boxes*/
  nb = 1 << (2*d) ;
  /*number of coefficients*/
  nc = afmm_coefficient_number(ns) ;
  sdist = 2*nc ;

  afmm_tree_source_data(t)    = Cs ;
  afmm_tree_source_stride(t)  = cstr ;
  afmm_tree_source_dist(t)    = cdist ;
  
  /*zero the coefficients before accumulating*/
  if ( zero_expansions )
    memset(t->Cs[d], 0, nb*nq*(N+1)*sdist*sizeof(AFMM_REAL)) ;

  boxes = t->boxes[d] ;

  for ( i = 0 ; i < nb ; i ++ ) {
    im = (guint64)i ;
    AFMM_FUNCTION_NAME(afmm_box_location_from_index)(im, d,
						     t->rmin, t->rmax,
						     t->zmin, t->zmax,
						     &r, &rb, &z, &zb) ;
    r += 0.5*rb ; z += 0.5*zb ;
    S = boxes[i].Cs ;
    for ( j = 0 ; j < boxes[i].n ; j ++ ) {
      idx = t->ip[boxes[i].i+j] ;
      rz = afmm_tree_point_index(t, idx) ;
      C = &(Cs[idx*cdist]) ;
      for ( s = 0 ; s < nq ; s ++ ) {
	AFMM_FUNCTION_NAME(afmm_source_moments)(rz[0]-r, rz[1]-z,
						&(C[s*cstr]), N, ns,
						&(S[s*sdist*(N+1)]), sdist) ;
      }
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_box_location_from_index)(guint64 idx,
						      guint level,
						      AFMM_REAL rmin,
						      AFMM_REAL rmax,
						      AFMM_REAL zmin,
						      AFMM_REAL zmax,
						      AFMM_REAL *r,
						      AFMM_REAL *rb,
						      AFMM_REAL *z,
						      AFMM_REAL *zb)

{
  guint nb ;
  guint32 i, j ;

  nb = 1 << level ;
  if ( !(idx < (1 << 2*level)) ) {
    g_error("%s: box %lu is not at level %u", __FUNCTION__, idx, level) ;
  }

  *rb = (rmax-rmin)/nb ;
  *zb = (zmax-zmin)/nb ;

  afmm_index_decode(idx, &i, &j) ;

  *r = rmin + i*(*rb) ;
  *z = zmin + j*(*zb) ;

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_expansion_eval_field)(afmm_tree_t *t,
							guint level,
							AFMM_REAL r,
							AFMM_REAL z,
							AFMM_REAL *f,
							AFMM_REAL *work)

{
  guint nb ;
  guint64 i ;
  gint nd, N, L, cdist, fdist, ns ;
  afmm_box_t *boxes ;
  AFMM_REAL *C, *dG, r1, z1, rb, zb ;
  
  g_assert(level <= afmm_tree_depth(t)) ;

  nb = 1 << (2*level) ;
  boxes = t->boxes[level] ;

  dG = work ;
  N  = afmm_tree_mode_number(t) ;
  ns = afmm_tree_source_size(t) ;
  L = t->order_s[level] ;
  cdist = 2*afmm_coefficient_number(L) ;
  fdist = 2*(N+2) ;
  memset(f, 0, fdist*ns*sizeof(AFMM_REAL)) ;
  
  for ( i = 0 ; i < nb ; i ++ ) {
    if ( boxes[i].n != 0 ) {
      AFMM_FUNCTION_NAME(afmm_box_location_from_index)(i, level,
						       afmm_tree_r_min(t),
						       afmm_tree_r_max(t),
						       afmm_tree_z_min(t),
						       afmm_tree_z_max(t),
						       &r1, &rb, &z1, &zb) ;
      r1 += 0.5*rb ; z1 += 0.5*zb ;
      nd = AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N+1, L,
							      r, r1, z-z1, dG) ;
      C = boxes[i].Cs ;
      AFMM_FUNCTION_NAME(afmm_laplace_field_eval)(N, L, dG, nd,
						  C, cdist,
						  afmm_tree_source_size(t),
						  f, fdist,
						  TRUE, TRUE, FALSE) ;
    }
  }

  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_local_field_eval)(afmm_tree_t *t,
						    AFMM_REAL r,
						    AFMM_REAL z,
						    AFMM_REAL *f, gint fdist)

{
  guint64 i, idx ;
  guint32 nb ;
  guint i0, j0, j, k ;
  AFMM_REAL rz[2] = {r, z}, *P, rc, zc, rb, zb ;
  AFMM_REAL work[256] ;
  gint Lp, pdist, ns, N ;
  afmm_box_t *boxes ;

  boxes = t->boxes[afmm_tree_depth(t)] ;
  Lp = t->order_f[afmm_tree_depth(t)] ;
  pdist = 2*afmm_derivative_offset_2(Lp+1) ;
  ns = afmm_tree_source_size(t) ;
  N  = afmm_tree_mode_number(t) ;
  
  i = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rz,
					      afmm_tree_r_min(t),
					      afmm_tree_r_max(t),
					      afmm_tree_z_min(t),
					      afmm_tree_z_max(t)) ;
  i = afmm_point_locate_box(i, afmm_tree_depth(t)) ;
  AFMM_FUNCTION_NAME(afmm_box_location_from_index)(i, afmm_tree_depth(t),
						   afmm_tree_r_min(t),
						   afmm_tree_r_max(t),
						   afmm_tree_z_min(t),
						   afmm_tree_z_max(t),
						   &rc, &rb, &zc, &zb) ;
  rc += 0.5*rb ; zc += 0.5*zb ;
  P = boxes[i].Cf ;

  AFMM_FUNCTION_NAME(afmm_expansion_eval)(r-rc, z-zc, N, Lp, P, pdist,
					  ns, f, fdist) ;

  nb = 1 << afmm_tree_depth(t) ;
  afmm_index_decode(i, &i0, &j0) ;

  for ( j = 0 ; j < nb ; j ++ ) {
    for ( k = 0 ; k < nb ; k ++ ) {
      if  ( !afmm_grid_boxes_separated(i0, j0, j, k,
				       afmm_tree_box_separation(t)) ) {
	idx = afmm_index_encode(j, k) ;
	AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[idx]),
						  rz, f, fdist, FALSE,
						  work) ;
      }
    }
  }
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_upward_pass)(afmm_tree_t *t, guint level,
					  AFMM_REAL *work)

{
  guint npb, N ;
  guint64 i, c ;
  AFMM_REAL *Sp, *Sc, r, z, dr, dz ;
  gint pdist, cdist, ns ;
  afmm_box_t *bp, *bc ;
  
  g_assert(level > 1) ;
  g_assert(level <= afmm_tree_depth(t)) ;

  N = afmm_tree_mode_number(t) ;
  ns = afmm_tree_source_size(t) ;
  /*number of parent boxes*/
  npb = 1 << (2*(level - 1)) ;
  Sp = t->Cs[level-1] ;
  /*zero parent box source data*/
  pdist = 2*afmm_coefficient_number(t->order_s[level-1]) ;
  memset(Sp, 0, npb*pdist*(N+1)*ns*sizeof(AFMM_REAL)) ;
  bp = t->boxes[level-1] ;
  
  Sc = t->Cs[level] ;
  cdist = 2*afmm_coefficient_number(t->order_s[level]) ;
  bc = t->boxes[level] ;

  /*get box size information*/
  AFMM_FUNCTION_NAME(afmm_box_location_from_index)(0, level-1,
						   t->rmin, t->rmax,
						   t->zmin, t->zmax,
						   &r, &dr, &z, &dz) ;
  /*moments are shifted by one quarter of the parent box size*/
  dr *= 0.25 ; dz *= 0.25 ;

  /*traverse parent level boxes*/
  for ( i = 0 ; i < npb ; i ++ ) {
    Sp = bp[i].Cs ;
    bp[i].n = 0 ;
    c = afmm_box_first_child(i) ;

    Sc = bc[c+0].Cs ;
    bp[i].n += bc[c+0].n ;
    AFMM_FUNCTION_NAME(afmm_moments_shift)(N,
					   t->order_s[level  ], Sc, cdist,
					   ns,  dr,  dz,
					   t->order_s[level-1], Sp, pdist) ;
    Sc = bc[c+1].Cs ;
    bp[i].n += bc[c+1].n ;
    AFMM_FUNCTION_NAME(afmm_moments_shift)(N,
					   t->order_s[level  ], Sc, cdist,
					   ns, -dr,  dz,
					   t->order_s[level-1], Sp, pdist) ;

    Sc = bc[c+2].Cs ;
    bp[i].n += bc[c+2].n ;
    AFMM_FUNCTION_NAME(afmm_moments_shift)(N,
					   t->order_s[level  ], Sc, cdist,
					   ns,  dr, -dz,
					   t->order_s[level-1], Sp, pdist) ;

    Sc = bc[c+3].Cs ;
    bp[i].n += bc[c+3].n ;
    AFMM_FUNCTION_NAME(afmm_moments_shift)(N,
					   t->order_s[level  ], Sc, cdist,
					   ns, -dr, -dz,
					   t->order_s[level-1], Sp, pdist) ;
    
  }
  
  return 0 ;
}

static void box_expansion_shift_s2l(afmm_box_t *boxes, gint N,
				    gint Ls, gint Lp, gint ns,
				    guint i, guint j, guint di, guint dj,
				    guint nb, AFMM_REAL etol,
				    AFMM_REAL *s2l1, AFMM_REAL *s2l2,
				    AFMM_REAL *s2l3, AFMM_REAL *s2l4)

{
  gint sdist, pdist ;
  guint64 is, ip ;
  AFMM_REAL *S, *P ;
  guint64 sp, pp ;
  
  if ( i + di >= nb ) return ;
  
  sdist = 2*afmm_derivative_offset_2(Ls+1) ;
  pdist = 2*afmm_derivative_offset_2(Lp+1) ;

  if ( j + dj < nb ) {
    is = afmm_index_encode(i, j) ; ip = afmm_index_encode(i+di, j+dj) ;
    sp = afmm_box_parent(is) ; 
    pp = afmm_box_parent(ip) ;
    if ( !afmm_boxes_separated(sp, pp, etol) ) {
      S = boxes[is].Cs ;  P = boxes[ip].Cf ;
      AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l1, S, sdist, Ls, ns,
					       P, pdist, Lp) ;
      S = boxes[ip].Cs ; P = boxes[is].Cf ;
      AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l4, S, sdist, Ls, ns,
					       P, pdist, Lp) ;
    }
  }

  if ( di == 0 || dj == 0 ) return ;
  
  if ( j >= dj ) {
    is = afmm_index_encode(i, j) ; ip = afmm_index_encode(i+di, j-dj) ;
    sp = afmm_box_parent(is) ; 
    pp = afmm_box_parent(ip) ;
    if ( !afmm_boxes_separated(sp, pp, etol) ) {
      S = boxes[is].Cs ;  P = boxes[ip].Cf ;
      AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l3, S, sdist, Ls, ns,
						 P, pdist, Lp) ;
      S = boxes[ip].Cs ; P = boxes[is].Cf ;
      AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l2, S, sdist, Ls, ns,
						 P, pdist, Lp) ;
    }
  }
  
  return ;
}

gint AFMM_FUNCTION_NAME(afmm_downward_pass)(afmm_tree_t *t, guint level,
					    AFMM_REAL *work,
					    gboolean downward)

{
  guint64 i, j, c, nb ;
  afmm_box_t *boxes, *cboxes ;
  AFMM_REAL r, r1, dz, *dG, *s2lfo, *s2lfi, *s2lbo, *s2lbi ;
  AFMM_REAL dr, z, *Pp, *Pc ;
  gint nd, Ls, Lp, N, k, ns, ni ;
  gint s2lr, s2lc, pdist, cdist, order_p, order_c ;
  guint ilist[3*512] ;
#ifdef AFMM_PRECOMPUTE_S2L
  AFMM_REAL *Gs2l ;
#endif /*AFMM_PRECOMPUTE_S2L*/
  
  g_assert(level < afmm_tree_depth(t)) ;

  /*no downward pass from level 0*/
  if ( level < 1 ) return 0 ;

  N  = afmm_tree_mode_number(t) ;
  ns = afmm_tree_source_size(t) ;
  /*loop on boxes at level and transfer local expansions downwards*/
  if ( downward ) {
    /*
     * this is the usual behaviour, but can be switched off for
     * testing local interactions in isolation
     */
    nb = 1 << 2*level ;
    boxes  = t->boxes[level  ] ;
    cboxes = t->boxes[level+1] ;
    order_p = afmm_tree_field_order(t,level  ) ;    
    order_c = afmm_tree_field_order(t,level+1) ;    
    pdist = 2*afmm_derivative_offset_2(order_p+1) ;
    cdist = 2*afmm_derivative_offset_2(order_c+1) ;

    nb = 1 << level ;
    dr = afmm_tree_delta_r(t)/nb ;
    dz = afmm_tree_delta_z(t)/nb ;
    dr *= -0.25 ; dz *= -0.25 ;
    nb *= nb ;
    for ( i = 0 ; i < nb ; i ++ ) {
      Pp =  boxes[i].Cf ;
      c = afmm_box_first_child(i) ;
      Pc = cboxes[c+0].Cf ;
      AFMM_FUNCTION_NAME(afmm_expansion_shift)(N,
					       order_p, Pp, pdist,
					       ns, -dr, -dz,
					       order_c, Pc, cdist) ;
      Pc = cboxes[c+1].Cf ;
      AFMM_FUNCTION_NAME(afmm_expansion_shift)(N,
					       order_p, Pp, pdist,
					       ns,  dr, -dz,
					       order_c, Pc, cdist) ;
      Pc = cboxes[c+2].Cf ;
      AFMM_FUNCTION_NAME(afmm_expansion_shift)(N,
					       order_p, Pp, pdist,
					       ns, -dr,  dz,
					       order_c, Pc, cdist) ;
      Pc = cboxes[c+3].Cf ;
      AFMM_FUNCTION_NAME(afmm_expansion_shift)(N,
					       order_p, Pp, pdist,
					       ns,  dr,  dz,
					       order_c, Pc, cdist) ;
    }
  }
  
  /*traverse boxes at level+1 and compute source-to-expansion shifts*/
  /*number of boxes per side at level+1*/
  nb = 1 << (level+1) ;

  dr = afmm_tree_delta_r(t)/nb ;
  dz = afmm_tree_delta_z(t)/nb ;
  g_assert(dr == dz) ;
  
  boxes = t->boxes[level+1] ;
  Ls = t->order_s[level+1] ;
  Lp = t->order_f[level+1] ;
  /*size of transfer matrices*/
  s2lr = afmm_derivative_offset_2(Lp+1) ;
  s2lc = afmm_derivative_offset_2(Ls+1) ;

  dG = work ;
  s2lfo = &(dG[(N+1)*afmm_derivative_offset(Ls+Lp+2)]) ;
  s2lfi = &(s2lfo[(N+1)*s2lr*s2lc]) ;
  s2lbo = &(s2lfi[(N+1)*s2lr*s2lc]) ;
  s2lbi = &(s2lbo[(N+1)*s2lr*s2lc]) ;

  ilist[3*0+0] = 0 ; ilist[3*0+1] = 2 ; 
  ilist[3*1+0] = 0 ; ilist[3*1+1] = 3 ; 
  ilist[3*2+0] = 1 ; ilist[3*2+1] = 2 ; 
  ilist[3*3+0] = 1 ; ilist[3*3+1] = 3 ;
  ilist[3*4+0] = 2 ; ilist[3*4+1] = 0 ; 
  ilist[3*5+0] = 2 ; ilist[3*5+1] = 1 ; 
  ilist[3*6+0] = 2 ; ilist[3*6+1] = 2 ; 
  ilist[3*7+0] = 2 ; ilist[3*7+1] = 3 ; 
  ilist[3*8+0] = 3 ; ilist[3*8+1] = 0 ; 
  ilist[3*9+0] = 3 ; ilist[3*9+1] = 1 ; 
  ilist[3*10+0] = 3 ; ilist[3*10+1] = 2 ; 
  ilist[3*11+0] = 3 ; ilist[3*11+1] = 3 ; 
  ni = 12 ;
  nd = afmm_derivative_offset(t->L+1) ;

  for ( i = 0 ; i < nb ; i ++ ) {
#ifndef AFMM_PRECOMPUTE_S2L
    afmm_radius_interaction_list(level+1, i,
				 afmm_tree_box_separation(t),
				 ilist, &ni, 512, TRUE) ;
    g_assert(ni < 13) ;
    r1 = t->rmin + dr*(i+0.5) ;
#else /*AFMM_PRECOMPUTE_S2L*/
    Gs2l = (AFMM_REAL *)(t->G) ;
    Gs2l = &(Gs2l[i*12*(t->ngd)]) ;
#endif /*AFMM_PRECOMPUTE_S2L*/   

#ifndef AFMM_PRECOMPUTE_S2L
    for ( k = 0 ; k < ni ; k ++ ) {
      r = t->rmin + dr*(ilist[3*k+0]+0.5) ;
      z = ilist[3*k+1]*dz ;
      nd = AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N+1, Lp+Ls,
							      r, r1, z,
							      dG) ;
      memset(s2lfo, 0, 4*(N+1)*s2lr*s2lc*sizeof(AFMM_REAL)) ;
      AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrices)(N, Lp+Ls, dG, nd,
						    Ls, Lp, 1.0,
						    s2lfo, s2lfi,
						    s2lbo, s2lbi) ;
#else /*AFMM_PRECOMPUTE_S2L*/      
    for ( k = 0 ; (k < ni) && (ilist[3*k+0] + i < nb) ; k ++ ) {
      gint idx ;
      idx = afmm_tree_s2l_index((ilist[3*k+0]), ilist[3*k+1]) ;
      g_assert(idx >= 0) ; g_assert(idx < 12) ;
      memset(s2lfo, 0, 4*(N+1)*s2lr*s2lc*sizeof(AFMM_REAL)) ;
      AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrices)(N, Lp+Ls,
						    &(Gs2l[idx*(t->ngd)]),
						    nd,
						    Ls, Lp, dr,
						    s2lfo, s2lfi,
						    s2lbo, s2lbi) ;
#endif /*AFMM_PRECOMPUTE_S2L*/   
      for ( j = 0 ; j < nb ; j ++ ) {
	box_expansion_shift_s2l(boxes, N, Ls, Lp, ns, i, j,
#ifndef AFMM_PRECOMPUTE_S2L
				ilist[3*k+0]-i, ilist[3*k+1]-0, nb,
#else /*AFMM_PRECOMPUTE_S2L*/   
				ilist[3*k+0], ilist[3*k+1], nb,
#endif /*AFMM_PRECOMPUTE_S2L*/   
				afmm_tree_box_separation(t),
				s2lfo, s2lfi, s2lbo, s2lbi) ;
      }
    }
  }
  
  return 0 ;
}
					    
gint AFMM_FUNCTION_NAME(afmm_box_field_direct)(afmm_tree_t *t,
					       afmm_box_t *b,
					       AFMM_REAL *rz,
					       AFMM_REAL *f, gint fdist,
					       gboolean zero,
					       AFMM_REAL *work)

{
  gint k, idx, cstr, cdist, ns, N, sstr ;
  AFMM_REAL *rzsrc, *src ;

  cstr = afmm_tree_source_stride(t) ;
  cdist = afmm_tree_source_dist(t) ;
  src = afmm_tree_source_data(t) ;
  ns = afmm_tree_source_size(t) ;
  N  = afmm_tree_mode_number(t) ;

  sstr = (t->pstr)/sizeof(AFMM_REAL) ;
  if ( afmm_tree_sources_sorted(t) ) {
    idx = t->ip[b->i] ;
    rzsrc = afmm_tree_point_index(t, idx) ;
    AFMM_FUNCTION_NAME(afmm_laplace_field_direct_vec)(rzsrc, sstr,
						      &(src[idx*cdist]), cdist,
						      cstr, ns, b->n, N,
						      rz, f, fdist, 1, FALSE,
						      8, work) ;
    return 0 ;
  }
  
  for ( k = 0 ; k < b->n ; k ++ ) {
    idx = t->ip[b->i+k] ;
    /* fprintf(stderr, "   %d\n", idx) ; */
    rzsrc = afmm_tree_point_index(t, idx) ;
    AFMM_FUNCTION_NAME(afmm_laplace_field_direct)(rzsrc,
						  &(src[idx*cdist]),
						  cstr, ns, 1, N,
						  rz, f, fdist, 1, FALSE,
						  work) ;
  }

  return 0 ;
}

static void box_field(afmm_tree_t *t,
		      afmm_box_t *boxes, guint i, guint j,
		      AFMM_REAL *f, gint fdist, AFMM_REAL *work)

{
  gint idx, k, ns ;
  AFMM_REAL *rz ;
  
  ns = afmm_tree_source_size(t) ;
  for ( k = 0 ; k < boxes[i].np ; k ++ ) {
    idx = t->ifld[boxes[i].ip+k] ;
    /* fprintf(stderr, "%d\n", idx) ; */
    rz = afmm_tree_field_index(t, idx) ;
    AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[j]),
					      rz, &(f[idx*ns*fdist]), fdist,
					      FALSE, work) ;
  }

  return ;
}

gint AFMM_FUNCTION_NAME(afmm_tree_field_eval)(afmm_tree_t *t,
					      AFMM_REAL *f,
					      gint fdist)

{
  guint nb, id, i, j, k, l, i0, j0, level ;
  AFMM_REAL work[16384], rc, zc, *P, *rz ;
  gint idx, Lp, pdist, ns, N ;
  afmm_box_t *boxes ;

  level = afmm_tree_depth(t) ;
  nb = 1 << level ;
  boxes = t->boxes[afmm_tree_depth(t)] ;
  Lp = t->order_f[afmm_tree_depth(t)] ;
  pdist = 2*afmm_derivative_offset_2(Lp+1) ;
  ns = afmm_tree_source_size(t) ;
  N  = afmm_tree_mode_number(t) ;

  for ( i0 = 0 ; i0 < nb ; i0 ++ ) {
    rc = afmm_tree_box_centre_r(t,i0,afmm_tree_depth(t)) ;
    for ( j0 = 0 ; j0 < nb ; j0 ++ ) {
      zc = afmm_tree_box_centre_z(t,j0,afmm_tree_depth(t)) ;
      i = afmm_index_encode(i0, j0) ;
      P = boxes[i].Cf ;
      for ( j = 0 ; j < boxes[i].np ; j ++ ) {
	idx = t->ifld[boxes[i].ip+j] ;
	rz = afmm_tree_field_index(t, idx) ;
	AFMM_FUNCTION_NAME(afmm_expansion_eval)(rz[0]-rc, rz[1]-zc, N, Lp, P,
						pdist,
						ns, &(f[idx*ns*fdist]),
						fdist) ;
      }

      for ( k = 0 ; k < i0 ; k ++ ) {
	for ( l = 0 ; l < nb ; l ++ ) {
	  if  ( !afmm_grid_boxes_separated(i0, j0, k, l,
					   afmm_tree_box_separation(t)) ) {
	    id = afmm_index_encode(k, l) ;
	    box_field(t, boxes, i , id, f, fdist, work) ;
	    box_field(t, boxes, id, i , f, fdist, work) ;
	  }
	}
      }

      k = i0 ; 
      for ( l = 0 ; l < j0 ; l ++ ) {
	if  ( !afmm_grid_boxes_separated(i0, j0, k, l,
					 afmm_tree_box_separation(t)) ) {
	  id = afmm_index_encode(k, l) ;
	  box_field(t, boxes, i , id, f, fdist, work) ;
	  box_field(t, boxes, id, i , f, fdist, work) ;
	}
      }
      box_field(t, boxes, i , i, f, fdist, work) ;      
    }
  }
  
  return 0 ;
}

static gint compare_morton_point(gconstpointer a, gconstpointer b,
				 gpointer data)

{
  AFMM_REAL *rza = (AFMM_REAL *)a, *rzb = (AFMM_REAL *)b, *lim = data ;
  guint64 ma, mb ;
  
  ma = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rza,
					       lim[0], lim[1], lim[2], lim[3]) ;
  mb = AFMM_FUNCTION_NAME(afmm_point_index_2d)(rzb,
					       lim[0], lim[1], lim[2], lim[3]) ;
  if ( ma < mb ) return -1 ;
  if ( ma > mb ) return  1 ;
  
  return 0 ;
}
  
gint AFMM_FUNCTION_NAME(afmm_sort_point_list)(AFMM_REAL *pts, gint psize,
					      gint npts,
					      AFMM_REAL rmin, AFMM_REAL rmax,
					      AFMM_REAL zmin, AFMM_REAL zmax)

{
  AFMM_REAL data[] = {rmin, rmax, zmin, zmax} ;

  g_qsort_with_data(pts, npts, psize, compare_morton_point, (gpointer)data) ;
  
  return 0 ;
}

gint AFMM_FUNCTION_NAME(afmm_box_field_indirect)(afmm_tree_t *t,
						 guint64 i,
						 guint level,
						 AFMM_REAL *rz,
						 AFMM_REAL *f, gint fdist,
						 gboolean zero,
						 AFMM_REAL *work)

{
  gint N, ns, L, nd, sdist ;
  AFMM_REAL *S, rb, dr, zb, dz, *dG ;
  afmm_box_t *boxes ;
  
  ns = afmm_tree_source_size(t) ;
  N  = afmm_tree_mode_number(t) ;
  L  = afmm_tree_source_order(t, level) ;
  sdist = 2*afmm_coefficient_number(L) ;
  dG = work ;
  
  boxes = t->boxes[level] ;
  
  S = boxes[i].Cs ;

  AFMM_FUNCTION_NAME(afmm_box_location_from_index)(i, level,
						   t->rmin, t->rmax,
						   t->zmin, t->zmax,
						   &rb, &dr, &zb, &dz) ;
  rb += 0.5*dr ; zb += 0.5*dz ;

  nd = AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N+1, L, rz[0], rb,
							  rz[1]-zb, dG) ;
  
  AFMM_FUNCTION_NAME(afmm_laplace_field_eval)(N, L, dG, nd, S, sdist,
					      ns,
					      f, fdist, TRUE, TRUE, FALSE) ;
  
  return 0 ;
}
