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

static gboolean boxes_separated(AFMM_REAL r, AFMM_REAL r1, AFMM_REAL z,
				AFMM_REAL dr, gdouble esep)

{
  AFMM_REAL chi ;

  g_assert(r >= r1) ;

  if ( r != r1 ) r1 += dr ;

  chi = 0.5*(r*r + r1*r1 + z*z)/r/r1 ;
  
  return ( chi > 1.0 + esep ) ;
}

afmm_tree_t *AFMM_FUNCTION_NAME(afmm_tree_new)(AFMM_REAL rmin, AFMM_REAL rmax,
					       AFMM_REAL zmin, AFMM_REAL zmax,
					       guint maxpoints)

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
  afmm_tree_point_number(t) = 0 ;
  afmm_tree_mode_number(t)  = 0 ;
  afmm_tree_source_size(t)  = 0 ;
  t->ip = (guint32 *)g_malloc0(maxpoints*sizeof(guint32)) ;
  t->points = NULL ;

  afmm_tree_box_separation(t) = 3.0/2.0/sqrt(2.0) - 1.0 ;
  
  for ( i = 0 ; i <= AFMM_TREE_MAX_DEPTH ; i ++ ) t->boxes[i] = NULL ;
  t->boxes[0] = (afmm_box_t *)g_malloc(1*sizeof(afmm_box_t)) ;

  afmm_tree_r_min(t) = rmin ; 
  afmm_tree_r_max(t) = rmax ; 
  afmm_tree_z_min(t) = zmin ; 
  afmm_tree_z_max(t) = zmax ; 

  t->size = sizeof(AFMM_REAL) ;

  afmm_tree_source_data(t) = NULL ;
  
  return t ;
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

gint AFMM_FUNCTION_NAME(afmm_tree_add_points)(afmm_tree_t *t,
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
  
  /* g_assert(t->problem == WBFMM_PROBLEM_HELMHOLTZ ) ; */
  /* /\* g_assert(t->ns == 1) ; *\/ */
  
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
						   gboolean zero_expansions)
{
  guint32 nb, nc, i, j, s, ns, d, idx, sdist ;
  guint64 im ;
  afmm_box_t *boxes ;
  AFMM_REAL *rz, r, rb, z, zb, *S, *C ;
  gint nq = afmm_tree_source_size(t) ;
  gint N = afmm_tree_mode_number(t) ;

  /* g_assert(nq == 1) ; */

  /*depth of leaves*/
  d = afmm_tree_depth(t) ;
  /*order of singular expansions*/
  ns = t->order_s[d] ;
  /*number of boxes*/
  nb = 1 << (2*d) ;
  /*number of coefficients*/
  nc = afmm_coefficient_number(ns) ;
  /* sdist = nc*2*(N+1) ; */
  sdist = 2*nc ;

  afmm_tree_source_data(t) = Cs ;
  
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
      C = &(Cs[idx*nq*cdist]) ;
      for ( s = 0 ; s < nq ; s ++ ) {
	AFMM_FUNCTION_NAME(afmm_source_moments)(rz[0]-r, rz[1]-z,
						&(C[s*cdist]), N, ns,
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
      /* fprintf(stderr, "%lu %g %g %g %g\n", i, r1, z1, rb, zb) ; */
      nd = AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, L,
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
  guint64 i, idx, nbrs[9] ;
  guint32 x, y, nb ;
  AFMM_REAL rz[2] = {r, z}, *P, rc, zc, rb, zb, *Cs, r1 ;
  AFMM_REAL work[256] ;
  gint Lp, pdist, ns, N, nn, j, sdist ;
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

  /*add contributions from neighbour boxes*/
  nn = afmm_box_neighbours(i, afmm_tree_depth(t), nbrs) ;
  sdist = 2*(N+2) ;
  Cs = afmm_tree_source_data(t) ;

  for ( j = 0 ; j < nn ; j ++ ) {
    /*field from sources in box nbrs[j] to field point*/
    /* if ( boxes[nbrs[j]].n != 0 ) { */
    /*   fprintf(stderr, "direct contribution to box %lu from box %lu\n", */
    /* 	      i, nbrs[j]) ; */
    /* } */
    AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[nbrs[j]]),
					      Cs, sdist, ns, N,
					      rz, f, fdist, FALSE, work) ;
  }

  nb = 1 << afmm_tree_depth(t) ;
  afmm_index_decode(i, &x, &y) ;
  /* rc = rb*x ;  */
  for ( j = 0 ; j < 12 ; j ++ ) {
    /* r1 = rb*(x+afmm_ilist_di[j]) ; */
    r1 = MIN(rb*(x+afmm_ilist_di[j]), rb*x) ;
    rc = MAX(rb*(x+afmm_ilist_di[j]), rb*x) ;
    zc = zb*   afmm_ilist_dj[j]  ;
    if ( !boxes_separated(rc, r1, zc, rb, afmm_tree_box_separation(t)) ) {
      if ( x + afmm_ilist_di[j] < nb && y + afmm_ilist_dj[j] < nb ) {
	idx = afmm_index_encode(x + afmm_ilist_di[j],
				y + afmm_ilist_dj[j]) ;
	AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[idx]),
						  Cs, sdist, ns, N,
						  rz, f, fdist, FALSE, work) ;
      }
      if ( x + afmm_ilist_di[j] < nb && y >= afmm_ilist_dj[j] &&
	   afmm_ilist_dj[j] != 0 ) {
	idx = afmm_index_encode(x + afmm_ilist_di[j],
				y - afmm_ilist_dj[j]) ;
	AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[idx]),
						  Cs, sdist, ns, N,
						  rz, f, fdist, FALSE, work) ;
      }
    }
      
    if ( afmm_ilist_di[j] != 0 && x >= afmm_ilist_di[j] ) {
      r1 = MIN(rb*(x-afmm_ilist_di[j]), rb*x) ;
      rc = MAX(rb*(x-afmm_ilist_di[j]), rb*x) ;
      /* r1 = rb*(x-afmm_ilist_di[j]) ; */
      if ( !boxes_separated(rc, r1, zc, rb, afmm_tree_box_separation(t)) ) {
	if ( x >= afmm_ilist_di[j] && y + afmm_ilist_dj[j] < nb ) {
	  idx = afmm_index_encode(x - afmm_ilist_di[j],
				  y + afmm_ilist_dj[j]) ;
	  AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[idx]),
						    Cs, sdist, ns, N,
						    rz, f, fdist, FALSE,
						    work) ;
	}
	if ( x >= afmm_ilist_di[j] && y >= afmm_ilist_dj[j] &&
	     afmm_ilist_dj[j] != 0 ) {
	  idx = afmm_index_encode(x - afmm_ilist_di[j],
				  y - afmm_ilist_dj[j]) ;
	  AFMM_FUNCTION_NAME(afmm_box_field_direct)(t, &(boxes[idx]),
						    Cs, sdist, ns, N,
						    rz, f, fdist, FALSE,
						    work) ;
	}
	
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

/* static void box_expansion_shift(afmm_box_t *boxes, gint N, */
/* 				gint Ls, gint Lp, gint ns, */
/* 				guint i, guint j, guint di, guint dj, */
/* 				guint nb, AFMM_REAL *dG, gint nd) */

/* { */
/*   gint sdist, pdist ; */
/*   guint64 is, ip ; */
/*   AFMM_REAL *S, *P ; */

/*   if ( i + di >= nb ) return ; */
  
/*   sdist = 2*afmm_derivative_offset_2(Ls+1) ; */
/*   pdist = 2*afmm_derivative_offset_2(Lp+1) ; */

/*   if ( j + dj < nb ) { */
/*     is = afmm_index_encode(i, j) ; ip = afmm_index_encode(i+di, j+dj) ; */
/*     S = boxes[is].Cs ;  P = boxes[ip].Cf ; */
/*     AFMM_FUNCTION_NAME(afmm_laplace_source_to_local)(N, Ls+Lp+1, dG, nd, */
/* 						     S, sdist, Ls, ns, */
/* 						     P, pdist, Lp, TRUE, TRUE) ; */
/*     is = afmm_index_encode(i+di, j+dj) ; ip = afmm_index_encode(i, j) ; */
/*     S = boxes[is].Cs ;  P = boxes[ip].Cf ; */
/*     AFMM_FUNCTION_NAME(afmm_laplace_source_to_local)(N, Ls+Lp+1, dG, nd, */
/* 						     S, sdist, Ls, ns, */
/* 						     P, pdist, Lp, */
/* 						     FALSE, FALSE) ; */
/*   } */

/*   /\* return ; *\/ */
/*   if ( di == 0 || dj == 0) return ; */
  
/*   if ( j >= dj ) { */
/*     is = afmm_index_encode(i, j) ; ip = afmm_index_encode(i+di, j-dj) ; */
/*     S = boxes[is].Cs ;  P = boxes[ip].Cf ; */
/*     AFMM_FUNCTION_NAME(afmm_laplace_source_to_local)(N, Ls+Lp+1, dG, nd, */
/* 						     S, sdist, Ls, ns, */
/* 						     P, pdist, Lp, */
/* 						     FALSE, TRUE) ; */
/*     is = afmm_index_encode(i+di, j-dj) ; ip = afmm_index_encode(i, j) ; */
/*     S = boxes[is].Cs ;  P = boxes[ip].Cf ; */
/*     AFMM_FUNCTION_NAME(afmm_laplace_source_to_local)(N, Ls+Lp+1, dG, nd, */
/* 						     S, sdist, Ls, ns, */
/* 						     P, pdist, Lp, */
/* 						     TRUE, FALSE) ; */
/*   } */
  
/*   return ; */
/* } */

static gboolean boxes_are_neighbours(guint64 i, guint64 j)

{
  guint32 ix, iy, jx, jy ;

  afmm_index_decode(i, &ix, &iy) ;
  afmm_index_decode(j, &jx, &jy) ;

  if ( ix > jx+1 ) return FALSE ;
  if ( jx > ix+1 ) return FALSE ;
  if ( iy > jy+1 ) return FALSE ;
  if ( jy > iy+1 ) return FALSE ;
  
  return TRUE ;
}

static void box_expansion_shift_s2l(afmm_box_t *boxes, gint N,
				    gint Ls, gint Lp, gint ns,
				    guint i, guint j, guint di, guint dj,
				    guint nb,
				    AFMM_REAL *s2l1, AFMM_REAL *s2l2,
				    AFMM_REAL *s2l3, AFMM_REAL *s2l4)

{
  gint sdist, pdist ;
  guint64 is, ip ;
  AFMM_REAL *S, *P ;
  /* guint64 icheck ; */
  guint64 sp, pp ;
  
  if ( i + di >= nb ) return ;

  /* if ( nb == 1  ) icheck =   0 ; */
  /* if ( nb == 2  ) icheck =   1 ; */
  /* if ( nb == 4  ) icheck =   7 ; */
  /* if ( nb == 8  ) icheck =  28 ; */
  /* if ( nb == 16 ) icheck = 113 ; */
  /* if ( nb == 32 ) icheck = 454 ; */
  /* if ( nb == 1  ) icheck =   0 ; */
  /* if ( nb == 2  ) icheck =   1 ; */
  /* if ( nb == 4  ) icheck =   5 ; */
  /* if ( nb == 8  ) icheck =  22 ; */
  /* if ( nb == 16 ) icheck =  91 ; */
  /* if ( nb == 32 ) icheck = 367 ; */
  /* if ( nb == 64 ) icheck =   8 ; */
  
  sdist = 2*afmm_derivative_offset_2(Ls+1) ;
  pdist = 2*afmm_derivative_offset_2(Lp+1) ;

  if ( j + dj < nb ) {
    is = afmm_index_encode(i, j) ; ip = afmm_index_encode(i+di, j+dj) ;
    sp = afmm_box_parent(is) ; 
    pp = afmm_box_parent(ip) ;
    if ( boxes_are_neighbours(sp, pp) ) {
    S = boxes[is].Cs ;  P = boxes[ip].Cf ;
    /* if ( ip == icheck && S[0] != 0 ) */
    /*   fprintf(stderr, "interaction: nb = %u; i=%u, j=%u, di=%u, dj=%u\n", */
    /* 	      nb, i, j, di, dj) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l1, S, sdist, Ls, ns,
					       P, pdist, Lp) ;
    S = boxes[ip].Cs ; P = boxes[is].Cf ;
    /* if ( is == icheck && S[0] != 0 ) */
    /*   fprintf(stderr, "interaction: nb = %u; i=%u, j=%u, di=%u, dj=%u\n", */
    /* 	      nb, i, j, di, dj) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l4, S, sdist, Ls, ns,
					       P, pdist, Lp) ;
    }
  }

  if ( di == 0 || dj == 0 ) return ;
  
  if ( j >= dj ) {
    is = afmm_index_encode(i, j) ; ip = afmm_index_encode(i+di, j-dj) ;
    sp = afmm_box_parent(is) ; 
    pp = afmm_box_parent(ip) ;
    if ( boxes_are_neighbours(sp, pp) ) {
    S = boxes[is].Cs ;  P = boxes[ip].Cf ;
    /* if ( ip == icheck && S[0] != 0 ) */
    /*   fprintf(stderr, "interaction: nb = %u; i=%u, j=%u, di=%u, dj=%u\n", */
    /* 	      nb, i, j, di, dj) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l3, S, sdist, Ls, ns,
					       P, pdist, Lp) ;
    S = boxes[ip].Cs ; P = boxes[is].Cf ;
    /* if ( is == icheck && S[0] != 0 ) */
    /*   fprintf(stderr, "interaction: nb = %u; i=%u, j=%u, di=%u, dj=%u\n", */
    /* 	      nb, i, j, di, dj) ; */
    AFMM_FUNCTION_NAME(afmm_laplace_shift_s2l)(N, s2l2, S, sdist, Ls, ns,
					       P, pdist, Lp) ;
    }
  }
  
  return ;
}

static void box_separation(afmm_tree_t *t, guint i, guint di, guint dj,
			   guint nb,
			   AFMM_REAL *r, AFMM_REAL *r1, AFMM_REAL *dz)

{
  *r1 = t->rmin + (t->rmax - t->rmin)/nb*(i+   0.5) ;
  *r  = t->rmin + (t->rmax - t->rmin)/nb*(i+di+0.5) ;
  *dz = (t->zmax - t->zmin)/nb*dj ;
    
  return ;
}

gint AFMM_FUNCTION_NAME(afmm_downward_pass)(afmm_tree_t *t, guint level,
					    AFMM_REAL *work,
					    gboolean downward)

{
  guint64 i, j, c, nb ;
  afmm_box_t *boxes, *cboxes ;
  AFMM_REAL r, r1, dz, *dG, *s2l1, *s2l2, *s2l3, *s2l4 ;
  AFMM_REAL dr, z, *Pp, *Pc ;
  gint nd, Ls, Lp, N, k, ns ;
  gint s2lr, s2lc, pdist, cdist, order_p, order_c ;
  AFMM_REAL rt, zt, drt, dzt ;
  
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
    
    for ( i = 0 ; i < nb ; i ++ ) {
      AFMM_FUNCTION_NAME(afmm_box_location_from_index)(i, level,
						       afmm_tree_r_min(t),
						       afmm_tree_r_max(t),
						       afmm_tree_z_min(t),
						       afmm_tree_z_max(t),
						       &r, &dr, &z, &dz) ;
      dr *= -0.25 ; dz *= -0.25 ;
      Pp =  boxes[i].Cf ;
      c = afmm_box_first_child(i) ;
      AFMM_FUNCTION_NAME(afmm_box_location_from_index)(c, level+1,
						       afmm_tree_r_min(t),
						       afmm_tree_r_max(t),
						       afmm_tree_z_min(t),
						       afmm_tree_z_max(t),
						       &rt, &drt, &zt, &dzt) ;
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
  /* nboxes = 1 << 2*(level+1) ; */
  nb = 1 << (level+1) ;

  dr = (afmm_tree_r_max(t) - afmm_tree_r_min(t))/nb ;
  /* dz = (afmm_tree_z_max(t) - afmm_tree_z_min(t))/nb ; */
  
  boxes = t->boxes[level+1] ;
  Ls = t->order_s[level+1] ;
  Lp = t->order_f[level+1] ;
  /*size of transfer matrices*/
  s2lr = afmm_derivative_offset_2(Lp+1) ;
  s2lc = afmm_derivative_offset_2(Ls+1) ;

  dG = work ;
  s2l1 = &(dG[(N+1)*afmm_derivative_offset(Ls+Lp+2)]) ;
  s2l2 = &(s2l1[(N+1)*s2lr*s2lc]) ;
  s2l3 = &(s2l2[(N+1)*s2lr*s2lc]) ;
  s2l4 = &(s2l3[(N+1)*s2lr*s2lc]) ;

  for ( i = 0 ; i < nb ; i ++ ) {
    for ( k = 0 ; (k < 12) && (i+afmm_ilist_di[k] < nb) ; k ++ ) {
      box_separation(t, i, afmm_ilist_di[k], afmm_ilist_dj[k],
		     nb, &r, &r1, &dz) ;
      if ( boxes_separated(r-0.5*dr, r1-0.5*dr, dz, dr,
			   afmm_tree_box_separation(t)) ) {
	nd = AFMM_FUNCTION_NAME(afmm_laplace_gfunc_derivatives)(N, Lp+Ls+1,
								r, r1, dz,
								dG) ;
	memset(s2l1, 0, 4*(N+1)*s2lr*s2lc*sizeof(AFMM_REAL)) ;
	AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrix)(N, Lp+Ls+1, dG, nd, Ls, Lp,
						    TRUE, TRUE, s2l1) ;
	AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrix)(N, Lp+Ls+1, dG, nd, Ls, Lp,
						    TRUE, FALSE, s2l2) ;
	AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrix)(N, Lp+Ls+1, dG, nd, Ls, Lp,
						    FALSE, TRUE, s2l3) ;
	AFMM_FUNCTION_NAME(afmm_laplace_s2l_matrix)(N, Lp+Ls+1, dG, nd, Ls, Lp,
						    FALSE, FALSE, s2l4) ;
	for ( j = 0 ; j < nb ; j ++ ) {
	  /* box_expansion_shift(boxes, N, Ls, Lp, ns, i, j, di[k], dj[k], nb, */
	  /* 		    work, nd) ; */
	  box_expansion_shift_s2l(boxes, N, Ls, Lp, ns, i, j,
				  afmm_ilist_di[k], afmm_ilist_dj[k], nb,
				  s2l1, s2l2, s2l3, s2l4) ;
	}
      } else {
	fprintf(stderr, "level = %u; i = %lu; di = %u; dj = %u;\n",
		level+1, i, afmm_ilist_di[k], afmm_ilist_dj[k]) ;
	fprintf(stderr, "r=%lg; r1=%lg; dr=%lg; dz=%lg;\n", r, r1, dr, dz) ;
      }
    }
  }
  
  return 0 ;
}
					    
gint AFMM_FUNCTION_NAME(afmm_box_field_direct)(afmm_tree_t *t,
					       afmm_box_t *b,
					       AFMM_REAL *src, gint sdist,
					       gint ns, gint N,
					       AFMM_REAL *rz,
					       AFMM_REAL *f, gint fdist,
					       gboolean zero,
					       AFMM_REAL *work)

{
  gint k, idx ;
  AFMM_REAL *rzsrc ;
  
  for ( k = 0 ; k < b->n ; k ++ ) {
    idx = t->ip[b->i+k] ;
    rzsrc = afmm_tree_point_index(t, idx) ;
    AFMM_FUNCTION_NAME(afmm_laplace_field_direct)(rzsrc,
						  &(src[idx*ns*sdist]),
						  sdist, ns, 1, N,
						  rz, f, fdist, 1, FALSE,
						  work) ;
  }

  return 0 ;
}
					       
