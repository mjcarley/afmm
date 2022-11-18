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

gint afmm_tree_add_level(afmm_tree_t *t)

{
  guint nb, d ;
  
  if ( (d = t->depth) >= AFMM_TREE_MAX_DEPTH ) 
    g_error("%s: tree already at maximum depth %u", 
	    __FUNCTION__, AFMM_TREE_MAX_DEPTH) ;

  t->depth ++ ;
  d = t->depth ;
  nb = 1 << d ;
  t->boxes[d] = (afmm_box_t *)g_malloc0(nb*nb*sizeof(afmm_box_t)) ;
  
  return 0 ;
}

guint64 afmm_box_parent(guint64 idx)

{
  guint64 p ;

  p = idx >> 2 ;

  return p ;
}

guint64 afmm_box_first_child(guint64 idx)

{
  guint64 c ;

  c = idx << 2 ;

  return c ;
}

guint64 afmm_point_locate_box(guint64 x, guint level)

{
  guint64 i ;

  i = x >> 2*(20 - level) ;

  return i ;
}

/* 
 * Morton encoding and decoding from:
 * https://stackoverflow.com/questions/30539347/2d-morton-code-encode-decode-64bits
 * (could use a bit-twiddling implementation)
 */

static guint32 morton_1(guint64 x)

{
    x = x & 0x5555555555555555;
    x = (x | (x >> 1))  & 0x3333333333333333;
    x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0F;
    x = (x | (x >> 4))  & 0x00FF00FF00FF00FF;
    x = (x | (x >> 8))  & 0x0000FFFF0000FFFF;
    x = (x | (x >> 16)) & 0x00000000FFFFFFFF;
    return (guint32)x ;
}

guint64 afmm_index_encode(guint32 x, guint32 y)
  
{
  guint64 z = 0 ;

  for ( gint i = 0 ; i < sizeof(x) * 8 ; i ++ ) {
    z |= (x & (guint64)1 << i) << i | (y & (guint64)1 << i) << (i + 1) ;
  }

  return z ;
}

gint afmm_index_decode(guint64 d, guint32 *x, guint32 *y)

{
  *x = morton_1(d) ;
  *y = morton_1(d >> 1) ;

  return 0 ;
}

static gint add_neighbour(guint64 b[], gint n, guint32 i, guint32 j)

{

  b[n] = afmm_index_encode(i, j) ;
  
  return n+1 ;
}

gint afmm_box_neighbours(guint64 i, guint level, guint64 b[])

{
  guint32 x, y, nb ;
  gint n = 1 ;

  nb = 1 << level ;

  b[0] = i ;

  afmm_index_decode(i, &x, &y) ;
  if ( x > 0 ) {
    n = add_neighbour(b, n, x-1, y) ;
    if ( y >    0 ) n = add_neighbour(b, n, x-1, y-1) ;
    if ( y < nb-1 ) n = add_neighbour(b, n, x-1, y+1) ;    
  }

  if ( x < nb-1 ) {
    n = add_neighbour(b, n, x+1, y) ;
    if ( y >    0 ) n = add_neighbour(b, n, x+1, y-1) ;
    if ( y < nb-1 ) n = add_neighbour(b, n, x+1, y+1) ;    
  }

  if ( y >    0 ) n = add_neighbour(b, n, x, y-1) ;    
  if ( y < nb-1 ) n = add_neighbour(b, n, x, y+1) ;    

  return n ;
}

gboolean afmm_grid_boxes_separated(guint ir, guint iz,
				   guint jr, guint jz,
				   gdouble etol)

{
  guint rmin, rmax,  zmin, zmax ;

  if ( ir == jr && iz == jz ) return FALSE ;
  
  if ( ir == jr ) {
    rmin = ir + 1 ; rmax = ir + 1 ;
    zmin = MIN(iz, jz) + 1 ;
    zmax = MAX(iz, jz) ;
    return (rmin*rmin + rmax*rmax + (zmax-zmin)*(zmax-zmin) >
	    2.0*rmin*rmax*(1.0 + etol)) ;
  }
  
  rmin = MIN(ir, jr) + 1 ;
  rmax = MAX(ir, jr)     ;
  if ( iz == jz ) {
    /* zmin = zmax = iz ; */
    zmin = iz ;
    zmax = iz ;
  } else {
    zmin = MIN(iz, jz) + 1 ;
    zmax = MAX(iz, jz) ;
  }
  
  return (rmin*rmin + rmax*rmax + (zmax-zmin)*(zmax-zmin) >
	  2.0*rmin*rmax*(1.0 + etol)) ;
}

gboolean afmm_boxes_separated(guint i, guint j, gdouble etol)

{
  guint ir, iz, jr, jz ;

  afmm_index_decode(i, &ir, &iz) ;
  afmm_index_decode(j, &jr, &jz) ;

  return afmm_grid_boxes_separated(ir, iz, jr, jz, etol) ;
}

gint afmm_radius_interaction_list(guint d, guint i0,
				  gdouble etol,
				  guint *ilist, gint *ni, gint nimax)

/*
 * ilist stride is 3: i j interaction_type
 */
  
{
  guint i, j, nb ;

  *ni = 0 ;
  nb = 1 << d ;

  /*
   * this looks for interaction boxes outboard of (greater radius
   * than) i0
   */
  for ( i = i0 ; i < nb ; i ++ ) {
    for ( j = 0 ; j < nb ; j ++ ) {
      if ( !afmm_grid_boxes_separated(i0, 0, i, j, etol) ) {
	ilist[(*ni)*3+0] = i ;
	ilist[(*ni)*3+1] = j ;
	ilist[(*ni)*3+2] = AFMM_INTERACTION_DIRECT ;
	g_assert((*ni) < nimax) ;
	(*ni) ++ ;
      } else {
	/*check separation of parent boxes*/
	if ( !afmm_grid_boxes_separated(i0/2, 0, i/2, j/2, etol) ) {
	  ilist[(*ni)*3+0] = i ;
	  ilist[(*ni)*3+1] = j ;
	  ilist[(*ni)*3+2] = AFMM_INTERACTION_S2L ;
	  g_assert((*ni) < nimax) ;
	  (*ni) ++ ;
	}
      }
    }
  }

  return 0 ;
}
				  
static void write_radius_box(FILE *f, guint *ilist)

{
  guint i, j, inter ;
  gdouble di ;
  
  i = ilist[0] ; j = ilist[1] ; inter = ilist[2] ;
  di = 1.0 ;
  
  fprintf(f, "%lg %lg  %lg %lg   %lg %lg   %lg %lg  %u\n",
	  di*i, di*j, di*(i+1), di*j, di*(i+1), di*(j+1), di*i, di*(j+1),
	  inter) ;
  
  return ;
}

gint afmm_write_radius_interaction_list(guint d, guint i,
					gdouble etol, FILE *f)

{
  guint ilist[512*3] ;
  gint ni, j ;
  
  afmm_radius_interaction_list(d, i, etol, ilist, &ni, 512) ;

  for ( j = 0 ; j < ni ; j ++ ) {
    write_radius_box(f, &(ilist[3*j])) ;
  }
  
  return 0 ;
}
