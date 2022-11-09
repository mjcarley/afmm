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
  guint nb ;

  if ( t->depth >= AFMM_TREE_MAX_DEPTH ) 
    g_error("%s: tree already at maximum depth %u", 
	    __FUNCTION__, AFMM_TREE_MAX_DEPTH) ;

  t->depth ++ ;
  nb = 1 << (2*(t->depth)) ;
  t->boxes[t->depth] = (afmm_box_t *)g_malloc0(nb*sizeof(afmm_box_t)) ;

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

gboolean afmm_boxes_are_separated(afmm_tree_t *t, guint level,
				  guint i, guint j)

{
  gdouble ri, zi, rj, zj, dr, dz, chi, r, r1, z ;
  
  afmm_box_location_from_index(i, level,
			       afmm_tree_r_min(t), afmm_tree_r_max(t),
			       afmm_tree_z_min(t), afmm_tree_z_max(t),
			       &ri, &dr, &zi, &dz) ;
  afmm_box_location_from_index(j, level,
			       afmm_tree_r_min(t), afmm_tree_r_max(t),
			       afmm_tree_z_min(t), afmm_tree_z_max(t),
			       &rj, &dr, &zj, &dz) ;

  z = MAX(zi, zj) - MIN(zi, zj) ;
  if ( z != 0 ) z -= dz ;

  r1 = MIN(ri, rj) ; r = MAX(ri, rj) ;
  g_assert(r >= r1) ;

  if ( r != r1 ) r1 += dr ;

  chi = 0.5*(r*r + r1*r1 + z*z)/r/r1 ;
  
  return (chi > 1 + afmm_tree_box_separation(t)) ;
}
