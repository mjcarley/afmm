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

static void add_interactions(guint32 x,  guint32 y,
			     guint32 xp, guint32 yp,
			     guint64 *ilist, gint *ni)

{
  guint64 c, p ;
  guint32 xc, yc ;
  gint i ;

  p = afmm_index_encode(xp, yp) ;
  c = afmm_box_first_child(p) ;
  for ( i = 0 ; i < 4 ; i ++ ) {
    afmm_index_decode(c, &xc, &yc) ;  
    /* if ( ((xc > x+1) || (yc > y+1) || (yc < y-1)) && (xc >= x) ) { */
    /*   ilist[(*ni)] = c ; (*ni) ++ ; */
    /* } */
    if ( ( yc  > y + 1) ||
	 (( yc == y + 1) && ( xc > x + 1)) ||
	 (( yc == y + 0) && ( xc > x + 1)) ||
	 (( yc == y - 1) && ( xc > x + 1)) ||
	 (( yc  < y - 1) && ( xc > x)) ) {
      ilist[(*ni)] = c ; (*ni) ++ ;
    }
    
    c ++ ;
  }
  
  return ;
}

gint afmm_interaction_list(guint64 b, guint depth,
			   guint64 *ilist, gint *ni)

{
  guint64 p ;
  guint32 xp, yp, lp, x, y ;

  *ni = 0 ;
  afmm_index_decode(b, &x, &y) ;

  /*identify parent*/
  p = afmm_box_parent(b) ;
  afmm_index_decode(p, &xp, &yp) ;

  /*lp is number of boxes per side at parent level*/
  lp = 1 << (depth - 1) ;

  /*consider three parent-level boxes at depth-1*/
  if ( xp + 1 < lp ) {
    add_interactions(x, y, xp + 1, yp + 0, ilist, ni) ;
    if ( yp + 1 < lp ) {
      add_interactions(x, y, xp + 1, yp + 1, ilist, ni) ;
    }
    if ( yp - 1 >= 0 ) {
      add_interactions(x, y, xp + 1, yp - 1, ilist, ni) ;
    }
  }
  if ( yp + 1 < lp ) {
    add_interactions(x, y, xp + 0, yp + 1, ilist, ni) ;
  }
  if ( yp - 1 >= 0 ) {
    add_interactions(x, y, xp + 0, yp - 1, ilist, ni) ;
  }
  
  return 0 ;
}

