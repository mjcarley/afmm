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

#ifndef AFMM_PRIVATE_H_INCLUDED 
#define AFMM_PRIVATE_H_INCLUDED

#include <stdio.h>

#ifdef AFMM_SINGLE_PRECISION

#define AFMM_REAL gfloat

#define AFMM_FUNCTION_NAME(_func) _func##_f

#define SQRT(_x) sqrtf((_x))
#define CBRT(_x) cbrtf((_x))
#define SIN(_x) sinf((_x))
#define COS(_x) cosf((_x))
#define ACOS(_x) acosf((_x))
#define ACOSH(_x) acoshf((_x))
#define ATAN(_x) atanf((_x))
#define ATAN2(_y,_x) atan2f((_y),(_x))
#define LOG(_x) logf((_x))
#define EXP(_x) expf((_x))

extern const gfloat AFMM_FACTORIALS_F[] ;
#define afmm_factorial(_n) (AFMM_FACTORIALS_F[(_n)])

extern const gfloat AFMM_N_SQUARED_F[] ;

#else

#define AFMM_REAL gdouble

#define AFMM_FUNCTION_NAME(_func) _func

#define SQRT(_x) sqrt((_x))
#define CBRT(_x) cbrt((_x))
#define SIN(_x) sin((_x))
#define COS(_x) cos((_x))
#define ACOS(_x) acos((_x))
#define ACOSH(_x) acosh((_x))
#define ATAN(_x) atan((_x))
#define ATAN2(_y,_x) atan2((_y),(_x))
#define LOG(_x) log((_x))
#define EXP(_x) exp((_x))

extern const gdouble AFMM_FACTORIALS[] ;
#define afmm_factorial(_n) (AFMM_FACTORIALS[(_n)])

extern const gdouble AFMM_N_SQUARED[] ;

#endif /*AFMM_SINGLE_PRECISION*/

extern const gint AFMM_BINOMIALS[] ;
/* #define afmm_binomial(_m,_k) (AFMM_BINOMIALS[(_m)*((_m)+1)/2+(_k)]) */

#define afmm_binomial(_m,_k)					\
  (exp(lgamma((_m)+1) - lgamma((_k)+1) - lgamma((_m)-(_k)+1)))

/* gdouble afmm_binomial(guint n, guint k) ; */

#define afmm_tree_point_index(_t,_i)			\
  ((AFMM_REAL *)(&((_t)->points[(_i)*((_t)->pstr)])))

extern const guint afmm_ilist_di[], afmm_ilist_dj[] ;

#endif /*AFMM_PRIVATE_H_INCLUDED*/
