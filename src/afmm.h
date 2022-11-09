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

#ifndef _AFMM_H_INCLUDED_
#define _AFMM_H_INCLUDED_

#include <fftw3.h>

#define AFMM_TREE_MAX_DEPTH 16
#define AFMM_INDEX_SCALE (1LU << 63)
#define AFMM_INDEX_SHIFT (1U << 20)

typedef struct _afmm_box_t afmm_box_t ;
struct _afmm_box_t {
  guint32 i, /**< index of first source point in box */
    n ;      /**< number of source points in box */
  gpointer Cs, Cf ; /*coefficients of source and field expansions*/
} ;

typedef struct _afmm_tree_t afmm_tree_t ;
struct  _afmm_tree_t {
  gdouble rmin, rmax, zmin, zmax, esep ;
  afmm_box_t *boxes[AFMM_TREE_MAX_DEPTH+1] ;
  guint
  N,
    maxpoints,
    npoints,
    *ip,
    ns,
    depth,
    order_s[AFMM_TREE_MAX_DEPTH+1],
    order_f[AFMM_TREE_MAX_DEPTH+1] ;
  gchar *points ;
  gpointer Cs[AFMM_TREE_MAX_DEPTH+1], Cf[AFMM_TREE_MAX_DEPTH+1], S ; 
  gsize
  size,
    pstr ;
} ;

#define afmm_tree_point_number_max(_t) ((_t)->maxpoints)
#define afmm_tree_point_number(_t)     ((_t)->npoints)
#define afmm_tree_mode_number(_t)      ((_t)->N)
#define afmm_tree_source_size(_t)      ((_t)->ns)
#define afmm_tree_r_min(_t)            ((_t)->rmin)
#define afmm_tree_r_max(_t)            ((_t)->rmax)
#define afmm_tree_delta_r(_t)          ((afmm_tree_r_max(_t)) - \
					(afmm_tree_r_min(_t)))
#define afmm_tree_z_min(_t)            ((_t)->zmin)
#define afmm_tree_z_max(_t)            ((_t)->zmax)
#define afmm_tree_delta_z(_t)          ((afmm_tree_z_max(_t)) - \
					(afmm_tree_z_min(_t)))
#define afmm_tree_depth(_t)            ((_t)->depth)
#define afmm_tree_source_order(_t,_l)  ((_t)->order_s[(_l)])
#define afmm_tree_field_order(_t,_l)   ((_t)->order_f[(_l)])
#define afmm_tree_source_data(_t)      ((_t)->S)
#define afmm_tree_box_separation(_t)   ((_t)->esep)
/* #define afmm_tree_ */

#define afmm_derivative_offset(_q) ((_q)*((_q)+1)*((_q)+2)/6)
#define afmm_derivative_index_ijk(_i,_j,_k)	\
  ((2*((_i)+(_j)+(_k))-(_i)+3)*(_i)/2+(_j))
#define afmm_derivative_index_ij(_i,_j)	(((_i)+(_j))*((_i)+(_j)+1)/2+(_i))
#define afmm_derivative_offset_2(_i) 	((_i)*((_i)+1)/2)
#define afmm_coefficient_number(_i)		\
  afmm_derivative_offset_2(((_i)+1))

gint afmm_elliptic_KE(gdouble k, gdouble *K, gdouble *E) ;
gint afmm_elliptic_KE_f(gfloat k, gfloat *K, gfloat *E) ;

gint afmm_laplace_gfunc(gint N,
			gdouble r, gdouble r1, gdouble z,
			gdouble *G, gint str) ;
gint afmm_laplace_gfunc_f(gint N,
			  gfloat r, gfloat r1, gfloat z,
			  gfloat *G, gint str) ;

gint afmm_laplace_gfunc_derivatives(gint N, gint L,
				    gdouble r, gdouble r1, gdouble z,
				    gdouble *dG) ;
gint afmm_laplace_gfunc_derivatives_f(gint N, gint L,
				      gfloat r, gfloat r1, gfloat z,
				      gfloat *dG) ;
gdouble afmm_gfunc_series_eval(gint L, gdouble *dG,
			       gdouble dr, gdouble dr1, gdouble dx) ;
gfloat afmm_gfunc_series_eval_f(gint L, gfloat *dG,
				gfloat dr, gfloat dr1, gfloat dx) ;
gint afmm_gfunc_coefficients_write(FILE *f, gint L, gdouble *dG) ;
gint afmm_gfunc_coefficients_write_f(FILE *f, gint L, gfloat *dG) ;

gint afmm_laplace_field_direct(gdouble *rzsrc, gdouble *src, gint sstr,
			       gint ns,
			       gint nsrc,
			       gint N,
			       gdouble *rzfld, gdouble *fld, gint fstr,
			       gint nfld,
			       gboolean zero,
			       gdouble *work) ;
gint afmm_laplace_field_direct_f(gfloat *rzsrc, gfloat *src, gint sstr,
				 gint ns,
				 gint nsrc,
				 gint N,
				 gfloat *rzfld, gfloat *fld, gint fstr,
				 gint nfld,
				 gboolean zero,
				 gfloat *work) ;

gint afmm_laplace_gfunc_comp(gint N,
			     gdouble r, gdouble r1, gdouble z,
			     gdouble *G) ;
gint afmm_laplace_gfunc_comp_f(gint N,
			       gfloat r, gfloat r1, gfloat z,
			       gfloat *G) ;

gint afmm_laplace_field_eval(gint N, gint L, gdouble *dG, gint nd,
			     gdouble *S, gint sdist, gint ns,
			     gdouble *f, gint fdist,
			     gboolean forward, gboolean outward,
			     gboolean zero) ;
gint afmm_laplace_field_eval_f(gint N, gint L, gfloat *dG, gint nd,
			       gfloat *S, gint sdist, gint ns,
			       gfloat *f, gint fdist,
			       gboolean forward, gboolean outward,
			       gboolean zero) ;
gint afmm_laplace_source_to_local(gint N, gint L,
				  gdouble *dG, gint nd,
				  gdouble *S, gint sdist, gint LS,
				  gint ns,
				  gdouble *P, gint pdist, gint LP,
				  gboolean forward, gboolean outward) ;
gint afmm_laplace_source_to_local_f(gint N, gint L,
				    gfloat *dG, gint nd,
				    gfloat *S, gint sdist, gint LS,
				    gint ns,
				    gfloat *P, gint pdist, gint LP,
				    gboolean forward, gboolean outward) ;
gint afmm_laplace_s2l_matrix(gint N, gint L, gdouble *dG, gint nd,
			     gint LS, gint LP,
			     gboolean forward, gboolean outward,
			     gdouble *S2L) ;
gint afmm_laplace_s2l_matrix_f(gint N, gint L, gfloat *dG, gint nd,
			       gint LS, gint LP,
			       gboolean forward, gboolean outward,
			       gfloat *S2L) ;
gint afmm_laplace_shift_s2l(gint N, gdouble *S2L,
			    gdouble *S, gint sdist, gint LS,
			    gint ns,
			    gdouble *P, gint pdist, gint LP) ;
gint afmm_laplace_shift_s2l_f(gint N, gfloat *S2L,
			      gfloat *S, gint sdist, gint LS,
			      gint ns,
			      gfloat *P, gint pdist, gint LP) ;

gint afmm_expansion_eval(gdouble dr, gdouble dz,
			 gint N, gint L,
			 gdouble *P, gint pdist,
			 gint ns,
			 gdouble *f, gint fdist) ;
gint afmm_expansion_eval_f(gfloat dr, gfloat dz,
			   gint N, gint L,
			   gfloat *P, gint pdist,
			   gint ns,
			   gfloat *f, gint fdist) ;

gint afmm_laplace_modes_to_field(gdouble *modes, gint ns, gint nm,
				 gint N, fftw_plan plan) ;
gint afmm_laplace_modes_to_field_f(gfloat *modes, gint ns, gint nm,
				   gint N, fftwf_plan plan) ;
fftw_plan afmm_laplace_modes_to_field_plan(gdouble *modes, gint ns, gint np,
					   gint dist, gint N) ;
fftwf_plan afmm_laplace_modes_to_field_plan_f(gfloat *modes, gint ns, gint np,
					      gint dist, gint N) ;

gdouble afmm_hypergeometric_2F1(gdouble a, gdouble b, gdouble c, gdouble z) ;
gfloat afmm_hypergeometric_2F1_f(gfloat a, gfloat b, gfloat c, gfloat z) ;

gint afmm_legendre_Q(gint N, gint M, gdouble chi, gdouble *Q, gint ldq) ;
gint afmm_legendre_Q_f(gint N, gint M, gfloat chi, gfloat *Q, gint ldq) ;

gint afmm_source_moments(gdouble r1, gdouble z1,
			 gdouble *C, gint N, gint L, gdouble *S, gint sdist) ;
gint afmm_source_moments_f(gfloat r1, gfloat z1,
			   gfloat *C, gint N, gint L, gfloat *S, gint sdist) ;
gint afmm_moments_shift(gint N,
			gint Mi,
			gdouble *Si, gint idist, gint ns,
			gdouble dr, gdouble dz,
			gint Mo,
			gdouble *So, gint odist) ;
gint afmm_moments_shift_f(gint N,
			  gint Mi,
			  gfloat *Si, gint idist, gint ns,
			  gfloat dr, gfloat dz,
			  gint Mo,
			  gfloat *So, gint odist) ;
gint afmm_expansion_shift(gint N,
			  gint Mi,
			  gdouble *Pi, gint idist, gint ns,
			  gdouble dr, gdouble dz,
			  gint Mo,
			  gdouble *Po, gint odist) ;
gint afmm_expansion_shift_f(gint N,
			    gint Mi,
			    gfloat *Pi, gint idist, gint ns,
			    gfloat dr, gfloat dz,
			    gint Mo,
			    gfloat *Po, gint odist) ;

afmm_tree_t *afmm_tree_new(gdouble rmin, gdouble rmax,
			   gdouble zmin, gdouble zmax,
			   guint maxpoints) ;
afmm_tree_t *afmm_tree_new_f(gfloat rmin, gfloat rmax,
			     gfloat zmin, gfloat zmax,
			     guint maxpoints) ;
guint64 afmm_point_index_2d(gdouble *rz,
			    gdouble rmin, gdouble rmax,
			    gdouble zmin, gdouble zmax) ;
guint64 afmm_point_index_2d_f(gfloat *rz,
			      gfloat rmin, gfloat rmax,
			      gfloat zmin, gfloat zmax) ;
gint afmm_tree_add_points(afmm_tree_t *t,
			  gpointer pts, gsize pstr,
			  guint npts, gboolean sorted) ;
gint afmm_tree_add_points_f(afmm_tree_t *t,
			    gpointer pts, gsize pstr,
			    guint npts, gboolean sorted) ;
gint afmm_tree_add_level(afmm_tree_t *t) ;
gint afmm_tree_refine(afmm_tree_t *t) ;
gint afmm_tree_refine_f(afmm_tree_t *t) ;
gint afmm_index_decode(guint64 d, guint32 *x, guint32 *y) ;
guint64 afmm_index_encode(guint32 x, guint32 y) ;
gboolean afmm_boxes_are_separated(afmm_tree_t *t, guint level,
				  guint i, guint j) ;

guint64 afmm_box_first_child(guint64 idx) ;
guint64 afmm_box_parent(guint64 idx) ;
gint afmm_box_neighbours(guint64 i, guint level, guint64 b[]) ;
guint64 afmm_point_locate_box(guint64 x, guint level) ;
gint afmm_box_location_from_index(guint64 im, guint level,
				  gdouble rmin, gdouble rmax,
				  gdouble zmin, gdouble zmax,
				  gdouble *r, gdouble *rb,
				  gdouble *z, gdouble *zb) ;
gint afmm_box_location_from_index_f(guint64 im, guint level,
				    gfloat rmin, gfloat rmax,
				    gfloat zmin, gfloat zmax,
				    gfloat *r, gfloat *rb,
				    gfloat *z, gfloat *zb) ;

gint afmm_tree_coefficient_init(afmm_tree_t *t, guint l,
				 guint nf, guint ns) ;
gint afmm_tree_coefficient_init_f(afmm_tree_t *t, guint l,
				  guint nf, guint ns) ;
gint afmm_tree_leaf_expansions(afmm_tree_t *t,
			       gdouble *C, gint cdist,
			       gboolean zero_expansions) ;
gint afmm_tree_leaf_expansions_f(afmm_tree_t *t,
				 gfloat *C, gint cdist,
				 gboolean zero_expansions) ;
gint afmm_tree_expansion_eval_field(afmm_tree_t *t, guint d,
				    gdouble r, gdouble z,
				    gdouble *f, gdouble *work) ;
gint afmm_tree_expansion_eval_field_f(afmm_tree_t *t, guint d,
				      gfloat r, gfloat z,
				      gfloat *f, gfloat *work) ;
gint afmm_tree_local_field_eval(afmm_tree_t *t, gdouble r, gdouble z,
				gdouble *f, gint fdist) ;
gint afmm_tree_local_field_eval_f(afmm_tree_t *t, gfloat r, gfloat z,
				  gfloat *f, gint fdist) ;
gint afmm_box_field_direct(afmm_tree_t *t, afmm_box_t *b,
			   gdouble *src, gint sdist, gint ns, gint N,
			   gdouble *rz, gdouble *f, gint fdist,
			   gboolean zero, gdouble *work) ;
gint afmm_box_field_direct_f(afmm_tree_t *t, afmm_box_t *b,
			     gfloat *src, gint sdist, gint ns, gint N,
			     gfloat *rz, gfloat *f, gint fdist,
			     gboolean zero, gfloat *work) ;
  
gint afmm_interaction_list(guint64 b, guint depth, guint64 *ilist, gint *ni) ;
gint afmm_upward_pass(afmm_tree_t *t, guint level, gdouble *work) ;
gint afmm_upward_pass_f(afmm_tree_t *t, guint level, gfloat *work) ;
gint afmm_downward_pass(afmm_tree_t *t, guint level, gdouble *work,
			gboolean downward) ;
gint afmm_downward_pass_f(afmm_tree_t *t, guint level, gfloat *work,
			gboolean downward) ;

#endif /*_AFMM_H_INCLUDED_*/
