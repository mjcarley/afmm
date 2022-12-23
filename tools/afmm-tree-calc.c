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

#include <afmm.h>

#include "afmm-private.h"

#define AFMM_INTERLEAVE_SOURCE_DATA

GTimer *timer ;
gchar *progname ;

static void points_read(FILE *f, gint *np, guint *N, guint *ns,
			gdouble **rz, gint *pstr,
			gdouble **src, gint *sdist, gint *sstr)

{
  gint i, j, k ;
  
  fscanf(f, "%d %d %d", np, N, ns) ;

  /*FFTW3 stride for C2R transforms*/
  (*sstr) = 2*((*N)+2) ;
  *sdist = (*sstr)*(*ns) ;
  *pstr = 2 ;
#ifdef AFMM_INTERLEAVE_SOURCE_DATA
  *sdist += 2 ;
  *pstr = *sdist ;
#endif /*AFMM_INTERLEAVE_SOURCE_DATA*/
  *rz = (gdouble *)fftw_malloc((*np)*(*pstr)*sizeof(gdouble)) ;

  if ( *ns != 0 ) {
#ifdef AFMM_INTERLEAVE_SOURCE_DATA
    *src = &((*rz)[2]) ;
#else /*AFMM_INTERLEAVE_SOURCE_DATA*/
    *src = (gdouble *)fftw_malloc((*np)*(*sdist)*sizeof(gdouble)) ;
#endif /*AFMM_INTERLEAVE_SOURCE_DATA*/    
  } else {
    src = NULL ;
  }

  for ( i = 0 ; i < *np ; i ++ ) {
    fscanf(f, "%lg %lg",
	   &((*rz)[i*(*pstr)+0]), &((*rz)[i*(*pstr)+1])) ;
    if ( src != NULL ) {
      for ( j = 0 ; j < *ns ; j ++ ) {
	for ( k = 0 ; k <= *N ; k ++ ) {
	  fscanf(f, "%lg", &((*src)[i*(*sdist) + j*(*sstr) + 2*k + 0])) ;
	  fscanf(f, "%lg", &((*src)[i*(*sdist) + j*(*sstr) + 2*k + 1])) ;
	}
      }
    }
  }
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gdouble rmin, rmax, zmin, zmax, *rz1, *rz, *C, *fld, *work ;
  gdouble etol ;
  guint ns, N, Nf, order[48] = {0}, order_s, order_f, order_max ;
  guint depth, order_s_max, order_f_max ;
  /* guint r_trace, z_trace ; */
  gint i, j, k, nrz1, cdist, nfld, testdepth, wsize, order_inc ;
  gint fpstr, fstr, pstr, cstr ;
  gchar ch, *sfile, *ffile ;
  afmm_tree_t *tree ;
  gboolean downward, sort_sources ;
  /* trace_field, trace_source ; */
  FILE *input ;
  
  rmin = 0 ; rmax = 1 ;
  zmin = 0 ; zmax = 1 ;

  order_f = 4 ; order_s = 4 ; depth = 5 ;
  testdepth = -1 ;
  sfile = ffile = NULL ;
  /*set to FALSE to suppress downward propagation of field
    expansions (so only local interactions are included)*/
  downward = TRUE ;
  sort_sources = FALSE ;
  /* r_trace = z_trace = 65536 ; */
  /* trace_field = FALSE ; trace_source = FALSE ; */
  /* trace_field = TRUE ; */
  /* trace_source = TRUE ; */
  etol = 0.0 ;
  order_inc = 2 ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  while ( (ch = getopt(argc, argv, "D:d:e:F:f:lOo:R:r:S:s:Z:z:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'D': depth = atoi(optarg) ; break ;
    case 'd': testdepth = atoi(optarg) ; break ;
    case 'e': etol = atof(optarg) ; break ;
    case 'F': order_f = atoi(optarg) ; break ;
    case 'f': ffile = g_strdup(optarg) ; break ;
    /* case 'I': r_trace = atoi(optarg) ; break ; */
    /* case 'J': z_trace = atoi(optarg) ; break ; */
    case 'l': downward = FALSE ; break ;
    case 'O': sort_sources = TRUE ; break ;
    case 'o': order_inc = atoi(optarg) ; break ;
    case 'R': rmax = atof(optarg) ; break ; 
    case 'r': rmin = atof(optarg) ; break ; 
    case 'S': order_s = atoi(optarg) ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    case 'Z': zmax = atof(optarg) ; break ; 
    case 'z': zmin = atof(optarg) ; break ; 
    }
  }

  if ( ffile == NULL ) {
    fprintf(stderr, "%s: field point file must be defined\n", progname) ;
    exit(1) ;    
  }
  input = fopen(ffile, "r") ;
  if ( input == NULL ) {
    fprintf(stderr, "%s: cannot open field point file %s\n", progname, ffile) ;
    exit(1) ;
  }

  points_read(input, &nfld, &Nf, &ns, &rz, &fpstr, &fld, &fstr, &cstr) ;

  fclose(input) ;
  fprintf(stderr, "%s: %d field points\n", progname, nfld) ;
  
  if ( sfile == NULL ) {
    fprintf(stderr, "%s: source file must be defined\n", progname) ;
    exit(1) ;    
  }
  input = fopen(sfile, "r") ;
  if ( input == NULL ) {
    fprintf(stderr, "%s: cannot open source file %s\n", progname, sfile) ;
    exit(1) ;
  }

  points_read(input, &nrz1, &N, &ns, &rz1, &pstr, &C, &cdist, &cstr) ;

  fclose(input) ;

  if ( sort_sources ) {
    afmm_sort_point_list(rz1, pstr*sizeof(gdouble), nrz1,
			       rmin, rmax, zmin, zmax) ;
  }
  
  fprintf(stderr, "%s: %d source points\n", progname, nrz1) ;
  fprintf(stderr, "%s: modal order %d\n", progname, N) ;
  fprintf(stderr, "%s: %d source components; %lg\n",
	  progname, ns, g_timer_elapsed(timer, NULL)) ;
  
  fld = (gdouble *)fftw_malloc(nfld*(2*N+4)*ns*sizeof(gdouble)) ;
  memset(fld, 0, nfld*(2*N+4)*ns*sizeof(gdouble)) ;

  tree = afmm_tree_new(rmin, rmax, zmin, zmax, nrz1, nfld) ;
  if ( etol != 0.0 ) afmm_tree_box_separation(tree) = etol ;
  afmm_tree_source_size(tree) = ns ;
  afmm_tree_mode_number(tree) = N ;

  /*set expansion orders at each level, for source and field*/
  order[2*depth+0] = order_s ;
  order[2*depth+1] = order_f ;
  order_max = order_s + order_f ;
  order_s_max = order_s ;
  order_f_max = order_f ;
  for ( i = depth-1 ; i > 0 ; i -- ) {
    order[2*i+0] = order[2*(i+1)+0] + order_inc ;
    order[2*i+1] = order[2*(i+1)+1] + order_inc ;
    order_max = MAX(order_max, order[2*i+0] + order[2*i+1]) ;
    order_s_max = MAX(order_s_max, order[2*i+0]) ;
    order_f_max = MAX(order_f_max, order[2*i+1]) ;
  }

  /*size workspace for generation of Green's function derivatives*/
  fprintf(stderr, "%s: maximum order %d\n", progname, order_max) ;
  wsize = (N+1)*afmm_derivative_offset(order_max+2) +
    4*(N+1)*
    afmm_derivative_offset_2(order_f_max+1)*
    afmm_derivative_offset_2(order_s_max+1) ;
    
  fprintf(stderr, "%s: workspace size: %d\n", progname, wsize) ;
  work = (gdouble *)g_malloc0(wsize*sizeof(gdouble)) ;

  /*add points to the tree*/
  afmm_tree_add_sources(tree, rz1, pstr*sizeof(gdouble), nrz1,
			      sort_sources) ;
  afmm_tree_add_field(tree, rz, fpstr*sizeof(gdouble), nfld, FALSE) ;
  
  for ( i = 1 ; i <= depth ; i ++ ) {
    afmm_tree_refine(tree) ;
  }
  
  /* if ( r_trace < 65536 && z_trace == 65536 ) { */

  /*   afmm_trace_interactions(tree, r_trace, stdout) ; */

  /*   return 0 ; */
  /* } */

  /* if ( r_trace < 65536 && z_trace < 65536 ) { */

  /*   afmm_write_interaction_list(tree, afmm_tree_depth(tree), */
  /* 				r_trace, z_trace, stdout) ; */

  /*   return 0 ; */
  /* } */
  
  /* if ( trace_source ) { */
  /*   i = 0 ; */
  /*   guint64 idx, jdx ; */
  /*   idx = afmm_point_index_2d(&(rz1[2*i]), */
  /* 				    afmm_tree_r_min(tree), */
  /* 				    afmm_tree_r_max(tree), */
  /* 				    afmm_tree_z_min(tree), */
  /* 				    afmm_tree_z_max(tree)) ; */
  /*   fprintf(stderr, "%d:", i) ; */
  /*   for ( j = 0 ; j <= afmm_tree_depth(tree) ; j ++ ) { */
  /*     jdx = afmm_point_locate_box(idx, j) ; */
  /*     fprintf(stderr, " %lu", jdx) ; */
  /*   } */
  /*   fprintf(stderr, "\n") ; */
  /*   return 0 ; */
  /* } */
  
  /*source data*/
  /* cdist = (2*N+4) ; */
  
  for ( i = 1 ; i <= depth ; i ++ ) {
    afmm_tree_coefficient_init(tree, i, order[2*i+1], order[2*i+0]) ;
  }

#ifdef AFMM_PRECOMPUTE_S2L
  fprintf(stderr, "%s: initializing shift matrices; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;  
  afmm_tree_init_s2l(tree, 0, 0) ;
  fprintf(stderr, "%s: shift matrices initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;  
#endif /*AFMM_PRECOMPUTE_S2L*/
  
  fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  afmm_tree_leaf_expansions(tree, C, cdist, cstr, TRUE) ;

  fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: upward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = depth ; i >= 2 ; i -- ) {
    afmm_upward_pass(tree, i, work) ;
  }

  fprintf(stderr, "%s: upward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  if ( testdepth != -1 ) {
    /*check field evaluation by summing box source expansion fields*/
    fprintf(stderr, "%s: evaluating field using level %d boxes\n",
	    progname, testdepth) ;
    for ( i = 0 ; i < nfld ; i ++ ) {
      /* fprintf(stderr, "%d\n", i) ; */
      afmm_tree_expansion_eval_field(tree, testdepth,
					   rz[2*i+0], rz[2*i+1],
					   &(fld[i*ns*(2*N+4)]), work) ;
    }

    fprintf(stderr, "%s: field evaluated; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    for ( i = 0 ; i < nfld ; i ++ ) {
      fprintf(stdout, "%lg %lg", rz[2*i+0], rz[2*i+1]) ;
      for ( j = 0 ; j < ns ; j ++ ) {
	for ( k = 0 ; k < 2*(N+1) ; k ++ ) {
	  fprintf(stdout, " %1.16e", fld[i*ns*2*(N+2) + j*2*(N+2) + k]) ;
	}
      }
      fprintf(stdout, "\n") ;
    }
    return 0 ;
  }  

  fprintf(stderr, "%s: downward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 1 ; i < depth ; i ++ ) {
    afmm_downward_pass(tree, i, work, downward) ;
  }

  fprintf(stderr, "%s: downward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  /* for ( i = 0 ; i < nfld ; i ++ ) { */
  /*   if ( trace_field ) { */
  /*     guint64 idx, jdx ; */
  /*     idx = afmm_point_index_2d(&(rz[2*i]), */
  /* 				      afmm_tree_r_min(tree), */
  /* 				      afmm_tree_r_max(tree), */
  /* 				      afmm_tree_z_min(tree), */
  /* 				      afmm_tree_z_max(tree)) ; */
  /*     fprintf(stderr, "%d:", i) ; */
  /*     for ( j = 0 ; j <= afmm_tree_depth(tree) ; j ++ ) { */
  /* 	jdx = afmm_point_locate_box(idx, j) ; */
  /* 	fprintf(stderr, " %lu", jdx) ; */
  /*     } */
  /*     fprintf(stderr, "\n") ; */
  /*   } */
  /*   afmm_tree_local_field_eval(tree,  */
  /* 				     rz[2*i+0], rz[2*i+1], */
  /* 				     &(fld[i*ns*(2*N+4)]), 2*N+4) ; */
  /* } */

  afmm_tree_field_eval(tree, fld, 2*N+4) ;
  
  fprintf(stderr, "%s: field evaluated; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( i = 0 ; i < nfld ; i ++ ) {
    fprintf(stdout, "%lg %lg", rz[2*i+0], rz[2*i+1]) ;
    for ( j = 0 ; j < ns ; j ++ ) {
      for ( k = 0 ; k < 2*(N+1) ; k ++ ) {
	fprintf(stdout, " %1.16e", fld[i*ns*2*(N+2) + j*2*(N+2) + k]) ;
      }
    }
    fprintf(stdout, "\n") ;
  }
  
  return 0 ;
}

