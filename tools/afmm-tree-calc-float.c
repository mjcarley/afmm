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

GTimer *timer ;
gchar *progname ;

static void points_read(FILE *f, gint *np, guint *N, guint *ns,
			gfloat **rz, gfloat **src)

{
  gint i, j, k, str ;
  
  fscanf(f, "%d %d %d", np, N, ns) ;

  /*FFTW3 stride for C2R transforms*/
  str = 2*((*N)+2) ;

  *rz = (gfloat *)g_malloc0((*np)*2*sizeof(gfloat)) ;
  if ( *ns != 0 ) {
    *src = (gfloat *)
      fftw_malloc((*np)*str*(*ns)*sizeof(gfloat)) ;
  } else {
    src = NULL ;
  }

  for ( i = 0 ; i < *np ; i ++ ) {
    fscanf(f, "%g %g", &((*rz)[2*i+0]), &((*rz)[2*i+1])) ;
    if ( src != NULL ) {
      for ( j = 0 ; j < *ns ; j ++ ) {
	for ( k = 0 ; k <= *N ; k ++ ) {
	  fscanf(f, "%g", &((*src)[i*(*ns)*str + j*str + 2*k + 0])) ;
	  fscanf(f, "%g", &((*src)[i*(*ns)*str + j*str + 2*k + 1])) ;
	}
      }
    }
  }
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gfloat rmin, rmax, zmin, zmax, *rz1, *rz, *C, *fld, *work ;
  guint npts, ns, N, Nf, order[48] = {0}, order_s, order_f, order_max ;
  guint depth, order_s_max, order_f_max, r_trace, z_trace ;
  gint i, j, k, nrz1, cdist, nfld, testdepth, wsize ;
  gchar ch, *sfile, *ffile ;
  afmm_tree_t *tree ;
  gboolean downward, trace_field, trace_source ;
  FILE *input ;
  
  rmin = 0 ; rmax = 1 ;
  zmin = 0 ; zmax = 1 ;

  order_f = 4 ; order_s = 4 ; depth = 5 ;
  testdepth = -1 ;
  sfile = ffile = NULL ;
  /*set to FALSE to suppress downward propagation of field
    expansions (so only local interactions are included)*/
  downward = TRUE ;
  r_trace = z_trace = 65536 ;
  trace_field = FALSE ; trace_source = FALSE ;
  /* trace_field = TRUE ; */
  /* trace_source = TRUE ; */
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  while ( (ch = getopt(argc, argv, "D:d:F:f:I:J:lR:r:S:s:Z:z:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'D': depth = atoi(optarg) ; break ;
    case 'd': testdepth = atoi(optarg) ; break ;
    case 'F': order_f = atoi(optarg) ; break ;
    case 'f': ffile = g_strdup(optarg) ; break ;
    case 'I': r_trace = atoi(optarg) ; break ;
    case 'J': z_trace = atoi(optarg) ; break ;
    case 'l': downward = FALSE ; break ;
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

  points_read(input, &nfld, &Nf, &ns, &rz, &fld) ;

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

  points_read(input, &nrz1, &N, &ns, &rz1, &C) ;

  fclose(input) ;
  
  fprintf(stderr, "%s: %d source points\n", progname, nrz1) ;
  fprintf(stderr, "%s: modal order %d\n", progname, N) ;
  fprintf(stderr, "%s: %d source components\n", progname, ns) ;
  
  fld = (gfloat *)fftw_malloc(nfld*(2*N+4)*ns*sizeof(gfloat)) ;
  memset(fld, 0, nfld*(2*N+4)*ns*sizeof(gfloat)) ;

  npts = nrz1 + nfld ;

  tree = afmm_tree_new_f(rmin, rmax, zmin, zmax, npts) ;

  /* afmm_tree_box_separation(tree) = 1e-2 ; */
  
  afmm_tree_source_size(tree) = ns ;
  afmm_tree_mode_number(tree) = N ;

  /*set expansion orders at each level, for source and field*/
  order[2*depth+0] = order_s ;
  order[2*depth+1] = order_f ;
  order_max = order_s + order_f ;
  order_s_max = order_s ;
  order_f_max = order_f ;
  for ( i = depth-1 ; i > 0 ; i -- ) {
    order[2*i+0] = order[2*(i+1)+0] + 4 ;
    order[2*i+1] = order[2*(i+1)+1] + 4 ;
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
  work = (gfloat *)g_malloc0(wsize*sizeof(gfloat)) ;
  
  /*add points to the tree*/
  afmm_tree_add_points_f(tree, rz1, 2*sizeof(gfloat), nrz1, FALSE) ;
  
  for ( i = 1 ; i <= depth ; i ++ ) {
    afmm_tree_refine_f(tree) ;
  }

  if ( r_trace < 65536 && z_trace == 65536 ) {

    afmm_trace_interactions(tree, r_trace, stdout) ;

    return 0 ;
  }

  /* if ( r_trace < 65536 && z_trace < 65536 ) { */

  /*   afmm_write_interaction_list(tree, afmm_tree_depth(tree), */
  /* 				r_trace, z_trace, stdout) ; */

  /*   return 0 ; */
  /* } */
  
  if ( trace_source ) {
    i = 0 ;
    guint64 idx, jdx ;
    idx = afmm_point_index_2d_f(&(rz1[2*i]),
				    afmm_tree_r_min(tree),
				    afmm_tree_r_max(tree),
				    afmm_tree_z_min(tree),
				    afmm_tree_z_max(tree)) ;
    fprintf(stderr, "%d:", i) ;
    for ( j = 0 ; j <= afmm_tree_depth(tree) ; j ++ ) {
      jdx = afmm_point_locate_box(idx, j) ;
      fprintf(stderr, " %lu", jdx) ;
    }
    fprintf(stderr, "\n") ;
    return 0 ;
  }
  
  /*source data*/
  cdist = (2*N+4) ;
  
  for ( i = 1 ; i <= depth ; i ++ ) {
    afmm_tree_coefficient_init(tree, i, order[2*i+1], order[2*i+0]) ;
  }

  fprintf(stderr, "%s: initializing leaf expansions; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  afmm_tree_leaf_expansions_f(tree, C, cdist, TRUE) ;

  fprintf(stderr, "%s: leaf expansions initialized; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  fprintf(stderr, "%s: upward pass; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = depth ; i >= 2 ; i -- ) {
    afmm_upward_pass_f(tree, i, work) ;
  }

  fprintf(stderr, "%s: upward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  if ( testdepth != -1 ) {
    /*check field evaluation by summing box source expansion fields*/
    fprintf(stderr, "%s: evaluating field using level %d boxes\n",
	    progname, testdepth) ;
    for ( i = 0 ; i < nfld ; i ++ ) {
      /* fprintf(stderr, "%d\n", i) ; */
      afmm_tree_expansion_eval_field_f(tree, testdepth,
					   rz[2*i+0], rz[2*i+1],
					   &(fld[i*ns*(2*N+4)]), work) ;
    }

    fprintf(stderr, "%s: field evaluated; %lg\n",
	    progname, g_timer_elapsed(timer, NULL)) ;
    for ( i = 0 ; i < nfld ; i ++ ) {
      fprintf(stdout, "%g %g", rz[2*i+0], rz[2*i+1]) ;
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
    afmm_downward_pass_f(tree, i, work, downward) ;
  }

  fprintf(stderr, "%s: downward pass completed; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  for ( i = 0 ; i < nfld ; i ++ ) {
    /* fprintf(stderr, "%d\n", i) ; */
    /* afmm_tree_expansion_eval_field_f(tree, testdepth, */
    /* 					 rz[2*i+0], rz[2*i+1], */
    /* 					 &(fld[i*ns*(2*N+4)]), work) ; */
    if ( trace_field ) {
      guint64 idx, jdx ;
      idx = afmm_point_index_2d_f(&(rz[2*i]),
				      afmm_tree_r_min(tree),
				      afmm_tree_r_max(tree),
				      afmm_tree_z_min(tree),
				      afmm_tree_z_max(tree)) ;
      fprintf(stderr, "%d:", i) ;
      for ( j = 0 ; j <= afmm_tree_depth(tree) ; j ++ ) {
	jdx = afmm_point_locate_box(idx, j) ;
	fprintf(stderr, " %lu", jdx) ;
      }
      fprintf(stderr, "\n") ;
    }
    afmm_tree_local_field_eval_f(tree, 
				     rz[2*i+0], rz[2*i+1],
				     &(fld[i*ns*(2*N+4)]), 2*N+4) ;
  }
  
  fprintf(stderr, "%s: field evaluated; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  for ( i = 0 ; i < nfld ; i ++ ) {
    fprintf(stdout, "%g %g", rz[2*i+0], rz[2*i+1]) ;
    for ( j = 0 ; j < ns ; j ++ ) {
      for ( k = 0 ; k < 2*(N+1) ; k ++ ) {
	fprintf(stdout, " %1.16e", fld[i*ns*2*(N+2) + j*2*(N+2) + k]) ;
      }
    }
    fprintf(stdout, "\n") ;
  }
  
  return 0 ;
}

