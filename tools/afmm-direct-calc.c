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

#include <fftw3.h>

#include <afmm.h>

#include "afmm-private.h"

GTimer *timer ;
gchar *progname ;

static void points_read(FILE *f, gint *np, gint *N, gint *ns,
			gdouble **rz, gdouble **src)

{
  gint i, j, k, str ;
  
  fscanf(f, "%d %d %d", np, N, ns) ;

  /*FFTW3 stride for C2R transforms*/
  str = 2*((*N)+2) ;

  *rz = (gdouble *)g_malloc0((*np)*2*sizeof(gdouble)) ;
  if ( *ns != 0 ) {
    *src = (gdouble *)
      fftw_malloc((*np)*str*(*ns)*sizeof(gdouble)) ;
  } else {
    src = NULL ;
  }

  for ( i = 0 ; i < *np ; i ++ ) {
    fscanf(f, "%lg %lg", &((*rz)[2*i+0]), &((*rz)[2*i+1])) ;
    if ( src != NULL ) {
      /* for ( j = 0 ; j < str*(*ns) ; j ++ ) { */
      for ( j = 0 ; j < *ns ; j ++ ) {
	for ( k = 0 ; k <= *N ; k ++ ) {
	  fscanf(f, "%lg", &((*src)[i*(*ns)*str + j*str + 2*k + 0])) ;
	  fscanf(f, "%lg", &((*src)[i*(*ns)*str + j*str + 2*k + 1])) ;
	}
      }
    }
  }
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gdouble *rz1, *rz, *src, *fld, *work ;
  gint nsrc, nfld, N, Nf, ns, i, j, n, sstr, fstr ;
  gchar ch, *sfile, *ffile ;
  FILE *input, *output ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  output = stdout ;
  sfile = NULL ; ffile = NULL ;
  while ( (ch = getopt(argc, argv, "f:s:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'f': ffile = g_strdup(optarg) ; break ;
    case 's': sfile = g_strdup(optarg) ; break ;
    }
  }

  if ( sfile == NULL ) {
    fprintf(stderr, "%s: source file must be defined\n", progname) ;
    exit(1) ;    
  }
  input = fopen(sfile, "r") ;
  if ( input == NULL ) {
    fprintf(stderr, "%s: cannot open source file %s\n", progname, sfile) ;
    exit(1) ;
  }

  points_read(input, &nsrc, &N, &ns, &rz1, &src) ;

  fclose(input) ;
  
  fprintf(stderr, "%s: %d source points\n", progname, nsrc) ;
  fprintf(stderr, "%s: modal order %d\n", progname, N) ;
  fprintf(stderr, "%s: %d source components\n", progname, ns) ;

  if ( ffile == NULL ) {
    fprintf(stderr, "%s: field point file must be defined\n", progname) ;
    exit(1) ;    
  }
  input = fopen(ffile, "r") ;
  if ( input == NULL ) {
    fprintf(stderr, "%s: cannot open field point file %s\n", progname, ffile) ;
    exit(1) ;
  }

  points_read(input, &nfld, &Nf, &i, &rz, &fld) ;

  fclose(input) ;
  
  fprintf(stderr, "%s: %d field points\n", progname, nfld) ;
  /* fprintf(stderr, "%s: modal order %d\n", progname, Nf) ; */

  /* g_assert(Nf == N) ; */
  fld = (gdouble *)
    fftw_malloc(nfld*2*(N+2)*ns*sizeof(gdouble)) ;
  
  work = (gdouble *)g_malloc0(32*(N+1)*sizeof(gdouble)) ;

  sstr = 2*(N+2) ;
  fstr = 2*(N+2) ;

  fprintf(stderr, "%s: evaluating field; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;
  
  /* afmm_laplace_field_direct(rz1, src, sstr, ns, nsrc, N, */
  /* 				  rz, fld, fstr, nfld, TRUE, work) ; */
  for ( i = 0 ; i < nfld ; i ++ ) {
    afmm_laplace_field_direct_vec(rz1, 2, src, ns*sstr, sstr, ns, nsrc, N,
					&(rz[2*i]), &(fld[i*ns*fstr]), fstr, 1,
					TRUE, 32, work) ;
  }
  
  fprintf(stderr, "%s: field evaluated; %lg\n",
	  progname, g_timer_elapsed(timer, NULL)) ;

  /* fprintf(output, "%d %d %d\n", nfld, N, ns) ; */
  for ( i = 0 ; i < nfld ; i ++ ) {
    fprintf(output, "%lg %lg", rz[2*i+0], rz[2*i+1]) ;
    for ( j = 0 ; j < ns ; j ++ ) {
      for ( n = 0 ; n <= N ; n ++ ) {
	fprintf(output, " %1.16e %1.16e",
		fld[i*ns*fstr + j*fstr + 2*n+0],
		fld[i*ns*fstr + j*fstr + 2*n+1]) ;
      }
    }
    fprintf(output, "\n") ;
  }    
  
  return 0 ;
}

