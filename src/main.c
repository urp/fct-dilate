/*
 * main.c
 * Copyright (C) Peter Urban 2009 <s9peurba@stud.uni-saarland.de>
 *
 * main.c is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * main.c is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

# include <stdio.h>
# include <time.h>
# include <omp.h>

# include "common.h"
# include "image.h"
# include "args.h"
# include "pgmio.h"

int     main( int argc, char** argv )
{
  printf( "starting dilation...\n" );
  printf( "using %d processes\n", omp_get_num_procs() );

  //initialize program parameters
  args_t args;

  //permanent storage of values referenced inside the args variable
  param_storage_t storage = initialize_storage();

  parse_args( &args, &storage, argc, argv );

  image_t image = load_pgm_image( args.inpath, args.boundary );

  if( ! image.data )
    return 1;

  printf( "main: starting %ld iterations\n", args.stepcount );
  time_t tick = time(0);

  for( args.step.t = 0; args.step.t < args.stepcount; ++args.step.t )
  {
    printf( "\rstep %ld", args.step.t+1 );

    // integrate
    args.step.integrate( image, args.step );

    //write output file
    if( args.write_step )
      write_pgm_image( args.outpath, image );

    fflush(stdout);
  }

  time_t tack = time(0);
  printf( "\rcompleted in %g seconds\n", difftime( tack, tick ) );

  write_pgm_image( args.outpath, image );
  disalloc_image( &image );

  return 0;
}
