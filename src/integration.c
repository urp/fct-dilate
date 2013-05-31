/***************************************************************************
 *            integration.c
 *
 *  Mon Aug 24 14:17:25 2009
 *  Copyright  2009  urp
 *  <urp@<host>>
 ****************************************************************************/

/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor Boston, MA 02110-1301,  USA
 */

#include "integration.h"

#include "alloc.h"
#include "stdio.h"

void	euler_step( image_t image,
                    step_t  step
	          )
{ 
  initialize_step(image,&step);
	
  reflecting_boundaries(image);
  image_t image_old = clone_image(image);
  
  long i,j;

  //parallel block
  #pragma omp parallel shared(image,step) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for(i=image.boundary;i<(image.nx+image.boundary);++i)
      for(j=image.boundary;j<(image.ny+image.boundary);++j)
      { 
        float ut = step.scheme(image_old,i,j,step);
        image.data[i][j] = image.data[i][j] + step.stepsize * ut;
      }
  }
  disalloc_image(&image_old);
}


void	heun_step ( image_t 	image,
                    step_t	step
	              )
{ 
  initialize_step(image,&step);
	
  reflecting_boundaries(image);
  image_t image_half = image;
  alloc_image(&image_half);
  
  long i,j;
  const long nx=image.nx+image.boundary;
  const long ny=image.ny+image.boundary;

  //parallel block
  #pragma omp parallel shared(image,step) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for(i=image.boundary;i<nx;++i)
      for(j=image.boundary;j<ny;++j)
      {
        float ut = step.scheme(image,i,j,step);
        image_half.data[i][j] = image.data[i][j] + step.stepsize * ut;
      }
  }
	
  reflecting_boundaries(image_half);
  initialize_step(image_half,&step);
	
  //parallel block
  #pragma omp parallel shared(image,image_half,step) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for(i=image.boundary;i<nx;++i)
      for(j=image.boundary;j<ny;++j)
      {
        float ut = step.scheme(image_half,i,j,step);
        image.data[i][j] = .5f*(   image.data[i][j] 
                                 + image_half.data[i][j] 
                                 + step.stepsize * ut    );
      }
  }
  disalloc_image(&image_half);
}

