/***************************************************************************
 *            fct.c
 *
 *  Thu Jul  9 12:31:57 2009
 *  Copyright  2009  Peter Urban
 *  <s9peurba@stud.uni-saarland.de>
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
 
 
#include "fct.h"
#include "rt.h"

#include <assert.h>
#include <stdio.h>
#include <tgmath.h>

//calculates the update of pixel (i,j) using the flux correxted transport scheme
float	fct_scheme( image_t  image,
                    long     i,
                    long     j,
	                step_t   step)
{ 
  const float	hstepsize = .5f*step.stepsize;

  //get pixel values from rotated grid
  const float x0y0 = image.data[i][j];
  const float p1y0 = step.rot(image,i,j, 1, 0,&step);
  const float m1y0 = step.rot(image,i,j,-1, 0,&step);
  const float x0p1 = step.rot(image,i,j, 0, 1,&step);
  const float x0m1 = step.rot(image,i,j, 0,-1,&step);

  //first order derivatives
  const float uxf = p1y0 - x0y0;
  const float uxb = x0y0 - m1y0;
  const float uyf = x0p1 - x0y0; 
  const float uyb = x0y0 - x0m1;

  //early termination when pixel is local maximum
  if( (uxf>=0.f && uyf>=0.f) && (uxb<=0.f && uyb<=0.f)  )
    return 0.f;//scheme.norm(uhxc,uhyc,scheme.norm_info);

  //get pixel values from rotated grid
  const float p2y0 = step.rot(image,i,j, 2, 0,&step);
  const float m2y0 = step.rot(image,i,j,-2, 0,&step);
  const float x0p2 = step.rot(image,i,j, 0, 2,&step);
  const float x0m2 = step.rot(image,i,j, 0,-2,&step);

  //shifted first order derivatives
  const float uxff= p2y0 - p1y0;
  const float uxbb= m1y0 - m2y0;
  const float uyff= x0p2 - x0p1;
  const float uybb= x0m1 - x0m2;

  //backward diffusion flows
  const float gxp = minmod3( uxb ,hstepsize * uxf,uxff );
  const float gxm = minmod3( uxbb,hstepsize * uxb,uxf  );
  const float gyp = minmod3( uyb ,hstepsize * uyf,uyff );
  const float gym = minmod3( uybb,hstepsize * uyb,uyf  );

  //higher order flow component
  const float_pair Ch = { .x = hstepsize*(p1y0-m1y0),
						  .y = hstepsize*(x0p1-x0m1) };

  //diffusive flow component
  const float_pair Cd = { .x = fabs(Ch.x) + gxp - gxm,
						  .y = fabs(Ch.y) + gyp - gym  };
  
  //applying the norm representing the structuring element 
  //and correct by substracting the diffusive "flow" component
  return  	step.norm(Ch.x,Ch.y,i,j,step.norm_info) 
		  - step.norm(Cd.x,Cd.y,i,j,step.norm_info);
}

//euler step with timestep size = 1 and without calling initialize_step
void	add_update_step( image_t image, // image to manipulate
                         step_t	 step   // step specification
                       )
{ 
  reflecting_boundaries(image);
  image_t image_old = clone_image(image);
  
  long i,j;
  #pragma omp parallel shared(image,step) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for(i=image.boundary;i<(image.nx+image.boundary);++i)
      for(j=image.boundary;j<(image.ny+image.boundary);++j)
	image.data[i][j] = image.data[i][j] + step.scheme(image_old,i,j,step);
  }
  disalloc_image(&image_old);
}

void	flux_corrected_step( image_t image,
                             step_t	 step
                           )
{ 
  //predictor step
  rouy_tourin_step (image,step);

  //prepare step information
  step.scheme = &fct_scheme;

  //corrector step
  add_update_step(image,step);
}

