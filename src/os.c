/***************************************************************************
 *            os.c
 *
 *  Wed Jul  8 16:39:45 2009
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

#include	"os.h"

#include        <math.h>
#include	<stdio.h>

//calculates the update of pixel (i,j) using the method of osher&sethian
float	os_scheme(  image_t image,
                    long    i,
                    long    j,
                    step_t  step)
{
  //get pixel values from rotated grid
  float x0y0 = image.data[i][j];
  float p1y0 = step.rot(image,i,j, 1, 0,&step);
  float m1y0 = step.rot(image,i,j,-1, 0,&step);
  float x0p1 = step.rot(image,i,j, 0, 1,&step);
  float x0m1 = step.rot(image,i,j, 0,-1,&step);

  //first order derivatives
  float xf = p1y0 - x0y0;
  float xb = x0y0 - m1y0;
  float yf = x0p1 - x0y0;
  float yb = x0y0 - x0m1;

  //early termination when pixel is local maximum
  if( (xf<=0.f && yf<=0.f) && (xb>=0.f && yb>=0.f) )
    return 0.f;

  //get pixel values from rotated grid
  float p2y0 = step.rot(image,i,j, 2, 0,&step);
  float m2y0 = step.rot(image,i,j,-2, 0,&step);
  float x0p2 = step.rot(image,i,j, 0, 2,&step);
  float x0m2 = step.rot(image,i,j, 0,-2,&step);

  //second order derivatives
  float xff= p2y0 - p1y0 - xf;
  float xcc= p1y0 - x0y0 - xb;
  float xbb= xb   - m1y0 + m2y0;
  float yff= x0p2 - x0p1 - yf;
  float ycc= x0p1 - x0y0 - yb;
  float ybb= yb   - x0m1 + x0m2;

  //add "inflow" from both sides
  float_pair du = { .x =        max(xf - .5f*minmod(xff,xcc),0.f) 
                        + fabs( min(xb + .5f*minmod(xcc,xbb),0.f) ),
                    .y =        max(yf - .5f*minmod(yff,ycc),0.f)
                        + fabs( min(yb + .5f*minmod(ycc,ybb),0.f) ) };

  //applying the norm representing the structuring element
  return step.norm(du.x,du.y,i,j,step.norm_info);
}

void osher_sethian_step( image_t image,
                         step_t	 step	
                       )
{
  //prepare step information
  step.scheme = &os_scheme;
	
  //integrate step
  heun_step(image,step);
}
