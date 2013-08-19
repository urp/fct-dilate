/***************************************************************************
 *            rt.c
 *
 *  Tue Jul  7 15:08:21 2009
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

#include <stdio.h>
#include <math.h>

#include "alloc.h"
#include "common.h"
#include "rt.h"
#include "integration.h"

float   rt_scheme(  image_t image,
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

  //simple upwind scheme
  float_pair du = { .x = max( max(p1y0-x0y0,0.f), max(m1y0-x0y0,0.f) ),
                    .y = max( max(x0p1-x0y0,0.f), max(x0m1-x0y0,0.f) ) };

  //applying the norm representing the structuring element
  return step.norm(du.x,du.y,i,j,step.norm_info);
}

void    rouy_tourin_step( image_t   image,
                      step_t    step
            )
{
  //prepare step information
  step.scheme = &rt_scheme;

  euler_step(image,step);
}
