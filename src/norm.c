/***************************************************************************
 *            norm.c
 *
 *  Mon Jul 20 18:12:45 2009
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

#include "norm.h"
#include "args.h"

#include <tgmath.h>
#include <assert.h>

float p_norm( float const	x,	// vector components
              float const	y,
			  long  const	i,	// current pixel coords
			  long  const	j,
              void* const	info// norm specific parameters
            )
{ 
  assert(info);
  float p = ((p_norm_storage_t*)info)->p;
  assert(p>0.f);
  //printf("p_norm: p %g\n",p);
  return  powf( powf(fabs(x),p)+powf(fabs(y),p) ,1.f/p);
}


float 	circle_norm( float const	x,		// vector components
                     float const	y,
	                 long  const	i,		// current pixel coords
	                 long  const	j,
                     void* const	info	// norm specific parameters
	               )
{ return sqrtf(x*x+y*y); }


float 	diamond_norm( float const	x,	// vector components
                      float const	y,
	                  long  const	i,	// current pixel coords
	                  long  const	j,
                      void* const	info// norm specific parameters
	                )
{ return max(fabs(x),fabs(y)); }


float elliptic_norm(float const	x,		// vector components
					float const	y,
					long  const	i,		// current pixel coords
					long  const	j,
					void* const	info	// norm specific parameters
				   )
{ assert(info);
  float_pair axes = ((elliptic_norm_storage_t*)info)->main_axes;
  return sqrtf(axes.x*axes.x*x*x + axes.y*axes.y*y*y);
} 


float adaptive_se_program_norm( float const	x,		// vector components
								float const	y,
								long  const	i,		// current pixel coords
								long  const	j,
								void* const	info	// norm specific parameters
							  )
{ assert(info);
  adaptive_se_program_storage_t* params = (adaptive_se_program_storage_t*)info;
  float_pair axes = { .x = 1.f-expf(-10.f*params->lambda1.data[i][j]/params->lambda1_max),
					  .y = 1.f-expf(-10.f*params->lambda2.data[i][j]/params->lambda2_max)
					};
  return sqrtf(axes.x*axes.x*x*x + axes.y*axes.y*y*y);
} 
