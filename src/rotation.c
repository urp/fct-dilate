/***************************************************************************
 *            rotation.c
 *
 *  Mon Aug 24 14:06:31 2009
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
#include "rotation.h"
#include "structure-tensor.h"

#include "alloc.h"
#include "pgmio.h"
#include "tgmath.h"
#include "float.h"
#include "assert.h"

float   no_rotation(image_t image,
                    long    center_i,
                    long    center_j,
                    long    i,
                    long    j,
                    step_t* step
                   )
{ 
  return image.data[center_i+i][center_j+j]; 
}

float fixed_rotation( image_t image,
                      long    center_i,
                      long    center_j,
                      long    i,
                      long    j,
			          step_t* step
		            )
{ 
  assert(step->rot_info);
  float angle =   -1.f
	            * ((fixed_rotation_storage_t*)(step->rot_info))->angle;
  float sina  = sinf(angle);
  float cosa  = //cos(angle);
                sqrtf(1.f-sina*sina);
  float_pair coords = { .x = center_i + i*cosa - j*sina,
                        .y = center_j + i*sina + j*cosa };
  return bilinear_interpolation(image,coords);
}


void auto_rotation_prepare( image_t	const	image,
                            step_t* const	step )
{  
  auto_rotation_storage_t* params = (auto_rotation_storage_t*) step->rot_info;

  // allocate storage for structure tensor entries
  image_t stxx = create_compatible_image(image);
  image_t styy = create_compatible_image(image);
  image_t stxy = create_compatible_image(image);  

  // compute structure tensor entries
  compute_structure_tensors(stxx,styy,stxy,image,
							params->sigma,params->sigma_pre);

  // get eigenvalue decomposition
  image_t lambda1 = create_compatible_image(image);
  //not required 
  image_t lambda2 = create_compatible_image(image);  
  image_t e1x 	  = create_compatible_image(image);
  image_t e1y 	  = create_compatible_image(image);  
  decompose_structure_tensors(lambda1,lambda2,e1x,e1y,stxx,styy,stxy);
  
  // free structure tensor storage
  disalloc_image(&stxx);
  disalloc_image(&styy);
  disalloc_image(&stxy);

  // save new roatation angles
  ensure_image_compatibility(&params->angles,image);   
  compute_rotation_angles(params->angles,
                          lambda1,
                          lambda2,
                          e1x,
                          e1y);

  disalloc_image(&lambda1);
  disalloc_image(&lambda2);
  disalloc_image(&e1x);
  disalloc_image(&e1y);
}

float auto_rotation( image_t image,
                     long    center_i,
                     long    center_j,
                     long    i,
                     long    j,
                     step_t* step
		           )
{
  auto_rotation_storage_t* params = (auto_rotation_storage_t*) step->rot_info;

  float angle = -params->angles.data[center_i][center_j];
  
  float sina    = sinf(angle);
  float cosa    = cosf(angle);
  
  float_pair coords = { .x = center_i + i*cosa - j*sina,
                        .y = center_j + i*sina + j*cosa };
  return bilinear_interpolation(image,coords);
}
 
void adaptive_se_program_prepare( image_t	const	image,
								  step_t*  const	step )
{  
  adaptive_se_program_storage_t* params;
  params = (adaptive_se_program_storage_t*) step->rot_info;
  
  // allocate storage for structure tensor entries
  image_t stxx = create_compatible_image(image);
  image_t styy = create_compatible_image(image);
  image_t stxy = create_compatible_image(image);  

  // compute structure tensor entries
  compute_structure_tensors(stxx,styy,stxy,image,
							params->sigma,params->sigma_pre);

  // get eigenvalue decomposition
  image_t e1x = create_compatible_image(image);
  image_t e1y = create_compatible_image(image);  
  ensure_image_compatibility(&params->lambda1,image);   
  ensure_image_compatibility(&params->lambda2,image);   
  decompose_structure_tensors(params->lambda1,params->lambda2,e1x,e1y,
                              stxx,styy,stxy);
  // free structure tensor storage
  disalloc_image(&stxx);
  disalloc_image(&styy);
  disalloc_image(&stxy);

  // save new roatation angles
  ensure_image_compatibility(&params->angles,image);   
  compute_rotation_angles(params->angles,
                          params->lambda1,
                          params->lambda2,
                          e1x,
                          e1y
                         );
  
  // save maximal eigenvalues
  const float_pair l1minmax = image_min_max(params->lambda1);
  const float_pair l2minmax = image_min_max(params->lambda2);
  params->lambda1_max = l1minmax.y;
  params->lambda2_max = l2minmax.y;

  disalloc_image(&e1x);
  disalloc_image(&e1y);
}

// align structure element along image edges 
float adaptive_se_program_rotation( image_t image,
									long    center_i,
									long    center_j,
									long    i,
									long    j,
									step_t* step
								   )
{
  adaptive_se_program_storage_t* params;
  params = (adaptive_se_program_storage_t*) step->rot_info;

  float angle = -params->angles.data[center_i][center_j]+.5*M_PI;
  
  float sina    = sinf(angle);
  float cosa    = cosf(angle);
  
  float_pair coords = { .x = center_i + i*cosa - j*sina,
                        .y = center_j + i*sina + j*cosa };
  return bilinear_interpolation(image,coords);
}
