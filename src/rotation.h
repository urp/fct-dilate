/***************************************************************************
 *            rotation.h
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

#ifndef DILATE_ROTATION_H
#define DILATE_ROTATION_H

#include "common.h"
#include "image.h"
#include "args.h"

float   no_rotation( image_t image,    // input image
                     long    center_i, // x-coord of structuring element center
                     long    center_j, // y-coord of structuring element center
                     long    i,        // unrotated x-coord of returned sample
                     long    j,        // unrotated y-coord of returned sample
                     step_t* step      // (not used) timestep-related data
                   );

//rotate the grid about the angle stored at rot_angle
float fixed_rotation( image_t image,	// input image
                      long    center_i, // x-coord of the SE center  
                      long    center_j, // y-coord of the SE center
                      long    i,		// unrotated x-coord of returned sample
                      long    j,        // unrotated y-coord of returned sample
                      step_t* step      // timestep-related data
                    );

float auto_rotation( image_t image,     // input image
                     long    center_i,  // x-coord of the SE center  
                     long    center_j,  // y-coord of the SE center
                     long    i,	        // unrotated x-coord of returned sample
                     long    j,         // unrotated y-coord of returned sample
                     step_t* step       // timestep-related data
                   );
void auto_rotation_prepare( image_t	const	image,
                            step_t* const	step );

void adaptive_se_program_prepare( image_t	const	image,
										   step_t*  const	step );

// align structure element along image edges
float adaptive_se_program_rotation( image_t image,
									long    center_i,
									long    center_j,
									long    i,
									long    j,
									step_t* step
								   );


#endif
