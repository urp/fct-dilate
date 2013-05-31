/***************************************************************************
 *            norm.h
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

#ifndef DILATE_NORM_H
#define DILATE_NORM_H

float p_norm( float const	x,	// vector components
              float const	y,
			  long  const	i,	// current pixel coords
			  long  const	j,
              void* const	info// norm specific parameters
            );

float 	circle_norm( float const	x, 	// vector components
                     float const	y,
	                 long  const	i, 	// current pixel coords
	                 long  const	j,
                     void* const	info// norm specific parameters
	               );

float 	diamond_norm( float const	x,	// vector components
                      float const	y,
	                  long  const	i,	// current pixel coords
	                  long  const	j,
                      void* const	info// norm specific parameters
	                );

float	elliptic_norm(float const	x,		// vector components
                      float const	y,
	                  long  const	i,		// current pixel coords
	                  long  const	j,
                      void* const	info	// norm specific parameters
	                 );

float adaptive_se_program_norm( float const	x,		// vector components
								float const	y,
								long  const	i,		// current pixel coords
								long  const	j,
								void* const	info	// norm specific parameters
							  );

#endif
