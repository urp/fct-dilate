/***************************************************************************
 *            integration.h
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

#ifndef DILATE_INTEGRATION_H
#define DILATE_INTEGRATION_H

#include "common.h"
#include "image.h"
#include "args.h"

 // simple integrate of first order
void	euler_step( image_t	image,	    // image to manipulate 
                    step_t	step        // step specification
                  );

// 2-step second-order integratorge-kutta method
void	heun_step ( image_t image,      // image to manipulate
                    step_t	step        // step specification
                  );

#endif
