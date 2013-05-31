/***************************************************************************
 *            os.h
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

#ifndef MORPH_OSHER_SETHIAN_H
#define MORPH_OSHER_SETHIAN_H

#include "common.h"
#include "image.h"
#include "args.h"

//calculates the update of pixel (i,j) using the method of osher&sethian
float	os_scheme(  image_t	image,
		    long	i,
		    long	j,
		    step_t	step);

void	osher_sethian_step( image_t image,
	                    step_t  step
			  );
#endif
