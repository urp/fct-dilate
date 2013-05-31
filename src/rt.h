/***************************************************************************
 *            rt.h
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

#ifndef	DILATE_ROUY_TOURIN_H
#define DILATE_ROUY_TOURIN_H

#include "common.h"
#include "args.h"
#include "image.h"

//calculates the update of pixel (i,j) using the method of rouy&tourin
float rt_scheme( image_t image,
                 long    i,
                 long    j,
                 step_t  step
               );




void rouy_tourin_step( image_t image,
                       step_t  step
                     );

#endif

