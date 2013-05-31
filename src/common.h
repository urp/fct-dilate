/***************************************************************************
 *            common.h
 *
 *  Wed Jul  8 16:25:53 2009
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

#ifndef MORPH_COMMON_H
#define MORPH_COMMON_H 

#include <stdlib.h>

float	max(float a,float b);

float	min(float a,float b);

float   minmod(float a,float b);

float   minmod3(float a,float b,float c);

typedef struct
{ 
  float x;
  float y; 
} float_pair;

#endif
