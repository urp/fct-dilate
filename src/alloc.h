/***************************************************************************
 *            alloc.h
 *
 *  Tue Jul  7 21:30:48 2009
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
 
#ifndef	MORPH_ALLOC_H
#define MORPH_ALLOC_H

#include "common.h"

// allocates storage for a vector of size n 
void alloc_vector(float **vector,long n);

// allocates storage for matrix of size nx * ny 
void alloc_matrix(float*** matrix,long nx,long ny);

// disallocates storage for a vector of size n
void disalloc_vector(float *vector);

// disallocates storage for matrix of size nx * ny
void disalloc_matrix(float **matrix,long nx,long ny);

// copies content from in to out 
void copy_matrix(float **in,float **out,long nx,long ny);

#endif
