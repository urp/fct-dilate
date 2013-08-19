/***************************************************************************
 *            image.h
 *
 *  Tue Aug 25 22:17:29 2009
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

#ifndef DILATE_IMAGE_H
#define DILATE_IMAGE_H

#include "common.h"

typedef struct
{
  float** data;     // the array containing the data
  long    nx;       // the x dimension of the image (without boundary)
  long    ny;       // the y dimension of the image (without boundary)
  long    boundary; // the width of the boundary surrounding the image data
} image_t;

// allocates storage for image
void alloc_image(image_t* image);


// returns an image with the same dimension and boundary as the given reference
// image
image_t create_compatible_image(image_t const ref);

// returns 1 if image dimensions and boundary size agree
// for the two given images and both have storage
// otherwise 0 is returned
int is_image_compatible(image_t const a,
                        image_t const b );

// reallocate storage if unallocated or if at least one image dimension is
// not the same size as in the reference image
void ensure_image_compatibility(image_t*        image,  // output
                                image_t const   ref     // reference image
                               );

// disallocates storage for image
void disalloc_image(image_t* image);

// allocates an image and copies content from 'source' to the returned image
image_t clone_image(image_t const source);

// copys the contents of image.data from one image to another
void copy_image_data(image_t const  from,
                     image_t        to
                    );

float_pair  image_min_max(image_t const image);

void print_image_info(char* title,image_t const image);

void rescale_image(image_t image,float nmin,float nmax);

// fills the boundary of the given image with values
// of the pixels on the other side of the boundary
void    reflecting_boundaries(image_t);

// use samples from image to interpolate value at location r
float   bilinear_interpolation(image_t image,float_pair r);


#endif
