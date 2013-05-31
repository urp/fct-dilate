/***************************************************************************
 *            image.c
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

#include <stdio.h>
#include <tgmath.h>
#include <assert.h>

#include "alloc.h"
#include "image.h"
#include "float.h"


// allocates storage for image
void alloc_image(image_t* image)
{
  const long bound = 2*image->boundary;
  alloc_matrix( &image->data,image->nx+bound,image->ny+bound);
}

// returns 1 if image dimensions and boundary size agree
// for the two given images and both have storage
// otherwise 0 is returned
int is_image_compatible(image_t const a,
                        image_t const b )
{
  return    a.data!=(void*)0
         && b.data!=(void*)0
         && a.nx==b.nx
         && a.ny==b.ny
         && a.boundary==b.boundary;
}

// returns an image with the same dimension and boundary as the given reference
// image
image_t create_compatible_image(image_t const   ref)
{ image_t out=ref;
  alloc_image(&out);
  return out;
}

// reallocate storage if unallocated or if at least one image dimension is
// not tthe same size as in the reference image
void ensure_image_compatibility(image_t*        image,  // output
                                image_t const   ref     // reference image
                               )
{
  // allocate storage for eigenalues and eigenvectors
  if( !is_image_compatible(*image,ref) )
  { if(image->data != 0)
      disalloc_image(image);
    *image=ref;
    alloc_image(image);
  }
}

// disallocates storage for image
void disalloc_image(image_t* image)
{
  const long bound = 2*image->boundary;
  disalloc_matrix(image->data,image->nx+bound,image->ny+bound);
  image->data = (void*)0;
}

// allocates an image and copies content from 'orig'
// to the new image having its own storage
image_t clone_image(image_t const orig)
{
  image_t out=orig;
  alloc_image(&out);
  const long bound = 2*orig.boundary;
  copy_matrix(orig.data,out.data,orig.nx+bound,orig.ny+bound);
  return out;
}

void copy_image_data(image_t const from,image_t to)
{
  assert(from.nx==to.nx);
  assert(from.ny==to.ny);
  assert(from.boundary==to.boundary);

  const long bound = 2*to.boundary;
  copy_matrix(from.data,to.data,to.nx+bound,to.ny+bound);
}

float_pair  image_min_max(image_t const image)
{
  float_pair minmax = { FLT_MAX,FLT_MIN };

  long i,j;
  const long nx = image.nx+image.boundary;
  const long ny = image.ny+image.boundary;

  for (j=image.boundary; j<ny; j++)
    for (i=image.boundary; i<nx; i++)
    { minmax.x = min(minmax.x,image.data[i][j]);
      minmax.y = max(minmax.y,image.data[i][j]);
    }
  return minmax;
}

void print_image_info(char* title,image_t const image)
{
  float_pair mm = image_min_max(image);
  printf("print_image_info\t|image %s:\tsize (%ld,%ld)\tboundary %ld\tmin %g\tmax %g\n",
         title,image.nx,image.ny,image.boundary,mm.x,mm.y);
}

void rescale_image(image_t image,float nmin,float nmax)
{
  assert(nmin<nmax);
  float_pair minmax=image_min_max(image);
  long i,j;
  const long nx = image.nx+image.boundary;
  const long ny = image.ny+image.boundary;

  float shift = nmin-minmax.x;
  float scale = (nmax-nmin)/(minmax.y-minmax.x);

  for (j=image.boundary; j<ny; j++)
    for (i=image.boundary; i<nx; i++)
    { image.data[i][j] += shift;
      image.data[i][j] *= scale;
    }
}

void reflecting_boundaries(image_t image)
{
  long    const bd   = image.boundary;
  long    const nx   = image.nx+bd;
  long    const ny   = image.ny+bd;
  float** const u    = image.data;

  long    i,j;
  #pragma omp parallel private(i,j)
  {
    // top & bottom
    #pragma omp for schedule(static) nowait
    for(i = bd; i < nx; i++)
    {
      for(j = 0; j < bd; j++)
      { u[i][bd-j-1] = u[i][bd+j];
        u[i][ny+j]   = u[i][ny-j-1];
      }
    }

    // left & right
    #pragma omp for schedule(static) nowait
    for(j = bd; j < ny; j++)
    {
      for(i = 0; i < bd; i++)
      { u[bd-i-1][j] = u[bd+i][j];
        u[nx+i]  [j] = u[nx-i-1][j];
      }
    }

    #pragma omp for schedule(static) nowait
    for(i = 0; i < bd; i++)
      for(j = 0; j < bd; j++)
      {
        //top left
        u[bd-i-1][bd-j-1]=u[bd+i][bd+j];
        //top right
        u[nx+i]  [bd-j-1]=u[nx-i-1][bd+j];
        //bottom left
        u[bd-i-1][ny+j]  =u[bd+i][ny-j-1];
        //bottom right
        u[nx+i]  [ny+j]  =u[nx-i-1][ny-j-1];
      }
  }
  return;
}

float   bilinear_interpolation(image_t image,float_pair r)
{
  float**    const   u  = image.data;
  long       const   i0 = floorf(r.x);
  long       const   j0 = floorf(r.y);
  float_pair         dr = { .x = r.x-i0,
                            .y = r.y-j0 };

  return    u[i0]  [j0]  *(1.f-dr.x)*(1.f-dr.y)
          + u[i0+1][j0]  *(    dr.x)*(1.f-dr.y)
          + u[i0]  [j0+1]*(1.f-dr.x)*(    dr.y)
          + u[i0+1][j0+1]*(    dr.x)*(    dr.y);
}
