/***************************************************************************
 *            convolution.c
 *
 *  Tue Jul 21 11:00:39 2009
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

#include "convolution.h"
#include "alloc.h"

#include <tgmath.h>
#include <assert.h>

void gaussian_convolution( image_t image, float const sigma )
{
  gaussian_convolution_explicit( image, sigma, 0.05 );
}

// semi-implicit gauss-seidel solver using the diffusion-reaction equation
void gaussian_convolution_gauss_seidel( image_t             image,
                                        float        const  sigma,
                                        unsigned int const  iterations )
{
  float alpha = sigma*sigma/2.;

  image_t old = clone_image(image);

  image_t B_lower = create_compatible_image(image);
  image_t B_upper = create_compatible_image(image);
  image_t B_left  = create_compatible_image(image);
  image_t B_right = create_compatible_image(image);
  image_t B_center = create_compatible_image(image);

  unsigned int t;
  for( t=0; t<iterations; ++t )
  {
    unsigned int i,j;

    const unsigned int nx=image.nx+image.boundary;
    const unsigned int ny=image.ny+image.boundary;

    reflecting_boundaries(image);

    // diffusivities

    for(i=image.boundary;i<nx;++i)
    {
      for(j=image.boundary;j<ny;++j)
      {
         // tikhonov regularizer
         B_lower.data[i][j] = - alpha * 2 * fabs( image.data[i][j] - image.data[i-1][j] );
         B_upper.data[i][j] = - alpha * 2 * fabs( image.data[i][j] - image.data[i+1][j] );
         B_left.data [i][j] = - alpha * 2 * fabs( image.data[i][j] - image.data[i][j-1] );
         B_right.data[i][j] = - alpha * 2 * fabs( image.data[i][j] - image.data[i][j+1] );
         B_center.data[i][j] = 1. - ( B_lower.data[i][j] + B_upper.data[i][j] + B_left.data[i][j] + B_right.data[i][j] );

      }
    }

    // Gauss-Seidel solver

    //#pragma omp parallel shared(image) private(i,j)
    {
      //#pragma omp for schedule(static) nowait
      for(i=image.boundary;i<nx;++i)
      {
        for(j=image.boundary;j<ny;++j)
        {
          //reflecting_boundaries(image);
          image.data[i][j] = 1./B_center.data[i][j] * ( old.data[i][j]
                                                 - B_lower.data[i][j] * image.data[i-1][j]
                                                 - B_upper.data[i][j] * image.data[i+1][j]
                                                 - B_left .data[i][j] * image.data[i][j-1]
                                                 - B_right.data[i][j] * image.data[i][j+1]
                                                 ) ;
        }
      }

    }

  }
  disalloc_image(&old);

}


void gaussian_convolution_explicit( image_t         image,
                                    float   const   sigma,
                                    float   const   timestep )
{
  assert(image.boundary>=1);

  const float           time    = sigma*sigma/2;
  const unsigned int    n       = ceilf(time/timestep);
  const float           stepsize = time/n;
  image_t  old   = image;
  alloc_image(&old);

  unsigned int i,j,t;
  for(t=0;t<n;++t)
  {
    reflecting_boundaries(image);
    copy_image_data(image,old);

    const long nx=image.nx+image.boundary;
    const long ny=image.ny+image.boundary;

    #pragma omp parallel shared(image,old) private(i,j)
    {
      #pragma omp for schedule(static) nowait
      for(i=image.boundary;i<nx;++i)
        for(j=image.boundary;j<ny;++j)
        {
          image.data[i][j] += (- 4.f * old.data[i][j]
                               +       old.data[i+1][j]
                               +       old.data[i-1][j]
                               +       old.data[i][j+1]
                               +       old.data[i][j-1] ) * stepsize;
        }
    }
  }
  disalloc_image(&old);
}
