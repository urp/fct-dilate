/***************************************************************************
 *            structure-tensor.c
 *
 *  Mon Aug 24 11:48:14 2009
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

#include "structure-tensor.h"

#include "convolution.h"
#include "tgmath.h"
#include "float.h"
#include "assert.h"
#include "pgmio.h"

//computes the eigenvalues and -vectors of a real symetric 2x2-matrix
void sym_eigenvalue_decomposition(float_pair*   const   lambda,
                                  float_pair*   const   ev1,
                                  float_pair*   const   ev2,
                                  float         const   xx,
                                  float         const   yy,
                                  float         const   xy
                                 )
{
  const float dxxyy = (xx-yy);

  const float d  = sqrtf(4.f*xy*xy + dxxyy*dxxyy );


  lambda->x  = .5f*(xx+yy + d);
  lambda->y  = .5f*(xx+yy - d);

  if(lambda->x == 0)
    return;

  if( xx>yy )
  { ev1->x     = d+dxxyy;
    ev1->y     = 2.f*xy;
  }else
  { ev1->x     = 2.f*xy;
    ev1->y     = d-dxxyy;
  }

  if(fabs(ev1->x)>FLT_EPSILON || fabs(ev1->y)>FLT_EPSILON) //debugging
  {
    const float len = sqrtf(  ev1->x*ev1->x + ev1->y*ev1->y );
    ev1->x/=len;
    ev1->y/=len;
  }

  ev2->x=-ev1->y;
  ev2->y= ev1->x;
};

void compute_structure_tensors(image_t  stxx,   // structure tensor elements J_11
                               image_t  styy,   // structure tensor elements J_22
                               image_t  stxy,   // structure tensor elements J_12
                               image_t  const   image,   // input image
                               float        const   sigma,   // averaging nbh of J
                               float        const   sigma_pre// presmoothing
                              )
{
   // make copy to working copy
  image_t cpy = clone_image(image);

  print_image_info("pre_presmoothed",cpy);
  write_pgm_image("pre_presmoothed.pgm",cpy);

  // presmooth input image to get rid of noise
  gaussian_convolution(cpy,sigma_pre);
  //gaussian_convolution_gauss_seidel(cpy,sigma_pre,100);

  print_image_info("presmoothed",cpy);
  write_pgm_image("presmoothed.pgm",cpy);


  //calculate structure tensors
  long i,j;

  const long nx=image.nx+image.boundary;
  const long ny=image.ny+image.boundary;

/*  #pragma omp parallel shared(cpy, \
                              structure_tensor_xx, \
                              structure_tensor_yy, \
                              structure_tensor_xy ) \
                       private(i,j)*/
  {
    //#pragma omp for schedule(static) nowait
    for(i=image.boundary;i<nx;++i)
      for(j=image.boundary;j<ny;++j)
      {
        float_pair du = { .5f*(cpy.data[i+1][j]-cpy.data[i-1][j]),
                          .5f*(cpy.data[i][j+1]-cpy.data[i][j-1]) };

        stxx.data[i][j] = du.x*du.x;
        styy.data[i][j] = du.y*du.y;
        stxy.data[i][j] = du.x*du.y;
      }
  }

  disalloc_image(&cpy);

  // smooth result to involve structure information from the neighborhood
  gaussian_convolution(stxx,sigma);
  gaussian_convolution(styy,sigma);
  gaussian_convolution(stxy,sigma);


  print_image_info("structxx",stxx);
  rescale_image(stxx,0.f,255.f);
  write_pgm_image("test/structxx.pgm",stxx);

  print_image_info("structyy",styy);
  rescale_image(styy,0.f,255.f);
  write_pgm_image("test/structyy.pgm",styy);

  print_image_info("structxy",stxy);
  rescale_image(stxy,0.f,255.f);
  write_pgm_image("test/structxy.pgm",stxy);

}

// decomposes the symmetric matrix of the structure tensor into its eigenvalues
// and orthonormal eigenvectors - the second eigenvector
void decompose_structure_tensors(image_t        lambda1, // larger  eigenvalue
                                 image_t        lambda2, // smaller eigenvalue
                                 image_t        e1x,     // 1st eigenvec(x-comp)
                                 image_t        e1y,     // 1st eigenvec(y-comp)
                                 image_t const  stxx,    // ST elements J_11
                                 image_t const  styy,    // ST elements J_22
                                 image_t const  stxy     // ST elements J_12
                                )
{
  assert(is_image_compatible(lambda1,lambda2));
  assert(is_image_compatible(lambda1,e1x));
  assert(is_image_compatible(lambda1,e1y));
  assert(is_image_compatible(lambda1,stxx));
  assert(is_image_compatible(lambda1,styy));
  assert(is_image_compatible(lambda1,stxy));

  //store eigenvalues and vectors
  long i,j;
  const long nx = stxx.nx+stxx.boundary;
  const long ny = stxx.ny+stxx.boundary;
  #pragma omp parallel shared(lambda1,lambda2,e1x,e1y) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for(i=stxx.boundary;i<nx;++i)
      for(j=stxx.boundary;j<ny;++j)
      {
        float_pair lambda,ev1,ev2;
        sym_eigenvalue_decomposition(&lambda,&ev1,&ev2,
                                     stxx.data[i][j],
                                     styy.data[i][j],
                                     stxy.data[i][j]  );
        lambda1.data[i][j]  = lambda.x;
        lambda2.data[i][j]  = lambda.y;
        e1x.data[i][j]      = ev1.x;
        e1y.data[i][j]      = ev1.y;
      }
  }
  //debugging
  image_t dblambda1 = clone_image(lambda1);
  print_image_info("lambda1",dblambda1);
  rescale_image(dblambda1,0.f,255.f);
  write_pgm_image("test/lambda1.pgm",dblambda1);
  disalloc_image(&dblambda1);

  image_t dblambda2 = clone_image(lambda2);
  print_image_info("lambda2",dblambda2);
  rescale_image(dblambda2,0.f,255.f);
  write_pgm_image("test/lambda2.pgm",dblambda2);
  disalloc_image(&dblambda2);
}

// computes the angle between the eigenvector direction and
// the positive x-axis for every image pixel.
void compute_rotation_angles(image_t        angles,  // output
                             image_t const  lambda1, // larger  eigenvalue
                             image_t const  lambda2, // smaller eigenvalue
                             image_t const  e1x,     // 1st eigenvector(x-comp.)
                             image_t const  e1y      // 1st eigenvector(y-comp.)
                            )
{
  assert(is_image_compatible(angles,lambda1));
  assert(is_image_compatible(angles,lambda2));
  assert(is_image_compatible(angles,e1x));
  assert(is_image_compatible(angles,e1y));

  long i,j;
  const long nx=angles.nx+angles.boundary;
  const long ny=angles.ny+angles.boundary;

  //calculate local rotation angles
  #pragma omp parallel shared(angles) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for(i=angles.boundary;i<nx;++i)
      for(j=angles.boundary;j<ny;++j)
      {
        angles.data[i][j] = lambda1.data[i][j]<FLT_EPSILON
                            ? 0.f
                    : atan2f(e1y.data[i][j],e1x.data[i][j]);
      }
  }

  //debugging
  image_t dbangles = clone_image(angles);
  print_image_info("angles",dbangles);
  rescale_image(dbangles,0.f,255.f);
  write_pgm_image("test/angles.pgm",dbangles);
  disalloc_image(&dbangles);
}
