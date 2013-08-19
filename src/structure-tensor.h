/***************************************************************************
 *            structure-tensor.h
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

#ifndef DILATE_STRUCTURE_TENSOR_H
#define DILATE_STRUCTURE_TENSOR_H

#include "common.h"
#include "image.h"

//computes the eigenvalues and -vectors of a real symetric 2x2-matrix
void sym_eigenvalue_decomposition(float_pair*   const   lambda,
                                  float_pair*   const   ev1,
                                  float_pair*   const   ev2,
                                  float         const   xx,
                                  float         const   yy,
                                  float         const   xy
                                 );

void compute_structure_tensors(image_t  stxx,   // structure tensor elements J11
                               image_t  styy,   // structure tensor elements J22
                               image_t  stxy,   // structure tensor elements J12
                               image_t  const   image,   // input image
                               float    const   sigma,   // averaging nb of J
                               float    const   sigma_pre// presmoothing
                              );

// decomposes the symmetric matrix of the structure tensor into its eigenvalues
// and orthonormal eigenvectors - the second eigenvector is (-e1y,e1x)
void decompose_structure_tensors(image_t        lambda1, // larger  eigenvalue
                                 image_t        lambda2, // smaller eigenvalue
                                 image_t        e1x,     // 1st eigenvec(x-comp)
                                 image_t        e1y,     // 1st eigenvec(y-comp)
                                 image_t const  stxx,    // ST elements J_11
                                 image_t const  styy,    // ST elements J_22
                                 image_t const  stxy     // ST elements J_12
                                );

// computes the angle between the eigenvector direction and
// the positive x-axis for every image pixel.
void compute_rotation_angles(image_t        angles,  // output
                             image_t const  lambda1, // larger  eigenvalue
                             image_t const  lambda2, // smaller eigenvalue
                             image_t const  e1x,     // 1st eigenvector(x-comp.)
                             image_t const  e1y      // 1st eigenvector(y-comp.)
                            );

#endif
