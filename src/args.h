/***************************************************************************
 *            initargs.h
 *
 *  Thu Jul  9 23:53:27 2009
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

#ifndef DILATE_ARGS_H
#define DILATE_ARGS_H

#include "common.h"
#include "image.h"

//encapsulates all information neccessary to do an integration step
struct step_t
{
          // ???remove??? the number of the current timestep 
  long    t;
          // the temporal stepsize
  float   stepsize;
          // integrators the integration method
  void  (*integrate)(image_t,struct step_t);
          // the scheme used by the integrate 
          // !!! this value is automatically assigned by the integrate function !!!
  float	(*scheme)(image_t,long,long,struct step_t);
          // scheme specific information
          // !!! this value is automatically assigned by the integrate function !!!
  void*   scheme_info;
		  // prepares data shared by other plugins which work together in a program
  void	(*program_prepare)(image_t const,struct step_t* const);
          // grid rotation function used to handle rotated structure elements
  float (*rot)(image_t,long,long,long,long,struct step_t*);
		  // prepares the data required in the next step
  void	(*rot_prepare)(image_t const,struct step_t* const);
          // rotation-function specific information
  void*	  rot_info;

		  // the norm determines the shape of the structuring element
  float (*norm)(float,float,long,long,void*);
		  // prepares the data required in the next step
  void	(*norm_prepare)(image_t,struct step_t*);
          // norm specific information
  void*   norm_info;
};

typedef struct step_t step_t;

//the record containing all information nesseccary to dilate an image
typedef struct
{
  //path of the input image
  char*	 inpath;
  //path of the output image
  char*	 outpath;
  //the boundary used to implement the boundary conditions
  long	 boundary;
  //the number of timesteps to perform
  long	 stepcount;
  //write output after every step
  int	 write_step;
  //all information regarding the steps
  step_t step;
} args_t;

// used to initialize data required for the next step
void initialize_step(image_t image,step_t* step);

// provides a permanent storage of values referenced in step_t
typedef struct 
{ 
  struct p_norm_storage_t
  {
    float      p;         	//parameter for p_norm
  } p_norm;

  struct elliptic_norm_storage_t
  {
    float_pair main_axes; 	//parameter for elliptic_norm
  } elliptic_norm;

  struct fixed_rotation_storage_t
  {	
    float      angle; 		//rotation angle of structure element
  } fixed_rotation;

  struct auto_rotation_storage_t
  {	
	float      sigma_pre;	//std. deviation of presmoothing    kernel 
	float      sigma;		//std. deviation of strcture tensor kernel 
	image_t    angles; 		//rotation angles in rad
  } auto_rotation;

  struct adaptive_se_program_storage_t
  {
	float      sigma_pre;	//std. deviation of presmoothing    kernel
	float      sigma;		//std. deviation of strcture tensor kernel
	float      lambda1_max;
	float 	   lambda2_max;
	image_t    angles; 		//rotation angles in rad
	image_t	   lambda1; 	//larger  eigenvalue
	image_t	   lambda2;		//smaller eigenvalue
  } adaptive_se_program;
} param_storage_t;

typedef struct p_norm_storage_t         		p_norm_storage_t; 
typedef struct elliptic_norm_storage_t  		elliptic_norm_storage_t;
typedef struct fixed_rotation_storage_t 		fixed_rotation_storage_t;
typedef struct auto_rotation_storage_t  		auto_rotation_storage_t;
typedef struct adaptive_se_program_storage_t	adaptive_se_program_storage_t;

param_storage_t initialize_storage();

void parse_args(args_t* args,param_storage_t* storage,int argc,char** argv);

#endif
