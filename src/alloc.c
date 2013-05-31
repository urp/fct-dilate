/***************************************************************************
 *            alloc.c
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

#include	<assert.h>
#include 	<stdlib.h>
#include 	<stdio.h>

#include	"alloc.h"

// allocates storage for a vector of size n 
void alloc_vector( float **vector, // vector 
		   long   n	)  // size 
{
  *vector = (float *) malloc (n * sizeof(float));
  if (*vector == 0)
  {
     printf("alloc_vector: not enough storage available\n");
     exit(1);
  }
  return;
}

// allocates storage for matrix of size nx * ny 
void alloc_matrix(	float*** matrix,  	// matrix
			long   	  nx,         	// size in x direction
			long      ny		// size in y direction
                 )         
{
  *matrix = (float **) malloc (nx * sizeof(float *));

  if (*matrix == NULL)
  {
    printf("alloc_matrix: not enough storage available\n");
    exit(1);
  }

  long i;
  for (i=0; i<nx; i++)
  {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
    {
      printf("alloc_matrix: not enough storage available\n");
      exit(1);
    }
  }
  return;
}

// disallocates storage for a vector of size n
void disalloc_vector(float *vector)
{ free(vector); }

// disallocates storage for matrix of size nx * ny
void disalloc_matrix(	float **matrix,	// matrix
			long   nx,     	// size in x direction
			long   ny      	// size in y direction
                    )
{
  long i;
  for (i=0; i<nx; i++)
    free(matrix[i]);
  free(matrix);
}

void copy_matrix( float**	from,
				  float**	to,
				  long 		nx,
				  long		ny
				)
{ 
  long i,j;
  #pragma omp parallel shared(from,to,nx,ny) private(i,j)
  {
    #pragma omp for schedule(static) nowait
    for (i=0; i<nx; i++)
      for (j=0; j<ny; j++)
		to[i][j] = from[i][j];
  }
}
