/***************************************************************************
 *            pgmio.c
 *
 *  Tue Jul  7 19:12:59 2009
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

#include "common.h"
#include "alloc.h"
#include "pgmio.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

#define 	PGMIO_MAX_ROW_SIZE 255

void  		readline(char* line,FILE* inimage)
{
  if(fgets(line,PGMIO_MAX_ROW_SIZE,inimage)==0) 
  { printf("load_pgm_image\t|ERROR: could not read from file\n"); exit(1); };
}

image_t load_pgm_image(char* in,long boundary)
{
  FILE*	  inimage;
  char    row[PGMIO_MAX_ROW_SIZE];
  image_t image = { .data=0,.nx=0,.ny=0,.boundary=boundary };
  
  // open pgm file and read header
  //printf("load_pgm_image: reading file %s ...\n",in);
  inimage = fopen(in,"r");
  
  if(!inimage) // file not open
  { printf("load_pgm_image\t|ERROR: Could not open file %s\n",in);
	return image; 
  }
  
  readline(row, inimage);
  readline(row, inimage);  
  while (row[0]=='#') 
  { printf("%s",row);
    readline(row, inimage);
  }
  // read size
  sscanf (row, "%ld %ld", &(image.nx), &(image.ny));
  readline(row, inimage);  
  printf("load_pgm_image\t|size (%ld,%ld)\n",image.nx,image.ny);

  // allocate storage
  alloc_image(&image);
  long i,j;
  // read image data
  for (j=boundary; j<(image.ny+boundary); j++)
   for (i=boundary; i<(image.nx+boundary); i++)
     image.data[i][j] = (float) getc (inimage);
  fclose(inimage);
  return image;
}


void write_pgm_image(char* out,image_t image)
{
  printf("write_pgm_image\t|saving image %s...",out);   
  FILE   *outimage;
  char rowSize[10];
  char colSize[10];
  char Color[10];
	
  sprintf(rowSize,"%ld",image.ny);
  sprintf(colSize,"%ld",image.nx);
  sprintf(Color,"%d",255);
  //create target file
  outimage = fopen(out,"wb");
  
  if(!outimage) // file not open
  { printf("write_pgm_image\t|ERROR: Could not open file %s\n",out);
	return; 
  }

  //format for pgm
  fputs("P5",outimage);					
  fputc('\n',outimage);
  //comments
  fputs("#this is a pgm file in binaray format",outimage);	
  fputc('\n',outimage); 
  //write number of rows and cols and max value of grey level
  fputs(colSize,outimage);				
  fputc(' ',outimage);
  
  fputs(rowSize,outimage);
  fputc('\n',outimage);
  
  fputs(Color,outimage);		
  fputc('\n',outimage);

  //clamp pixel values to [0,255]
  long i,j;
  for (j=image.boundary; j<(image.ny+image.boundary); j++)				
    for (i=image.boundary; i<(image.nx+image.boundary); i++)
      fputc( (unsigned char) max(min(image.data[i][j],255.f),0.f) , outimage);
		  
  fclose(outimage); 
  printf("complete\n"); 
}
