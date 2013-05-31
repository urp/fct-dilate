/***************************************************************************
 *            common.c
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


#include "common.h"
#include "pgmio.h"
#include "alloc.h"

#include <float.h>
#include <stdio.h>
#include <tgmath.h>
#include <assert.h>

float	max(float a,float b)
{ return a>b?a:b; };

float	min(float a,float b)
{ return a<b?a:b; };

float	sgn(float a)
{ return a>0.f ? 1.f : ( a<0.f ? -1.f : 0.f ); }

float   minmod(float a,float b)
{
  if(a>0.f && b>0.f)
    return min(a,b);
  if(a<0.f && b<0.f)
    return max(a,b);
  return 0.f;
}

float   minmod3(float a,float b,float c)
{ 
  if(a>0.f && b>0.f && c>0.f)
    return min(a,min(b,c));
  if(a<0.f && b<0.f && c<0.f)
    return max(a,max(b,c));
  return 0.f;

}
