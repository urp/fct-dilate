/***************************************************************************
 *            initargs.c
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

#include "args.h"
#include "rt.h"
#include "os.h"
#include "fct.h"
#include "norm.h"
#include "rotation.h"
#include "integration.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <tgmath.h>
#include <float.h>
typedef enum { NO_KEY=0,KEYCHAR=1,KEYWORD=2} argkey_t;

argkey_t is_key(const char* str)
{
  const size_t len=strlen(str);
  if( len==2 && str[0]=='-' && str[1]!='-') return KEYCHAR;
  if( len>2  && str[0]=='-' && str[1]=='-') return KEYWORD;
  return NO_KEY;
}

int find_arg(const char   keychar,
             const char*  keyword,
             int          argc,
             char** argv)
{ size_t i;
  for(i=0;i<argc;++i)
  {
    argkey_t token = is_key(argv[i]);
    if(token==KEYWORD && strcmp(keyword,argv[i]+2)==0)
    return i;
    if(token==KEYCHAR && argv[i][1]==keychar)
        return i;
  }
  return -1;
}

void initialize_step(image_t image,step_t* step)
{
  if(step->program_prepare)
    step->program_prepare(image,step);

  if(step->rot_prepare)
    step->rot_prepare(image,step);

  if(step->norm_prepare)
    step->norm_prepare(image,step);
  /*
  step->init_state=1;
  step->rot(image,1,1,0,0,step);
  step->init_state=0;
  */
}

param_storage_t initialize_storage()
{
  param_storage_t storage   =
  { .p_norm                 = { .p          = 100.f },

    .elliptic_norm          = { .main_axes  = { .x = .01f,
                                                .y = 2.f
                                              }
                              },

    .fixed_rotation         = { .angle      = .1f*M_PI },

    .auto_rotation          = { .sigma_pre  = .5f,
                                .sigma      = 1.5f,
                                .angles     = { .data     = 0,
                                                .nx       = 0,
                                                .ny       = 0,
                                                .boundary = 0
                                              }
                              },
    .adaptive_se_program    = { .sigma_pre  = .5f,
                                .sigma      = 1.5f,
                                .angles     = { .data     = 0,
                                                .nx       = 0,
                                                .ny       = 0,
                                                .boundary = 0
                                              },
                                .lambda1    = { .data     = 0,
                                                .nx       = 0,
                                                .ny       = 0,
                                                .boundary = 0
                                              },
                                .lambda2    = { .data     = 0,
                                                .nx       = 0,
                                                .ny       = 0,
                                                .boundary = 0
                                              },
                                .lambda1_max= FLT_MIN,
                                .lambda2_max= FLT_MIN
                              }
  };
  return storage;
}

void print_usage_description()
{
  printf("Usage: dilate -i INPUT -o OUTPUT \n");
  printf("\t   [--method {rt|os|fct}]\n");
  printf("\t   [--steps STEPCOUNT] [-s STEPSIZE]\n");
  printf("\t   [--norm {circle|diamond|p [P]|elliptic [X_AXIS Y_AXIS]}]\n");
  printf("\t   [--rot {none|fixed [ANGLE]|auto [SIGMA [SIGMA_PRE]]}]\n");
}

void key_parse_exception(char* keyword)
{ printf("key_parse_exception: argument '%s' required but not specified\n",keyword);
  print_usage_description();
  exit(1);
}

void argument_parse_exception(char* keyword)
{ printf("parse_args: value argument '%s' has illegal format\n",keyword);
  print_usage_description();
  exit(1);
}

void parse_args(args_t* args,param_storage_t* params,int argc,char** argv)
{
  printf("parse_args: initializing arguments...\n");

  //boundary is set to minimal required value
  //implicitly by other options
  args->boundary = 0;


  // input file
  char* keyword = "input";
  char  keychar = 'i';
  int   pos = find_arg(keychar,keyword,argc,argv);
  if(pos!=-1 && ++pos<argc)
  { args->inpath=argv[pos];
    printf("setting '%s'\t\t= %s\n",keyword,args->inpath);
  }else
    key_parse_exception(keyword);


  // output file
  keyword = "output";
  keychar = 'o';
  pos = find_arg(keychar,keyword,argc,argv);
  if(pos!=-1 && ++pos<argc)
  { args->outpath=argv[pos];
    printf("setting '%s'\t= %s\n",keyword,args->outpath);
  }else
    key_parse_exception(keyword);


  // stepcount
  keyword = "steps";
  keychar = 'N';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t\t= ",keyword);
  if(pos!=-1 && ++pos<argc)
  { long count=strtol(argv[pos], (char **)NULL, 10);
    if(count==0)
      argument_parse_exception(keyword);
    args->stepcount=count;
    printf("%ld\n",args->stepcount);
  }else
  { args->stepcount=100;
    printf("%ld (default)\n",args->stepcount);
  }

  // stepsize
  keyword = "stepsize";
  keychar = 's';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t= ",keyword);
  if(pos!=-1 && ++pos<argc)
  { float size= (float) strtod(argv[pos], (char **)NULL);
    if(size==0.f)
      argument_parse_exception(keyword);
    args->step.stepsize=size;
    printf("%g\n",args->step.stepsize);
  }else
  { args->step.stepsize=.1f;
    printf("%g (default)\n",args->step.stepsize);
  }

  // write-step
  keyword = "write-steps";
  keychar = 'w';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t= ",keyword);
  if(pos!=-1)
  { printf("true\n");
    args->write_step=1;
  }else
  { printf("false\n");
    args->write_step=0;
  }


  // method
  keyword = "method";
  keychar = 'm';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t= ",keyword);
  if(pos!=-1 && ++pos<argc)
  {
    if(strcmp(argv[pos],"rt")==0)
    { args->step.integrate=&rouy_tourin_step;
      if(args->boundary<1)
        args->boundary=1;
      printf("rouy tourin\n");
    }
    if(strcmp(argv[pos],"os")==0)
    { args->step.integrate=&osher_sethian_step;
      if(args->boundary<2)
        args->boundary=2;
      printf("osher sethian\n");
    }
    if(strcmp(argv[pos],"fct")==0)
    { args->step.integrate=&flux_corrected_step;
      if(args->boundary<2)
        args->boundary=2;
      printf("flux corrected\n");
    }
  }else
  { args->step.integrate=&flux_corrected_step;
    if(args->boundary<2)
      args->boundary=2;
    printf("flux_corrected (default)\n");
  }

  // program
  keyword = "program";
  keychar = 'p';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t\t= ",keyword);
  //default parameters
  args->step.program_prepare=(void*)0;
  args->step.rot_prepare =(void*)0;
  args->step.rot_info =(void*)0;
  args->step.norm_prepare=(void*)0;
  args->step.norm_info=(void*)0;
  if(pos!=-1 && ++pos<argc)
  {
    //auto_rotation
    if(strcmp(argv[pos],"adaptive_se")==0)
    {

      args->step.program_prepare=&adaptive_se_program_prepare;

      args->step.rot        =&adaptive_se_program_rotation;
      args->step.norm       =&adaptive_se_program_norm;

      args->step.rot_info   =&params->adaptive_se_program;
      args->step.norm_info  =&params->adaptive_se_program;
      args->boundary++;
      printf("adaptive_se");

      if( (pos+1)<argc && !is_key(argv[pos+1]))
      { // read sigma
        const float sigma = (float) strtod(argv[pos+1], (char **)NULL);
        if(sigma<0.f)
          argument_parse_exception(keyword);
        params->adaptive_se_program.sigma=sigma;
        printf(" (sigma=%g",params->adaptive_se_program.sigma);
      }else
      { // default sigma
        params->adaptive_se_program.sigma=1.5f;
        printf(" (sigma=%g (default)",params->adaptive_se_program.sigma);
      }

      if( (pos+2)<argc && !is_key(argv[pos+1]) && !is_key(argv[pos+2]))
      { // read sigma_pre
        const float sigma_pre = (float) strtod(argv[pos+2], (char **)NULL);
        if(sigma_pre<0.f)
          argument_parse_exception(keyword);
        params->adaptive_se_program.sigma_pre=sigma_pre;
        printf(",sigma_pre=%g)\n",params->adaptive_se_program.sigma_pre);
      }else
      { // default sigma_pre
        params->adaptive_se_program.sigma_pre=.5f;
        printf(",sigma_pre=%g (default))\n",params->adaptive_se_program.sigma_pre);
      }
      return;
    }
  }

  // norm
  keyword = "norm";
  keychar = 'n';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t\t= ",keyword);
  //parameter default
  args->step.norm_prepare=(void*)0;
  args->step.norm_info=(void*)0;
  if(pos!=-1 && ++pos<argc)
  {
    //circle_norm
    if(strcmp(argv[pos],"circle")==0)
    { args->step.norm     =&circle_norm;
      printf("circle_norm\n");
    }
    //diamond_norm
    if(strcmp(argv[pos],"diamond")==0)
    { args->step.norm     =&diamond_norm;
      printf("diamond_norm\n");
    }
    //p_norm
    if(strcmp(argv[pos],"p")==0)
    {
      if(++pos<argc && !is_key(argv[pos]))
      { // read parameter
        float p= (float) strtod(argv[pos], (char **)NULL);
        if(p==0.f)
          argument_parse_exception(keyword);
        params->p_norm.p = p;
        printf("p_norm (p=%g)\n",params->p_norm.p);
      }else
      { // default parameter
        params->p_norm.p = 100.f;
        printf("p_norm (default p=%g)\n",params->p_norm.p);
      }
      args->step.norm=&p_norm;
      args->step.norm_info=&params->p_norm;
    }
    //elliptic_norm
    if(strcmp(argv[pos],"elliptic")==0)
    {
      printf("elliptic_norm");
      if( (pos+2)<argc && !is_key(argv[pos+1]) && !is_key(argv[pos+2]) )
      { // read parameters
        float ax= (float) strtod(argv[pos+1], (char **)NULL);
        float ay= (float) strtod(argv[pos+2], (char **)NULL);
        if(ax<=0.f||ay<=0.f)
          argument_parse_exception(keyword);
        params->elliptic_norm.main_axes.x=ax;
        params->elliptic_norm.main_axes.y=ay;
        printf(" (a=%g,b=%g)\n",params->elliptic_norm.main_axes.x,
                                params->elliptic_norm.main_axes.y);
      }else
      { // default parameters
        params->elliptic_norm.main_axes.x=1.f;
        params->elliptic_norm.main_axes.y=.1f;
        printf(" (a=%g (default),b=%g (default))\n",
               params->elliptic_norm.main_axes.x,
               params->elliptic_norm.main_axes.y    );
      }
      args->step.norm=&elliptic_norm;
      args->step.norm_info=&params->elliptic_norm;
    }
  }else
  { args->step.norm=&circle_norm;
    printf("circle_norm (default)\n");
  }

  // grid rotation
  keyword = "rotation";
  keychar = 'r';
  pos = find_arg(keychar,keyword,argc,argv);
  printf("setting '%s'\t= ",keyword);
  //parameter default
  args->step.rot_prepare=(void*)0;
  args->step.rot_info=(void*)0;
  if(pos!=-1 && ++pos<argc)
  {
    //no_rotation
    if(strcmp(argv[pos],"none")==0)
    { args->step.rot     =&no_rotation;
      printf("no_rotation\n");
    }

    //fixed_rotation
    if(strcmp(argv[pos],"fixed")==0)
    {
      printf("fixed_rotation");
      if(++pos<argc && !is_key(argv[pos]))
      { // read parameters
        const float angle= (float) strtod(argv[pos], (char **)NULL);
        if(angle==0.f)
          argument_parse_exception(keyword);
        params->fixed_rotation.angle=angle;
        printf(" (angle=%g)\n",params->fixed_rotation.angle);
      }else
      { // default parameters
        params->fixed_rotation.angle = 0.f;
        printf(" (default angle=%g)\n",params->fixed_rotation.angle);
      }
      args->step.rot=&fixed_rotation;
      args->step.rot_info=&params->fixed_rotation;
      args->boundary++;
    }

    //auto_rotation
    if(strcmp(argv[pos],"auto")==0)
    { args->step.rot        =&auto_rotation;
      args->step.rot_prepare=&auto_rotation_prepare;
      args->step.rot_info   =&params->auto_rotation;
      args->boundary++;
      printf("auto_rotation");

      if( (pos+1)<argc && !is_key(argv[pos+1]))
      { // read sigma
        const float sigma = (float) strtod(argv[pos+1], (char **)NULL);
        if(sigma<0.f)
          argument_parse_exception(keyword);
        params->auto_rotation.sigma=sigma;
        printf(" (sigma=%g",params->auto_rotation.sigma);
      }else
      { // default sigma
        params->auto_rotation.sigma=1.5f;
        printf(" (sigma=%g (default)",params->auto_rotation.sigma);
      }

      if( (pos+2)<argc && !is_key(argv[pos+1]) && !is_key(argv[pos+2]))
      { // read sigma_pre
        const float sigma_pre = (float) strtod(argv[pos+2], (char **)NULL);
        if(sigma_pre<0.f)
          argument_parse_exception(keyword);
        params->auto_rotation.sigma_pre=sigma_pre;
        printf(",sigma_pre=%g)\n",params->auto_rotation.sigma_pre);
      }else
      { // default sigma_pre
        params->auto_rotation.sigma_pre=.5f;
        printf(",sigma_pre=%g (default))\n",params->auto_rotation.sigma_pre);
      }

    }

  }else
  { args->step.rot=&no_rotation;
    printf("no_rotation (default)\n");
  }

}
