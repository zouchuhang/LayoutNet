/*----------------------------------------------------------------------------

  LSD - Line Segment Detector on digital images

  Copyright (c) 2007-2011 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** @file lsd_cmd.c
    Command line interface for LSD module (Line Segment Detector).
    @author rafael grompone von gioi <grompone@gmail.com>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Definition of the command line interface.                                 */
#define USE "                                                                  \
#name: lsd                                                                     \
#author: rafael grompone von gioi <grompone@gmail.com>                         \
#version: 1.6 of November 11, 2011                                             \
#year: 2007-2011                                                               \
#desc: Line Segment Detector                                                   \
#opt: scale | s | double | 0.8 | 0.0 | |                                       \
      Scale image by Gaussian filter before processing.                        \
#opt: sigma_coef | c | double | 0.6 | 0.0 | |                                  \
      Sigma for Gaussian filter is computed as sigma_coef/scale.               \
#opt: quant | q | double | 2.0 | 0.0 | |                                       \
      Bound to quantization error on the gradient norm.                        \
#opt: ang_th | a | double | 22.5 | 0.0 | 180.0 |                               \
      Gradient angle tolerance in degrees.                                     \
#opt: log_eps | e | double | 0.0 | | | Detection threshold, -log10(max. NFA)   \
#opt: density_th | d | double | 0.7 | 0.0 | 1.0 |                              \
      Minimal density of region points in a rectangle to be accepted.          \
#opt: n_bins | b | int | 1024 | 1 | |                                          \
      Number of bins in 'ordering' of gradient modulus.                        \
#opt: reg | R | str | | | |                                                    \
      Output image: owner LS number at each pixel. Scaled size. (PGM)          \
#opt: epsfile | P | str | | | | Output line segments into EPS file 'epsfile'.  \
#opt: svgfile | S | str | | | | Output line segments into SVG file 'svgfile'.  \
#opt: width | W | double | 1.5 | | |                                           \
      LS width used in EPS and SVG files. If <=0, use detected values.         \
#req: in  | | str | | | | Input image (PGM)                                    \
#req: out | | str | | | |                                                      \
      Line Segment output (each ascii line: x1,y1,x2,y2,width,p,-log10(NFA) )  \
"
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "lsd.h"

#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/** Fatal error, print a message to standard-error output and exit.
 */
static void error(char * msg)
{
  fprintf(stderr,"%s\n",msg);
  exit(EXIT_FAILURE);
}


/*----------------------------------------------------------------------------*/
/*--------------------- Command Line interface handling ----------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#define FIELD_LENGTH 160
#define HELP_OPTION    "--help"
#define VERSION_OPTION "--version"

/*----------------------------------------------------------------------------*/
/** Structure to store one argument definition and read value.
 */
struct argument
{
  char name[FIELD_LENGTH];     /* name to internally identify the argument */
  char desc[FIELD_LENGTH];     /* description */
  char id;                     /* letter used with '-' to use the option */
  char type;                   /* i=int, d=double, s=str, b=bool */
  int required;
  int assigned;
  int def_value;               /* true or false, a default value is assigned? */
  char d_value[FIELD_LENGTH];  /* default value */
  char s_value[FIELD_LENGTH];  /* string found, also the value if 'str' */
  int  i_value;
  double f_value;
  int    min_set;              /* true or false, is minimal value set? */
  double min;
  int    max_set;              /* true or false, is maximal value set? */
  double max;
};

/*----------------------------------------------------------------------------*/
/** Structure to store the full set of argument definitions and its values.
 */
struct arguments
{
  char name[FIELD_LENGTH];
  char author[FIELD_LENGTH];
  char version[FIELD_LENGTH];
  char desc[FIELD_LENGTH];
  char compiled[FIELD_LENGTH];
  char year[FIELD_LENGTH];
  int  arg_num;
  int  arg_allocated;
  struct argument * args;
};

/*----------------------------------------------------------------------------*/
/** Free an 'arguments' structure.
 */
static void free_arguments(struct arguments * arg)
{
  if( arg == NULL || arg->args == NULL )
    error("Error: NULL pointer at 'free_arguments'.");
  free( (void *) arg->args );
  free( (void *) arg );
}

/*----------------------------------------------------------------------------*/
/** Accepted characters in field identifier: numbers, letters and '_'.
 */
static int is_id_char(int c)
{
  return c=='_' || isalpha(c) || isdigit(c);
}

/*----------------------------------------------------------------------------*/
/** Read next field definition in an argument definition.
 */
static char * get_next_field(char * p, char * id, char * value)
{
  int n;

  /* check input */
  if( p == NULL || id == NULL || value == NULL )
    error("Error: invalid input to 'get_next_field'.");

  /* search for field id */
  while( isspace(*p) ) ++p; /* skip spaces */
  if( *p != '#' ) error("Error: missing '#' in 'use description'.");
  ++p;
  for( n=0; is_id_char(*p) && n<FIELD_LENGTH; n++ ) id[n] = *(p++);
  if( n >= FIELD_LENGTH ) error("Error: field too long in 'use description'.");
  id[n] = '\0';
  if( *(p++) != ':' ) error("Error: missing ':' in 'use description'.");

  /* search for field value */
  while( isspace(*p) ) ++p; /* skip spaces */
  for( n=0; *p != '#' && *p != '\0' && n<FIELD_LENGTH; n++ ) value[n] = *(p++);
  if( n >= FIELD_LENGTH ) error("Error: field too long in 'use description'.");
  value[n] = '\0';

  /* remove spaces at the end of the field */
  while( --n >= 0 && isspace(value[n]) ) value[n] = '\0';

  return p;
}

/*----------------------------------------------------------------------------*/
/** Read next token in an argument definition.
 */
static char * get_next_token(char * p, char div, char * value)
{
  int n;

  /* check input */
  if( p == NULL || value == NULL )
    error("Error: invalid input to 'get_next_token'.");

  if( *p == '\0' )
    error("Error: argument token expected in 'use description'.");

  while( isspace(*p) ) ++p; /* skip spaces */
  for( n=0; *p!=div && *p!='\0' && n<FIELD_LENGTH; n++) value[n] = *(p++);
  if( n >= FIELD_LENGTH ) error("Error: field too long in 'use description'.");
  value[n] = '\0';
  while( --n >= 0 && isspace(value[n]) ) value[n] = '\0';

  /* remove 'div' at the end of the token, if present */
  if( *p == div ) ++p;

  return p;
}

/*----------------------------------------------------------------------------*/
/** Process one argument description.
 */
static void process_new_argument(char * id, char * value,struct arguments * arg)
{
  char token[FIELD_LENGTH];
  char * p;
  int i;

  /* check input */
  if( id == NULL || value == NULL || arg == NULL )
    error("Error: invalid input to 'process_new_argument'.");

  /* allocate memory if needed */
  if( arg->arg_num >= arg->arg_allocated )
    {
      arg->arg_allocated *= 2;
      arg->args = (struct argument *) realloc( (void *) arg->args,
                                 arg->arg_allocated * sizeof(struct argument) );
      if( arg->args == NULL ) error("Error: not enough memory.");
    }

  /* argument name */
  p = get_next_token(value,'|',arg->args[arg->arg_num].name);
  for( i=0; i<arg->arg_num; i++ )
    if( strcmp(arg->args[i].name,arg->args[arg->arg_num].name) == 0 )
      error("Error: argument name used twice in 'use description'.");

  /* 'option' letter - to be used with '-' to identify option */
  p = get_next_token(p,'|',token);
  if( strcmp(id,"opt") == 0 )
    {
      arg->args[arg->arg_num].required = FALSE;
      if( strlen(token) != 1 )
        error("Error: invalid option letter in 'use description'.");
      arg->args[arg->arg_num].id = token[0];
      if( !isalpha(arg->args[arg->arg_num].id) )
        error("Error: option id must be a letter in 'use description'.");
      for( i=0; i<arg->arg_num; i++ )
        if( !(arg->args[i].required) &&
            arg->args[i].id == arg->args[arg->arg_num].id )
          error("Error: option letter used twice in 'use description'.");
    }
  else /* must be 'req' - required argument */
    {
      arg->args[arg->arg_num].required = TRUE;
      if( strlen(token) > 0 )
        error("Error: unused option letter in 'use description'.");
      arg->args[arg->arg_num].id = 0;
    }

  /* argument type */
  p = get_next_token(p,'|',token);
  if( strcmp(token,"int") == 0 )         arg->args[arg->arg_num].type = 'i';
  else if( strcmp(token,"double") == 0 ) arg->args[arg->arg_num].type = 'd';
  else if( strcmp(token,"str") == 0 )    arg->args[arg->arg_num].type = 's';
  else if( strcmp(token,"bool") == 0 )   arg->args[arg->arg_num].type = 'b';
  else error("Error: unknown argument type in 'use description'.");

  /* required arguments can't be boolean */
  if( arg->args[arg->arg_num].required && arg->args[arg->arg_num].type == 'b' )
    error("Error: required arguments can't be boolean in 'use description'.");

  /* default value */
  p = get_next_token(p,'|',token);
  if( strlen(token) > 0 )
    {
      if( arg->args[arg->arg_num].required )
       error("Error: default value in required argument in 'use description'.");
      arg->args[arg->arg_num].def_value = TRUE;
      arg->args[arg->arg_num].assigned = TRUE;
      strcpy(arg->args[arg->arg_num].d_value,token);
      strcpy(arg->args[arg->arg_num].s_value,token);
      if( arg->args[arg->arg_num].type == 'i' )
        arg->args[arg->arg_num].i_value = atoi(token);
      if( arg->args[arg->arg_num].type == 'd' )
        arg->args[arg->arg_num].f_value = atof(token);
    }
  else
    {
      arg->args[arg->arg_num].def_value = FALSE;
      arg->args[arg->arg_num].s_value[0] = '\0';
      arg->args[arg->arg_num].assigned = FALSE;
    }

  /* required arguments can't have default value */
  if( arg->args[arg->arg_num].required && arg->args[arg->arg_num].def_value )
   error("Error: required args can't have default value in 'use description'.");

  /* min value */
  p = get_next_token(p,'|',token);
  if( strlen(token) > 0 )
    {
      arg->args[arg->arg_num].min_set = TRUE;
      arg->args[arg->arg_num].min = atof(token);
    }
  else
    {
      arg->args[arg->arg_num].min_set = FALSE;
    }

  /* max value */
  p = get_next_token(p,'|',token);
  if( strlen(token) > 0 )
    {
      arg->args[arg->arg_num].max_set = TRUE;
      arg->args[arg->arg_num].max = atof(token);
    }
  else
    {
      arg->args[arg->arg_num].max_set = FALSE;
    }

  /* argument description */
  p = get_next_token(p,'|',arg->args[arg->arg_num].desc);

  /* the field should end there */
  if( *p != '\0' )
    error("Error: too many tokens in one argument in 'use description'.");

  arg->arg_num++;
}

/*----------------------------------------------------------------------------*/
/** Process an argument definition.
 */
static void process_argument_description( char * desc, struct arguments * arg )
{
  char id[FIELD_LENGTH];
  char value[FIELD_LENGTH];

  /* check input */
  if( desc == NULL || arg == NULL )
    error("Error: invalid input to 'process_argument_description'.");

  /* initialize 'arg' */
  arg->name[0] = '\0';
  arg->author[0] = '\0';
  arg->version[0] = '\0';
  arg->year[0] = '\0';
  arg->desc[0] = '\0';
  arg->compiled[0] = '\0';
  arg->arg_num = 0;
  arg->arg_allocated = 2;
  arg->args = (struct argument *)
                      malloc( arg->arg_allocated * sizeof(struct argument) );
  if( arg->args == NULL ) error("Error: not enough memory.");

  /* assign compilation date and time */
  strcat(arg->compiled,__DATE__);
  strcat(arg->compiled," ");
  strcat(arg->compiled,__TIME__);

  /* process description */
  while( *desc != '\0' )
    {
      desc = get_next_field(desc,id,value);

      if( strcmp(id,"name") == 0 )
        {
          if( arg->name[0] != '\0' )
            error("Error: multiple 'name' fields in 'use description'.");
          strcpy(arg->name,value);
        }
      else if( strcmp(id,"author") == 0 )
        {
          if( arg->author[0] != '\0' )
            error("Error: multiple 'author' fields in 'use description'.");
          strcpy(arg->author,value);
        }
      else if( strcmp(id,"version") == 0 )
        {
          if( arg->version[0] != '\0' )
            error("Error: multiple 'version' fields in 'use description'.");
          strcpy(arg->version,value);
        }
      else if( strcmp(id,"year") == 0 )
        {
          if( arg->year[0] != '\0' )
            error("Error: multiple 'year' fields in 'use description'.");
          strcpy(arg->year,value);
        }
      else if( strcmp(id,"desc") == 0 )
        {
          if( arg->desc[0] != '\0' )
            error("Error: multiple 'desc' fields in 'use description'.");
          strcpy(arg->desc,value);
        }
      else if( strcmp(id,"opt") == 0 || strcmp(id,"req") == 0 )
        {
          process_new_argument(id,value,arg);
        }
      else
        {
          error("Error: unknown token in 'use description'.");
        }
    }

  /* verify required arguments */
  if( arg->name[0] == '\0' )
    error("Error: program name is required in 'use description'.");
  if( arg->author[0] == '\0' )
    error("Error: author name is required in 'use description'.");
  if( arg->version[0] == '\0' )
    error("Error: version is required in 'use description'.");
  if( arg->desc[0] == '\0' )
    error("Error: program description is required in 'use description'.");
  if( arg->year[0] == '\0' )
    error("Error: year is required in 'use description'.");
}

/*----------------------------------------------------------------------------*/
/** Print version.
 */
static void print_version(struct arguments * arg, FILE * f)
{
  if( arg == NULL || f == NULL )
    error("Error: invalid input to 'print_version'.");
  fprintf(f,"Version %s, compiled %s\n",arg->version,arg->compiled);
}

/*----------------------------------------------------------------------------*/
/** Print command line interface help and exit.
 */
static void use(struct arguments * arg)
{
  int i;

  if( arg == NULL ) error("Error: invalid input to 'use'.");

  fprintf(stderr,"----------------------------------------");
  fprintf(stderr,"----------------------------------------\n");
  fprintf(stderr,"This is %s, ",arg->name);
  print_version(arg,stderr);
  fprintf(stderr,"%s\n",arg->desc);
  fprintf(stderr,"Copyright (c) %s %s\n",arg->year,arg->author);
  fprintf(stderr,"----------------------------------------");
  fprintf(stderr,"----------------------------------------\n");
  fprintf(stderr,"\nUsage: %s",arg->name);

  /* always present 'help' option */
  fprintf(stderr," [%s]",HELP_OPTION);

  /* always present version option */
  fprintf(stderr," [%s]",VERSION_OPTION);

  for(i=0;i<arg->arg_num;i++)
    if( !(arg->args[i].required) )
      {
        fprintf(stderr," [-%c",arg->args[i].id);
        if( arg->args[i].type != 'b' ) fprintf(stderr," %s",arg->args[i].name);
        fprintf(stderr,"]");
      }

  for(i=0;i<arg->arg_num;i++)
    if( arg->args[i].required )
      fprintf(stderr," %s",arg->args[i].name);

  fprintf(stderr,"\n\n");

  /* option description */
  fprintf(stderr,"  %s\tPrint this help message and exit.\n",
          HELP_OPTION);
  fprintf(stderr,"  %s\tPrint version and compilation date/time and exit.\n",
          VERSION_OPTION);
  for(i=0;i<arg->arg_num;i++)
    if( !(arg->args[i].required) )
      {
        fprintf(stderr,"  -%c",arg->args[i].id);
        if( arg->args[i].type != 'b' )
          {
            fprintf(stderr," %s",arg->args[i].name);
          }
        fprintf(stderr,"\t%s\n",arg->args[i].desc);
        if( arg->args[i].type == 'i' )
          {
            fprintf(stderr,"\t\t'%s' is integer",arg->args[i].name);
            fprintf(stderr,", range [");
            if( arg->args[i].min_set )
              fprintf(stderr,"%d,",(int)arg->args[i].min);
            else fprintf(stderr,"-inf,");
            if( arg->args[i].max_set )
              fprintf(stderr,"%d]",(int)arg->args[i].max);
            else fprintf(stderr,"inf]");
            if( arg->args[i].def_value )
              fprintf(stderr,", default value %d",atoi(arg->args[i].d_value));
            fprintf(stderr,"\n");
          }
        if( arg->args[i].type == 'd' )
          {
            fprintf(stderr,"\t\t'%s' is double",arg->args[i].name);
            fprintf(stderr,", range [");
            if( arg->args[i].min_set ) fprintf(stderr,"%g,",arg->args[i].min);
            else fprintf(stderr,"-inf,");
            if( arg->args[i].max_set ) fprintf(stderr,"%g]",arg->args[i].max);
            else fprintf(stderr,"inf]");
            if( arg->args[i].def_value )
              fprintf(stderr,", default value %g",atof(arg->args[i].d_value));
            fprintf(stderr,"\n");
          }
      }

  for(i=0;i<arg->arg_num;i++)
    if( arg->args[i].required )
      {
        fprintf(stderr,"  %s",arg->args[i].name);
        fprintf(stderr,"\t%s\n",arg->args[i].desc);
        if( arg->args[i].type == 'i' )
          {
            fprintf(stderr,"\t\t'%s' is integer",arg->args[i].name);
            fprintf(stderr,", range [");
            if( arg->args[i].min_set )
              fprintf(stderr,"%d,",(int)arg->args[i].min);
            else fprintf(stderr,"-inf,");
            if( arg->args[i].max_set )
              fprintf(stderr,"%d]",(int)arg->args[i].max);
            else fprintf(stderr,"inf]");
            fprintf(stderr,"\n");
          }
        if( arg->args[i].type == 'd' )
          {
            fprintf(stderr,"\t\t'%s' is double",arg->args[i].name);
            fprintf(stderr,", range [");
            if( arg->args[i].min_set ) fprintf(stderr,"%f,",arg->args[i].min);
            else fprintf(stderr,"-inf,");
            if( arg->args[i].max_set ) fprintf(stderr,"%f]",arg->args[i].max);
            else fprintf(stderr,"inf]");
            fprintf(stderr,"\n");
          }
      }

  fprintf(stderr,"\n");

  free_arguments(arg);

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/** Evaluate arguments.
 */
static void evaluate_arguments(int argc, char ** argv, struct arguments * arg)
{
  int in_required_args = FALSE;
  int n,i;

  /* check input */
  if( argc <= 0 ) error("Error: unexpected command line: missing command.");
  if( argv == NULL || arg == NULL )
    error("Error: invalid input to 'evaluate_arguments'.");

  for( n=1; !in_required_args && n < argc; n++ )
    {
      /* when an argument do not start with "-" it is not optional.
         but, if the argument is just "-", then is a non optional
         argument with value "-", and will be analyzed later.  */
      if( argv[n][0] != '-' || (argv[n][0]=='-' && strlen(argv[n])== 1) )
        {
          in_required_args = TRUE;
          --n;
          continue;
        }

      if( strlen(argv[n]) != 2 )
        {
          /* check if it is the especial option 'help' */
          if( strcmp(argv[n],HELP_OPTION) == 0 ) use(arg);

          /* check if it is the special option 'version' */
          if( strcmp(argv[n],VERSION_OPTION) == 0 )
            {
              print_version(arg,stdout);
              free_arguments(arg);
              exit(EXIT_SUCCESS);
            }

          /* otherwise is a bad option */
          fprintf(stderr,"Error: %s ",argv[n]);
          error("unrecognized option.");
        }

      for( i=0; i<arg->arg_num; i++ )
        if( !(arg->args[i].required) && arg->args[i].id == argv[n][1] )
          {
            arg->args[i].assigned = TRUE;
            if( arg->args[i].type != 'b' )
              {
                /* go for the value */
                ++n;

                /* a value is expected */
                if( n >= argc )
                  {
                    fprintf(stderr,"Error: in '%s': ",argv[n-1]);
                    error("a value was expected.");
                  }
                if( strlen(argv[n]) > FIELD_LENGTH )
                  {
                    fprintf(stderr,"Error: in '%s': ",argv[n-1]);
                    error("value too long.");
                  }
                strcpy(arg->args[i].s_value,argv[n]);
                if( arg->args[i].type == 'i' )
                  {
                    arg->args[i].i_value = atoi(argv[n]);
                    if( arg->args[i].min_set &&
                        arg->args[i].i_value < (int) arg->args[i].min )
                      {
                        fprintf(stderr,"Error: in '%s': ",argv[n-1]);
                        error("value out of range.");
                      }
                    if( arg->args[i].max_set &&
                        arg->args[i].i_value > (int) arg->args[i].max )
                      {
                        fprintf(stderr,"Error: in '%s': ",argv[n-1]);
                        error("value out of range.");
                      }
                  }
                if( arg->args[i].type == 'd' )
                  {
                    arg->args[i].f_value = atof(argv[n]);
                    if( arg->args[i].min_set &&
                        arg->args[i].f_value < arg->args[i].min )
                      {
                        fprintf(stderr,"Error: in '%s': ",argv[n-1]);
                        error("value out of range.");
                      }
                    if( arg->args[i].max_set &&
                        arg->args[i].f_value > arg->args[i].max )
                      {
                        fprintf(stderr,"Error: in '%s': ",argv[n-1]);
                        error("value out of range.");
                      }
                  }
              }
            i = arg->arg_num; /* argument found, stop search */
          }
    }

  for( i=0; n<argc && i<arg->arg_num; i++ )
    if( arg->args[i].required )
      {
        arg->args[i].assigned = TRUE;
        strcpy(arg->args[i].s_value,argv[n]);
        if( arg->args[i].type == 'i' )
          {
            arg->args[i].i_value = atoi(argv[n]);
            if( arg->args[i].min_set &&
                arg->args[i].i_value < (int) arg->args[i].min )
              {
                fprintf(stderr,"Error: in '%s': ",arg->args[i].name);
                error("value out of range.");
              }
            if( arg->args[i].max_set &&
                arg->args[i].i_value > (int) arg->args[i].max )
              {
                fprintf(stderr,"Error: in '%s': ",arg->args[i].name);
                error("value out of range.");
              }
          }
        if( arg->args[i].type == 'd' )
          {
            arg->args[i].f_value = atof(argv[n]);
            if( arg->args[i].min_set &&
                arg->args[i].f_value < arg->args[i].min )
              {
                fprintf(stderr,"Error: in '%s': ",arg->args[i].name);
                error("value out of range.");
              }
            if( arg->args[i].max_set &&
                arg->args[i].f_value > arg->args[i].max )
              {
                fprintf(stderr,"Error: in '%s': ",arg->args[i].name);
                error("value out of range.");
              }
          }
        ++n;
      }
}

/*----------------------------------------------------------------------------*/
/** Process and evaluate a program arguments.
 */
static struct arguments * process_arguments(char * desc, int argc, char ** argv)
{
  struct arguments * arg;
  int i;

  /* check input */
  if( desc == NULL || argv == NULL )
    error("Error: invalid input to 'process_arguments'.");

  /* get memory */
  arg = (struct arguments *) malloc(sizeof(struct arguments));
  if( arg == NULL ) error("Error: not enough memory.");

  process_argument_description(desc,arg);
  evaluate_arguments(argc,argv,arg);

  /* if there are missing arguments print the 'use' information */
  for(i=0; i<arg->arg_num; i++)
    if( arg->args[i].required && !(arg->args[i].assigned) ) use(arg);

  return arg;
}

/*----------------------------------------------------------------------------*/
/** Test if an argument has a defined value.
 */
static int is_assigned(struct arguments * arg, char * name)
{
  int i;
  if( arg == NULL || name == NULL )
    error("Error: invalid input to 'is_assigned'.");
  for(i=0; i<arg->arg_num; i++)
    if( strcmp(name,arg->args[i].name) == 0 ) return arg->args[i].assigned;
  error("Error: is_assigned: unknown argument.");
  return -1; /* useless, just to prevent warning in strict compilers */
}

/*----------------------------------------------------------------------------*/
/** Get the value of a string argument.
 */
static char * get_str(struct arguments * arg, char * name)
{
  int i;

  if( arg == NULL || name == NULL )
    error("Error: invalid input to 'get_str'.");

  for(i=0; i<arg->arg_num; i++)
    if( strcmp(name,arg->args[i].name) == 0 )
      {
        if( arg->args[i].type == 's' )
          {
            if( !(arg->args[i].assigned) ) return NULL;
            return arg->args[i].s_value;
          }
        else error("Error: get_str: the parameter is not a double.");
      }
  error("Error: get_str: unknown argument.");
  return NULL; /* useless, just to prevent warning in strict compilers */
}

/*----------------------------------------------------------------------------*/
/** Get the value of an integer argument.
 */
static int get_int(struct arguments * arg, char * name)
{
  int i;

  if( arg == NULL || name == NULL )
    error("Error: invalid input to 'get_int'.");

  for(i=0; i<arg->arg_num; i++)
    if( strcmp(name,arg->args[i].name) == 0 )
      {
        if( !(arg->args[i].assigned) )
          error("Error: get_int: parameter not assigned.");
        if( arg->args[i].type == 'i' ) return arg->args[i].i_value;
        else error("Error: get_int: the parameter is not an integer.");
      }
  error("Error: get_int: unknown argument.");
  return -1; /* useless, just to prevent warning in strict compilers */
}

/*----------------------------------------------------------------------------*/
/** Get the value of a double argument.
 */
static double get_double(struct arguments * arg, char * name)
{
  int i;

  if( arg == NULL || name == NULL )
    error("Error: invalid input to 'get_double'.");

  for(i=0; i<arg->arg_num; i++)
    if( strcmp(name,arg->args[i].name) == 0 )
      {
        if( !(arg->args[i].assigned) )
          error("Error: get_double: parameter not assigned.");
        if( arg->args[i].type == 'd' ) return arg->args[i].f_value;
        else error("Error: get_double: the parameter is not a double.");
      }
  error("Error: get_double: unknown argument.");
  return -1.0; /* useless, just to prevent warning in strict compilers */
}


/*----------------------------------------------------------------------------*/
/*------------------------------ PGM image I/O -------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Skip white characters and comments in a PGM file.
 */
static void skip_whites_and_comments(FILE * f)
{
  int c;
  do
    {
      while(isspace(c=getc(f))); /* skip spaces */
      if(c=='#') /* skip comments */
        while( c!='\n' && c!='\r' && c!=EOF )
          c=getc(f);
    }
  while( c == '#' || isspace(c) );
  if( c != EOF && ungetc(c,f) == EOF )
    error("Error: unable to 'ungetc' while reading PGM file.");
}

/*----------------------------------------------------------------------------*/
/** Read a ASCII number from a PGM file.
 */
static int get_num(FILE * f)
{
  int num,c;

  while(isspace(c=getc(f)));
  if(!isdigit(c)) error("Error: corrupted PGM file.");
  num = c - '0';
  while( isdigit(c=getc(f)) ) num = 10 * num + c - '0';
  if( c != EOF && ungetc(c,f) == EOF )
    error("Error: unable to 'ungetc' while reading PGM file.");

  return num;
}

/*----------------------------------------------------------------------------*/
/** Read a PGM file into an double image.
    If the name is "-" the file is read from standard input.
 */
static double * read_pgm_image_double(int * X, int * Y, char * name)
{
  FILE * f;
  int c,bin;
  int xsize,ysize,depth,x,y;
  double * image;

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdin;
  else f = fopen(name,"rb");
  if( f == NULL ) error("Error: unable to open input image file.");

  /* read header */
  if( getc(f) != 'P' ) error("Error: not a PGM file!");
  if( (c=getc(f)) == '2' ) bin = FALSE;
  else if( c == '5' ) bin = TRUE;
  else error("Error: not a PGM file!");
  skip_whites_and_comments(f);
  xsize = get_num(f);            /* X size */
  if(xsize<=0) error("Error: X size <=0, invalid PGM file\n");
  skip_whites_and_comments(f);
  ysize = get_num(f);            /* Y size */
  if(ysize<=0) error("Error: Y size <=0, invalid PGM file\n");
  skip_whites_and_comments(f);
  depth = get_num(f);            /* depth */
  if(depth<=0) fprintf(stderr,"Warning: depth<=0, probably invalid PGM file\n");
  /* white before data */
  if(!isspace(c=getc(f))) error("Error: corrupted PGM file.");

  /* get memory */
  image = (double *) calloc( (size_t) (xsize*ysize), sizeof(double) );
  if( image == NULL ) error("Error: not enough memory.");

  /* read data */
  for(y=0;y<ysize;y++)
    for(x=0;x<xsize;x++)
      image[ x + y * xsize ] = bin ? (double) getc(f)
                                   : (double) get_num(f);

  /* close file if needed */
  if( f != stdin && fclose(f) == EOF )
    error("Error: unable to close file while reading PGM file.");

  /* return image */
  *X = xsize;
  *Y = ysize;
  return image;
}

/*----------------------------------------------------------------------------*/
/** Write an int image into a PGM file.
    If the name is "-" the file is written to standard output.
 */
static void write_pgm_image_int(int * image, int xsize, int ysize, char * name)
{
  FILE * f;
  int x,y,n,v,max,min;

  /* check input */
  if( image == NULL || xsize <= 0 || ysize <= 0 )
    error("Error: invalid input image to write_pgm_image_int.");

  /* check min and max values */
  max = min = 0;
  for(y=0; y<ysize; y++)
    for(x=0; x<xsize; x++)
      {
        v = image[ x + y * xsize ];
        if( v > max ) max = v;
        if( v < min ) min = v;
      }
  if( min < 0 ) fprintf(stderr,
        "Warning: write_pgm_image_int: negative values in '%s'.\n",name);
  if( max > 65535 ) fprintf(stderr,
        "Warning: write_pgm_image_int: values exceeding 65535 in '%s'.\n",name);

  /* open file */
  if( strcmp(name,"-") == 0 ) f = stdout;
  else f = fopen(name,"w");
  if( f == NULL ) error("Error: unable to open output image file.");

  /* write header */
  fprintf(f,"P2\n");
  fprintf(f,"%d %d\n",xsize,ysize);
  fprintf(f,"%d\n",max);

  /* write data */
  for(n=0,y=0; y<ysize; y++)
    for(x=0; x<xsize; x++)
      {
        fprintf(f,"%d ",image[ x + y * xsize ]);
        if(++n==8)  /* lines should not be longer than 70 characters  */
          {
            fprintf(f,"\n");
            n = 0;
          }
      }

  /* close file if needed */
  if( f != stdout && fclose(f) == EOF )
    error("Error: unable to close file while writing PGM file.");
}


/*----------------------------------------------------------------------------*/
/*----------------------------- Write EPS File -------------------------------*/
/*----------------------------------------------------------------------------*/
/** Write line segments into an EPS file.
    If the name is "-" the file is written to standard output.

    According to

      Adobe "Encapsulated PostScript File Format Specification",
      Version 3.0, 1 May 1992,

    and

      Adobe "PostScript(R) LANGUAGE REFERENCE", third edition, 1999.
 */
static void write_eps( double * segs, int n, int dim,
                       char * filename, int xsize, int ysize, double width )
{
  FILE * eps;
  int i;

  /* check input */
  if( segs == NULL || n < 0 || dim <= 0 )
    error("Error: invalid line segment list in write_eps.");
  if( xsize <= 0 || ysize <= 0 )
    error("Error: invalid image size in write_eps.");

  /* open file */
  if( strcmp(filename,"-") == 0 ) eps = stdout;
  else eps = fopen(filename,"w");
  if( eps == NULL ) error("Error: unable to open EPS output file.");

  /* write EPS header */
  fprintf(eps,"%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(eps,"%%%%BoundingBox: 0 0 %d %d\n",xsize,ysize);
  fprintf(eps,"%%%%Creator: LSD, Line Segment Detector\n");
  fprintf(eps,"%%%%Title: (%s)\n",filename);
  fprintf(eps,"%%%%EndComments\n");

  /* write line segments */
  for(i=0;i<n;i++)
    {
      fprintf( eps,"newpath %f %f moveto %f %f lineto %f setlinewidth stroke\n",
               segs[i*dim+0],
               (double) ysize - segs[i*dim+1],
               segs[i*dim+2],
               (double) ysize - segs[i*dim+3],
               width <= 0.0 ? segs[i*dim+4] : width );
    }

  /* close EPS file */
  fprintf(eps,"showpage\n");
  fprintf(eps,"%%%%EOF\n");
  if( eps != stdout && fclose(eps) == EOF )
    error("Error: unable to close file while writing EPS file.");
}


/*----------------------------------------------------------------------------*/
/*----------------------------- Write SVG File -------------------------------*/
/*----------------------------------------------------------------------------*/
/** Write line segments into a SVG file.
    If the name is "-" the file is written to standard output.
*/
static void write_svg( double * segs, int n, int dim,
                       char * filename, int xsize, int ysize, double width )
{
  FILE * svg;
  int i;

  /* check input */
  if( segs == NULL || n < 0 || dim <= 0 )
    error("Error: invalid line segment list in write_svg.");
  if( xsize <= 0 || ysize <= 0 )
    error("Error: invalid image size in write_svg.");

  /* open file */
  if( strcmp(filename,"-") == 0 ) svg = stdout;
  else svg = fopen(filename,"w");
  if( svg == NULL ) error("Error: unable to open SVG output file.");

  /* write SVG header */
  fprintf(svg,"<?xml version=\"1.0\" standalone=\"no\"?>\n");
  fprintf(svg,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
  fprintf(svg," \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
  fprintf(svg,"<svg width=\"%dpx\" height=\"%dpx\" ",xsize,ysize);
  fprintf(svg,"version=\"1.1\"\n xmlns=\"http://www.w3.org/2000/svg\" ");
  fprintf(svg,"xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n");

  /* write line segments */
  for(i=0;i<n;i++)
    {
      fprintf(svg,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" ",
              segs[i*dim+0],segs[i*dim+1],
              segs[i*dim+2],segs[i*dim+3]);
      fprintf(svg,"stroke-width=\"%f\" stroke=\"black\" />\n",
              width <= 0.0 ? segs[i*dim+4] : width);
    }

  /* close SVG file */
  fprintf(svg,"</svg>\n");
  if( svg != stdout && fclose(svg) == EOF )
    error("Error: unable to close file while writing SVG file.");
}


/*----------------------------------------------------------------------------*/
/*                                    Main                                    */
/*----------------------------------------------------------------------------*/
/** Main function call
 */
int main(int argc, char ** argv)
{
  struct arguments * arg = process_arguments(USE,argc,argv);
  FILE * output;
  double * image;
  int X,Y;
  double * segs;
  int n;
  int dim = 7;
  int * region;
  int regX,regY;
  int i,j;

  /* read input file */
  image = read_pgm_image_double(&X,&Y,get_str(arg,"in"));

  /* execute LSD */
  segs = LineSegmentDetection( &n, image, X, Y,
                               get_double(arg,"scale"),
                               get_double(arg,"sigma_coef"),
                               get_double(arg,"quant"),
                               get_double(arg,"ang_th"),
                               get_double(arg,"log_eps"),
                               get_double(arg,"density_th"),
                               get_int(arg,"n_bins"),
                               is_assigned(arg,"reg") ? &region : NULL,
                               &regX, &regY );

  /* output */
  if( strcmp(get_str(arg,"out"),"-") == 0 ) output = stdout;
  else output = fopen(get_str(arg,"out"),"w");
  if( output == NULL ) error("Error: unable to open ASCII output file.");
  for(i=0;i<n;i++)
    {
      for(j=0;j<dim;j++)
        fprintf(output,"%f ",segs[i*dim+j]);
      fprintf(output,"\n");
    }
  if( output != stdout && fclose(output) == EOF ) /* close file if needed */
    error("Error: unable to close file while output file.");

  /* store region output if needed */
  if(is_assigned(arg,"reg"))
    {
      write_pgm_image_int(region,regX,regY,get_str(arg,"reg"));
      free( (void *) region );
    }

  /* create EPS output if needed */
  if(is_assigned(arg,"epsfile"))
    write_eps(segs,n,dim,get_str(arg,"epsfile"),X,Y,get_double(arg,"width"));

  /* create SVG output if needed */
  if(is_assigned(arg,"svgfile"))
    write_svg(segs,n,dim,get_str(arg,"svgfile"),X,Y,get_double(arg,"width"));

  /* free memory */
  free( (void *) image );
  free( (void *) segs );
  free_arguments(arg);

  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
