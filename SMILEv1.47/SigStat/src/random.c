/******************************************************************************/
/* SMILE v1.47 - Extraction of structured motifs common to several sequences  */
/* Copyright (C) 2004 L.Marsan (lama -AT- prism.uvsq.fr)                      */
/*                                                                            */
/* This program is free software; you can redistribute it and/or              */
/* modify it under the terms of the GNU General Public License                */
/* as published by the Free Software Foundation; either version 2             */
/* of the License, or (at your option) any later version.                     */
/*                                                                            */
/* This program is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with this program; if not, write to the Free Software                */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
/******************************************************************************/

/* (C) by "coward" <coward -AT- ii.uib.no> from 1996-1999 */
/*   modified and extended by Laurent Marsan <lama@prism.uvsq.fr> 1999-2004.  */


#include <stdio.h>

#define MULTIPLIER 69069
#define SHIFT          1
#define MODULUS    256*256*256*128
#define INVMOD     ( (double) 1 / ((double) MODULUS)) / ((double) 2)

#define DEFAULTSEED 1000001


static unsigned int seed;

/* changeseed **************************************************************************************************************

DESCRIPTION: initiate the static variable seed, used as seed value in functions generating random numbers

SIDE EFFECTS: the value of the static variable seed is changed

RETURN VALUE: the old value of the static variable seed

***************************************************************************************************************************/
unsigned int changeseed(unsigned int Iseed)
{
  unsigned int oldseed;
  
  oldseed = seed;
  seed = Iseed;
  
  return oldseed;
}


/* unif01 ******************************************************************************************************************

DESCRIPTION: Draw a random number from a uniform distribution on the interval [0,1]

SIDE EFFECTS: the value of the static variable seed is changed

RETURN VALUE: the random number

***************************************************************************************************************************/
double unif01(void)
{
  double random;
  
  seed = MULTIPLIER * seed + SHIFT;
  random = ((double) seed) * INVMOD;

  return random;
}


/* readseed  **************************************************************************************************************

DESCRIPTION: read new seed value from file, use standard value if nonexistant

SIDE EFFECTS: the value of the static variable seed is changed

RETURN VALUE: the seed value

***************************************************************************************************************************/
unsigned int readseed(char filename[])
{
  FILE *fp;

  seed = DEFAULTSEED;
  fp = fopen(filename, "r");
  if (fp != NULL)
  { fscanf(fp, "%u", &seed);
    fclose(fp);
  }
  return seed;
}


/* writeseed  *************************************************************************************************************

DESCRIPTION: write seed value to seed file

SIDE EFFECTS: none

RETURN VALUE: 1 if OK (return value from fprintf)

*****************************************************************************************************************************/
int writeseed(char filename[])
{
  FILE *fp;
  int r;
  
  fp = fopen(filename, "w");
  if (fp != NULL)
    r = fprintf(fp, "%u\n", seed);
  else
    r = -1;
  fclose(fp);
  return r;
}





