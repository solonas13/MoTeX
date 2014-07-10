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

#include <global.h>

int symb2int[256]={-1};

/******************************************************************************/
/* fatalError                                                                 */
/******************************************************************************/
/* Gestion des erreurs FATALES!                                               */
/******************************************************************************/
void fatalError(char *msg)
{
fprintf(stderr,"Error: %s\n",msg);
exit(1);
}

/******************************************************************************/
/* warning                                                                    */
/******************************************************************************/
void warning(char *msg)
{
fprintf(stderr,"Warning: %s\n",msg);
}

/******************************************************************************/
/* initEntiers                                                                */
/******************************************************************************/
void initEntiers(void)
{
int i;

for(i=0;i<256;i++)
    symb2int[i] = -1;

symb2int[(int)'A']  = 0;
symb2int[(int)'C']  = 1;
symb2int[(int)'G']  = 2;
symb2int[(int)'T']  = 3;
symb2int[(int)'$']  = 4;
}

/******************************************************************************/
/* entiers                                                                    */
/******************************************************************************/
/* Conversion des symboles vers les entiers                                   */
/******************************************************************************/
int entiers(char c)
{
int ret;

ret = symb2int[(int)c];
if(ret == -1)
    {
    fprintf(stderr,"-- %c --",c);
    fatalError("Unknown character, non A,C,G,T or $");
    }

return ret;
}


/******************************************************************************/
/* entree                                                                     */
/******************************************************************************/
void entree(void)
{
printf("\n-- Type ENTER\n");
fflush(stdin);
getchar();
}
