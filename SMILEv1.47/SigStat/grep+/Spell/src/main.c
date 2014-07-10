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

#include <grep+.h>



/******************************************************************************/
/******************************************************************************/
/********************************** MAIN **************************************/
/******************************************************************************/
/******************************************************************************/
int main(int argc, char **argv)
{
Criteres    criteres;
FILE        *f=NULL;
unsigned int nbocc, nboccex;
Arbre       b;
char        prout[100];

initCriteres(&criteres);

/* setBloc(&criteres,2); */
/* setLongueurBloc(&criteres,0,5); */
/* setLongueurBloc(&criteres,1,5); */
/* setErreurGlobal(&criteres,0); */
/* setErreurBloc(&criteres,0,0); */
/* setErreurBloc(&criteres,1,0); */
/* setSaut(&criteres,0,14,16); */

setBloc(&criteres,1);
setLongueurBloc(&criteres,0,11);
setErreurGlobal(&criteres,1);
if(valideCriteres(&criteres) == FAUX)
    return 1;

f = fopen("pos","w");
if(f==NULL)
    printf("CHIIIIIIIER\n");

creeArbreSuffixeFromFile(&b, "seq2",maxLongMod(criteres));
chercheMot(b,&criteres, "AAAAA",&nbocc, &nboccex, f);
printf("%d %d!\n",nbocc, nboccex);
libereArbreSuffixeFromFile(b);

/* creeArbreSuffixeFromFile(&b, "ficseq",maxLongMod(criteres)); */
/* chercheMot(b,&criteres, "AAAAAC_AGTGTT",&nbocc, &nboccex, f); */
/* printf("%d %d!\n",nbocc, nboccex); */
/* libereArbreSuffixeFromFile(b); */

return(0); 
}
