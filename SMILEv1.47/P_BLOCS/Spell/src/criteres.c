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

#include <criteres.h>

/* From alphabet.c                                                            */
extern  int     nbSymbMod;

int             **code2Sauts;

/******************************************************************************/
/* FONCTIONS PRIVEES                                                          */
/******************************************************************************/
int     recFillTab(int bloc, P_Criteres cr, int *nbcodes, int **code2Sauts);



/******************************************************************************/
/* setCompoPal                                                                */
/******************************************************************************/
void    setCompoPal(P_Criteres cr, char **argv, int argc)
{
int i,j,bloc,nbsymb;

if( (cr->compo = malloc(nbSymbMod*sizeof(LongSeq))) == NULL)
        fatalError("initCriteres: cannot allocate 'cr->compo'\n");

for(i=0;i<cr->bloc;i++)
    {
    cr->compobloc[i]    = (LongSeq *) malloc(nbSymbMod*sizeof(LongSeq));
    if(!(cr->compobloc[i]))
        fatalError("setCompoPal: cannot allocate 'cr->compobloc[i]'\n");
    }

/* Initialisations                                                            */
cr->flag_compo      = FAUX;
for (i = 0; i != nbSymbMod; i++)
    cr->compo[i] = -1;
for (i = 0; i != cr->bloc; i++)
    {
    cr->flag_compobloc[i] = FAUX;
    for (j = 0; j != nbSymbMod; j++)
        cr->compobloc[i][j] = -1;
    }


while(argc>0 && (**argv!='p'))
    {
    bloc    = atoi(*argv);
    argc--; argv++;
    nbsymb  = atoi(*argv);
    argc--; argv++;

    for(i=0; i!=nbsymb; i++)
        {
        j   = str2nummod(*argv);
        argc--; argv++;

        if( j == -1)
            {
            fprintf(stderr,
                    "> Warning: composition in '%s' ignored, symbol is not in the models alphabet.\n",
                    *(argv-1));
            argc--; argv++;
            continue;
            }

        if(bloc == 0)
            {
            cr->flag_compo  = VRAI;
            cr->compo[j]    = atoi(*argv);
/*             printf("compo glob %s %d\n",*(argv-1),cr->compo[j]); */
            }
        else
            {
            cr->flag_compobloc[bloc-1]= VRAI;
            cr->compobloc[bloc-1][j]  = atoi(*argv);
/*             printf("compo bloc %s %d\n",*(argv-1),cr->compobloc[j]); */
            }

        argc--; argv++;
        }
    }

/* S'il n'y a pas de palindromes                                              */
if(argc==0)
    return;

cr->flag_palindrom = VRAI;

while(argc>0)
    {
    sscanf(*argv, "p%d/%d",&i, &j);
    cr->palindrom[i-1]   = j-1;
    argc--;
    argv++;
    }
}

/******************************************************************************/
/* addSaut2Code                                                               */
/******************************************************************************/
int addSaut2Code(int oldcode, LongSeq saut, LongSeq curbloc,  P_Criteres cr)
{
#if DEBUG_SAUT
printf("addSaut2Code: Je recois %d",oldcode);
#endif
if (curbloc == 0)
    {
#if DEBUG_SAUT
    printf(" et renvoie %d\n",saut -  cr->saut[0].min);
#endif
    return (saut -  cr->saut[0].min);
    }

oldcode *= cr->saut[curbloc].max - cr->saut[curbloc].min + 1;
oldcode += saut - cr->saut[curbloc].min;

#if DEBUG_SAUT
printf(" et renvoie %d\n",oldcode);
#endif
return oldcode;
}

/******************************************************************************/
/* initTabSauts                                                               */
/******************************************************************************/
void initTabSauts(P_Criteres cr)
{
int bloc = cr->bloc, **i,j, nbcodes;

nbcodes = 1;
for(j=0; j != bloc-1; j++)
    nbcodes *= cr->saut[j].max - cr->saut[j].min +1;


if ( (code2Sauts    = (int **) malloc(nbcodes * sizeof(int *)) ) == NULL)
        fatalError("criteres.c: initTabSauts: cannot allocate 'code2Sauts'\n");

for(j=0, i=code2Sauts; j != nbcodes; j++,i++)
    if ( (*i = (int *) malloc((bloc-1) * sizeof(int)) ) == NULL )
        fatalError("criteres.c: initTabSauts: cannot allocate 'code2Sauts[j]'\n");

recFillTab(0, cr, &nbcodes, code2Sauts);
}

/******************************************************************************/
/* recFillTab                                                                 */
/******************************************************************************/
int recFillTab(int bloc, P_Criteres cr, int *nbcodes, int **code2Sauts)
{
int i,j,k,pos=0,a;

if(bloc != cr->bloc-2)
    a = recFillTab(bloc+1, cr, nbcodes, code2Sauts);
else
    a = 1;

*nbcodes /= cr->saut[bloc].max - cr->saut[bloc].min +1;
for(i=0; i!=*nbcodes; i++)
    for(j=cr->saut[bloc].min; j!=cr->saut[bloc].max+1; j++)
        for(k=0; k!=a; k++)
            {
            code2Sauts[pos][bloc] = j;
            pos++;
            }
printf("\n");
return( a*(cr->saut[bloc].max - cr->saut[bloc].min +1));
}

/******************************************************************************/
/* allocBloc                                                                  */
/******************************************************************************/
void allocBloc(P_Criteres cr, int bloc)
{
int i;

cr->maxerrblocs = (LongSeq *)       malloc(bloc*sizeof(LongSeq));
cr->longbloc    = (Fourchette *)    malloc(bloc*sizeof(Fourchette));
cr->saut        = (Fourchette *)    malloc(bloc*sizeof(Fourchette));
cr->flag_compobloc = (Flag *)       malloc(bloc*sizeof(Flag));
cr->compobloc   = (LongSeq **)      malloc(bloc*sizeof(LongSeq *));
cr->palindrom   = (LongSeq *)       malloc(bloc*sizeof(LongSeq));
if(!cr->maxerrblocs || !cr->longbloc || !cr->saut || !cr->flag_compobloc
        || !cr->compobloc || !cr->palindrom)
    fatalError("criteres.h: allocBloc: allocation error\n");

/* Initialisations                                                            */
for (i = 0; i != bloc; i++)
    {
    cr->saut[i].min     = cr->saut[i].max = 0;
    cr->maxerrblocs[i]  = -1;
    cr->longbloc[i].min = cr->longbloc[i].max = -1;
    cr->palindrom[i]    = -1;
    }

cr->flag_palindrom  = FAUX;
}
