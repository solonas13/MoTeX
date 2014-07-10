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
FILE**  recInitFiles(FILE **f, int bloc, P_Criteres cr, char *buf, int *tab,
            FILE *namefile);

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
        fatalError("setCompo: cannot allocate 'cr->compobloc'\n");
    }

/* Initialisations                                                            */
cr->flag_compo  = FAUX;
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
/* initFiles                                                                  */
/******************************************************************************/
int    initFiles(FILE **f, char *nom, P_Criteres cr)
{
char    buf[500];
int     tab[100];
FILE *  namefile;

strcpy(buf, nom);

namefile    = fopen(NAMEFILE, "w");
if(namefile == NULL)
    return FAUX;

if(recInitFiles(f, 0, cr, buf,tab, namefile) == NULL)
    {
    fclose(namefile);
    return FAUX;
    }

fclose(namefile);
return VRAI;
}

/******************************************************************************/
/* recInitFiles                                                               */
/******************************************************************************/
FILE**  recInitFiles(FILE **f, int bloc, P_Criteres cr, char *buf, int tab[100],
            FILE *namefile)
{
int i;
/* int tmpi,tab2[100]; */
char tmp[20];
char *posfin;

posfin = buf + strlen(buf);
for(i=cr->saut[bloc].min+cr->delta[bloc]; i+cr->delta[bloc]<=cr->saut[bloc].max; i++)
    {
    if(bloc < cr->bloc-1)
        {
        tab[bloc]=i;
        sprintf(tmp,"[%d-%d]",i-cr->delta[bloc], i+cr->delta[bloc]);
        strcat(buf,tmp);

        f = recInitFiles(f, bloc+1, cr, buf, tab, namefile);

        if(f==NULL)
            return(NULL);
        }
    else
        {
/*          printf("\n%s =%d\n",buf,tmpi =delta2File(tab, cr)); */
/*      if ( ! file2Delta(tmpi, tab2, cr) ) */
/*          printf("bug?CDS<A>CVSD<A\n"); */
/*      for(i=0; i != cr->bloc-1; i++) */
/*          printf("%d ", tab2[i]); */
/*      printf("\n"); */
        *f  = fopen(buf,"w");

        if(*f == NULL)
            return NULL;

        fprintf(namefile, "%s\n",buf);

/* Insertion de l'espace necessaire en tete de fichier pour y mettre la       */
/* ligne d'informations apres extraction                                      */
/* J'ecris en tout 80 * 3 = 240 espaces et 80 '='  ==> 320 caracteres         */
        for(i=0; i!=3; i++)
            {
            fprintf(*f,"                                        ");
            fprintf(*f,"                                       \n");
            }
        fprintf(*f,"========================================");
        fprintf(*f,"=======================================\n");
        
        return (++f);
        }
    *posfin = '\0';
    }
return(f);
}


/******************************************************************************/
/* delta2File                                                                 */
/******************************************************************************/
int delta2File(LongSeq *deltatab, P_Criteres cr)
{
int pos=0, i;

pos =  deltatab[0] - cr->saut[0].min - cr->delta[0];

for(i=1; i != cr->bloc-1; i++)
    {
    pos     *= cr->saut[i].max - cr->saut[i].min + 1-cr->delta[i]*2;
    pos     += deltatab[i] - cr->saut[i].min - cr->delta[i];
    }

return pos;
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
/* file2Delta                                                                 */
/******************************************************************************/
/* A REFAIRE!! Vaut mieux avoir un tableau ou chaque colonne pointe vers
 *  le deltatab associe. */
int file2Delta(int pos, LongSeq *deltatab, P_Criteres cr)
{
int i, tmp;

for(i=0; i != cr->bloc-1; i++)
    {
    tmp     = cr->saut[i].max-cr->saut[i].min+1-cr->delta[i]*2;
    deltatab[cr->bloc-2-i]  = pos % tmp;
    pos     /= tmp;
    }

if(pos != 0)
    {
    warning("file2Delta: conversion error");
    return(0);
    }
return (1);
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
cr->delta       = (LongSeq *)       malloc(bloc*sizeof(LongSeq));
cr->palindrom   = (LongSeq *)       malloc(bloc*sizeof(LongSeq));
if(!cr->maxerrblocs || !cr->longbloc || !cr->saut || !cr->flag_compobloc
        || !cr->compobloc || !cr->delta || !cr->palindrom)
    fatalError("criteres.h: allocBloc: allocation error\n");

/* Initialisations                                                            */
for (i = 0; i != bloc; i++)
    {
    cr->saut[i].min     = cr->saut[i].max = 0;
    cr->maxerrblocs[i]  = -1;
    cr->longbloc[i].min = cr->longbloc[i].max = -1;
    cr->delta[i]        = 0;
    cr->palindrom[i]    = -1;
    }

cr->flag_palindrom  = FAUX;
}
