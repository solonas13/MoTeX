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

#include "criteres.h"

/******************************************************************************/
/* FONCTIONS PRIVEES                                                          */
/******************************************************************************/
void    libereTabSauts(P_Criteres);
int     recFillTab(int bloc, Criteres cr, int *nbcodes, int **code2Sauts);
Flag    allocBloc(P_Criteres cr, NbBlocs bloc);
/******************************************************************************/



/******************************************************************************/
/* Fonctions de positionnement/lecture des variables de classe                */
/******************************************************************************/
/* setBloc                                                                    */
Flag    setBloc(P_Criteres cr, NbBlocs bloc)
{
if( !allocBloc(cr, bloc))
    return FAUX;

cr->bloc    = bloc;

return VRAI;
}

/* setLongueurBloc                                                         */
Flag    setLongueurBloc(P_Criteres cr, NbBlocs num_bloc, LongSeq lon)
{
if (num_bloc >= cr->bloc || num_bloc < 0)
    return FAUX;

cr->longbloc[(int)num_bloc].max  = lon;

return VRAI;
}

/* setErreurGlobal                                                            */
void setErreurGlobal(P_Criteres cr, LongSeq err)
{
cr->maxerr    = err;
}

/* setErreurBloc                                                              */
Flag    setErreurBloc(P_Criteres cr, NbBlocs num_bloc, LongSeq err)
{
if (num_bloc >= cr->bloc || num_bloc < 0)
    return FAUX;

cr->maxerrblocs[(int)num_bloc]  = err;

return VRAI;
}

/* setSaut                                                                    */
Flag    setSaut(P_Criteres cr, NbBlocs num_bloc, LongSeq min, LongSeq max)
{
if (num_bloc >= cr->bloc-1 || num_bloc < 0)
    return FAUX;

cr->saut[(int)num_bloc].min  = min;
cr->saut[(int)num_bloc].max  = max;

return VRAI;
}

/* getBloc                                                                    */
NbBlocs getBloc(Criteres cr)
{return cr.bloc;}

/* getLongueurBloc                                                         */
LongSeq getLongueurBloc(Criteres cr,NbBlocs num_bloc)
{
if (num_bloc >= cr.bloc || num_bloc < 0)
    return -1;

return cr.longbloc[(int)num_bloc].max;
}

/* getErreurGlobal                                                            */
LongSeq getErreurGlobal(Criteres cr)
{return cr.maxerr;}

/* getErreurBloc                                                              */
LongSeq getErreur(Criteres cr, NbBlocs num_bloc)
{
if (num_bloc >= cr.bloc || num_bloc < 0)
    return -1;

return cr.maxerrblocs[(int)num_bloc];
}

/* getSaut                                                                    */
Fourchette getSaut(Criteres cr, NbBlocs num_bloc)
{
return cr.saut[(int)num_bloc];
}

/* maxLongMod                                                                 */
int     maxLongMod(Criteres cr)
{
int i, max;

max = cr.longbloc[0].max;

for(i=0;i<cr.bloc-1;i++)
    {
    max += cr.saut[i].max;
    max += cr.longbloc[i+1].max;
    }

return max;
}

/******************************************************************************/
/* afficheCriteres                                                            */
/******************************************************************************/
void    afficheCriteres(Criteres cr, FILE *f)
{
int i;

fprintf(f, "Blocs               %d\n",cr.bloc);
fprintf(f, "Longueur min totale %d\n",cr.longmod.min);
fprintf(f, "Longueur max totale %d\n",cr.longmod.max);
fprintf(f, "Erreurs globales    %d\n",cr.maxerr);
if(cr.bloc<=1)
    {
    fprintf(f, "\n\n");
    return;
    }

for(i=0;i<cr.bloc;i++)
    {
    fprintf(f, "\nBLOC %d\n",i);
    fprintf(f, "Longueur min    %d\n", cr.longbloc[i].min);
    fprintf(f, "Longueur max    %d\n", cr.longbloc[i].max);
    fprintf(f, "Erreurs         %d\n", cr.maxerrblocs[i]);

    if(i!=cr.bloc-1)
        {
        fprintf(f, "Saut min        %d\n",cr.saut[i].min);
        fprintf(f, "Saut max        %d\n",cr.saut[i].max);
        }
    }
fprintf(f, "\n\n");
}

/******************************************************************************/
/* verifCriteres                                                              */
/******************************************************************************/
Flag    verifCriteres(Criteres cr)
{
int i;

if (cr.maxerr==-1)
    {
    fprintf(stderr,
    "verifCriteres: global errors number needs to be fixed (setErreurGlobal\n");
    return FAUX;
    }

if(cr.bloc <= 0)
    {
    fprintf(stderr,
        "verifCriteres: boxes number has to be fixed (setBloc)\n");
    return FAUX;
    }

if(cr.longbloc[0].max <0 )
    {
    fprintf(stderr,
    "verifCriteres: first box's max length must be fixed (setLongueurBloc)\n");
    return FAUX;
    }

if(cr.bloc == 1)
    return VRAI;

for(i=0;i<cr.bloc;i++)
    {
    if(cr.maxerrblocs[i]==-1)
        {
        fprintf(stderr,
        "verifCriteres: max errors for box %d has to be fixed (setErreurBloc)\n", i+1);
        return FAUX;
        }

    if(cr.longbloc[i].max<0)
        {
        fprintf(stderr,
        "verifCriteres: max length for box %d has to be fixed (setLongueurBloc)\n", i+1);
        return FAUX;
        }

    if(i!=cr.bloc-1)
        {
        if(cr.saut[i].min==-1)
            {
            fprintf(stderr,
            "verifCriteres: min spacer length after box %d has to be fixed (setSaut)\n", i+1);
            return FAUX;
            }

        if(cr.saut[i].max<=0)
            {
            fprintf(stderr,
            "verifCriteres: max spacer length after box %d has to be fixed (setSaut)\n", i+1);
            return FAUX;
            }

        if(cr.saut[i].min>cr.saut[i].max)
            {
            fprintf(stderr,
            "verifCriteres: max spacer length has to be greater than min spacer length!\n");
            return FAUX;
            }
        }
    }
return VRAI;
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

oldcode *= cr->saut[(int)curbloc].max - cr->saut[(int)curbloc].min + 1;
oldcode += saut - cr->saut[(int)curbloc].min;

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

if(cr->code2Sauts != NULL)
    libereTabSauts(cr);

nbcodes = 1;
for(j=0; j != bloc-1; j++)
    nbcodes *= cr->saut[j].max - cr->saut[j].min +1;

if ( (cr->code2Sauts    = (int **) malloc(nbcodes * sizeof(int *)) ) == NULL)
        fatalError("criteres.c: initTabSauts: cannot allocate 'code2Sauts'\n");

for(j=0, i=cr->code2Sauts; j != nbcodes; j++,i++)
    if ( (*i = (int *) malloc((bloc-1) * sizeof(int)) ) == NULL )
        fatalError("criteres.c: initTabSauts: cannot allocate 'code2Sauts[j]'\n");

recFillTab(0, *cr, &nbcodes, cr->code2Sauts);
}

/******************************************************************************/
/* libereTabSauts                                                             */
/******************************************************************************/
void libereTabSauts(P_Criteres cr)
{
int bloc = cr->bloc, **i,j, nbcodes;

nbcodes = 1;
for(j=0; j != bloc-1; j++)
    nbcodes *= cr->saut[j].max - cr->saut[j].min +1;


for(j=0, i=cr->code2Sauts; j != nbcodes; j++,i++)
    free(*i);

free(cr->code2Sauts);

cr->code2Sauts  = NULL;
}


/******************************************************************************/
/* recFillTab                                                                 */
/******************************************************************************/
int recFillTab(int bloc, Criteres cr, int *nbcodes, int **code2Sauts)
{
int i,j,k,pos=0,a;

if(bloc != cr.bloc-2)
    a = recFillTab(bloc+1, cr, nbcodes, code2Sauts);
else
    a = 1;

*nbcodes /= cr.saut[bloc].max - cr.saut[bloc].min +1;
for(i=0; i!=*nbcodes; i++)
    for(j=cr.saut[bloc].min; j!=cr.saut[bloc].max+1; j++)
        for(k=0; k!=a; k++)
            {
            code2Sauts[pos][bloc] = j;
            pos++;
            }
printf("\n");
return( a*(cr.saut[bloc].max - cr.saut[bloc].min +1));
}

/******************************************************************************/
/* allocBloc                                                                  */
/******************************************************************************/
Flag allocBloc(P_Criteres cr, NbBlocs bloc)
{
int i;

if(cr->maxerrblocs != NULL)
    {
    libereTabSauts(cr);
    free(cr->longbloc);
    if(cr->bloc > 1)
        {
        free(cr->maxerrblocs);
        free(cr->saut);
        }
    }

cr->longbloc    = (Fourchette *)       malloc(bloc*sizeof(Fourchette));
if(!cr->longbloc)
    return FAUX;

cr->maxerrblocs = (LongSeq *)       malloc(bloc*sizeof(LongSeq));
if( !cr->maxerrblocs )
    return FAUX;

if(bloc > 1)
    {
    cr->saut        = (Fourchette *)    malloc((bloc-1)*sizeof(Fourchette));

    if( !cr->saut )
        return FAUX;
    }

/* Initialisations blocs                                                      */
for (i = 0; i != bloc; i++)
    {
    cr->maxerrblocs[i]  = -1;
    cr->longbloc[i].max = -1;
    cr->longbloc[i].min = -1;
    if(i!=bloc-1)
        cr->saut[i].min = cr->saut[i].max = -1;
    }

/* Initialisations generales                                                  */
cr->maxerr = cr->bloc = cr->multiblocs = -1;
cr->code2Sauts = NULL;
cr->nbtotseq   = -1;
cr->longmod.max = -1;
cr->longmod.min = -1;

return VRAI;
}

/******************************************************************************/
/* initCriteres                                                               */
/******************************************************************************/
void    initCriteres(P_Criteres cr)
{
cr->maxerr      = -1;
cr->bloc        = -1;
cr->maxerrblocs = NULL;
cr->longbloc    = NULL;
cr->saut        = NULL;
cr->multiblocs  = 0;
cr->code2Sauts  = NULL;
cr->nbtotseq    = 0;
cr->longmod.max  = 0;
cr->longmod.min  = 0;
}

/******************************************************************************/
/* chargeCriteres                                                             */
/******************************************************************************/
Flag    chargeCriteres(P_Criteres cr, char *line)
{
int     i,tmp, tmp2;

tmp = atoi(strtok(line,"% ")); /* Nb blocs */
if(tmp < 1)
    return FAUX;
setBloc(cr, tmp);

strtok(NULL,"%/ "); /* Quorum */
cr->nbtotseq    = atoi(strtok(NULL," /")); /* Nb total sequences */
cr->nbsymb      = atoi(strtok(NULL," ")); /* Nb symboles */
cr->longmod.min = atoi(strtok(NULL," ")); /* L min */

tmp=atoi(strtok(NULL," ")); /* L max */
cr->longmod.max = tmp;
if (cr->bloc == 1)
    setLongueurBloc(cr, 0, tmp);

setErreurGlobal(cr, atoi(strtok(NULL," "))); /* Err glob */

if(cr->bloc != 1)
    {
    for(i=0; i!=cr->bloc; i++)
        {
        strtok(NULL," "); /* L min bloc i */
        tmp=atoi(strtok(NULL," ")); /* L max bloc i */
        setLongueurBloc(cr, i, tmp);
    
        setErreurBloc(cr, i, atoi(strtok(NULL," "))); /* Err bloc i */
    
        if(i != cr->bloc-1)
            {
            tmp =atoi(strtok(NULL," ")); /* Saut min bloc i */
            tmp2=atoi(strtok(NULL," ")); /* Saut max bloc i */
    
            setSaut(cr, i, tmp, tmp2);
            }
        }
    }

strcpy(cr->ficalph, strtok(NULL," "));  /* Fichier alphabet des modeles */
strcpy(cr->alphaseq, strtok(NULL," ")); /* Alphabet des sequences */
/* afficheCriteres(*cr); */
return VRAI;
}

