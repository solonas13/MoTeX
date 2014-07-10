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

#ifndef _CRITERES_H
#define _CRITERES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"

/******************************************************************************/
/* STRUCTURE DE STOCKAGE DES CRITERES DE RECHERCHE                            */
/******************************************************************************/
typedef struct  struct_fourchette
    {
    LongSeq     min;
    LongSeq     max;
    }   Fourchette, *P_Fourchette;

/******************************************************************************/
/* Seuls sont necessaires pour la recherche: maxerr, bloc, maxerrbloc,        */
/* longbloc.max, saut, multiblocs, code2Sauts. Le reste est superflu          */
/* (utilise par sigstat ou converter) et n'a donc pas de methodes get/set.    */
/******************************************************************************/
typedef struct  struct_criteres
    {
    Fourchette  *longbloc;
    Fourchette  *saut;
    LongSeq     *maxerrblocs;
    int         **code2Sauts;
    char        alphaseq[128];
    char        ficalph[128];
    long int    nbsymb;
    Fourchette  longmod;
    LongSeq     maxerr;
    NbSeq       nbtotseq;
    NbBlocs     bloc;
    Flag        multiblocs;
    }   Criteres, *P_Criteres;


/******************************************************************************/
/* FONCTIONS PUBLIQUES                                                        */
/******************************************************************************/

Flag    setBloc(P_Criteres, NbBlocs);
Flag    setLongueurBloc(P_Criteres, NbBlocs num_bloc, LongSeq lon);
Flag    setErreurBloc(P_Criteres, NbBlocs num_bloc, LongSeq erreur);
void    setErreurGlobal(P_Criteres, LongSeq erreur);
Flag    setSaut(P_Criteres, NbBlocs num_bloc, LongSeq min, LongSeq max);
NbBlocs getBloc(Criteres);
LongSeq getLongueurBloc(Criteres,NbBlocs num_bloc);
LongSeq getErreurBloc(Criteres, NbBlocs num_bloc);
LongSeq getErreurGlobale(Criteres);
LongSeq getJoker(Criteres);
Fourchette getSaut(Criteres,NbBlocs num_bloc);
int     maxLongMod(Criteres);
void    afficheCriteres(Criteres, FILE *);
Flag    verifCriteres(Criteres);
void    initCriteres(Criteres *);
Flag    chargeCriteres(Criteres *, char *);


int     addSaut2Code(int oldcode, LongSeq saut, LongSeq curbloc, P_Criteres cr);

void    initTabSauts(P_Criteres);

Flag    allocBloc(P_Criteres cr, NbBlocs bloc);

#endif
