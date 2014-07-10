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
#include <global.h>
#include <alphabet.h>

/******************************************************************************/
/* STRUCTURE DE STOCKAGE DES CRITERES DE RECHERCHE                            */
/******************************************************************************/
typedef struct  struct_fourchette
    {
    LongSeq     min;
    LongSeq     max;
    }   Fourchette, *P_Fourchette;

typedef struct  struct_criteres
    {
    LongSeq     **compobloc;
    LongSeq     *maxerrblocs;
    LongSeq     *compo;
    Fourchette  *longbloc;
    Fourchette  *saut;
    Flag        *flag_compobloc;
    LongSeq     *palindrom;
    long int    nbsymb;
    Fourchette  longueur;
    LongSeq     maxerr;
    NbSeq       quorum;
    char        bloc;
    Flag        flag_compo;
    Flag        multiblocs;
    Flag        flag_palindrom;
    }   Criteres, *P_Criteres;


/******************************************************************************/
/* FONCTIONS PUBLIQUES                                                        */
/******************************************************************************/


void    setCompoPal(P_Criteres cr, char **argv, int argc);

int     addSaut2Code(int oldcode, LongSeq saut, LongSeq curbloc, P_Criteres cr);

void    initTabSauts(P_Criteres);

void    allocBloc(P_Criteres cr, int bloc);

#endif
