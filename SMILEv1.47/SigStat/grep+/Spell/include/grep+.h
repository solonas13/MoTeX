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

#ifndef _SPELL_H
#define _SPELL_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <time.h>
#include <math.h>

#include <global.h>
#include <criteres.h>
#include <occ.h>
#include <libfasta.h>
#include <pile_occ.h>
#include <model.h>
#include <struct_tab.h>
#include <bit_tab.h>
#include <liste_pos.h>
#include <construction.h>
#include <allocateurs.h>
/*#include <barre.h>*/
#include <alphabet.h>

/* Grain d'allocation des sequences Fasta */
#define GRAINSEQ 500

typedef struct struct_arbre
    {
    Noeud   *arbre;
    signed char **text;
    int     nbtxt;
    } Arbre;


/******************************************************************************/
/* 							PROTOTYPES PUBLICS                                */
/******************************************************************************/
/* explore les modeles */
void  chercheMot ( Arbre a, P_Criteres  cr, signed char *mot,
        unsigned int *nbocc, unsigned int *nboccex, FILE *f);

/* creation de l'arbre                                                        */
Flag    creeArbreSuffixeFromFile(Arbre *, char *fic, int maxlongmod,
            char *alphaseq);
Flag    creeArbreSuffixeFromArray(Arbre *, char **seq, int nbseq,
            int maxlongmod, char *alphaseq);
void    libereArbreSuffixeFromFile(Arbre);
void    libereArbreSuffixeFromArray(Arbre);

/* Initialise la recherche                                                    */
Flag    valideCriteres(P_Criteres cr);

#endif
