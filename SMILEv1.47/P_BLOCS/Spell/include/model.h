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

#ifndef _MODEL_INCLUDE
#define _MODEL_INCLUDE

#include <stdlib.h>
#include <global.h>

/* Structure de stockage d'un modele                          */

typedef struct model {
    int     *name;      /* Sequence du modele (codes alphabets) */
    LongSeq lon;        /* Taille du modele stocke dans name */
    LongSeq taille;     /* Taille de name */
    } Mod, *P_mod;


/* FONCTIONS */
/* allouer un modele */
P_mod   allocModel(void);

/* diminuer la taille d'un modele de 1 */
void    decrModel(P_mod);

/* ajoute un symbole au modele */
void    changeModel(P_mod , int);

/* Affiche le modele                                                           */
void AfficheModel(P_mod);


#endif

