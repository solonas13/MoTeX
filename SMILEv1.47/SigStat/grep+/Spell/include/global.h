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

#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <strings.h>
#include <malloc.h>
#include <symb.h>


#define FAUX    0
#define VRAI    1



typedef char Flag;

/******************************************************************************/
/* Flags                                                                      */
/******************************************************************************/
/* DEBUGGING                                                                  */
#define DEBUG_BASE  0       /* Debug base */
#define DEBUG_BT    0       /* Tableaux de bits */
#define DEBUG_SAUT  0       /* Procedures de saut */
#define DEBUG_PILE  0       /* Pile d'occurrences */
#define DEBUG_TREE  0       /* Arbre suffixe : HS bicoz Julien */


/******************************************************************************/
/* Define dependants du jeu de donnees                                        */
/******************************************************************************/
/* Grain d'allocation de la taille du modele                                  */
#define GRAIN_SIZMOD    1000

/******************************************************************************/
/* Caracteres speciaux                                                        */
/******************************************************************************/
/* => dans symb.h */

/******************************************************************************/
/* Types                                                                      */
/******************************************************************************/
/* Nombre de sequences                                                        */
#define NbSeq   int
/* Longueur de sequence                                                       */
#define LongSeq int
/* Nombre de blocs */
#define NbBlocs signed char


/******************************************************************************/
/* Active DEBUG_BASE si l'un des DEBUGs est active                            */
/******************************************************************************/
#if DEBUG_BT || DEBUG_SAUT || DEBUG_PILE || DEBUG_TREE
    #undef DEBUG_BASE
    #define DEBUG_BASE 1
#endif

/******************************************************************************/
/* Fonctions basiques                                                         */
/******************************************************************************/
void    fatalError(char *msg);

int     entiers(char);

void    entree(void);

void    initEntiers(void);

#endif
