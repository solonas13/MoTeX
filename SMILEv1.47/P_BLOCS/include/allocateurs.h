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

#ifndef ALLOCATEUR_H
#define ALLOCATEUR_H

#include<structures.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <bit_tab.h>


/*
                                              SINGLE_TAB_SIZE
			                   __________|___________
                                         /                       \
                         /  ----            
                         |  |  |   ----> [ | | | | | | | | | ....]
                         |  ----
                         |  |  |   ----> [ | | | | | | | | | ....]   
ALLOC_NOEUD_TAB_STEP -   ----
                         |  |  |   ----> [ | | | | | | | | | ....]
			 |  ----
			 |  |  |   ----> [ | | | | | | | | | ....]
			 \  ----
			    ....
 */

Noeud *Alloc_Noeud(void );
Feuille *Alloc_Feuille(void );
Liste *Alloc_Liste(void);
void Free_Liste(Liste *l);
void Free_All_Liste_Cell(void);
void Init_Allocateurs(void );


void Free_Arbre(Noeud *Racine);

typedef struct _cell
{
  struct _cell *suivant;
  unsigned char *data;
  int current;
  int max;
}Alloc_Cell;

typedef struct _allocateur
{
  Alloc_Cell *first;
  Alloc_Cell *last;
}Allocateur;

#endif
