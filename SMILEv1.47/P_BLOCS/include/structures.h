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

#ifndef _STRUCTURES_H
#define _STRUCTURES_H

#include<define.h>
#include<struct_tab.h>
#include<global.h>

typedef struct _Liste_Pos
{
  int *tab[2];
  int last_cell;
  int tab_size;
} ListePositions;

typedef struct _Noeud /* Ne pas changer l'ordre des 3 premiers champs! */
{
  LongSeq debut;
  NbSeq sequence_number;
  LongSeq fin;
  struct _Noeud *suffixe_link;
  struct _Noeud  **fils;
  Bit_Tab *sequences;
  int nb_element_bt;
}Noeud;

typedef struct _Feuille /* Ne pas changer l'ordre des 3 premiers champs! */
{
  LongSeq debut;
  NbSeq sequence_number;
  LongSeq fin_deb;
  Bit_Tab *sequences;
}Feuille;

typedef struct _Liste
{
  struct _Liste *suiv;
  Feuille *feuille;
} Liste;



extern unsigned char Translation_Table[255];

extern unsigned char **Sequence;
extern ListePositions *Liste_positions_fin;

extern int global_indice;
extern int current_sequence;

#endif
