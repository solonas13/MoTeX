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

#ifndef _LIST_POS_H
#define _LIST_POS_H

#include<define.h>
#include<structures.h>
#include<stdio.h>
#include<stdlib.h>
#include<criteres.h>

ListePositions *Alloc_ListePositions(int size);

int Ajoute_Position_Liste(ListePositions *lpos,int *deb_liste,int position,int change_seq);
int getValue(ListePositions *lpos,int i);
void setListeValue(ListePositions *lpos,int i,int value);
int getIndiceSuivant(ListePositions *lpos,int i);

void Free_ListePositions(ListePositions *lpos);
int	Print_Positions(FILE *f, Feuille *n, LongSeq longway, P_Criteres cr, int code);

#endif
