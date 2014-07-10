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

#include<liste_pos.h>


int *Positions = NULL;
int current_pos=0;

ListePositions *Alloc_ListePositions(int size)
{
  Positions =(int *)malloc(sizeof(int)*1000000);
  current_pos=0;
  return NULL;
}

int Ajoute_Position_Liste(ListePositions *lpos,int *deb_liste,int position,int change_seq)
{
  if (*deb_liste == -1)
    {
      *deb_liste = current_pos;
      Positions[current_pos] = position;
      current_pos++;
    }
  else
    {
      Positions[*deb_liste] = position;
    }
 
  return 1;
}

int getValue(ListePositions *lpos,int i)
{
  return Positions[i];
}

void setListeValue(ListePositions *lpos,int i,int value)
{
  Positions[i]=value;
}

int getIndiceSuivant(ListePositions *lpos,int i)
{
  return LISTE_END;
}

void Free_ListePositions(ListePositions *lpos)
{
}

