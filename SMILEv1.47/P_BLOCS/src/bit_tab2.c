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

#include<struct_tab.h>
#include<define.h>

void initBitTab(int nb_seq)
{}
Bit_Tab *AllocBitTab()
{
  return NULL;
}

void addBitTabValue(Bit_Tab **tab,int value)
{}

void fusionneBitTab(Bit_Tab **tab1,Bit_Tab *tab2)
{} // 1 <- 1 & 2

int  nbSequenceInBitTab(Bit_Tab *tab)
{
  return 0;
}

void printBitTab(Bit_Tab *tab)
{}
