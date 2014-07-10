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

#ifndef _BIT_TAB_H
#define _BIT_TAB_H

#include<string.h>
#include<struct_tab.h>
#include<define.h>

void initBitTab(int nb_seq);
Bit_Tab *AllocBitTab(void);

void addBitTabValue(Bit_Tab **tab,int value);
void fusionneBitTab(Bit_Tab **tab1,Bit_Tab *tab2); /* 1 <- 1 & 2 */
int  nbSequenceInBitTab(Bit_Tab *tab);

void printBitTab(Bit_Tab *tab);
void CopyBitTab(Bit_Tab **dest,Bit_Tab *src);
void ReinitBitTab(Bit_Tab **bt);
/*---------------------------------------------------------*/

void convertBitTab(Bit_Tab **tab); /* transforme tab de dyn a static... */

void addBitTabValueStatic(Bit_Tab **tab,int value);
void addBitTabValueDynamic(Bit_Tab **tab,int value);

int  nbSequenceInBitTabStatic(Bit_Tab *tab);
int  nbSequenceInBitTabDynamic(Bit_Tab *tab);

void fusionneBitTabStatic(Bit_Tab **tab1,Bit_Tab *tab2);
void fusionneBitTabDynamic(Bit_Tab **tab1,Bit_Tab *tab2);
Bit_Tab *AllocBitTabStatic(void);

#if DEBUG_JTREE
extern int nb_alloc_tab;
#endif
#endif 

