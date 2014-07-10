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

#ifndef _DEFINE_H
#define _DEFINE_H

int ALPHA_CARD;

#define LEAF_BIT             0x80000000
#define LEAF_BIT_INV         0x7FFFFFFF

#define POS_ALLOC_STEP       50
#define LISTE_CHANGE_BIT     0x80000000
#define LISTE_CHANGE_BIT_INV 0x7FFFFFFF
#define LISTE_END            0x0FFFFFFF

#define FEUILLE_ALLOC_STEP   100
#define NOEUD_ALLOC_STEP     100

#define DEBUG_JTREE  0


#define DEFAULT_WINDOW_SIZE_OPTION 4

#endif
