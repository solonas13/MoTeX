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

/******************************************************************************/
/* ALPHABET - Gestion des alphabets de SMILEv1.4                              */
/******************************************************************************/

#ifndef _ALPHABET
#define _ALPHABET

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <global.h>

#define equiv(i,j)  TabSymb[i][j]
#define MAXSYMBMOD  127             /* Nb max de symboles des modeles */

/******************************************************************************/
/* TYPES ABSTRAITS                                                            */
/******************************************************************************/
typedef unsigned char Symbole;

/******************************************************************************/
/* STRUCTURES                                                                 */
/******************************************************************************/




/******************************************************************************/
/* PROTOTYPES                                                                 */
/******************************************************************************/
void        initAlphabet(void);
Symbole *   chargeAlphabet(FILE *f, Symbole **seq, NbSeq nbseq);
int         str2nummod(char *str);
void        transAlphMod(Flag);

#endif


