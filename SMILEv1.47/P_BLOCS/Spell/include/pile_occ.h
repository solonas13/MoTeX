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

#ifndef _PILE_OCC_H
#define _PILE_OCC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <global.h>
#include <occ.h>
#include <liste_pos.h>

/* Grain du tableau des occurrences                                           */
#define GRAIN 2000

typedef struct struct_pile_occ
	{
	Occ				*occ;
/* Carte des positions des 'dummy'                                            */
	unsigned int	*carte;
/* Position courante dans carte                                               */
	unsigned int	pos_carte;
	unsigned int	size_carte;
	unsigned int	size;
	int				pos;
	} PileOcc, *P_PileOcc;


/******************************************************************************/
/* FONCTIONS PUBLIQUES                                                        */
/******************************************************************************/
P_PileOcc	creePileOcc(void);
void		ajouteDummy(P_PileOcc);
LongSeq		getPrecDummy(P_PileOcc);
void		ajouteInitOcc2Pile(P_PileOcc, Noeud *); 
void		ajouteOcc2Pile(P_PileOcc, Noeud *, int,LongSeq,LongSeq ,LongSeq,
				LongSeq, int);
void		transferePile2Pile(P_PileOcc, P_PileOcc);
int			copieLastOcc(P_PileOcc, P_PileOcc);
void		depileRec(P_PileOcc);
void		videPile(P_PileOcc);
void		liberePileOcc(P_PileOcc);

#if OCC
void		afficheLastOcc(FILE *f, P_PileOcc, LongSeq l, P_Criteres cr);
#endif


#if DEBUG_BASE
int			afficheOcc(FILE *f, P_occ o, LongSeq longmod, P_Criteres cr);
void		afficheOldOcc(P_PileOcc p, LongSeq l);
void		affichePileOcc(P_PileOcc);
#endif

#endif
