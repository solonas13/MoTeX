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

#ifndef _OCC_H
#define _OCC_H

#include<structures.h>

/* Structure de maillon d'une liste d'occurrences             */
typedef struct liste_occ {
    Noeud *             x;      /* noeud de depart de l'etat virtuel */
    LongSeq             lon;    /* longueur sur la transition        */
    LongSeq             xerr;   /* nb total d'erreurs de l'occurrence      */
    LongSeq             saut;   /* longueur totale de saut */
    LongSeq             blocerr;/* erreurs du bloc courant */
    int                 codesaut; /* code indiquant les != sauts empruntes */
    int                 num;    /* numero de la transition           */
    } Occ, *P_occ;

void    initOcc(P_occ);

#endif
