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

#ifndef _CONSTRUCTION_H
#define _CONSTRUCTION_H

#include<structures.h>
#include<global_fonctions.h>


Noeud *Construction_Arbre(unsigned char *S,int taille_fenetre);

void Premiere_Phase(Noeud *Racine,
		    int taille_fenetre,
		    Liste **debut_liste,
		    Liste **fin_liste,
		    Noeud **fin_liste_pere,
		    int *nb_element_liste,
		    int start_indice);

int  Deuxieme_Phase(Noeud *Racine,
		    int taille_fenetre,
		    Liste **debut_liste,
		    Liste **fin_liste,
		    Noeud **fin_liste_pere,
		    int *nb_element_liste,
		    int start_indice,
		    int fict,
		    int ini_res_type);

void Troisieme_Phase(Noeud *Racine,
		     int taille_fenetre,
		     Liste **debut_liste,
		     Liste **fin_liste,
		     Noeud **fin_liste_pere,
		     int *nb_element_liste,
		     int fictive);

Noeud *AjouteSequence(Noeud *Arbre,unsigned char *S,int taille_fenetre);

void CloseTheFirstPhase( Liste **debut_liste,
			 Liste **fin_liste,
			 Noeud **fin_liste_pere);


Noeud *CaseOneAddSequence(Noeud *Arbre,int taille_fenetre);

Noeud *CaseTwoAddSequence(Noeud *Arbre,Noeud *resultat, Noeud *pere,int taille_fenetre);

Noeud *CaseTreeAddSequence(Noeud *Arbre,Noeud *resultat,int deb,int fin,int taille_fenetre);

Noeud *CaseFourAddSequence(Noeud *Arbre,Noeud *resultat,Noeud *pere,int res_type,int position_arc,int i,int taille_fenetre);




#if DEBUG_JTREE
extern int nb_alloc_noeud;
extern int nb_alloc_feuille;
extern int nb_alloc_liste;
#endif

#endif
