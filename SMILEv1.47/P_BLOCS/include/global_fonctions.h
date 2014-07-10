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

#ifndef _GLOBAL_FONCTIONS_H
#define _GLOBAL_FONCTIONS_H

#include<structures.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<liste_pos.h>
#include<allocateurs.h>
#include<bit_tab.h>

void Init_All(unsigned char *Alphabet,int Joker,int nb_sequence);

void Ajoute_Fils_Au_Noeud(Noeud *N,Noeud *F);
Noeud *Get_Child_Start_Letter(Noeud *N,int indice);
int seg_taille(Noeud *N);

Noeud *Add_Fast_String(Noeud *N,int deb,int fin,int *type,Noeud **pere);
            /*
	      --> Retourne:
	        * Si ajout d'une feuille a l'arbre : 
		                      @ de la feuille.
				      type = 1.
	        * Si extension d'une feuille a l'arbre:
		                      @ de la feuille.
				      type = 2.
		* Si decoupe d'un arc avec creation d'une feuille:
		                      @ de la feuille cree.
				      type = 3.
	        * Si rien (chaine deja ds l'arbre) :
		                      @ du dernier noeud en amont.
				      type = - lg du dernier seg.
	     */
int compare_string(int d1,int f1,
		   int d2,int f2);
void Print_Tree(Noeud *N,int affichage,int stat);
void Print_Liste(Liste *liste);
Noeud *FindString(Noeud *N,int deb,int fin,Noeud **pere,int *restant,int *pos_in_edge);
				/*
				  cherche la chaine deb fin de la sequence
				  courante à partir de N.
				  retourne le sommet au bout de l'arc
				  contenant la chaine cherchée.
				  pere = pere du sommet retouné.
				  restant :
				      < 0 : coupure au milieu de l'arc
				            pos_in_edge : nb de car. commun sur l'arc
				      1   : la chaine est dans l'arbre et elle
				            aboutit à un sommet.
				      2   : la recherche aboutit à une feuille.
				            et elle est plus grande que l'arc.
				      3   : la chaine est plus courte que l'arc.
				            mais elle est contenue dans celui-ci.
				 */

void UpdateBit_TabForAllTree(Noeud *N);
void Print_BTTree_Debug(Noeud *N,int *cpt);

#endif
