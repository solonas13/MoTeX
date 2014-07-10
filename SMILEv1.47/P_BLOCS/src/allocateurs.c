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

#include<allocateurs.h>

Allocateur Noeud_Alloc;
Allocateur Feuille_Alloc;
Liste *FREE_LISTE=NULL;
Liste *ALLOC_LISTE=NULL;

#ifdef DEBUG_J_ALLOC
int nb_alloc_alloc_cell=0;
int nb_liste_cell=0;
#endif
Alloc_Cell *Alloc_Alloc_Cell(int obj_size)
{
  Alloc_Cell *cel;
#ifdef DEBUG_J_ALLOC
  nb_alloc_alloc_cell++;
#endif
  cel = (Alloc_Cell *) malloc(sizeof(Alloc_Cell));
  if (cel == NULL)
    {
      fprintf(stderr,"No Enougth space\nProgram Abord\n");
      exit(-1);
    }
  cel->suivant = NULL;
  cel->data = (unsigned char *)malloc(getpagesize());
  if (cel->data == NULL)
    {
      fprintf(stderr,"No Enougth space\nProgram Abord\n");
      exit(-1);
    }
  cel->current=0;
  cel->max = getpagesize()/obj_size;
  return cel;
}

Noeud *Alloc_Noeud(void )
{
  int i;
  Noeud *tmp;

  if (Noeud_Alloc.last->max==Noeud_Alloc.last->current)
    {
      Noeud_Alloc.last->suivant = Alloc_Alloc_Cell(sizeof(Noeud));
      Noeud_Alloc.last = Noeud_Alloc.last->suivant ;
      tmp = (Noeud *)(Noeud_Alloc.last->data );
      Noeud_Alloc.last->current+=1;
    }
  else
    {
      tmp = (Noeud *)(Noeud_Alloc.last->data + sizeof(Noeud)*Noeud_Alloc.last->current);
      Noeud_Alloc.last->current+=1;
    }
  
  tmp->debut        = 0   ;
  tmp->fin          = 0   ;
  tmp->suffixe_link = NULL;
  tmp->fils = (Noeud **)malloc(sizeof(Noeud *)*ALPHA_CARD);
  tmp->sequence_number =  current_sequence;
  if (!tmp->fils)
    {
      fprintf(stderr,"No Enougth space\nProgram Abord\n");
      exit(-1);
    }
  for(i=0;i<ALPHA_CARD;i++)
    tmp->fils[i]=NULL;
  
  tmp->sequences = AllocBitTab();
  
  return tmp;
}

Feuille *Alloc_Feuille(void )
{
  Feuille *tmp;
  
  if (Feuille_Alloc.last->max==Feuille_Alloc.last->current)
    {
      Feuille_Alloc.last->suivant = Alloc_Alloc_Cell(sizeof(Feuille));
      Feuille_Alloc.last = Feuille_Alloc.last->suivant ;
      tmp = (Feuille *)(Feuille_Alloc.last->data );
      Feuille_Alloc.last->current+=1;
    }
  else
    {
      tmp = (Feuille *)(Feuille_Alloc.last->data + sizeof(Feuille)*Feuille_Alloc.last->current);
      Feuille_Alloc.last->current+=1;
    }
  
  tmp->debut      = LEAF_BIT;
  tmp->fin_deb    = -1;
  tmp->sequence_number =  current_sequence;
  
  tmp->sequences = AllocBitTab();
  addBitTabValue(&(tmp->sequences),current_sequence);
  return tmp;
}


Liste *Alloc_Liste(void)
{
  Liste *l;
  if (FREE_LISTE==NULL)
    {
      l =(Liste *)malloc(sizeof(Liste));
#ifdef DEBUG_J_ALLOC
      nb_liste_cell++;
#endif
      if (!l)
	{
	  fprintf(stderr,"No Enougth space\nProgram Abord\n");
	  exit(-1);
	}
      l->suiv = NULL;
      l->feuille = NULL;
    }
  else
    {
      l = FREE_LISTE;
      FREE_LISTE = FREE_LISTE->suiv;
      l->suiv = NULL;
      l->feuille = NULL;
    }
  return l;
}

void Free_Liste(Liste *l)
{
  l->suiv = FREE_LISTE;
  FREE_LISTE = l;
}

void Free_All_Liste_Cell(void)
{
  Liste *tmp;
#ifdef DEBUG_J_ALLOC
  int nb=0;
#endif
  while(FREE_LISTE != NULL)
    {
      tmp = FREE_LISTE;
      FREE_LISTE = FREE_LISTE->suiv;
      free(tmp);
#ifdef DEBUG_J_ALLOC
      nb++;
#endif
    }
#ifdef DEBUG_J_ALLOC
  printf("%d/%d Liste Cell désalouées...\n",nb,nb_liste_cell);
#endif
  FREE_LISTE=NULL;
}

void Init_Allocateurs(void )
{
#ifdef DEBUG_J_ALLOC
  printf("Init Allocateur : %d alloc cell, %d liste cell\n",nb_alloc_alloc_cell,nb_liste_cell);
#endif
  Noeud_Alloc.first=Alloc_Alloc_Cell(sizeof(Noeud));
  Noeud_Alloc.last=Noeud_Alloc.first;
  Feuille_Alloc.first=Alloc_Alloc_Cell(sizeof(Feuille));
  Feuille_Alloc.last=Feuille_Alloc.first;
  FREE_LISTE=NULL;
}

void Free_Arbre(Noeud *Racine)
{
  Alloc_Cell *tmp,*tmp2;
  Noeud      *n;
  Feuille    *f;
  int        i;

#ifdef DEBUG_J_ALLOC
  int nb=0;
#endif
  tmp=Noeud_Alloc.first;
  while(tmp!=NULL)
    {
#ifdef DEBUG_J_ALLOC
      nb++;
#endif
      for(i=0; i<tmp->current;i++)
          {
          n = (Noeud *) (tmp->data + i*sizeof(Noeud));
          free(n->fils);
          free(n->sequences);
          }

      tmp2=tmp->suivant;
      free(tmp->data);
      free(tmp);
      tmp=tmp2;
    }
  tmp=Feuille_Alloc.first;
  while(tmp!=NULL)
    {
#ifdef DEBUG_J_ALLOC
      nb++;
#endif
      for(i=0; i<tmp->current;i++)
          {
          f = (Feuille *) (tmp->data + i*sizeof(Feuille));
          free(f->sequences);
          }

      tmp2=tmp->suivant;
      free(tmp->data);
      free(tmp);
      tmp=tmp2;
    }
  Noeud_Alloc.first=Noeud_Alloc.last=Feuille_Alloc.first=Feuille_Alloc.last=NULL;
#ifdef DEBUG_J_ALLOC
  printf("%d/%d alloc cell  OK \n",nb,nb_alloc_alloc_cell);
#endif
}



