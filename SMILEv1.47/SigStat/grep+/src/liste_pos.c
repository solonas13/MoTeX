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

ListePositions *Alloc_ListePositions(int size)
{
  ListePositions * tmp = (ListePositions *)malloc(sizeof(ListePositions));
  if (!tmp)
    {
      fprintf(stderr,"No Enough space\nProgram Abort\n");
      exit(-1);
    }
  tmp->tab[0] = (int *)malloc(sizeof(int) * size);
  tmp->tab[1] = (int *)malloc(sizeof(int) * size);

#if DEBUG_JTREE
  printf("Alloc Liste.... %d \n",sizeof(int) * size);
#endif

  if ((!tmp->tab[0]) || (!tmp->tab[1]))
    {
      fprintf(stderr,"No Enough space\nProgram Abort\n");
      exit(-1);
    }
  memset(tmp->tab[0],0,sizeof(int) * size);
  memset(tmp->tab[1],0,sizeof(int) * size);
  tmp->last_cell = 0;
  tmp->tab_size = size;
  return tmp;
}

int ChercheDerniereCelluleDansListe(ListePositions *lpos,int deb_liste)
{
  if (lpos == NULL)
    return -1;
  while(lpos->tab[1][deb_liste] != LISTE_END)
    deb_liste = lpos->tab[1][deb_liste] & LISTE_CHANGE_BIT_INV;
  return deb_liste;
}

int Ajoute_Position_Liste(ListePositions *lpos,int *deb_liste,int position,int change_seq)
{
  if (!lpos)
    return -2;
  
  if (lpos->last_cell == lpos->tab_size)
    {

#if DEBUG_JTREE
      printf("realloc .... LPOS\n");
#endif

      lpos->tab[0] = (int *)realloc(lpos->tab[0],sizeof(int)*(lpos->tab_size + POS_ALLOC_STEP));
      lpos->tab[1] = (int *)realloc(lpos->tab[1],sizeof(int)*(lpos->tab_size + POS_ALLOC_STEP));
      lpos->tab_size+=POS_ALLOC_STEP;
    }
  lpos->tab[0][lpos->last_cell] = position ;
  lpos->tab[1][lpos->last_cell] = LISTE_END ;
  
  if (*deb_liste!=-1)
    { 
      if (change_seq)
    lpos->tab[1][lpos->last_cell] = *deb_liste | LISTE_CHANGE_BIT;
      else
    lpos->tab[1][lpos->last_cell] = *deb_liste;
    }
  *deb_liste = lpos->last_cell;
  (lpos->last_cell)++;
  return 1;
}

int getValue(ListePositions *lpos,int i)
{
  if (lpos==NULL)
    return -2;
  if ((i<0) || (i>lpos->last_cell))
    return -3;
  return lpos->tab[0][i] ;/*& LISTE_CHANGE_BIT_INV; */
}

void setListeValue(ListePositions *lpos,int i,int value)
{
  if ((lpos==NULL) || (i<0) || (i>lpos->last_cell))
    return ;
  lpos->tab[0][i] = value;
}

int getIndiceSuivant(ListePositions *lpos,int i)
{
  if (lpos==NULL)
    return -2;
  if ((i<0) || (i>lpos->last_cell))
    return -3;
  return lpos->tab[1][i];
}


void Free_ListePositions(ListePositions *lpos)
{
  if (lpos == NULL)
    return ;
  free(lpos->tab[0]);
  free(lpos->tab[1]);
  free(lpos); 
}

int Print_Positions_Dynamique(FILE *f, Feuille *n, LongSeq longway,
        P_Criteres cr, int code)
{
int indice      = n->fin_deb,
    occurrence  = 0;

#if AFF_OCC
int longueur, *i,
    nb_element  = ((n->sequences[0] & 0x7F) << 8) | n->sequences[1];

nb_element--; 
longueur    = getValue(Liste_positions_fin,indice)-(n->debut&LEAF_BIT_INV)+longway;
#endif

while(indice != LISTE_END)
    {
#if AFF_OCC
    if(f!=NULL)
        {
        fprintf(f,"Seq %5d   Pos %5d",
            ((unsigned short int *)(n->sequences + 2))[nb_element],
            getValue(Liste_positions_fin,indice)-longueur);

        if(cr && cr->bloc != 1)
            {
            fprintf(f, "\tSaut ");
            for(i=cr->code2Sauts[code]; i!=cr->code2Sauts[code]+cr->bloc-1; i++)
                fprintf(f, "%d ",*i);
            }
        
        fprintf(f, "\n");
        }
#endif

    indice = getIndiceSuivant(Liste_positions_fin,indice);
#if AFF_OCC
    if (indice & LISTE_CHANGE_BIT)
        nb_element--;
#endif
    indice = indice & LISTE_CHANGE_BIT_INV;
    occurrence++;
    }
return occurrence;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int Print_Positions_Statique(FILE *f, Feuille *n, LongSeq longway,
        P_Criteres cr, int code)
{
unsigned char   mask = 0x01;
int             compteur = SIZE_STATIC_BIT_TAB-1,
                sequence = 8 * compteur + 6,
                indice = n->fin_deb,
                occurrence  = 0;

#if AFF_OCC
int             longueur, *i;
            
longueur    = getValue(Liste_positions_fin,indice)-(n->debut&LEAF_BIT_INV)+longway;
#endif

while ( (n->sequences[compteur] & mask) == 0 )
    {
    mask <<= 1;
    sequence--;

    if ( mask == 0 )
        {
        mask    = 0x01;
        compteur--;
        }
    }

while(indice != LISTE_END)
    {
#if AFF_OCC
    if(f!=NULL)
        {
        fprintf(f,"Seq %5d   Pos %5d",sequence, getValue(Liste_positions_fin,indice)-longueur);

        if(cr && cr->bloc != 1)
            {
            fprintf(f, "\tSaut ");
            for(i=cr->code2Sauts[code]; i!=cr->code2Sauts[code]+cr->bloc-1; i++)
                fprintf(f, "%d ",*i);
            }
        
        fprintf(f, "\n");
        }
#endif
    indice = getIndiceSuivant(Liste_positions_fin,indice);
    if ( (indice & LISTE_CHANGE_BIT) && ( indice != LISTE_END ) )
        {
        do
            {
            mask <<= 1;
            sequence--;
            
            if ( mask == 0 )
                {
                mask    = 0x01;
                compteur--;
                }
            }
        while ( (n->sequences[compteur] & mask) == 0 );
        }
    occurrence++;
    indice = indice & LISTE_CHANGE_BIT_INV;
    }
return occurrence;
}

int Print_Positions(FILE *f, Feuille *n, LongSeq longway, P_Criteres cr, int code)
{
/* #if DEBUG */
/* printf("J'arrive dans Print avec long %d\n",longway); */
/* #endif */

if (n->sequences[0] & 0x80)
    return (Print_Positions_Dynamique(f, n, longway, cr, code));

return (Print_Positions_Statique(f, n, longway, cr, code));
}


