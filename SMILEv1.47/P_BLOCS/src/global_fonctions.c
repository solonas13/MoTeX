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

#include<global_fonctions.h>

void Init_All(unsigned char *Alphabet,int Joker,int nb_sequence)
{
  int i;
  int l;
  
  initBitTab(nb_sequence);
  
  Init_Allocateurs();
  
  l=strlen((const char *)Alphabet);
  
  if (Joker)
    {
      for(i=0;i<255;i++)
	Translation_Table[i]=0;
      for(i=0;i<l;i++)
	Translation_Table[Alphabet[i]]=i+1;
      ALPHA_CARD = l + 1;
    }
  else
    {
      for(i=0;i<255;i++)
	Translation_Table[i]=0xFF;
      for(i=0;i<l;i++)
	Translation_Table[Alphabet[i]]=i;
      ALPHA_CARD = l;
    }

  Sequence =(unsigned char **) malloc (nb_sequence * sizeof(unsigned char *));
  current_sequence = 0;
}

void Ajoute_Fils_Au_Noeud(Noeud *N,Noeud *F)
{
  if ((!N) || (!F)
      || (N->debut & LEAF_BIT))
    return ;
  N->fils[Translation_Table[Sequence[F->sequence_number][F->debut & LEAF_BIT_INV]]] = F;
}


Noeud *Get_Child_Start_Letter(Noeud *N,int indice)
{
  if ((Translation_Table[Sequence[current_sequence][indice]] == 255)
      || (!N)
      || (N->debut & LEAF_BIT))
    return NULL;
  return N->fils[Translation_Table[Sequence[current_sequence][indice]]];
}

int seg_taille(Noeud *N)
{
  
  if ( N->debut & LEAF_BIT)
    {
      if ( getValue(Liste_positions_fin,((Feuille *)N)->fin_deb) == -1)
	return (global_indice-
		((((Feuille *)N)->debut) & LEAF_BIT_INV));
      else
	return ( getValue(Liste_positions_fin,((Feuille *)N)->fin_deb) - 
		((((Feuille *)N)->debut) & LEAF_BIT_INV));
    }
  else
    {
      if (N->fin == -1)
	return global_indice - N->debut;
      else
	return (N->fin - N->debut);
    }
}

Noeud *Add_Fast_String(Noeud *N,int deb,int fin,int *type,Noeud **pere)
{
  Noeud * res, *tmp_n;
  Feuille *tmp_f;
  int tmp_i,res_d;
  
  if (!N)
    return NULL;
  
  if (N->debut & LEAF_BIT)
    {
      *type = 2; /* Extension d'un feuille. */
      if (seg_taille(N)<(fin -deb))
	{
	  N->debut = (deb - (getValue(Liste_positions_fin,N->fin) - (N->debut & LEAF_BIT_INV)))| LEAF_BIT;
	  setListeValue(Liste_positions_fin,((Feuille *)N)->fin_deb,fin);
	}
      else
	if ((seg_taille(N)==(fin -deb)) &&((getValue(Liste_positions_fin,((Feuille *)N)->fin_deb))!=-1))
	  {
	    if ((N->sequence_number!=current_sequence) 
		||
		(getValue(Liste_positions_fin,((Feuille *)N)->fin_deb) != fin))
	      {
		N->debut = deb | LEAF_BIT;
		if (N->sequence_number==current_sequence)
		  {
		    Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)N)->fin_deb),fin,0);
		  }
		else
		  {
		    Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)N)->fin_deb),fin,1);
		    addBitTabValue(&(((Feuille *)N)->sequences),current_sequence);
		    N->sequence_number = current_sequence;
		  }
	      }
	  }
      return N;
    }
  
  res = Get_Child_Start_Letter(N,deb);
  
  if (res == NULL)
    {
      *type = 1; /* Creation d'une feuille. */
      tmp_f = Alloc_Feuille();
      tmp_f->debut = deb | LEAF_BIT;
      Ajoute_Position_Liste(Liste_positions_fin,&(tmp_f->fin_deb),fin,0);
      *pere = N;
      Ajoute_Fils_Au_Noeud(N,(Noeud *)tmp_f);
      return (Noeud *)tmp_f;
    }
  
  tmp_i = seg_taille(res);
  res_d = res->debut & LEAF_BIT_INV;
  
  if (deb + tmp_i >= fin)
    {
      if (Translation_Table[Sequence[current_sequence][fin-1]] 
	  == 
	  Translation_Table[Sequence[res->sequence_number][res_d + (fin - deb) - 1]])
	{
	  *type = deb-fin;
	  if (res->debut & LEAF_BIT)
	    {
	      if ((getValue(Liste_positions_fin,((Feuille *)res)->fin_deb)!=-1)
		  && ((fin - deb) == seg_taille(res)))
		{
		  if ((res->sequence_number != current_sequence)
		      || (getValue(Liste_positions_fin,((Feuille *)res)->fin_deb) != fin))
		    {
		      ((Feuille *)res)->debut = deb | LEAF_BIT;
		      if (res->sequence_number != current_sequence)
			{
			  Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)res)->fin_deb),fin,1);
			  addBitTabValue(&(((Feuille *)res)->sequences),current_sequence);
			  res->sequence_number = current_sequence;
			}
		      else
			Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)res)->fin_deb),fin,0);
		    }
		}
	    }
	  return N;
	}
      else
	{
	  tmp_n = Alloc_Noeud();
	  tmp_n->debut = res_d;
	  tmp_n->fin = res_d + fin - deb - 1;
	  tmp_n->sequence_number = res->sequence_number;
	  tmp_n->suffixe_link = N;
	  tmp_n->fils[Translation_Table[Sequence[res->sequence_number][tmp_n->fin]]]=res;
	  
	  if (res->debut & LEAF_BIT)
	    {
	      res->debut =  tmp_n->fin  | LEAF_BIT;
	    }
	  else
	    res->debut = tmp_n->fin;
	  
	  tmp_f = Alloc_Feuille();
	  *pere = tmp_n;
	  tmp_f->debut = (fin - 1) | LEAF_BIT;
	  Ajoute_Position_Liste(Liste_positions_fin,&(tmp_f->fin_deb),fin,0);
	  N->fils[Translation_Table[Sequence[tmp_n->sequence_number][tmp_n->debut]]]=tmp_n;
	  tmp_n->fils[Translation_Table[Sequence[tmp_f->sequence_number][tmp_f->debut&LEAF_BIT_INV]]]=(Noeud *)tmp_f;
	  *type = 3;
	  return (Noeud *)tmp_f;
	}
    }
  *pere = N;
  return Add_Fast_String(res,deb+tmp_i,fin,type,pere);
}

int compare_string(int d1,int f1,
		   int d2,int f2)
{
  int i,j;
  for(i=d1,j=d2;(i<f1) && (j<f2);i++,j++)
    if (Sequence[i]!=Sequence[j])
      break;
  if ((i==f1)
      &&
      (j==f2))
    return -1; /* A == B */
  if (i==f1)
    return -2; /* A incl B */
  if (j==f2)
    return -3; /* B incl A */
  
  return i-d1; /* n avec A[d1+k]=B[d2+k] pour k<n */
}

void Print_Tree_Indent(Noeud *N,int *nb_noeud,int *nb_feuilles,int *nb_fils,char *indent,int affichage)
{
  int i,tmp;
  int d,f,j;
  if (!N)
    return;
  
  if (N->debut & LEAF_BIT)
    {
      *nb_feuilles = *nb_feuilles + 1;
      if (affichage){
	printf("Feuille (%p): %d->(",N,N->debut & LEAF_BIT_INV);
      
	tmp = ((Feuille *)N)->fin_deb;
	while(tmp != LISTE_END)
	  {
	    printf("%d,",(getValue(Liste_positions_fin,tmp)!=-1)?getValue(Liste_positions_fin,tmp):-global_indice);
	    tmp = getIndiceSuivant(Liste_positions_fin,tmp);
	    if (tmp & LISTE_CHANGE_BIT)
	      printf("-,");
	    tmp = tmp & LISTE_CHANGE_BIT_INV;
	  }
	printf(") ["); 
	printBitTab(((Feuille *)N)->sequences);
	tmp = getValue(Liste_positions_fin,((Feuille *)N)->fin_deb);
	d=N->debut & LEAF_BIT_INV;
	f=(tmp==-1)?global_indice-1:tmp-1;
	printf("] {%d}= (",nbSequenceInBitTab(((Feuille *)N)->sequences));
	for(j=d;j<=f;j++)
	  putchar(Sequence[N->sequence_number][j]);
	printf(")\n");
      }
    }
  else
    {
      *nb_noeud = *nb_noeud + 1;
      if (affichage)
	{
	  printf("Noeud (%p): %d->%d (",N,N->debut,N->fin);
	  d=N->debut;
	  f=(N->fin-1)>=0?N->fin-1:0;
	  for(j=d;j<=f;j++)
	    putchar(Sequence[N->sequence_number][j]);
	  printf(") [");
	  printBitTab(N->sequences);
	  printf("]{%d}\n%s       sl : %p\n%s       fils :\n",nbSequenceInBitTab(N->sequences),indent,N->suffixe_link,indent);
	  strcat(indent,"\t");
	}
      for(i=0;i<ALPHA_CARD;i++)
	{
	  if (N->fils[i])
	    {
	      *nb_fils=*nb_fils + 1;
	      if (affichage)
		printf("%s--> %c :",indent,Sequence[N->fils[i]->sequence_number][N->fils[i]->debut & LEAF_BIT_INV]);
	      Print_Tree_Indent(N->fils[i],nb_noeud,nb_feuilles,nb_fils,indent,affichage);
	    }
	}
      indent[strlen(indent)-1]='\0';
      if (affichage)
	printf("%s       fin fils\n",indent);
    }
}

void Print_Tree(Noeud *N,int affichage,int stat)
{
  char *indent = (char *)malloc(1000);
  int nb_noeud=0,nb_feuilles=0;
  int nb_fils=0;
  indent[0]='\0';
  Print_Tree_Indent(N,&nb_noeud,&nb_feuilles,&nb_fils,indent, affichage);
  if (stat)
    {
      printf("nombre de Noeud   (%d + %d) : %d \n"
	     "nombre de Feuille (%d)      : %d \n"
	     "nombre de fils              : %d \n"
	     "nombre de fils/Noeud        : %f \n"
	     "cardinal de l'alphabet      : %d \n",
	     (int)sizeof(Noeud),(int)sizeof(Noeud *)*ALPHA_CARD,nb_noeud,
	     (int)sizeof(Feuille),nb_feuilles,
	     nb_fils,(double )((double )nb_fils/(double )nb_noeud),
	     ALPHA_CARD);
    }
  free(indent);
}

void Print_Liste(Liste *liste)
{
  Liste *tmp = liste;
  while(tmp)
    {
      printf("( %p )->",tmp->feuille);
      tmp = tmp->suiv;
    }
  printf("\n");
}

 /*  Recursif!!!! */
Noeud *FindString(Noeud *N,int deb,int fin,Noeud **pere,int *restant,int *pos_in_edge)
{
  Noeud *res;
  int start,end;
  int i,j;

  
  if (N->debut & LEAF_BIT)   /* Si N est un feuille. */
    {
      start = N->debut & LEAF_BIT_INV;
      end = getValue(Liste_positions_fin,((Feuille *) N)->fin_deb);
      for(i = start,j=deb; ( (i<end) && (j<fin)) ; i++,j++)
	if (Sequence[N->sequence_number][i] != Sequence[current_sequence][j] )
	  {
	    *restant = j - fin ;  /* <0 */
	    *pos_in_edge = i - start ;
	    return N;
	  }
      if ((i==end) && (j==fin))
	{
	  /* Pile poil... C'est la */
	  *restant = 1;
	  *pos_in_edge = 0;
	  return N;
	}
      if (i==end)
	{
	  *restant = 2;
	  *pos_in_edge = 0;
	  return N;
	}
      *restant = 3;
      *pos_in_edge = 0;
      return (N);
    }
  res = Get_Child_Start_Letter(N,deb);
  if (res == NULL)
    {
      *restant = deb - fin ;  /* <0 */
      *pos_in_edge = -1;
      return N;
    }
  
  if (res->debut & LEAF_BIT)
    {
      *pere = N;
      return FindString(res,deb,fin,pere,restant,pos_in_edge);
    }
  
  start = res->debut;
  end =  res->fin;
  *pere = N;
  for(i = start,j=deb; ( (i<end) && (j<fin)) ; i++,j++)
    if (Sequence[res->sequence_number][i] != Sequence[current_sequence][j] )
      {
	*restant = j - fin;  /* <0 */
	*pos_in_edge = i - start;
	return res;
      }
  if ((i==end) && (j==fin))
    {
      *restant = 1;
      *pos_in_edge = 0;
      return res;
    }
  
  if (i==end)
    {
      return FindString(res,deb + end - start,fin,pere,restant,pos_in_edge);
    }
  *restant = 3;
  *pos_in_edge = 0;
  return (res);
}



void UpdateBit_TabForAllTree(Noeud *N)
{
  int i,j;
  if (N==NULL)
    return;
  if (N->debut & LEAF_BIT)
    return;
  
  for(i=0;i<ALPHA_CARD;i++)
    UpdateBit_TabForAllTree(N->fils[i]);
  
  j=0;

  while(N->fils[j]==NULL) j++;
  if (N->fils[j]->debut & LEAF_BIT)
    CopyBitTab(&(N->sequences),((Feuille *)N->fils[j])->sequences);
  else
    CopyBitTab(&(N->sequences),N->fils[j]->sequences);
  
  for(i=j+1;i<ALPHA_CARD;i++)
    {
      if (N->fils[i])
	{
	  if (N->fils[i]->debut & LEAF_BIT)
	    fusionneBitTab(&(N->sequences),((Feuille *)N->fils[i])->sequences);
	  else
	    fusionneBitTab(&(N->sequences),N->fils[i]->sequences);
	}
    }
  N->nb_element_bt = nbSequenceInBitTab(N->sequences);
}


void Print_BTTree_Debug(Noeud *N,int *nb_noeud)
{
  int i,j,max;
  if (N==NULL)
    return;
  if (N->debut & LEAF_BIT)
    {
      *nb_noeud = *nb_noeud  + 1;
      printf("%d\t",*nb_noeud);
      max = getValue(Liste_positions_fin,((Feuille *)N)->fin_deb);
      for(j=N->debut & LEAF_BIT_INV;j<max;j++)
	printf("%c",Sequence[N->sequence_number][j]);
      printf("\t");
      printf("%d",nbSequenceInBitTab(((Feuille *)N)->sequences));
      printf("\t");
      printBitTab(((Feuille *)N)->sequences);
      printf("\n");
      return;
    }
  *nb_noeud = *nb_noeud + 1; 
  printf("%d\t",*nb_noeud);
  for(j=N->debut ;j<N->fin;j++)
    printf("%c",Sequence[N->sequence_number][j]);
  printf("\t");
  printf("%d",nbSequenceInBitTab(N->sequences));
  printf("\t");
  printBitTab(N->sequences);
  printf("\n");
  for(i=0;i<ALPHA_CARD;i++)
    Print_BTTree_Debug(N->fils[i],nb_noeud);
}





