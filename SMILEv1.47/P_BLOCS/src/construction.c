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

#include<construction.h>

Noeud *Construction_Arbre(unsigned char *S,int taille_fenetre)
{
  Liste *Debut_liste = Alloc_Liste();
  Liste *Fin_liste =NULL; 
  Noeud *Fin_liste_pere;
  
  int nb_element;
  int fictive;
  
  Noeud *Racine = Alloc_Noeud();
  Feuille *feuille = Alloc_Feuille();
  
  if (!S)
    {
      fprintf(stderr,"Construction_Arbre : Invalide String\nProgram Abord\n");
      exit(-1);
    }

  if (taille_fenetre <=0)
      {
      fprintf(stderr,"Construction_Arbre: taille fenetre = 0!\n");
      exit(-1);
      }
  
  Sequence[0] = S;
  
  Liste_positions_fin = Alloc_ListePositions(strlen((const char *) S)*20);
  
  feuille->debut = 0 | LEAF_BIT;
  Ajoute_Position_Liste(Liste_positions_fin,&(feuille->fin_deb),-1,0);
  
  Fin_liste_pere = Racine;
  
  nb_element = 1;
  Debut_liste->feuille = feuille;
  Fin_liste = Debut_liste;
  
  Ajoute_Fils_Au_Noeud(Racine,(Noeud *)feuille);

  if (taille_fenetre >= strlen((const char *) S))
    taille_fenetre=strlen((const char *) S);
  
  Premiere_Phase(Racine,taille_fenetre,&Debut_liste,&Fin_liste,&Fin_liste_pere,&nb_element,2);
  
  /*----------------------------------------*/
  
  fictive = Deuxieme_Phase(Racine,taille_fenetre,&Debut_liste,&Fin_liste,&Fin_liste_pere,&nb_element,taille_fenetre,0,0);
  
 
  /*----------------------------------------*/
  
  Troisieme_Phase(Racine,taille_fenetre,&Debut_liste,&Fin_liste,&Fin_liste_pere,&nb_element,fictive);
  
  /*----------------------------------------*/
  
#if DEBUG_JTREE
  printf("nb_alloc_noeud = %d ; nb_alloc_feuille = %d; nb_alloc_liste = %d; nb_alloc_tab = %d\n",nb_alloc_noeud,nb_alloc_feuille,nb_alloc_liste,nb_alloc_tab);
  printf("TAILLE TOTALE ALLOUEE : %d\n",nb_alloc_noeud*(sizeof(Noeud)+sizeof(Noeud *)*ALPHA_CARD) +
	 nb_alloc_feuille * (sizeof(Feuille))+
	 nb_alloc_liste * sizeof(Liste)+
	 Liste_positions_fin->tab_size * sizeof(int) * 2+
	 nb_alloc_tab * 2);
#endif

  return Racine;
}

void Premiere_Phase(Noeud *Racine,
		    int taille_fenetre,
		    Liste **debut_liste,
		    Liste **fin_liste,
		    Noeud **fin_liste_pere,
		    int *nb_element_liste,
		    int start_indice)
{
  int i,j;
  int result_type;
  
  Noeud *result;
  Noeud *result_pere;
  Noeud *last_created_node;
  
  /* Construction de l'arbre pour S[0..taille_fenetre] */
  for(i=start_indice;i<taille_fenetre;i++)
    {
/*       Reprise sur le père de la dernière feuille créee. */
      result = *fin_liste_pere;
/*       Initialisation pour le suffix_link */
      last_created_node = NULL;
/*       Valeur de l'indice globale */
      global_indice = i;
      
/* on remet la longeur du dernier segment Noeud - feuille  */
      result_type = -(seg_taille((Noeud *)((*fin_liste)->feuille)));
/* on insere S[*nb_element_liste ... i]. */
      for(j=*nb_element_liste;j<i;j++)
	{
/* si le dernier resultat est nul ou que c'est la Racine de  */
/* l'arbre ou que son lien suffixe est nul ... */
	  if ((!result) || (result == Racine) || (!result->suffixe_link))
	    result = Add_Fast_String(Racine,j,i,&result_type,&result_pere);
	  else
/* sinon on ajoute par lien suffixe... */
	    result = Add_Fast_String(result->suffixe_link,i+result_type,i,&result_type,&result_pere);
	  
/* Si il a eu création d'une feuille  */
/* 1: à un noeud déjà existant */
/* 3: coupure d'un arc avec création d'un noeud. */
	  if ((result_type == 1) || (result_type == 3))
	    {
/* On met l'indice de fin de la feuille crée à -1 (indice global) */
	      setListeValue(Liste_positions_fin,((Feuille *)result)->fin_deb,-1);
/* On ajoute la feuille créee à la liste des feuilles */
	      (*fin_liste)->suiv = Alloc_Liste();
	      (*fin_liste) = (*fin_liste)->suiv;
	      (*fin_liste)->feuille = (Feuille *)result;
/* On positionne la variable fin_liste_pere au pere */
/* de la feuille créee. */
	      *fin_liste_pere = result_pere;
/* On incremente le nombre de feuille dans la liste */
	      *nb_element_liste = *nb_element_liste + 1;
/* Si à l'étape precedente on a crée un noeud (lien suffixe) */
	      if (last_created_node)
		last_created_node->suffixe_link = result_pere;
/* Alors on positionne le lien suffixe sur le pere de la feuille. */
	      
/* on reinitialise la variable last_created_node a NULL */
	      last_created_node = NULL;
/*    Si on a crée un Noeud */
	      if (result_type == 3)
		last_created_node = result_pere;
	    }
/* la chaine est deja dans l'arbre.... */
	  else
	    if (result_type<0)
	      {
/* si il y a un lien suffixe pendent .... */
		if (last_created_node)
/* on le positionne */
		  last_created_node->suffixe_link = result;
/* on le reinitialise. */
		last_created_node = NULL;
/* on s'arrete */
		break;
	      }
	    else
/* sinon car 2 : rien */
	      last_created_node = NULL;
	  
	  if (result_type > 0)
/* Si le resultat est la creation d'une feuille */
/* alors on va repositionner de maniere a reprendre sur un noeud */
/* afin de faire l'ajout suivant par lien suffixe. */
	    { 
	      if (result_type==3)
/* Si il y a creation d'un noeud : on remonte */
/* au pere de celui-ci : il est mis dans le lien suffixe */
/* du noeud crée. et on recalcul la longeur a parcourir */
/* a partir du lien suffixe. */
		{
		  result_type =  - (seg_taille(result_pere) + 1);
		  result = result_pere->suffixe_link;
		  result_pere->suffixe_link = NULL;
		}
	      else
/* sinon : il y a juste eu ajout d'une feuille */
		{
/* la feuille vient d'être créée : elle mesure 1. */
		  result_type = -1;
/* on remet result sur le pere de la feuille. */
		  result = result_pere;
		}
	    }
	}
    }
}

int  Deuxieme_Phase(Noeud *Racine,
		    int taille_fenetre,
		    Liste **debut_liste,
		    Liste **fin_liste,
		    Noeud **fin_liste_pere,
		    int *nb_element_liste,
		    int start_indice,
		    int fict,
		    int ini_res_type)
{
  int i,j;
  int result_type=ini_res_type;
  int taille_sequence;
  int fictive = fict;
  
  Liste *tmp_liste;
  
  Noeud *result;
  Noeud *result_pere;
  Noeud *last_created_node;
  
  Feuille *tmp_feuille;
  Noeud   *tmp_noeud;
  
  taille_sequence = strlen((const char *) Sequence[current_sequence]);
  
/* Construction de l'arbre pour S[k...m-1] */
  for(i=start_indice;i<taille_sequence;i++)
    {
/* Reprise sur le père de la dernière feuille crée. */
      result = *fin_liste_pere;
/* Initialisation pour le suffix_link */
      last_created_node = NULL;
/* Valeur de l'indice globale */
      global_indice = i;
      
/* on remet la longeur du dernier segment Noeud - feuille  */
      if (result_type>0)
	result_type=-result_type;
      else
	result_type = -(seg_taille((Noeud *)((*fin_liste)->feuille))) ;
      
      if (fictive)
/* On reprend sur debut liste qui est une cellule fictive :  */
/* elle ne correspond pas à une feuille crée. mais l'extension  */
/* se fait dans cette feuille... */
	{
/* On verifie que l'extension en cours est bien dans la cellule. */
	  tmp_noeud = result;
	  
	  result = Add_Fast_String(result,i+result_type,i,&result_type,&result_pere);
	  if (result_type == 1)
	    {
	      result_type = -1;
	      result = result_pere;
	    }
	  else
/* On a cree un une feuille et un Noeud */
	    if (result_type==3)
	      {
/* Alors result est la feuille cree. */
/* result_pere est le noeud cree. */
/* result_pere->suffixe_link est le pere noeud cree. */
	     
		result_type = -(seg_taille(result_pere)+1);
		result = result_pere->suffixe_link;
		result_pere->suffixe_link = NULL;
		last_created_node = result_pere;
	      }
	    else
	      if (result_type == 2)
		{
		  result_type = -(seg_taille(result));
		  result      = tmp_noeud;
		}
/* Sinon : On a cherché à inserer une chaine de longueur k */
/* si elle est deja dans l'arbre, elle aboutie a une feuille */
/* on ajoute alors la position i dans la feuille: ce qui est  */
/* effectué par la fonction fast_string. */
/* ????? else result = result->suffixe_link ????? */
	}
      else
/* On fixe l'indice de fin de la premiere cellule de la liste(vrais feuille) à i: */
	{
	  setListeValue(Liste_positions_fin,(*debut_liste)->feuille->fin_deb,i);
	  if (result_type>0)
	    result_type=-result_type;
	  else
	    result_type = -(seg_taille((Noeud *)((*fin_liste)->feuille))) ;
	}
	
      fictive = 0;
/* on insere S[*nb_element_liste ... i]. */
      for(j=*nb_element_liste;j<i;j++)
	{
/* si le dernier résultat est nul ou que c'est la Racine de  */
/* l'arbre ou que son lien suffixe est nul ... */
	   result_pere = NULL;
       
	  if ((!result) || (result == Racine) || (!result->suffixe_link))
	    result = Add_Fast_String(Racine,j,i,&result_type,&result_pere);
	  else
/* sinon on ajoute par lien suffixe... */
	    result = Add_Fast_String(result->suffixe_link,i+result_type,i,&result_type,&result_pere);
/* Si il a eu création d'une feuille  */
/* 1: à un noeud déjà existant */
/* 3: coupure d'un arc avec création d'un noeud. */
	  if ((result_type == 1) || (result_type == 3))
	    {
/* On met l'indice de fin de la feuille crée à -1 (indice global) */
	      setListeValue(Liste_positions_fin,((Feuille *)result)->fin_deb,-1);
/* On ajoute la feuille créee à la liste des feuilles */
	      (*fin_liste)->suiv = Alloc_Liste();
	      (*fin_liste) = (*fin_liste)->suiv;
	      (*fin_liste)->feuille = (Feuille *)result;
/* On positionne la variable fin_liste_pere au pere */
/* de la feuille créee. */
	      *fin_liste_pere = result_pere;
/* On incremente le nombre de feuille dans la liste */
	      *nb_element_liste = *nb_element_liste + 1;
/* Si à l'étape precedente on a crée un noeud (lien suffixe) */
	      if (last_created_node)
		last_created_node->suffixe_link = result_pere;
/* Alors on positionne le lien suffixe sur le pere de la feuille. */
	      
/* on réinitialise la variable last_created_node a NULL */
	      last_created_node = NULL;
/* Si on a crée un Noeud */
	      if (result_type == 3)
		last_created_node = result_pere;
	    }
/* la chaine est deja dans l'arbre.... */
	  else
	    if (result_type<0)
	      {
/* si il y a un lien suffixe pendent .... */
		if (last_created_node)
		  last_created_node->suffixe_link = result; /* on le positionne */
		last_created_node = NULL; /* on le reinitialise. */
		
/* On doit verifier que l'on n'a suffisement avancé dans la liste */
/* .a.d. que *nb_element_liste >= (i+1) - k. */
		if (*nb_element_liste==(i-taille_fenetre+1))
		  {
/* 		    // On sait que le resultat est le pere d'une feuille */
/* 		    // car si *nb_element=i-taille_sequence alors la taille */
/* 		    // de la derniere chaine cherchée est k. comme l'arbre */
/* 		    // se coupe a la hauteur k... */
		    
/* 		    // De plus on sait que qu'il n'y a plus d'element dans la liste */
/* 		    // car on l'indente d'une fois au fur et à mesure que  */
/* 		    // l'on progresse.... */
		    tmp_feuille = (Feuille *)Get_Child_Start_Letter(result,i+result_type);

		    if (-seg_taille((Noeud *)tmp_feuille)==result_type)
		      result_type = 1;
		    if (tmp_feuille->debut & LEAF_BIT)
		      {
			(*fin_liste_pere) = result;
			(*fin_liste)->suiv = Alloc_Liste();
			(*fin_liste) = (*fin_liste)->suiv;
			(*fin_liste)->feuille = tmp_feuille;
			*nb_element_liste = *nb_element_liste + 1;
		      }
		    else
		      {
			(*fin_liste_pere) = (Noeud *)tmp_feuille;
			(*fin_liste)->suiv = Alloc_Liste();
			(*fin_liste) = (*fin_liste)->suiv;
			(*fin_liste)->feuille = tmp_feuille; /* FAUX */
			*nb_element_liste = *nb_element_liste + 1;
		      }
		    fictive = 1; /* On indique pour la prochaine reprise... */
/* 		    // Probleme : lors de la reprise au i suivant on va  */
/* 		    // initialise result a ce noeud et on va compter  */
/* 		    // la longeur N-Feuille pour la longeur du suffixe link. */
		  }
		
		break; /* on s'arrete */
	      }
	    else
	      last_created_node = NULL; /* sinon car 2 : rien */
	  
	  if (result_type > 0) /* Si le resultat est la creation d'une feuille */
/* 	    // alors on va repositionner de maniere a reprendre sur un noeud */
/* 	    // afin de faire l'ajout suivant par lien suffixe. */
	    { 
	      if (result_type==3) /* Si il y a creation d'un noeud : on remonte */
/* 		// au pere de celui-ci : il est mis dans le lien suffixe */
/* 		// du noeud crée. et on recalcul la longeur a parcourir */
/* 		// a partir du lien suffixe. */
		{
		  result_type =  - (seg_taille(result_pere) + 1);
		  result = result_pere->suffixe_link;
		  result_pere->suffixe_link = NULL;
		}
	      else /* sinon : il y a juste eu ajout d'une feuille */
		if (result_type==1)
		  {
		    result_type = -1; /* la feuille vient d'être créée : elle mesure 1*/
		    result = result_pere; /* on remet result sur le pere de la feuille*/
		  }
		else /* cas 2. */
		  {
		    result_type = - seg_taille(result);
		    result = result_pere;
		  }
	    }
	}
/*       // On avance dans la liste chainée de cellules: */
      {
	tmp_liste = (*debut_liste);
	(*debut_liste) = (*debut_liste)->suiv;
	Free_Liste(tmp_liste);
      }
    }
  return fictive;
}


void Troisieme_Phase(Noeud *Racine,
		     int taille_fenetre,
		     Liste **debut_liste,
		     Liste **fin_liste,
		     Noeud **fin_liste_pere,
		     int *nb_element_liste,
		     int fictive)
{
  int i,j;/*,lm; */
  int result_type;
  int taille_sequence;
  
  Noeud *result;
  Noeud *result_pere=NULL;
  Noeud *last_created_node;
  
  Noeud *tmp_noeud;
  Noeud *tmp_noeud2;

  Liste * tmp_liste;
  
  taille_sequence = strlen((const char *) Sequence[current_sequence]);
  
  i=taille_sequence; /* Construction de l'arbre pour S[taille_fenetre] */
  result = *fin_liste_pere; /* Reprise sur le père de la dernière feuille créee. */
  last_created_node = NULL; /* Initialisation pour le suffix_link */
  global_indice = i;        /* Valeur de l'indice globale */
  
  
  /* on remet la longeur du dernier segment Noeud - feuille  */
  if ((*fin_liste)->feuille->debut & LEAF_BIT)
    result_type = -(seg_taille((Noeud *)((*fin_liste)->feuille)));
  else
    result_type = -1;
  
  if (fictive)
    {
      result_pere = NULL;
      result = Add_Fast_String(result,i+result_type,i,&result_type,&result_pere);
      if (result_type==3) /* On a cree un une feuille et un Noeud */
	{
	  last_created_node = result_pere;
	  result_type = -(seg_taille(result_pere)+1);
	  result = result_pere->suffixe_link;
	  result_pere->suffixe_link = NULL;
	}
      if (result_type==1)/* on a juste cree une feuille */
	{
	  result_type = -1;
	  result = result_pere;
	}
    }
  else
    {
      while((*debut_liste))
	{
	  setListeValue(Liste_positions_fin,(*debut_liste)->feuille->fin_deb,i);
	  if ((*debut_liste)==(*fin_liste))
	    break;
	  tmp_liste = *debut_liste;
	  *debut_liste = (*debut_liste)->suiv;
	  Free_Liste(tmp_liste);
	}
      if ((*fin_liste)->feuille->debut & LEAF_BIT)
	result_type = -(seg_taille((Noeud *)((*fin_liste)->feuille)));
      else
	result_type = -1;
      Free_Liste((*fin_liste));
    }
  /* on insere S[*nb_element_liste ... i]. */
  for(j=*nb_element_liste;j<i;j++)
    {
      /* si le dernier resultat est nul ou que c'est la Racine de  */
      /* l'arbre ou que son lien suffixe est nul ... */
      result_pere=NULL;
      if ((!result) || (result == Racine) || (!result->suffixe_link))
	result = Add_Fast_String(Racine,j,i,&result_type,&result_pere);
      else
	/* sinon on ajoute par lien suffixe... */
	result = Add_Fast_String(result->suffixe_link,i+result_type,i,&result_type,&result_pere);
      /* Si il a eu création d'une feuille  */
      /* 1: à un noeud déjà existant */
      /* 3: coupure d'un arc avec création d'un noeud. */
      if ((result_type == 1) || (result_type == 3))
	{
	  /* On met l'indice de fin de la feuille crée à -1 (indice global) */
	  /* On ajoute la feuille créee à la liste des feuilles */
	  /* ATTENTION : A VERIFIER...
	  (*fin_liste)->suiv = Alloc_Liste();
	  (*fin_liste) = (*fin_liste)->suiv;
	  (*fin_liste)->feuille = (Feuille *)result;
	  */
	  /* On positionne la variable fin_liste_pere au pere */
	  /* de la feuille créee. */
	  *fin_liste_pere = result_pere;
	  /* On incremente le nombre de feuille dans la liste */
	  *nb_element_liste = *nb_element_liste + 1;
	  /* Si à l'étape precedente on a crée un noeud (lien suffixe) */
	  if (last_created_node)
	    last_created_node->suffixe_link = result_pere;
	  /* Alors on positionne le lien suffixe sur le pere de la feuille. */
	  
	  /* on reinitialise la variable last_created_node a NULL */
	  last_created_node = NULL;
	  if (result_type == 3) /* Si on a crée un Noeud */
	    last_created_node = result_pere;
	}
      else /* la chaine est deja dans l'arbre.... */
	if (result_type<0)
	  { 
	    /* on ne s'arrete pas : */
	    /* On coupe l'arc ... et on "emule" l'insertion d'une feuille... */
	    tmp_noeud2 = Get_Child_Start_Letter(result,i+result_type);
	    if (seg_taille(tmp_noeud2)>-result_type)
	      {
		/* si il y a un lien suffixe pendent .... */
		tmp_noeud = Alloc_Noeud();
		if (last_created_node)
		  last_created_node->suffixe_link = tmp_noeud; /* on le positionne */
		last_created_node = NULL;/* on le reinitialise. */
		
		tmp_noeud->debut = i + result_type;
		tmp_noeud->fin = i;
		if (tmp_noeud2->debut & LEAF_BIT) /* result est un feuille */
		  tmp_noeud2->debut = (tmp_noeud2->debut - result_type) | LEAF_BIT;
		else
		  tmp_noeud2->debut = (tmp_noeud2->debut - result_type);
		Ajoute_Fils_Au_Noeud(result,tmp_noeud);
		Ajoute_Fils_Au_Noeud(tmp_noeud,tmp_noeud2);
		last_created_node = tmp_noeud;
		result = result_pere;
	      }
	    else
	      {
		/* si il y a un lien suffixe pendent .... */
		if (last_created_node)
		  last_created_node->suffixe_link = result; /* on le positionne */
		last_created_node = NULL;/* on le reinitialise. */
		/*    ATTENTION : PAS SUR
		  if (tmp_noeud2->debut & LEAF_BIT)
		  &&()
		  {
		  tmp_noeud2->debut = (i + result_type) | LEAF_BIT;
		  if ((getValue(Liste_positions_fin,((Feuille *)tmp_noeud2)->fin_deb) != i )
		  || (((Feuille *)tmp_noeud2)->sequence_number != current_sequence))
		  Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)tmp_noeud2)->fin_deb),i,0);
		  }
		  else
		  {
		  tmp_noeud2->debut = i + result_type;
		  tmp_noeud2->fin = i;
		  }
		*/
	      }
	  }
	else
	  last_created_node = NULL; /* sinon car 2 : rien */
      
      if (result_type > 0) /* Si le resultat est la creation d'une feuille */
	/* alors on va repositionner de maniere a reprendre sur un noeud */
	/* afin de faire l'ajout suivant par lien suffixe. */
	{ 
	  if (result_type==3) /* Si il y a creation d'un noeud : on remonte */
	    /* au pere de celui-ci : il est mis dans le lien suffixe */
	    /* du noeud crée. et on recalcul la longeur a parcourir */
	    /* a partir du lien suffixe. */
	    {
	      result_type =  - (seg_taille(result_pere) + 1);
	      result = result_pere->suffixe_link;
	      result_pere->suffixe_link = NULL;
	    }
	  else /* sinon : il y a juste eu ajout d'une feuille */
	    {
	      result_type = -1; /* la feuille vient d'être créée : elle mesure 1. */
	      result = result_pere; /* on remet result sur le pere de la feuille. */
	    }
	}
    }
}

Noeud *AjouteSequence(Noeud *Arbre,unsigned char *S,int taille_fenetre)
{
  Noeud *pere = NULL;
  Noeud *resultat;

  int position_arc;
  int seg_lg;
  int res_type;
  int i=taille_fenetre,j;
  int lg_sequence;
  int decalage = taille_fenetre;
  
  current_sequence++;
  
  Sequence[current_sequence] = S;
  lg_sequence = strlen((const char *) S);
  if (taille_fenetre>=lg_sequence)
    taille_fenetre = lg_sequence;

  i=taille_fenetre;
  decalage = taille_fenetre;
  

  /* Recherché la chaine S[i..i+taille_fenetre] dans l'arbre  */
  /* pour i allant de 0 à lg(S)-taille_fenetre */
  resultat = Arbre;
  position_arc = 0;
  while(i<=lg_sequence)
    {
      pere=NULL;
      resultat = FindString(resultat,i-decalage,i,&pere,&res_type,&position_arc);
      if (position_arc == -1)
	{
	  /*resultat  est un noeud et n'a pas de fils pour la chaine i+res_type ... i */
	  return CaseTreeAddSequence(Arbre,resultat,i+res_type,i+res_type+1,taille_fenetre);
	}
      if (res_type == 1)
	{
	  seg_lg = seg_taille(resultat);
	  if (resultat->debut & LEAF_BIT)
	    {
	      resultat->debut = (i-seg_lg) | LEAF_BIT;
	      Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)resultat)->fin_deb),i,(resultat->sequence_number==current_sequence)?0:1);
	      resultat->sequence_number = current_sequence;
	      addBitTabValue(&(((Feuille *)resultat)->sequences),current_sequence);
	    }
	  i++;
	  decalage = seg_taille(resultat) + 1;
	  if (pere!=Arbre)
	    resultat = pere->suffixe_link;
	  else
	    {
	      decalage = taille_fenetre;
	      resultat = Arbre;
	    }
	}
      else
	if (res_type == 2) 
	  {
	    /* En principe impossible. */
	    printf("ERREUR    CAS 2 for sequence %d and string %d .. %d \n",current_sequence,i-decalage,i);
	    for(j=i-decalage;j<i;j++)
	      printf("%c",Sequence[current_sequence][j]);
	    printf("\n");
	    exit(0);
	  }
	else
	  if (res_type == 3)
	    {
	      /* En principe impossible. */
	    }
	  else
	    return CaseTreeAddSequence(Arbre,pere,i+res_type-position_arc,i+res_type+1,taille_fenetre);
    }
  
  i = lg_sequence;
  for(j=lg_sequence-taille_fenetre; j<lg_sequence;j++)
    {
      resultat = FindString(Arbre,j,i,&pere,&res_type,&position_arc);
      if (res_type == 1)
	{
	  seg_lg = seg_taille(resultat);
	  if (resultat->debut & LEAF_BIT)
	    {
	      if ((resultat->sequence_number != current_sequence) || 
		  (
		   (resultat->sequence_number == current_sequence)
		   && (resultat->debut&LEAF_BIT_INV) !=(i-seg_lg)
		   )
		  )
		{
		  resultat->debut = (i-seg_lg) | LEAF_BIT;
		  Ajoute_Position_Liste(Liste_positions_fin,&(((Feuille *)resultat)->fin_deb),i,(resultat->sequence_number==current_sequence)?0:1);
		  resultat->sequence_number = current_sequence;
		  addBitTabValue(&(((Feuille *)resultat)->sequences),current_sequence);
		}
	    }
	}
      else 
	printf("YOUPI II I I I %d \n",res_type);
    }
  return Arbre;
}


/* cas creation d'une feuille à la racine de l'arbre. */
Noeud *CaseOneAddSequence(Noeud *Arbre,int taille_fenetre)
{
  Liste *Debut_liste = Alloc_Liste(); 
  Liste *Fin_liste;
  Noeud *Fin_liste_pere;
  int nb_el_liste=0;
  int fictive = 0;
  
  Feuille *tmp_f;
  
  fprintf(stderr,"CASE ONE\n");
  
  tmp_f = Alloc_Feuille();
  tmp_f->debut = 0 | LEAF_BIT;
  Ajoute_Position_Liste(Liste_positions_fin,
			&(tmp_f->fin_deb),-1,0);
  Fin_liste_pere = Arbre;
  nb_el_liste = 1;
  Debut_liste->feuille = tmp_f;
  Fin_liste = Debut_liste;
  Ajoute_Fils_Au_Noeud(Arbre,(Noeud *)tmp_f);
  
  Premiere_Phase(Arbre,taille_fenetre,
		 &Debut_liste,&Fin_liste,&Fin_liste_pere,
		 &nb_el_liste,1);
  fictive = Deuxieme_Phase(Arbre,taille_fenetre,
			   &Debut_liste,&Fin_liste,&Fin_liste_pere,
			   &nb_el_liste,taille_fenetre,
			   0,0);
  Troisieme_Phase(Arbre,taille_fenetre,
		  &Debut_liste,&Fin_liste,&Fin_liste_pere,
		  &nb_el_liste,fictive);
  return Arbre;
}

/* cas de reprise au debut de 2eme phase */
Noeud *CaseTwoAddSequence(Noeud *Arbre,Noeud *resultat, Noeud *pere,int taille_fenetre)
{
  
  Liste *Debut_liste = Alloc_Liste(); 
  Liste *Fin_liste;
  Noeud *Fin_liste_pere;
  int nb_el_liste=0;
  int fictive = 0;
  
  fprintf(stderr,"CASE TWO\n");
  
  Debut_liste->feuille = (Feuille *)pere; 
  Fin_liste = Debut_liste;
  Fin_liste_pere = resultat;
  nb_el_liste = 1;
  
  fictive = Deuxieme_Phase(Arbre,taille_fenetre,&Debut_liste,
			   &Fin_liste,&Fin_liste_pere,&nb_el_liste,
			   taille_fenetre,1,0);
  
  Troisieme_Phase(Arbre,taille_fenetre,
		  &Debut_liste,&Fin_liste,&Fin_liste_pere,
		  &nb_el_liste,fictive);
  
  return Arbre;
}

/* cas de reprise en milieu de 1er phase */
Noeud *CaseTreeAddSequence(Noeud *Arbre,Noeud *resultat,int deb,int fin,int taille_fenetre)
{  
  
  Liste *Debut_liste = Alloc_Liste(); 
  Liste *Fin_liste=NULL;
  Noeud *Fin_liste_pere;
  Noeud *pere=NULL;
  
  int nb_el_liste=0;
  int fictive = 0;
  int j,lim,res_type;
  int start;
  
  Noeud *last_created = NULL;

  
  /* Dans  ce cas on doit faire un tour d'algo à la main..... */

  lim = fin;
  res_type = deb - fin;
  global_indice = lim;
  
  start = lim-taille_fenetre;
  if (start<0)
    start=0;
  for (j=start;j<lim;j++)
    {
      pere=NULL;
      resultat = Add_Fast_String(resultat,lim + res_type,lim,&res_type,&pere);
      
      if ((res_type == 1) || (res_type ==3))
	{
	  setListeValue(Liste_positions_fin,(((Feuille *)resultat)->fin_deb),-1);
	  if (Fin_liste)
	    {
	      Fin_liste->suiv = Alloc_Liste();
	      Fin_liste = Fin_liste->suiv;
	      Fin_liste->feuille = (Feuille *)resultat;
	    }
	  else 
	    {
	      Fin_liste = Debut_liste;
	      Debut_liste->feuille = (Feuille *)resultat;
	    }
	  
	  Fin_liste_pere = pere;
	  nb_el_liste++;
	  if (last_created)
	    last_created->suffixe_link = pere;
	  last_created = NULL;
	  if (res_type == 3)
	    last_created = pere;
	}
      else
	if (res_type<0)
	  {
	    if (last_created)
	      last_created->suffixe_link = resultat;
	    last_created = NULL;
	    break;
	  }
	else
	  last_created = NULL;
      if (res_type>0)
	{
	  if (res_type==3)
	    {
	      res_type = - seg_taille(pere) - 1;
	      resultat = pere->suffixe_link;
	      if (resultat!=Arbre)
		resultat = resultat->suffixe_link;
	      else
		res_type += 1;
	      pere->suffixe_link = NULL;
	    }
	  else
	    {
	      res_type = -1;
	      resultat = pere;
	      if (resultat == NULL)
		{
		  resultat = Arbre;
		  res_type = -(j+1);
		}
	      else
		if (resultat!=Arbre)
		  resultat = resultat->suffixe_link;
		else
		  res_type += 1;
	    }
	}
    }
  

  if (fin <= taille_fenetre)
    {
      Premiere_Phase(Arbre,taille_fenetre,
		     &Debut_liste,&Fin_liste,&Fin_liste_pere,
		     &nb_el_liste,lim+1);
      fictive = Deuxieme_Phase(Arbre,taille_fenetre,
			       &Debut_liste,&Fin_liste,&Fin_liste_pere,
			       &nb_el_liste,taille_fenetre,0,0);
      Troisieme_Phase(Arbre,taille_fenetre,
		      &Debut_liste,&Fin_liste,&Fin_liste_pere,
		      &nb_el_liste,fictive);
    }
  else
    {
      setListeValue(Liste_positions_fin,Debut_liste->feuille->fin_deb,lim);
      if (nb_el_liste == 1)
	fictive = 1;
      else
	{
	  res_type = -1;
	  Debut_liste = Debut_liste->suiv;
	  nb_el_liste--;
	}
      
      nb_el_liste = ((fictive)?lim:lim+1)-taille_fenetre+nb_el_liste;
      fictive = Deuxieme_Phase(Arbre,taille_fenetre,
			       &Debut_liste,&Fin_liste,&Fin_liste_pere,
			       &nb_el_liste,(fictive)?lim:lim+1,fictive,0);

      Troisieme_Phase(Arbre,taille_fenetre,
		      &Debut_liste,&Fin_liste,&Fin_liste_pere,
		      &nb_el_liste,fictive);
    }
  
  return Arbre;
}


Noeud *CaseFourAddSequence(Noeud *Arbre,Noeud *resultat,Noeud *pere,int res_type,int position_arc,int i,int taille_fenetre)
{
  
  Liste *Debut_liste = Alloc_Liste(); 
  Liste *Fin_liste;
  Noeud *Fin_liste_pere;

  int nb_el_liste=0;
  int fictive = 0;
  
  fprintf(stderr,"\nCASE FOUR------------------------------------------------------------------\n\n");
  /* On doit faire un tour d'algo avant de */
  /* Rentrer dans les fonctions standard... */

  
  global_indice = i;
  
  
  
  if (res_type != -1)
    {
      Debut_liste->feuille = (Feuille *)resultat;
      Fin_liste = Debut_liste;
      Fin_liste_pere = pere;
      nb_el_liste = i-taille_fenetre+1;
    }
  else
    {
      Debut_liste->feuille = (Feuille *)resultat;
      Fin_liste = Debut_liste;
      Fin_liste_pere = resultat;
      nb_el_liste = i-taille_fenetre+1;
    }
  
    printf("Arbre AVANT la Deuxieme phase\n");
    printf("Arbre = %p \n taille_fenetre = %d \n"
    "Debut_liste = %p ->feuille %p \n"
    "Fin_liste = %p ->feuille %p\n"
    "Fin_liste_pere = %p\n"
    "nb_element_liste=%d\n"
    "i = %d\n"
    "Fictive = 0\n", Arbre,taille_fenetre,Debut_liste,Debut_liste->feuille,Fin_liste,Fin_liste->feuille,Fin_liste_pere,nb_el_liste,i);
  
  
  fictive = Deuxieme_Phase(Arbre,taille_fenetre,&Debut_liste,
			   &Fin_liste,&Fin_liste_pere,&nb_el_liste,
			   i,1,0); 
  Troisieme_Phase(Arbre,taille_fenetre,&Debut_liste,
		  &Fin_liste,&Fin_liste_pere,&nb_el_liste,fictive);
  
  return Arbre;
  
}

void CloseTheFirstPhase( Liste **debut_liste,
			   Liste **fin_liste,
			   Noeud **fin_liste_pere)
{
  Liste *tmp = *debut_liste;

  printf("CLOSE THE FIRST CASE\n");
  while(tmp!=NULL)
    {
      if ((tmp->feuille) && (tmp->feuille->debut & LEAF_BIT))
	setListeValue(Liste_positions_fin,tmp->feuille->fin_deb,global_indice);
      tmp=tmp->suiv;
    }
  
}
