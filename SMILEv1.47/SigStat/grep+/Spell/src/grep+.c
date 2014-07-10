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

#include <grep+.h>


/******************************************************************************/
/*                          PROTOTYPES PRIVES                                 */
/******************************************************************************/
/* Gestion des modeles acceptes                                               */
void    keepModel(P_PileOcc, LongSeq, P_Criteres cr, unsigned int *nboccex,
            FILE *f);

/* essaie d'avancer d'une lettre dans un arc, et renvoie le noeud image */
Flag avanceBranche(P_occ, P_occ, int, int, Flag, P_Criteres, LongSeq, Flag);

/* Lancement du saut */
NbSeq gestionSaut(P_PileOcc pocc, P_Criteres, LongSeq curbloc,
        signed char **text);

/* Charge la sequence dans l'arbre a partir d'un fichier FASTA                */
Flag    chargeSequence(Arbre *a, char *fic);

/* Construction de l'arbre                                                    */
Flag    creeArbreSuffixe(Arbre *a, int maxlongmod, char *alphaseq);


/* EXTERNES from alphabet.c                                                   */
extern  int     nbSymbMod;
extern  int     nbSymbSeq;
extern  char    *nummod2str[127];
extern  int     carseq2num[127];
extern  Flag    TabSymb[127][127];
extern  int     numJOKER;
extern  int     numSAUT;


/******************************************************************************/
/******************************************************************************/
/************************ FONCTIONS DE BASE ***********************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/******************************************************************************/
/*********************** GESTION DES LISTES D'OCCURRENCES *********************/
/*********************************ET DES MODELES*******************************/
/******************************************************************************/

/******************************************************************************/
/* KeepModel                                                                  */
/******************************************************************************/
/* Affiche (ou stocke si necessaire) les modeles trouves                      */
/******************************************************************************/
void keepModel(P_PileOcc pocc, LongSeq l, P_Criteres cr, unsigned int *nboccex,
        FILE *f)
{
#if OCC
    afficheLastOcc(f, pocc, l, cr, nboccex);
#endif
}


/******************************************************************************/
/******************************************************************************/
/************************* RECHERCHE DES MODELES ******************************/
/******************************************************************************/
/******************************************************************************/


/******************************************************************************/
/* avanceBranche                                                              */
/******************************************************************************/
/* Essaie d'avancer d'une lettre dans un arc.                                 */
/* Renvoie 1 si reussi, 0 sinon.                                              */
/* La variable 'flag' indique si on est sur un noeud(1) ou une branche(0)     */
/******************************************************************************/
Flag avanceBranche( P_occ   next,       P_occ       tmp,    int symbol,
                    int     trans,      Flag        flag_noeud,
                    P_Criteres  cr,     LongSeq curbloc,
                    Flag    multiblocs)
{
/* Dans cette fonction, le code est duplique dans un souci de rapidite:       */
/* j'essaie de faire un max de tests eliminatoires avant affectations         */

/* Si la branche courante n'est pas epuisee... */
if (flag_noeud == FAUX)
    {
    if ( equiv(symbol, trans) )
        {
        next->xerr    = tmp->xerr;
        next->blocerr = tmp->blocerr;
        }
    else
        {
        next->xerr = tmp->xerr+1;
        
        if (next->xerr == cr->maxerr+1) /* si maxerr global atteint */
            return 0;
        
        if(multiblocs == VRAI)
            {
            next->blocerr = tmp->blocerr+1; /* si maxerr local atteint */
            if (next->blocerr == cr->maxerrblocs[curbloc]+1)
                return 0;
            }
        }
    
    next->x   = tmp->x;
    next->num = tmp->num;
    next->lon = tmp->lon+1;
    }
else
/* Si la branche courante est epuisee, on est sur une nouvelle branche */
    {
    next->x = tmp->x->fils[tmp->num];
    
    if (next->x->fils[trans] == NULL)
        return(0);
    
    if ( equiv(symbol, trans) )
        {
        next->xerr    = tmp->xerr;
        next->blocerr = tmp->blocerr;
        }
    else
        {
        next->xerr = tmp->xerr+1;
        
        if (next->xerr == cr->maxerr+1) /* si maxerr global atteint */
            return 0;
        
        if(multiblocs == VRAI)
            {
            next->blocerr = tmp->blocerr+1; /* si maxerr local atteint */
            if (next->blocerr == cr->maxerrblocs[curbloc]+1)
                return 0;
            }
        }
        
    next->num = trans;
    next->lon = 1;
    }


if(multiblocs == VRAI)
    {
    next->saut= tmp->saut;
    next->codesaut= tmp->codesaut;
    }

return(1);
}


/******************************************************************************/
/* sauteSymbole                                                               */
/******************************************************************************/
int sauteSymbole(Occ curocc, P_PileOcc pocc, P_Criteres cr, LongSeq curbloc,
        LongSeq longsaut, signed char **text)
{
LongSeq	lmaxbr;
Noeud	*tmpnoeud;
Occ		tmpocc;
int		res = 0;
int     trans;
char    carseq;


tmpnoeud = curocc.x->fils[curocc.num];

if (tmpnoeud->debut & LEAF_BIT)
   lmaxbr  = getValue(Liste_positions_fin,((Feuille *)tmpnoeud)->fin_deb)
       - (((Feuille *)tmpnoeud)->debut & LEAF_BIT_INV); 
else
   lmaxbr = tmpnoeud->fin - tmpnoeud->debut;

#if DEBUG_SAUT
printf("SauteSymbole: j'ai gere le saut pour %d, noeud %d, etat: %d/%d branche %d\n",longsaut,curocc.x,curocc.lon,lmaxbr,curocc.num);
printf("saut %d longsaut %d\n",curocc.saut,longsaut);
#endif

ajouteOcc2Pile(pocc, curocc.x, curocc.num, curocc.lon, curocc.xerr, 0, 
        curocc.saut+longsaut, addSaut2Code(curocc.codesaut, longsaut, curbloc,
        cr));
res++;
longsaut++;

if (curocc.lon != lmaxbr)  /* on est au milieu d'une branche */
	{
	curocc.lon++;

	carseq  = text[tmpnoeud->sequence_number]
            [(tmpnoeud->debut & LEAF_BIT_INV)+curocc.lon-1];

	if(carseq==FINAL)	/* si on rencontre un FINAL c'est fini */
		return res;

	if(longsaut<=cr->saut[curbloc].max)
		res += sauteSymbole(curocc, pocc, cr, curbloc, longsaut, text);
	}
else	/* sinon on est a un noeud, plusieurs trans sont possibles */
	{
	if(longsaut<=cr->saut[curbloc].max)
		{
		tmpocc.x    = tmpnoeud;
		tmpocc.lon  = 1;
		tmpocc.xerr = curocc.xerr;
        tmpocc.codesaut = curocc.codesaut;
        tmpocc.saut = curocc.saut;

		if ((tmpnoeud->debut & LEAF_BIT) == 0)
		  for (trans = 0; trans != nbSymbSeq; trans++)
			{
			if (tmpnoeud->fils[trans] != NULL)
				{
				tmpocc.num  = trans;

				res += sauteSymbole(tmpocc, pocc, cr, curbloc, longsaut,
                                    text);
				}
			}
		}
	}
#if DEBUG_SAUT
    printf("SauteSymbole: J'ai trouve %d occ\n",res);
#endif
return res;
}

/******************************************************************************/
/* sauteBranche                                                               */
/******************************************************************************/
int sauteBranche(Occ curocc, P_PileOcc pocc, P_Criteres cr, LongSeq curbloc,
        LongSeq longsaut, signed char **text)
{
LongSeq	lmaxbr;
Noeud *	tmpnoeud, *newtmpnoeud;
Occ		tmpocc;
int		res = 0, newlongsaut;
int	    trans;
char    carseq;



tmpnoeud = curocc.x->fils[curocc.num];

if (tmpnoeud->debut & LEAF_BIT)
   lmaxbr  = getValue(Liste_positions_fin,((Feuille *)tmpnoeud)->fin_deb) - (((Feuille *)tmpnoeud)->debut & LEAF_BIT_INV); 
else
   lmaxbr = tmpnoeud->fin - tmpnoeud->debut;

#if DEBUG_SAUT
printf("SauteBranche: j'ai gere le saut pour %d, noeud %d, etat: %d/%d branche %d\n",longsaut,curocc.x,curocc.lon,lmaxbr,curocc.num);
#endif

if (curocc.lon != lmaxbr)  /* on est au milieu d'une branche */
	{
	if ( lmaxbr-curocc.lon <= cr->saut[curbloc].min-longsaut )
		{	
		longsaut+=lmaxbr-curocc.lon;
		curocc.lon=lmaxbr;

#if DEBUG_SAUT
		printf("SauteBranche: milieuBr, fast, je vais au bout %d/%d br %d et lgsaut %d\n",curocc.lon,lmaxbr,curocc.num,longsaut);
#endif

		carseq  = text[tmpnoeud->sequence_number]
                [(tmpnoeud->debut & LEAF_BIT_INV)+lmaxbr-1];

		if(carseq!=FINAL)	/* si on rencontre un FINAL c'est fini */
			{
#if DEBUG_SAUT
			printf("SauteBranche: finBr=$, c'est fini\n");
#endif
			res += sauteBranche(curocc, pocc, cr, curbloc, longsaut, text);
			}
		}
	else
		{
		curocc.lon+=cr->saut[curbloc].min-longsaut;
		longsaut=cr->saut[curbloc].min;

#if DEBUG_SAUT
		printf("SauteBranche: milieuBr, minsaut ds Br, je m'arrete a %d/%d num %d et lgsaut %d\n",curocc.lon,lmaxbr,curocc.num,longsaut);
#endif

		res += sauteSymbole(curocc, pocc, cr, curbloc, longsaut, text);
		}
	}
else	/* sinon on est a un noeud, plusieurs trans sont possibles */
	{
	tmpocc.x    = tmpnoeud;
	tmpocc.xerr = curocc.xerr;
    tmpocc.codesaut = curocc.codesaut;
    tmpocc.saut = curocc.saut;
    
	if ((tmpnoeud->debut & LEAF_BIT) == 0)
	  for (trans = 0; trans != nbSymbSeq; trans++)
		{
		tmpocc.num  = trans;
		newlongsaut = longsaut;

		if (tmpnoeud->fils[trans] != NULL)
			{
			newtmpnoeud	= tmpnoeud->fils[trans];

			if (newtmpnoeud->debut & LEAF_BIT)
    		    lmaxbr  = getValue(Liste_positions_fin,
                        ((Feuille *)newtmpnoeud)->fin_deb)
                         - (newtmpnoeud->debut & LEAF_BIT_INV); 
			else
			    lmaxbr = newtmpnoeud->fin - newtmpnoeud->debut;
			
			if ( lmaxbr <= cr->saut[curbloc].min-longsaut )
				{
				newlongsaut+=lmaxbr;
				tmpocc.lon=lmaxbr;

#if DEBUG_SAUT
				printf("SauteBranche: noeud, fast, %d/%d, br %d, lgsaut %d\n",tmpocc.lon,lmaxbr,tmpocc.num,newlongsaut);
#endif

				carseq = text[newtmpnoeud->sequence_number]
                        [(newtmpnoeud->debut & LEAF_BIT_INV)+lmaxbr-1];

				if(carseq!=FINAL)	/* si on rencontre un FINAL c'est fini */
					{
#if DEBUG_SAUT
					printf("SauteBranche: finBr=$, c'est fini\n");
#endif
					res += sauteBranche(tmpocc, pocc, cr, curbloc,
                            newlongsaut, text);
					}
				}
			else
				{
				tmpocc.lon=cr->saut[curbloc].min-newlongsaut;
				newlongsaut=cr->saut[curbloc].min;

#if DEBUG_SAUT
				printf("SauteBranche2: noeud %d, minsaut ds Br, %d/%d, br %d, lgsaut %d\n",tmpocc.x, tmpocc.lon,lmaxbr,tmpocc.num,newlongsaut);
#endif

				res += sauteSymbole(tmpocc, pocc, cr, curbloc,
                        newlongsaut, text);
				}
			}
		}
	}
#if DEBUG_SAUT
    printf("SauteBranche: J'ai trouve %d occ\n",res);
#endif
return res;
}

/******************************************************************************/
/* gestionSaut                                                                */
/******************************************************************************/
NbSeq gestionSaut(P_PileOcc pocc, P_Criteres cr, NbSeq curbloc,
        signed char **text)
{
LongSeq	pos, precdummy;
Occ		curocc;
P_occ   tmpocc;
int		res = 0;

pos		= pocc->pos-1;
precdummy = getPrecDummy(pocc);
tmpocc  = pocc->occ+pos;

ajouteDummy(pocc);

while ( (pos != precdummy) && (tmpocc->x != NULL) )
	{
	curocc = *tmpocc;

	if (cr->saut[curbloc].min == 0)
		res+=sauteSymbole(curocc, pocc, cr, curbloc, 0, text);
	else
		res += sauteBranche(curocc, pocc, cr, curbloc, 0, text);
	pos--;
    tmpocc  = pocc->occ+pos;
	}

/* if(res==0) */
/* 	depileRec(pocc); */

return (res);
}

/******************************************************************************/
/* sommeBTOcc                                                                 */
/******************************************************************************/
/* Fait l'union des sequences d'une liste d'occurrence et renvoie le nombre   */
/* de ces sequences.                                                          */
/******************************************************************************/
NbSeq sommeBTOcc(P_PileOcc p, Bit_Tab ** bt)
{
LongSeq pos, precdummy;
P_occ   po;

ReinitBitTab(bt);
pos = p->pos-1;
if(pos < 0)
    fatalError("grep+.c: sommeBTOcc: wrong stack position\n");
po  = p->occ+pos;
precdummy = getPrecDummy(p);

while ((pos != precdummy) && (po->x != NULL))
    {
#if  DEBUG_BT
    printf("Fusion avec :    ");
#endif

    if(po->x->fils[po->num]->debut & LEAF_BIT)
        {
        fusionneBitTab(bt,((Feuille *)po->x->fils[po->num])->sequences);
#if  DEBUG_BT
        printBitTab(((Feuille *)po->x->fils[po->num])->sequences);
#endif
        }
    else
        {
        fusionneBitTab(bt,po->x->fils[po->num]->sequences);
#if  DEBUG_BT
        printBitTab(po->x->fils[po->num]->sequences);
#endif
        }
    pos--;
    po--;
    }

#if  DEBUG_BT
printf("Somme BT : \n");
printBitTab(*bt);
printf(" -> %d values\n", nbSequenceInBitTab(*bt));
#endif

return nbSequenceInBitTab(*bt);
}

/******************************************************************************/
/* chercheMot                                                                 */
/******************************************************************************/
/* Explore les modeles recursivement.                                         */
/******************************************************************************/
void  chercheMot (Arbre a, P_Criteres  cr, signed char *mot,
        unsigned int *nb_occ, unsigned int *nbocc_ex, FILE *f)
{
int         symbol,
            trans;
LongSeq     lmaxbr,
            pos, precdummy,
            longmod = 0,
            curbloc = 0;
P_occ       tmpocc;
int         nbocc = 0;
char        carseq;

static P_mod        model = NULL;
static P_PileOcc    pocc = NULL; 
static P_occ        next = NULL;
static Bit_Tab      *colors_model;
static char         flag=0;


/* PREPARATIFS...                                                             */
if(flag == 0)
    {
/* Allocation du modele utilise lors de la recherche                          */
    model   = allocModel();

/* Allocation du tableau de bits courant                                      */
    colors_model = AllocBitTab();
    ReinitBitTab(&colors_model);

/* Initialisation de l'occurrence courante                                    */
    next = (P_occ) calloc (1,sizeof(Occ));
    if (next == NULL)
        fatalError("doSpell: cannot allocate 'next'\n");

/* Allocation des piles d'occurrences                                         */
    pocc     = creePileOcc();
    flag     = 1;
    }
else
    {
    model->name[0] = '\0';
    model->lon     = 0;

    videPile(pocc);
    }

initOcc(next);

/* Ajout de l'occurrence nulle dans la pile d'occurrence                      */
ajouteInitOcc2Pile(pocc, a.arbre);



/* CONDITION D'EXTENSION                                                      */
for(; *mot!=-1; mot++) 
    {
    symbol  = (int) *mot;
#if DEBUG_BASE
    printf("LONGMOD %d: j'etends %s vers %s%c\n",longmod,model->name,
                model->name,lettres[symbol]);
#endif

    pos     = pocc->pos-1;
    tmpocc = pocc->occ+pos;
    precdummy = getPrecDummy(pocc);
    ajouteDummy(pocc);

    nbocc   = 0;

#if DEBUG_BASE
    printf("J'ENTRE (l=%d  symbol=%d model=%s)\n",longmod,symbol,
            model->name); 
#endif
    while ((pos != precdummy) && (tmpocc->x != NULL))
        {
        lmaxbr = ((tmpocc->x->fils)[tmpocc->num]->debut & LEAF_BIT)?
                getValue(Liste_positions_fin,
                ((Feuille *)tmpocc->x->fils[tmpocc->num])->fin_deb) 
                - (tmpocc->x->fils[tmpocc->num]->debut & LEAF_BIT_INV) :
                tmpocc->x->fils[tmpocc->num]->fin
                - tmpocc->x->fils[tmpocc->num]->debut;

#if DEBUG_BASE
        if(longmod!=0)
            {
            printf("Je traite l'occ:%p  num %d  lon %d  saut %d  codesaut %d (longmod=%d)\n",
                    tmpocc->x,tmpocc->num,tmpocc->lon,tmpocc->saut,
                    tmpocc->codesaut, longmod);
            afficheOcc(stdout, tmpocc, longmod,0);
            printf("...et je trouve:\n");
            }
#endif

/* on est au milieu d'une branche - une transition possible */
        if (tmpocc->lon != lmaxbr)
            {
            carseq = a.text[tmpocc->x->fils[tmpocc->num]->sequence_number]
                [ (tmpocc->x->fils[tmpocc->num]->debut & LEAF_BIT_INV)
                + tmpocc->lon];

            if ( (carseq != FINAL)
                && (avanceBranche(next, tmpocc, symbol, carseq2num[(int)carseq],
                        0, cr, curbloc, cr->multiblocs) ) )
                {
                ajouteOcc2Pile(pocc, next->x, next->num, next->lon,
                    next->xerr,next->blocerr, next->saut, next->codesaut);
#if DEBUG_BASE
                printf("occ:%p  num %d  lon %d  saut %d  codesaut %d (longmod=%d)\n",
                        next->x,next->num,next->lon,next->saut,
                        next->codesaut, longmod);
                afficheOcc(stdout, next, longmod+1,0);
#endif

                nbocc++;
                }
            }
/* sinon on est a un noeud, plusieurs trans sont eventuellement possibles */
        else
            {
            for (trans = 0; trans != nbSymbSeq; trans++)
                {
                tmpocc=pocc->occ+pos;
                if (avanceBranche(next, tmpocc, symbol, trans, 1, cr,
                    curbloc, cr->multiblocs))
                    {
                    ajouteOcc2Pile(pocc, next->x, next->num, next->lon,
                        next->xerr, next->blocerr, next->saut,
                        next->codesaut);

#if DEBUG_BASE
                    printf("occ:%p  num %d  lon %d  saut %d  codesaut %d (longmod=%d)\n",
                            next->x,next->num,next->lon,next->saut,
                            next->codesaut, longmod);
                    afficheOcc(stdout, next, longmod+1, 0);
#endif
                    nbocc++;
                    }
                }
            }         

/* Si on n'a plus d'occurrences dans la pile                                  */
        if(pos == 0)
            {
#if DEBUG_BASE
            printf("break avec %d occ\n",nbocc);
#endif
            break;
            }

        pos--;
        tmpocc=pocc->occ+pos;

#if DEBUG_PILE
        printf("pos pile %d (adresse %p), len mod %d, nbocc %d\n",pos,
                tmpocc,longmod,nbocc);
        printf("x %p\n",tmpocc->x);
#endif
        }

#if DEBUG_BASE
    printf("J'ai trouve %d occ\n",nbocc);
    afficheOldOcc(pocc, longmod+1);
#endif


    if (nbocc == 0)
        {
        if(nb_occ != NULL)
            *nb_occ = 0;
        if(nbocc_ex!=NULL)
            *nbocc_ex = 0;
        return;
        }

/***************/
/* CAS DU SAUT */
/***************/
    if ( cr->multiblocs && ( *(mot+1)==numSAUT ))
        {
        if  ( gestionSaut(pocc, cr, curbloc, a.text) == 0 ) 
            {
            if(nb_occ != NULL)
                *nb_occ = 0;
            if(nbocc_ex!=NULL)
                *nbocc_ex = 0;
            return;
            }
            
        changeModel(model, symbol);
        changeModel(model, numSAUT);
                
#if DEBUG_SAUT
        afficheOldOcc(pocc,longmod+1);
#endif
        curbloc++;
        mot++;
        }
    else
        { 
#if DEBUG_BASE
        printf("nbocc = %d\n",nbocc);
#endif
    
        changeModel(model, symbol);
        }
    longmod++;
    }

*nb_occ  = sommeBTOcc(pocc, &colors_model);
if(nbocc_ex == NULL && f == NULL)
    return;

keepModel(pocc, longmod, cr, nbocc_ex, f);

return;
}





/******************************************************************************/
/* valideCriteres                                                             */
/******************************************************************************/
Flag    valideCriteres(P_Criteres cr)
{
/* Verification de la coherence des criteres                                  */
if(!verifCriteres(*cr))
    return FAUX;

initEntiers();

if(cr->bloc > 1)
    {
    initTabSauts(cr);
    cr->multiblocs = VRAI;
    }
else
    cr->multiblocs = FAUX;
 
return VRAI;
}

/******************************************************************************/
/* termineRecherche                                                           */
/******************************************************************************/
/* void    termineRecherche(void) */
/* { */
/* free(next); */
/* liberePileOcc(pocc); */
/* } */


/******************************************************************************/
/* chargeSequence                                                             */
/* Lit les sequences dans un fichier FASTA, rajoute un caractere FINAL.       */
/******************************************************************************/
Flag    chargeSequence(Arbre *a, char *fic)
{
FastaSequence   **seq;
Flag            readok;
int             taille, siztxt, i;
FILE            *fasta;

/* Allocations                                                                */
seq     = (FastaSequence **) malloc(GRAINSEQ * sizeof(FastaSequence *));
a->text = (signed char **)   malloc(GRAINSEQ * sizeof(signed char *));
if(!seq || !a->text)
    fatalError("charSequence: cannot allocate 'seq/text'\n");
siztxt  = GRAINSEQ;

/* Ouverture du fichier contenant les sequences                               */
fasta      = fopen (fic,"r"); 
if(fasta == NULL)
    {
    fprintf(stderr,"charSequence: cannot open fasta file '%s'\n",fic);
    return FAUX;
    }

readok  = 1;
a->nbtxt= 0;

/* Stockage des sequences en memoire                                          */
do
    {
    if(a->nbtxt == siztxt)
        {
        siztxt  *= 2;
        seq     = (FastaSequence **)
                    realloc(seq,siztxt * sizeof(FastaSequence *));
        a->text = (signed char **)
                    realloc(a->text, siztxt * sizeof(signed char *));
        if(!seq || !a->text)
            fatalError("chargeSequence: cannot reallocate 'seq/text'\n");
        }

    seq[a->nbtxt] = NewFastaSequence();
    readok      = ReadFastaSequence(fasta, seq[a->nbtxt]);
    if (readok)
        {
        taille            = seq[a->nbtxt]->length+1;
        a->text[a->nbtxt] = (signed char *) malloc ((taille+2) * sizeof(signed char)); 
        if (a->text[a->nbtxt] == NULL)
            fatalError("chargeSequence: cannot allocate 'text'\n");

        strcpy((char *) a->text[a->nbtxt],seq[a->nbtxt]->seq);
        a->text[a->nbtxt][taille-1] = FINAL;
        a->text[a->nbtxt][taille]   = '\0';
        (a->nbtxt)++;
        }
    }
while (readok);

fclose(fasta); 

/* Liberation de la structure Fasta                                           */
for(i=0;i != a->nbtxt;i++)
    FreeFastaSequence(seq[i]);
free(seq);

return VRAI;
}

/******************************************************************************/
/* creeArbreSuffixe                                                           */
/******************************************************************************/
Flag creeArbreSuffixe(Arbre *a, int maxlongmod, char *alphaseq)
{
int i;
Noeud   *root_pere;


if(alphaseq==NULL)
    alphaseq    = chargeAlphaSeq((Symbole **) a->text, a->nbtxt, NULL);


/* Construction de l'arbre compact generalise                                 */
/* fprintf(stderr, "** Constructing suffix tree **\n"); */
/* barre(a->nbtxt); */
Init_All((unsigned char *) alphaseq, 0, a->nbtxt);
a->arbre = Construction_Arbre((unsigned char *)a->text[0], maxlongmod);
/* barre(0); */

for (i = 1; i != a->nbtxt; i++)
    {
    a->arbre    = AjouteSequence(a->arbre,(unsigned char *)a->text[i],
        maxlongmod);
/*     barre(0); */
    }

UpdateBit_TabForAllTree(a->arbre);

/* Creation du faux pere du pere (pour faciliter la recursion)                */
root_pere = Alloc_Noeud();

root_pere->fils[Translation_Table[FINAL]]        = a->arbre;
root_pere->sequence_number                       = 0;
a->arbre->debut                                  = 0;
a->arbre->fin                                    = 1;
a->arbre->sequence_number                        = 0;

a->arbre    = root_pere;

return VRAI;
}

/******************************************************************************/
/* creeArbreSuffixeFromFile                                                   */
/******************************************************************************/
Flag creeArbreSuffixeFromFile(Arbre *a, char *fic, int maxlongmod,
        char *alphaseq)
{
if( !chargeSequence(a, fic) ||  a->nbtxt < 1 )
    {
    fprintf(stderr,"Not enough sequences (<1)\n");
    return FAUX;
    }

creeArbreSuffixe(a, maxlongmod, alphaseq);

return VRAI;
}

/******************************************************************************/
/* creeArbreSuffixeFromArray                                                  */
/******************************************************************************/
Flag creeArbreSuffixeFromArray(Arbre *a, char **seq, int nbseq, int maxlongmod,
        char *alphaseq)
{
if( nbseq < 1 )
    {
    fprintf(stderr,"Not enough sequences (<1)\n");
    return FAUX;
    }
a->nbtxt    = nbseq;
a->text     = (signed char **) seq;

creeArbreSuffixe(a, maxlongmod, alphaseq);

return VRAI;
}

/******************************************************************************/
/* libereArbreSuffixeFromFile                                                 */
/******************************************************************************/
void    libereArbreSuffixeFromFile(Arbre a)
{
int i;

libereArbreSuffixeFromArray(a);

for(i=0; i!=a.nbtxt; i++)
   free(a.text[i]);
free(a.text);
}


/******************************************************************************/
/* libereArbreSuffixeFromArray                                                */
/******************************************************************************/
void    libereArbreSuffixeFromArray(Arbre a)
{
Free_Arbre(a.arbre);
Free_All_Liste_Cell();
Free_ListePositions(Liste_positions_fin);
}





/******************************************************************************/
/******************************************************************************/
/********************************** MAIN **************************************/
/******************************************************************************/
/******************************************************************************/
/*
int main(int argc, char **argv)
{
Criteres    criteres, cr;
FILE        *f=NULL;
unsigned int nbocc, nboccex;
Arbre       a,b;

initCriteres(&criteres);
initCriteres(&cr);

setBloc(&criteres,2);
setLongueurBloc(&criteres,0,8);
setLongueurBloc(&criteres,1,4);
setErreurGlobal(&criteres,1);
setErreurBloc(&criteres,0,0);
setErreurBloc(&criteres,1,1);
setSaut(&criteres,0,8,10);
if(valideCriteres(&criteres) == FAUX)
    return 1;

setBloc(&cr,3);
setLongueurBloc(&cr,0,3);
setLongueurBloc(&cr,1,3);
setLongueurBloc(&cr,2,3);
setErreurGlobal(&cr,1);
setErreurBloc(&cr,0,1);
setErreurBloc(&cr,1,1);
setErreurBloc(&cr,2,1);
setSaut(&cr,0,2,5);
setSaut(&cr,1,0,3);
if(valideCriteres(&cr) == FAUX)
    return 1;

f = fopen("out","w");

printf("Je construit l'arbre de taille %d\n",maxLongMod(cr));
creeArbreSuffixeFromFile(&b, "seq", maxLongMod(cr));
chercheMot(b,&cr, "AA_CA_AA",&nbocc, &nboccex, f);
printf("%d %d!\n",nbocc, nboccex);
libereArbreSuffixe(b);

creeArbreSuffixe(&a, "ficseq",maxLongMod(criteres));
chercheMot(a,&criteres, "GCGACATA_GATG",&nbocc, &nboccex, f);
printf("%d %d!\n",nbocc, nboccex);
libereArbreSuffixe(a);

return(0); 
}
*/
