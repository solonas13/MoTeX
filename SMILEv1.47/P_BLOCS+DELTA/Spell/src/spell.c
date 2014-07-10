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

#include <spell.h>


/******************************************************************************/
/*                          PROTOTYPES PRIVES                                 */
/******************************************************************************/
/* Gestion des modeles acceptes                                               */
void    keepModel(P_mod , P_PileOcc, NbSeq, NbSeq, LongSeq, LongSeq *,
        int *cursaut, P_Criteres cr);

/* essaie d'avancer d'une lettre dans un arc, et renvoie le noeud image */
Flag avanceBranche(P_occ, P_occ, int, int, Flag, P_Criteres, LongSeq);

/* Lancement du saut */
NbSeq gestionSaut(P_mod model, P_PileOcc pocc, P_PileOcc poccnew, P_Criteres,
        int *nbocc_saut, int maxinter, LongSeq curbloc);

/* explore les modeles */
Flag  spellModels ( P_PileOcc   pocc,
                    P_PileOcc   poccnew,        P_PileOcc   poccsaut,
                    LongSeq     longmod,        LongSeq     longcurbloc,
                    LongSeq     curbloc,
                    P_mod       model,          P_occ       next,
                    Bit_Tab     **colors_model, NbSeq       nbseq,
                    NbSeq       tmp_quorum,     P_Criteres  cr,
                    int         **nbocc_saut,   Bit_Tab     ***bt_intermediaire,
                    int         deb_occ,        int         fin_occ,
                    int         *cursaut,       LongSeq     *longbloc,
                    LongSeq     *posdebbloc);

/* Calcule le BT union de tous les BT occurrences                             */
NbSeq   sommeBTOcc(P_PileOcc, Bit_Tab **);

/* Compute CPU time                                                           */
static  float PrintCpuTime(char);



/******************************************************************************/
/* VARIABLES GLOBALES                                                         */
/******************************************************************************/
int     *nbmod;
LongSeq *maxlongmod, **maxlongbloc;
signed char    ** text;
FILE    ** f;

/* EXTERNES from alphabet.c                                                   */
extern  int     nbSymbMod;
extern  int     nbSymbSeq;
extern  char    *nummod2str[127];
extern  int     carseq2num[127];
extern  int     comp[127];
extern  Flag    TabSymb[127][127];
extern  int     numJOKER;
extern  int     numSAUT;




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
void keepModel(P_mod model, P_PileOcc pocc, NbSeq nbseq, NbSeq quorum,
        LongSeq l, LongSeq *longbloc, int *cursaut, P_Criteres cr)
{
int     numfile, i, j;
FILE    *g;
LongSeq *lb, *mb;


numfile = delta2File(cursaut, cr);
nbmod[numfile]++;
g = f[numfile];

#if  DEBUG_BASE
        printf("MODELE %s valide!\n",model->name);
#endif

j   = model->lon;
for(i=0; i!=j; i++) 
        fprintf(g,"%s", nummod2str[model->name[i]]);
fprintf(g," ");
for(i=0; i!=j; i++) 
    {
    if(model->name[i]==numJOKER)
        fprintf(g,"%c",JOKERinterne);
    else if(model->name[i]==numSAUT)
        fprintf(g,"%c",SAUTinterne);
    else
        fprintf(g,"%c", model->name[i]+SHIFTALPHA);
    }
fprintf(g," %d", quorum);

if(l > maxlongmod[numfile])
        maxlongmod[numfile]  = l;

for(i=0, lb=longbloc, mb=maxlongbloc[numfile]; i!=cr->bloc; i++, lb++, mb++)
    if(*lb > *mb)
        *mb  = *lb;
    
#if OCC
    #if AFF_OCC
        fprintf(g,"\n");
    #else
        fprintf(g,"\t");
    #endif
    afficheLastOcc(g, pocc, l, cr);
#else
    fprintf(g,"\n");
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
Flag avanceBranche( P_occ   next,       P_occ       tmp,        int symbol,
                    int     trans,      Flag        flag_noeud,
                    P_Criteres  cr,     LongSeq curbloc)
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
        
        next->blocerr = tmp->blocerr+1; /* si maxerr local atteint */
        if (next->blocerr == cr->maxerrblocs[curbloc]+1)
            return 0;
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
        
        next->blocerr = tmp->blocerr+1; /* si maxerr local atteint */
        if (next->blocerr == cr->maxerrblocs[curbloc]+1)
            return 0;
        }
    
    next->num = trans;
    next->lon = 1;
    }


next->saut= tmp->saut;
next->codesaut= tmp->codesaut;

return(1);
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
if (pos < 0)
    fatalError("spell.c: sommeBTOcc: wrong position in stack!\n");
po  = p->occ+pos;

precdummy   = getPrecDummy(p);

while ((pos != precdummy) && (po->x != NULL))
    {
#if DEBUG_BT
    printf("Fusion avec :    ");
#endif
      
    if (po->x->fils[po->num]->debut & LEAF_BIT)
        {
        fusionneBitTab(bt,((Feuille *)po->x->fils[po->num])->sequences);
#if DEBUG_BT
        printBitTab(((Feuille *)po->x->fils[po->num])->sequences);
#endif
        }
    else
        {
        fusionneBitTab(bt,po->x->fils[po->num]->sequences);
#if DEBUG_BT
        printBitTab(po->x->fils[po->num]->sequences);
#endif
        }
    pos--;
    po--;
    }

#if DEBUG_BT
printf("Somme BT : \n");
printBitTab(*bt);
printf(" -> %d values\n", nbSequenceInBitTab(*bt));
#endif

return nbSequenceInBitTab(*bt);
}

/******************************************************************************/
/* sommeBTOccPartielle                                                        */
/******************************************************************************/
/* Fait l'union des BR d'une partie de la pile d'occurrences                  */
/******************************************************************************/
void sommeBTOccPartielle(   P_PileOcc p, Bit_Tab **bt, int nb_tranches,
                            int *nboccsaut)
{
LongSeq pos = p->pos-1;
P_occ   po  = p->occ+pos;
int     i, *ptr_saut;
Bit_Tab **curbt;
  
if (pos < 0)
    fatalError("spell.c: sommeBTOccPartielle: wrong position in stack!\n");

curbt       = bt+nb_tranches-1;
ptr_saut    = nboccsaut+nb_tranches-1;
while ( nb_tranches != 0 )
    {
    ReinitBitTab(curbt);
    
    for ( i=*ptr_saut; i!=0; i--)
        {
        if (po->x->fils[po->num]->debut & LEAF_BIT)
            fusionneBitTab(curbt,((Feuille *)po->x->fils[po->num])->sequences);
        else
            fusionneBitTab(curbt,po->x->fils[po->num]->sequences);

        po--;
        }
    nb_tranches--;
    ptr_saut--;
    curbt--;
    }

#if DEBUG_BT
printf("Somme BT : \n");
printBitTab(*curbt);
printf(" -> %d values\n", nbSequenceInBitTab(*curbt));
#endif
}

/******************************************************************************/
/* sauteSymbole                                                               */
/******************************************************************************/
int sauteSymbole(Occ curocc, P_mod model, P_PileOcc pocc, Fourchette *range,
        P_Criteres cr, LongSeq longsaut, int longmod, LongSeq curbloc)
{
LongSeq lmaxbr;
Noeud   *tmpnoeud;
int     res,
        trans;
char    carseq;


tmpnoeud = curocc.x->fils[curocc.num];

 if (tmpnoeud->debut & LEAF_BIT)
   lmaxbr  = getValue(Liste_positions_fin,((Feuille *)tmpnoeud)->fin_deb)
           - (((Feuille *)tmpnoeud)->debut & LEAF_BIT_INV); 
 else
   lmaxbr = tmpnoeud->fin - tmpnoeud->debut;

#if DEBUG_SAUT
printf("SauteSymbole: je gere le saut pour %d, noeud %p, etat: %d/%d branche %d\n",
        longsaut,curocc.x,curocc.lon,lmaxbr,curocc.num);
#endif


if (curocc.lon != lmaxbr)  /* on est au milieu d'une branche */
    {
    curocc.lon++;

    carseq = text[tmpnoeud->sequence_number]
        [(tmpnoeud->debut & LEAF_BIT_INV)+curocc.lon-1];

    if(carseq==FINAL)  /* si on rencontre un $ c'est fini */
        return 0;

    if(longsaut<=(*range).max)
        {
        ajouteOcc2Pile(pocc, curocc.x, curocc.num, curocc.lon, curocc.xerr,
                curocc.blocerr, curocc.saut+1,
                addSaut2Code(curocc.blocerr, longsaut, curbloc, cr));
#if DEBUG_SAUT
        printf("Occ (br) num %d  lon %d  saut %d\n",curocc.num,curocc.lon,
                longsaut);
        afficheOcc(stdout, &curocc, longmod+1, 0);
#endif
        return 1;
        }

    return 0;
    }
else    /* sinon on est a un noeud, plusieurs trans sont possibles */
    {
    res = 0;
    if(longsaut<=(*range).max)
        {
        curocc.x    = tmpnoeud;
        curocc.lon  = 1;

        if ((tmpnoeud->debut & LEAF_BIT) == 0)
          for (trans = 0; trans != nbSymbSeq; trans++)
            {
            if (tmpnoeud->fils[trans] != NULL)
                {
                curocc.num  = trans;

                res++;
                ajouteOcc2Pile(pocc, curocc.x, curocc.num, curocc.lon,
                        curocc.xerr, curocc.blocerr, curocc.saut+1,
                        addSaut2Code(curocc.blocerr, longsaut, curbloc, cr));
#if DEBUG_SAUT
                printf("Occ (nd) num %d  lon %d  saut %d\n", curocc.num,
                        curocc.lon,longsaut);
                afficheOcc(stdout, &curocc, longmod, 0);
#endif
                }
            }
        }
    return res;
    }
}


/******************************************************************************/
/* saute2MinSaut                                                              */
/******************************************************************************/
/* Avance les occ jusqu'a saut min. Les occurrences trouvees sont ajoutees
 * a la pile, avec cela de specifique que le champ 'blocerr' sert temporaire-
 * ment de stockage pour 'codesaut' */
/******************************************************************************/
int saute2MinSaut(Occ curocc, P_mod model, P_PileOcc pocc, Fourchette *range,
        P_Criteres cr, LongSeq longsaut, int longmod, LongSeq curbloc)
{
LongSeq lmaxbr;
Noeud * tmpnoeud, *newtmpnoeud;
int     res = 0, newlongsaut,
        trans;
char    carseq;


tmpnoeud = curocc.x->fils[curocc.num];

if (tmpnoeud->debut & LEAF_BIT)
    lmaxbr  = getValue(Liste_positions_fin,((Feuille *)tmpnoeud)->fin_deb)
            - (((Feuille *)tmpnoeud)->debut & LEAF_BIT_INV); 
else
    lmaxbr = tmpnoeud->fin - tmpnoeud->debut;

#if DEBUG_SAUT
printf("Saute2MinSaut: j'ai gere le saut pour %d, noeud %p, etat: %d/%d branche %d\n",
        longsaut,curocc.x,curocc.lon,lmaxbr,curocc.num);
#endif

if (curocc.lon != lmaxbr)  /* on est au milieu d'une branche */
    {
    if ( lmaxbr-curocc.lon <= (*range).min-longsaut )
        /* si on peut aller au bout de la branche */
        {   
        longsaut    += lmaxbr-curocc.lon;
        curocc.lon  = lmaxbr;

#if DEBUG_SAUT
     printf("Saute2MinSaut: milieuBr, fast, je vais au bout %d/%d br %d et lgsaut %d\n",
             curocc.lon,lmaxbr,curocc.num,longsaut);
#endif

        carseq = text[tmpnoeud->sequence_number]
            [(tmpnoeud->debut & LEAF_BIT_INV)+lmaxbr-1];

        if(carseq != FINAL)
         /* si on rencontre un $ c'est fini */
            res += saute2MinSaut(curocc, model, pocc, range, cr, longsaut,
                                longmod, curbloc);
#if DEBUG_SAUT
        else
            printf("Une qui tombe avant minima\n");
#endif
        }
    else
        {
        curocc.lon  += (*range).min-longsaut;
        longsaut    = (*range).min;

#if DEBUG_SAUT
        printf("Saute2MinSaut: milieuBr, minsaut ds Br, je m'arrete a %d/%d num %d et lgsaut %d\n",
                curocc.lon,lmaxbr,curocc.num,longsaut);
#endif

        res++;
        ajouteOcc2Pile(pocc, curocc.x, curocc.num, curocc.lon, curocc.xerr,
                curocc.codesaut, curocc.saut+longsaut,
                addSaut2Code(curocc.codesaut, longsaut, curbloc, cr) );
#if DEBUG_SAUT
        printf("Occ (br) num %d  lon %d  saut %d\n", curocc.num,curocc.lon,
                longsaut);
        afficheOcc(stdout, &curocc, longmod, 0);
#endif
        }
    }
else    /* sinon on est a un noeud, plusieurs trans sont possibles */
    {
    curocc.x    = tmpnoeud;

    if ((tmpnoeud->debut & LEAF_BIT) == 0)
      {
      for (trans = 0; trans != nbSymbSeq; trans++)
        {
        curocc.num  = trans;
        newlongsaut = longsaut;

        if (tmpnoeud->fils[trans] != NULL)
            {
            newtmpnoeud = tmpnoeud->fils[trans];

            if (newtmpnoeud->debut & LEAF_BIT)
                lmaxbr  = getValue(Liste_positions_fin,
                        ((Feuille *)newtmpnoeud)->fin_deb) 
                        - (newtmpnoeud->debut & LEAF_BIT_INV); 
            else
                lmaxbr = newtmpnoeud->fin - newtmpnoeud->debut;
            
            if ( lmaxbr <= (*range).min-longsaut )
                {
                newlongsaut += lmaxbr;
                curocc.lon  = lmaxbr;

#if DEBUG_SAUT
                printf("Saute2MinSaut: noeud, fast, %d/%d, br %d, lgsaut %d\n",
                        curocc.lon,lmaxbr,curocc.num,newlongsaut);
#endif

                carseq = text[newtmpnoeud->sequence_number]
                    [(newtmpnoeud->debut & LEAF_BIT_INV)+lmaxbr-1];
                if(carseq != FINAL)  /* si on rencontre un $ c'est fini */
                    {
#if DEBUG_SAUT
                    printf("Saute2MinSaut: finBr=$, c'est fini\n");
#endif
                    res += saute2MinSaut(curocc, model, pocc, range, cr,
                            newlongsaut, longmod, curbloc);
                    }
                }
            else
                {
                curocc.lon  = (*range).min-newlongsaut;
                newlongsaut = (*range).min;

#if DEBUG_SAUT
                printf("Saute2MinSaut2: noeud %p, minsaut ds Br, %d/%d, br %d, lgsaut %d\n",
                        curocc.x, curocc.lon,lmaxbr,curocc.num,newlongsaut);
#endif

                res++;
                ajouteOcc2Pile(pocc, curocc.x, curocc.num, curocc.lon,
                        curocc.xerr, curocc.codesaut, curocc.saut+newlongsaut,
                        addSaut2Code(curocc.codesaut, newlongsaut, curbloc,
                        cr));
#if DEBUG_SAUT
                printf("Occ (nd) num %d  lon %d  saut %d\n", curocc.num,
                        curocc.lon,longsaut);
                afficheOcc(stdout, &curocc, longmod, 0);
#endif
                }
            }
        }
      }
    }
return res;
}

/******************************************************************************/
/* initBlocerrPile                                                            */
/******************************************************************************/
/* Sert a reinitialiser le champ 'blocerr' des occurrences trouvees pour
 * le saut, et qui etaient utilisess pour 'codesaut'. */
/******************************************************************************/
void    initBlocerrPile(P_PileOcc p)
{
LongSeq pos     = p->pos-1,
        precdummy;
P_occ   tmpocc  = p->occ+pos;

precdummy   = getPrecDummy(p);

while ((pos != precdummy) && (tmpocc->x != NULL) )
    {  
    tmpocc->blocerr = 0;
    pos--;
    tmpocc--;
    }
}

/******************************************************************************/
/* gestionSaut                                                                */
/******************************************************************************/
/* Gere le saut en classant les occurrences obtenues par longueur de saut     */
/* Renvoie le nombre de tranches de saut                                      */
/******************************************************************************/
NbSeq gestionSaut( P_mod      model, P_PileOcc    pocc, P_PileOcc     poccnew,
                   P_Criteres cr,   int    *nbocc_saut,
                   int        maxinter, LongSeq curbloc)
{
LongSeq     pos, precdummy;
P_occ       tmpocc;
Fourchette  range = cr->saut[curbloc];

int     i, oldpos, nb_tranches = 0, tmpres = 0, newoldpos, *ptr_saut;

#if DEBUG_SAUT
printf("Et je saute (modele %s)...\n",model->name);
#endif


ptr_saut    = nbocc_saut;
/* Initialisation du tableau de stockage des nb d'occ par tranche de saut     */
for (i=0; i != maxinter; i++, ptr_saut++)
    *ptr_saut   = 0;
ptr_saut=nbocc_saut;

/* Releve de la fin de pile d'occurrences                                     */
ajouteDummy(pocc);
oldpos  = pocc->pos-1;

#if DEBUG_SAUT
printf("Le dummy de la pile des new occ est en %d\n",oldpos);
#endif

/* Extension au seuil minimal des occurrences courantes                       */
if (range.min != 0)
    {
/* Releve des positions dans la pile des nouvelles occurrences                */
    pos     = poccnew->pos-1;
    precdummy = getPrecDummy(poccnew);
    tmpocc  = poccnew->occ+pos;
    while ((pos != precdummy) && (tmpocc->x != NULL) )
        {
        tmpres += saute2MinSaut(*tmpocc, model, pocc, &range, cr, 0,
                        model->lon+1, curbloc);

        pos--;
        tmpocc  = poccnew->occ+pos; /** ATTENTION: c'est obligatoire
                                      cause realloc!! */
        }
    
    if ( tmpres == 0 )
        {
        depileRec(pocc);
        return 0;
        }
    }
else
    {
    tmpres  = copieLastOcc(pocc, poccnew);

#if DEBUG_SAUT
    for(i=tmpres, tmpocc=pocc->occ+(pocc->pos-1); i!=0; i--, tmpocc--)
        printf("Occ num %d  lon %d  saut %d\n", tmpocc->num,tmpocc->lon,
                tmpocc->saut);
#endif

    if (tmpres == 0)
        {
        depileRec(pocc);
        return 0;
        }

/* Mise a jour des codesauts!                                                 */
    if(curbloc != 0 || range.max!=0)
        {
        pos     = pocc->pos-1;
        precdummy = getPrecDummy(pocc);
        tmpocc  = pocc->occ+pos;
        while ((pos != precdummy) && (tmpocc->x != NULL))
            {
            if(range.max!=0)
                tmpocc->blocerr = tmpocc->codesaut;
            if(curbloc !=0)
                tmpocc->codesaut= addSaut2Code(tmpocc->codesaut, 0, curbloc, cr);
            pos--;
            tmpocc--;
            }
        }
    }

*ptr_saut   = tmpres;
ptr_saut++;
nb_tranches++;

#if DEBUG_SAUT
printf("Etape de saut %d, j'en trouve %d\n",range.min,tmpres);
#endif


if ( range.min == range.max )
    {
    if(range.max!=0)
        initBlocerrPile(pocc);
    return 1;
    }

/* Extension du seuil minimal au seuil maximal                                */
for ( i = range.min; i != range.max; i++ )
    {
    pos     = newoldpos = pocc->pos-1;
    tmpocc  = pocc->occ+pos;
    tmpres  = 0;

#if DEBUG_SAUT
    printf("Je m'arreterai en %d\n",oldpos+1);
#endif

    while (pos != oldpos)
        {
#if DEBUG_SAUT
        printf("saut: ext occ No %d\n",pos);
#endif

        tmpres  += sauteSymbole(*tmpocc, model, pocc, &range, cr, i+1,
                model->lon+1, curbloc);
        pos--;
        tmpocc  = pocc->occ+pos; /** ATTENTION: c'est obligatoire cause
                                   realloc!! */
        }
    oldpos  = newoldpos;

    *ptr_saut   = tmpres;

#if DEBUG_SAUT
    printf("Etape de saut %d, j'en trouve %d\n",i+1,tmpres);
#endif

    if (tmpres == 0)
        {
#if DEBUG_SAUT
        break;
#endif
        initBlocerrPile(pocc);
        return nb_tranches;
        }

    ptr_saut++;
    nb_tranches++;
    }

#if DEBUG_SAUT
printf("PILE APRES SAUT:\n");
pos = pocc->pos-1;
tmpocc = pocc->occ+pos;

while(tmpocc && tmpocc->lon != -1)
    {
    printf("*** Occ num %d  lon %d  saut %d  codesaut %d\n",
            tmpocc->num,tmpocc->lon, tmpocc->saut, tmpocc->codesaut);
    afficheOcc(stdout, tmpocc, model->lon, 0);
    tmpocc--;
    }

tmpres=0;
printf("\nTABLEAU NB_TRANCHES EN SORTIE DE SAUT: (nbtranches=%d)\n",nb_tranches);

for(i=0; i<nb_tranches;i++)
    {
    printf("%d\t",nbocc_saut[i]);
    tmpres+=nbocc_saut[i];
    }
printf("\nAu total %d occ dans la pile de saut\n",tmpres);
#endif

initBlocerrPile(pocc);

return nb_tranches;
}




/******************************************************************************/
/* spellModels                                                                */
/******************************************************************************/
/* Explore les modeles recursivement.                                         */
/******************************************************************************/
Flag  spellModels ( P_PileOcc   pocc,
                    P_PileOcc   poccnew,        P_PileOcc   poccsaut,
                    LongSeq     longmod,        LongSeq     longcurbloc,
                    LongSeq     curbloc,
                    P_mod       model,          P_occ       next,
                    Bit_Tab     **colors_model, NbSeq       nbseq,
                    NbSeq       tmp_quorum,     P_Criteres  cr,
                    int         **nbocc_saut,   Bit_Tab     ***bt_intermediaire,
                    int         deb_occ,        int         fin_occ,
                    int         *cursaut,       LongSeq     *longbloc,
                    LongSeq     *posdebbloc)
{
Flag        zarb_back = 0,
            zarb_ext  = 0;
int         symbol,
            trans;
LongSeq     lmaxbr,
            pos, precdummy,
            minsaut, maxsaut, delta,
            palbloc;
long int    maxseq;
P_occ       tmpocc;
int         tmpint, *tmpintptr,
            nbnewmod = 0,
            nbocc, i, j, k, nb_tranches, new_deb_occ, new_fin_occ;
NbSeq       tmp_quorum2;
char        carseq;


#if DEBUG_SAUT
int     tmpint2;
#endif

if(longmod==3)
    barre(0);

/* CONDITION D'EXTENSION                                                      */
if ( ( (cr->longueur.max == 0) || (longmod < cr->longueur.max))
    && ( 
        ( (cr->longbloc[curbloc].max == 0)
        || (longcurbloc < cr->longbloc[curbloc].max) ) )
    &&
       ( (cr->flag_palindrom == FAUX)
       || cr->palindrom[curbloc] == -1
       || longcurbloc != longbloc[(int)(cr->palindrom[curbloc])]
     )
    )
    {
/* Boucle sur les symboles de l'alphabet pourl'extension du modele ************/
    for (symbol = 0; symbol != nbSymbMod; symbol++)
        {
#if DEBUG_BASE
        printf("LONGMOD %d: j'etends %s vers %s%c\n",longmod,model->name,
            model->name,lettres[symbol]);
#endif

/* Pas de JOKER en premiere position                                          */
        if(longmod == 0 && symbol == numJOKER)
            continue;


/* Gestion de la composition des modeles **************************************/
        if (cr->flag_compo == VRAI || cr->flag_compobloc[curbloc] == VRAI )
            {
            if ( (cr->compobloc[curbloc][symbol] == 0)
                || (cr->compo[symbol] == 0) )
                continue;
            }

/* Gestion des palindromes                                                    */
        if (cr->flag_palindrom)
            {
            if (longcurbloc == 0)
                posdebbloc[curbloc] = model->lon;

            if (cr->palindrom[curbloc] != -1)
                {
                palbloc = cr->palindrom[curbloc];
    
                if (symbol!=
                   comp[model->name[posdebbloc[palbloc]+longbloc[palbloc]-1-longcurbloc]])
                    continue;
                }
            }

/* Init variables de pile d'occs                                              */
        if ( fin_occ == 0 )
            pos     = pocc->pos-1;
        else
            pos     = fin_occ;
        precdummy = getPrecDummy(pocc);
        tmpocc = pocc->occ+pos;
        videPile(poccnew);

        maxseq  = 0;
        nbocc   = 0;
#if DEBUG_BASE
        printf("J'ENTRE (l=%d  symbol=%d model=%s quorum=%d)\n",longmod,symbol,
                model->name,tmp_quorum); 
#endif
        while ( (pos!=precdummy) && (tmpocc->x != NULL) && (fin_occ == 0  || pos != deb_occ))
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
                carseq = text[tmpocc->x->fils[tmpocc->num]->sequence_number]
                    [ (tmpocc->x->fils[tmpocc->num]->debut & LEAF_BIT_INV)
                    + tmpocc->lon];

                if ( (carseq != FINAL)
                    && (avanceBranche(next, tmpocc, symbol,
                        carseq2num[(int) carseq],0, cr, curbloc) ) )
                    {
                    ajouteOcc2Pile(poccnew, next->x, next->num, next->lon,
                        next->xerr,next->blocerr, next->saut, next->codesaut);
#if DEBUG_BASE
                    printf("occ:%p  num %d  lon %d  saut %d  codesaut %d (longmod=%d)\n",
                            next->x,next->num,next->lon,next->saut,
                            next->codesaut, longmod);
                    afficheOcc(stdout, next, longmod+1,0);
#endif

                    nbocc++;

                    if (next->x->fils[next->num]->debut & LEAF_BIT)
                        {
                        maxseq += nbSequenceInBitTab(
                            ((Feuille *)next->x->fils[next->num])->sequences);
#if DEBUG_BT
                        printf("nb seq in bt (br): %d \n",
                            nbSequenceInBitTab(((Feuille *)
                            next->x->fils[next->num])->sequences));
#endif
                        }
                    else
                        {
                        maxseq += next->x->fils[next->num]->nb_element_bt;
#if DEBUG_BT
                        printf("nb seq in bt (br): %d \n",
                                nbSequenceInBitTab(
                                next->x->fils[next->num]->sequences));
#endif
                        }
                    }
                }
/* sinon on est a un noeud, plusieurs trans sont eventuellement possibles */
            else
                {
                for (trans = 0; trans != nbSymbSeq; trans++)
                    {
                    tmpocc = pocc->occ+pos;
                    if (avanceBranche(next, tmpocc, symbol, trans, 1, cr,
                        curbloc))
                        {
                        ajouteOcc2Pile(poccnew, next->x, next->num, next->lon,
                            next->xerr, next->blocerr, next->saut,
                            next->codesaut);

#if DEBUG_BASE
                        printf("occ:%p  num %d  lon %d  saut %d  codesaut %d (longmod=%d)\n",
                                next->x,next->num,next->lon,next->saut,
                                next->codesaut, longmod);
                        afficheOcc(stdout, next, longmod+1, 0);
#endif
                        nbocc++;
                        if (next->x->fils[next->num]->debut & LEAF_BIT)
                            {
                            maxseq += nbSequenceInBitTab(((Feuille *)
                                next->x->fils[next->num])->sequences);
#if DEBUG_BT
                            printf("nb seq in bt (nd): %d \n",
                                    nbSequenceInBitTab(((Feuille *)
                                    next->x->fils[next->num])->sequences));
#endif
                            }
                        else
                            {
                            maxseq += next->x->fils[next->num]->nb_element_bt;
#if DEBUG_BT
                            printf("nb seq in bt (nd): %d \n",
                                    nbSequenceInBitTab(
                                    next->x->fils[next->num]->sequences));
#endif
                            }
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
            tmpocc = pocc->occ+pos;

#if DEBUG_PILE
            printf("pos pile %d (adresse %p), len mod %d, nbocc %d\n",pos,
                    tmpocc,longmod,nbocc);
            printf("x %p\n",tmpocc->x);
#endif
            }

#if DEBUG_BASE
printf("J'ai trouve %d occ\n",nbocc);
afficheOldOcc(poccnew, longmod+1);
#endif


        if (nbocc == 0)
            continue;

/***************/
/* CAS DU SAUT */
/***************/
        tmp_quorum2 = -1;

        if (    (curbloc != cr->bloc-1)
            &&  (longcurbloc+1 >= cr->longbloc[curbloc].min )
            &&  (maxseq >= cr->quorum)
            &&  ( (tmp_quorum2 = sommeBTOcc(poccnew, colors_model) )
                     >= cr->quorum) 
            &&  ( (nb_tranches = gestionSaut(model, pocc, poccnew, cr,
                *nbocc_saut, cr->maxinter, curbloc)) != 0 ) )
            {
            minsaut = cr->saut[curbloc].min;
            maxsaut = cr->saut[curbloc].max;
            delta   = cr->delta[curbloc];

/* Calcul des BT pour chaque tranche de saut                                  */
            sommeBTOccPartielle(pocc, *bt_intermediaire, nb_tranches,
                    *nbocc_saut);
            
/* deb et fin stockent les positions des occ interessantes dans la pile       */
            new_deb_occ = pocc->carte[pocc->pos_carte-1];
            new_fin_occ = new_deb_occ;
            tmpintptr   = *nbocc_saut;
            for ( i=0; i<2*delta+1; i++, tmpintptr++)
                new_fin_occ += *tmpintptr;

            changeModel(model, symbol);
            changeModel(model, numSAUT);
                
            ajouteDummy(poccsaut);

            tmpint  = copieLastOcc(poccsaut,poccnew);
#if DEBUG_SAUT
            printf("J'ai copie %d occ de Pnew->Psaut\n",tmpint);
            afficheOldOcc(poccnew,longmod+1);
#endif
            videPile(poccnew);
            zarb_ext = 1;
            
            if ( cr->flag_compo == VRAI || cr->flag_compobloc[curbloc] == VRAI )
                {
                cr->compo[symbol]--;
                cr->compobloc[curbloc][symbol]--;
                }

            longbloc[curbloc]   = longcurbloc+1;


            k       = 0;
            tmpint  = nb_tranches;
            for(i=minsaut+delta; (i+delta<=maxsaut)
                    && (tmpint != 0); i++, tmpint--, k++) 
                {
/* Mise a jour des positions des occurrences interessantes dans la pile       */
                if(k!=0)
                    {
                    new_deb_occ += (*nbocc_saut)[k-1];
                    new_fin_occ += (*nbocc_saut)[k+2*delta];
                    }

                cursaut[curbloc]    = i;
                
                ReinitBitTab(colors_model);

                for(j=i-delta; j!=i+delta+1; j++)
                    fusionneBitTab(colors_model,
                            (*bt_intermediaire)[j-minsaut]);
                    
                if ( nbSequenceInBitTab(*colors_model) >= cr->quorum )
                    {   
#if DEBUG_SAUT
                    printf("Modele %s : La tranche [%d-%d] est acceptee avec\n",
                            model->name,i-delta, i+delta);
                    printf("fin occ: %d  deb occ: %d\n",new_fin_occ,new_deb_occ);
                    for(tmpint2=new_fin_occ; tmpint2!=new_deb_occ;tmpint2--)
                        {
                        printf("Pos pile %d:\n",tmpint2);
                        afficheOcc(stdout,pocc->occ+tmpint2, longmod+1, 0);
                        }
#endif
    
                    zarb_back += spellModels(pocc, poccnew,poccsaut,
                        longmod+1, 0, curbloc+1, model, next, colors_model,
                        nbseq, tmp_quorum2, cr, nbocc_saut+1,
                        bt_intermediaire+1, new_deb_occ, new_fin_occ, cursaut,
                        longbloc, posdebbloc);
                    }
                }


            if ( cr->flag_compo == VRAI || cr->flag_compobloc[curbloc] == VRAI )
                {
                cr->compo[symbol]++;
                cr->compobloc[curbloc][symbol]++;
                }
                
            decrModel(model); /* vire la premiere lettre du nouveau bloc */
            decrModel(model); /* vire le symbole de saut */

            videPile(poccnew);
            tmpint  = copieLastOcc(poccnew,poccsaut);
#if DEBUG_SAUT
            printf("J'ai copie %d occ de Psaut->Pnew\n",tmpint);
            afficheOldOcc(poccnew,longmod+1);
#endif
            depileRec(poccsaut);
            depileRec(pocc);
            }
        
#if DEBUG_BASE
        printf("nbocc = %d\n",nbocc);
        if (nbocc<=0)
            printf("Sortie nbocc\n");
        else if (maxseq < cr->quorum)
            printf("Sortie maxseq %ld\n",maxseq);
        else 
            printf("Calcul de quorum (maxseq = %ld) tmp_quorum2=%d\n",
                    maxseq,tmp_quorum2);
#endif

    

        if ( (maxseq >= cr->quorum)
            &&  ( tmp_quorum2!=-1 ? tmp_quorum2 >= cr->quorum:
            (tmp_quorum2 = sommeBTOcc(poccnew, colors_model) ) >= cr->quorum))
            {
#if DEBUG_BASE
            printf("Accepte (res quorum=%d)\n", tmp_quorum2);
#endif

            if(symbol == numJOKER)
                zarb_ext = 1;
            else
                nbnewmod++;

            changeModel(model,symbol);

            ajouteDummy(pocc);

            transferePile2Pile(pocc, poccnew);

            if ( cr->flag_compo == VRAI || cr->flag_compobloc[curbloc] == VRAI )
                {
                cr->compo[symbol]--;
                cr->compobloc[curbloc][symbol]--;
                }

            zarb_back += spellModels(pocc, poccnew,poccsaut, longmod+1,
                longcurbloc+1, curbloc, model, next, colors_model, nbseq,
                tmp_quorum2, cr, nbocc_saut, bt_intermediaire, 0, 0, cursaut,
                longbloc, posdebbloc);

            if ( cr->flag_compo == VRAI || cr->flag_compobloc[curbloc] == VRAI )
                {
                cr->compo[symbol]++;
                cr->compobloc[curbloc][symbol]++;
                }

            depileRec(pocc);

            /* on decremente la longueur du modele de 1 */
            decrModel(model);
            } 
#if DEBUG_BASE
        else
                printf("Refuse curquorum=%d quorum=%d\n",tmp_quorum2,cr->quorum);
#endif
        }

/* Si: il n'y pas eu d'extension REGULIERE, la longueur courante est valide, */
/* ET [il n'y a pas eu d'extension bizarre (joker, saut) OU ces extensions   */
/* ont pose un probleme (modele se terminant par jokers)] */
    if ( (curbloc == cr->bloc-1)
        /*&& (nbnewmod == 0)*/ && (longmod >= cr->longueur.min)
        && ( longcurbloc >= cr->longbloc[curbloc].min )
        && ( cr->flag_palindrom == FAUX || cr->palindrom[curbloc] == -1
            || longcurbloc == longbloc[(int)(cr->palindrom[curbloc])] )
        && ( (zarb_back!=0) || (zarb_ext==0) ) )
        {
        /* A VIRER? ce test est il inutile? */
        if ( (model->name[model->lon-1] != numJOKER)
            && (model->name[model->lon-1] != numSAUT) )
            {
            longbloc[curbloc]    = longcurbloc;

            keepModel(model, pocc, nbseq, tmp_quorum, longmod, longbloc,
                    cursaut, cr);
            return(0);
            }
        else
            return(1);
        }
    }
else if ( (curbloc == cr->bloc-1))
/*         && ( ( (cr->longueur.max != 0) && (longmod == cr->longueur.max) ) */
/*         || ( (cr->longbloc[curbloc].max != 0) */
/*         && (longcurbloc <= cr->longbloc[curbloc].max) ) ) ) */
    {
    /* A VIRER? ce test est il inutile? */
    if ( (model->name[model->lon-1] != numJOKER)
        && (model->name[model->lon-1] != numSAUT))
        {
        longbloc[curbloc]    = longcurbloc;
        keepModel(model, pocc, nbseq, tmp_quorum, longmod, longbloc, cursaut,
                cr);
        return(0);
        }
    else
        return(1);
    }
return(0);
}


/******************************************************************************/
/* doSpell                                                                   */
/******************************************************************************/
/* Lance la recursion sur les modeles.                                        */
/******************************************************************************/
void doSpell(P_Criteres cr, NbSeq nbseq, Noeud  *root)
{
/* bloc est un indicateur du bloc en cours de construction */
P_mod       model;
P_PileOcc   pocc, poccnew, poccsaut;
P_occ       next;
Bit_Tab     *colors_model, ***bt_intermediaire;
int         **nbocc_saut, i, j, *cursaut;
LongSeq     *longbloc;
LongSeq     *posdebbloc=NULL;

/* Creation du faux pere de la racine (pour faciliter la recursion)                */
Noeud *root_pere = Alloc_Noeud();
root_pere->fils[Translation_Table[FINAL]] = root;
root_pere->sequence_number              = 0;
root->debut                             = 0;
root->fin                               = 1;
root->sequence_number                   = 0;
 
/* Allocation du modele                                                       */
model      = allocModel();
model->lon = 0;

/* Allocation du tableau de longueur de bloc courant                          */
longbloc    = (LongSeq *) malloc(cr->bloc * sizeof(LongSeq));
if(longbloc==NULL)
    fatalError("doSpell: cannot allocate 'longbloc'");

/* Allocation de l'occurrence courante                                        */
next = (P_occ) calloc (1,sizeof(Occ));
if (next == NULL)
    fatalError("doSpell: cannot allocate 'next'\n");
initOcc(next);

/* Allocation du tableau de bits courant et BT intermediaires                 */
bt_intermediaire    = (Bit_Tab ***) malloc((cr->bloc-1) * sizeof(Bit_Tab **));
if (bt_intermediaire == NULL)
    fatalError("doSpell: cannot allocate 'bt_intermediaire'\n");
for (i=0; i!=cr->bloc-1; i++)
    {
    bt_intermediaire[i] = (Bit_Tab **) malloc(cr->maxsaut * sizeof(Bit_Tab *));
    if (bt_intermediaire[i] == NULL)
        fatalError("doSpell: cannot allocate 'bt_intermediaire[i]'\n");
    for (j=0; j!=cr->maxsaut; j++)
        bt_intermediaire[i][j]  = AllocBitTab();
    }
colors_model = AllocBitTab();
ReinitBitTab(&colors_model);
    
    
/* Allocation du tableau de comptage du nb d'occ par tranche de saut          */
nbocc_saut  = (int **) malloc( (cr->bloc-1) * sizeof(int *));
if (nbocc_saut == NULL)
    fatalError("doSpell: cannot allocate 'nbocc_saut'\n");
for(i=0; i!=cr->bloc-1; i++)
    {
    nbocc_saut[i]   = (int *) malloc(cr->maxsaut * sizeof(int));
    if (nbocc_saut[i] == NULL)
        fatalError("doSpell: cannot allocate 'nbocc_saut[i]'\n");
    }

/* Allocation du tableau d'indication des saut effectues                      */
cursaut     = (int *) malloc((cr->bloc-1) * sizeof(int) );
if (cursaut == NULL)
    fatalError("doSpell: cannot allocate 'cursaut'\n");


/* Allocation des piles d'occurrences                                         */
poccsaut = creePileOcc();
poccnew  = creePileOcc();
pocc     = creePileOcc();

/* Allocation de la structure de palindromes                                  */
if ( cr->flag_palindrom )
    {
    posdebbloc = (LongSeq *) malloc(cr->bloc *  sizeof(LongSeq));
    if (posdebbloc == NULL)
        fatalError("doSpell: cannot allocate 'posdebbloc'\n");
    }

/* Ajout de l'occurrence nulle dans la pile d'occurrence                      */
ajouteInitOcc2Pile(pocc, root_pere);

 
fprintf(stderr,"** Models extraction **\n");
barre((int)pow((double)nbSymbMod, 3.0));
/* ...et lancement de la recursion                                            */
spellModels(pocc, poccnew, poccsaut, 0, 0, 0, model, next, &colors_model,
        nbseq, 0, cr, nbocc_saut, bt_intermediaire, 0, 0, cursaut, longbloc,
        posdebbloc);


/* Liberation des structures                                                  */
free(next);
free(model->name);
free(model);
for (i=0; i!=cr->bloc-1; i++)
    {
    for (j=0; j!=cr->maxsaut; j++)
        free(bt_intermediaire[i][j]);
    free(bt_intermediaire[i]);
    }
free(bt_intermediaire);
free(colors_model);
free(cursaut);
for(i=0; i<cr->bloc-1; i++)
    free(nbocc_saut[i]);
free(nbocc_saut);
liberePileOcc(pocc);
liberePileOcc(poccnew);
liberePileOcc(poccsaut);
}  



/******************************************************************************/
/******************************************************************************/
/********************************** MAIN **************************************/
/******************************************************************************/
/******************************************************************************/
int main(int argc, char **argv)
{
FILE            *g; 
FastaSequence   **seq;
Flag            readok;
char            infini = 0;
NbSeq           nbtxt;
float           quorum = 0.0, user_time;
int             i, j, taille, siztxt, posarg, nbmodtot, nbfiles;
long int        nbsymb;
LongSeq         maxlongsaut = 0;
Noeud           *arbre_suffixe; 
Criteres        criteres;
LongSeq         maxlongmodel, *deltatab;
Symbole         *alphaseq;

posarg  = 4;

/* QUORUM                                                                     */
quorum = atof(argv[posarg++]);

/* BLOCS                                                                      */
criteres.bloc = atoi(argv[posarg++]);
if(criteres.bloc < 2)
    fatalError("Incorrect boxes number\n");
allocBloc(&criteres, criteres.bloc);

/* LONGUEUR MIN                                                               */
criteres.longueur.min = (LongSeq)atoi(argv[posarg++]);

/* LONGUEUR MAX                                                               */
criteres.longueur.max = (LongSeq)atoi(argv[posarg++]);
 if ( criteres.longueur.max == 0 )
   infini = 1;

/* ERREURS GLOBALES                                                           */
criteres.maxerr = (LongSeq)atoi(argv[posarg++]);

maxlongmodel    = 0;

/* PARAMETRES BLOCS ***********************************************************/
if(criteres.bloc > 1)
    {
    for(i = 0; i != criteres.bloc; i++ )
        {

/* LONGUEUR MIN BLOC                                                          */
        criteres.longbloc[i].min = (LongSeq)atoi(argv[posarg++]);

/* LONGUEUR MAX BLOC                                                          */
        criteres.longbloc[i].max = (LongSeq)atoi(argv[posarg++]);
        if ( criteres.longbloc[i].max == 0 )
          infini = 1;
        maxlongmodel += criteres.longbloc[i].max;

/* ERREURS BLOC                                                               */
        criteres.maxerrblocs[i] = atoi(argv[posarg++]);

        if(i != criteres.bloc-1 )
            {
/* SAUT MIN BLOC                                                              */
            criteres.saut[i].min = (LongSeq)atoi(argv[posarg++]);

/* SAUT MAX BLOC                                                              */
            criteres.saut[i].max = (LongSeq)atoi(argv[posarg++]);
            
            maxlongsaut += criteres.saut[i].max;

/* DELTA BLOC                                                                 */
            criteres.delta[i] = (LongSeq) atoi(argv[posarg++]);
            }
        }
    }

 if ( infini == 0 )
   {
     if (maxlongmodel < criteres.longueur.max)
       maxlongmodel = criteres.longueur.max;
     
     maxlongmodel += maxlongsaut; 
   }
 else
  maxlongmodel = INT_MAX;
 


/******************************************************************************/
/* DEBUT DU TRAITEMENT                                                        */
/******************************************************************************/

/* Allocations                                                                */
seq     = (FastaSequence **) malloc(GRAINSEQ * sizeof(FastaSequence *));
text    = (signed char **) malloc(GRAINSEQ * sizeof(signed char *));
if(!seq || !text)
    fatalError("main: seq/text: cannot allocate\n");
siztxt  = GRAINSEQ;

/* Ouverture du fichier contenant les sequences                               */
g       = fopen (argv[2],"r"); 
if(g==NULL)
    fatalError("main: cannot open FASTA file");
readok  = 1;
nbtxt   = 0;
nbsymb  = 0;

/* Stockage des sequences en memoire                                          */
do
    {
    if(nbtxt == siztxt)
        {
        siztxt  *= 2;
        seq  = (FastaSequence **) realloc(seq,siztxt * sizeof(FastaSequence *));
        text = (signed char **)  realloc(text, siztxt * sizeof(signed char *));
        if(!seq || !text)
            fatalError("main: seq/text: cannot reallocate\n");
        }

    seq[nbtxt] = NewFastaSequence();
    readok     = ReadFastaSequence(g, seq[nbtxt]);
    if (readok)
        {
        nbsymb          += seq[nbtxt]->length;
        taille          = seq[nbtxt]->length+1;
        text[nbtxt]     = (signed char*) malloc ((taille+2)*sizeof(signed char)); 
        if (text[nbtxt] == NULL)
            fatalError("main: cannot allocate 'text'\n");

        strcpy((char *) text[nbtxt],seq[nbtxt]->seq);
        text[nbtxt][taille-1] = FINAL;
        text[nbtxt][taille]   = '\0';
        nbtxt++;
        }
    }
while (readok);

fclose(g); 



if (nbtxt == 0)
    fatalError("main: no sequences in FASTA file");

if (nbtxt == 1)
    fatalError("main: one sequence only in FASTA file");
        
criteres.nbsymb = nbsymb;

if (quorum == 0.0)
    criteres.quorum = (NbSeq) ceil( (double) (70*nbtxt)/100.0);
else
    criteres.quorum = (NbSeq) ceil( (double) (quorum*nbtxt)/100.0);

if(criteres.quorum==1)
    warning("quorum value is 1 sequence");
else if(criteres.quorum<1)
    fatalError("quorum value is lower than 1 sequence");
/******************************************************************************/
/* Chargement alphabet sequences et modeles                                   */
if(!(g=fopen(argv[1],"r")))
        fatalError("main: cannot open alphabet file\n");

initAlphabet();

if(!(alphaseq    = chargeAlphabet(g, (Symbole **) text, nbtxt)))
        fatalError("main: wrong alphabet file format\n");

fclose(g);

/******************************************************************************/
/* COMPOSITION (traitee apres car besoin alphabet modeles)                    */
/* S'il reste des arguments, c'est la composition                             */
setCompoPal(&criteres, argv+posarg, argc-posarg);

/* Transformation de l'alphabet (ex: AG => R)                                 */
transAlphMod(criteres.flag_palindrom);

/******************************************************************************/
/* Construction de l'arbre compact generalise                                 */
fprintf(stderr, "** Suffix tree construction **\n");
barre(nbtxt);
Init_All(alphaseq,0,nbtxt);
arbre_suffixe = Construction_Arbre((unsigned char *)text[0], maxlongmodel);
barre(0);

for (i = 1; i != nbtxt; i++)
    {
    arbre_suffixe=AjouteSequence(arbre_suffixe,(unsigned char *)text[i],maxlongmodel);
    barre(0);
    }
fprintf(stderr,"\n");

/******************************************************************************/
/* Liberation de la structure Fasta                                           */
for(i=0;i != nbtxt;i++)
    FreeFastaSequence(seq[i]);
free(seq);


UpdateBit_TabForAllTree(arbre_suffixe);

/* if (flag_tree == VRAI) */
/* Print_Tree(arbre_suffixe,1,0); */



/************************ enumeration des resultats ***************************/
printf("Extraction is going to be made with the following parameters:\n");
printf("FASTA file:                    %s\n",argv[2]);
printf("Alphabet file:                 %s\n",argv[1]);
printf("Output file:                   %s\n",argv[3]);
printf("Total min length:              %d\n",criteres.longueur.min);
if (criteres.longueur.max == 0)
    printf("Total max length:              MAX\n");
else
    printf("Total max length:              %d\n",criteres.longueur.max);
printf("Boxes:                         %d\n",criteres.bloc);
printf("Total number of subst.:        %d\n",criteres.maxerr);
printf("Quorum:                        %f%% (%d sequences in %d)\n\n",
        quorum,criteres.quorum,nbtxt);

if (criteres.flag_compo)
    {
    for (i = 0; i != nbSymbMod; i++)
        {
        if (criteres.compo[i] != -1)
            printf("Total max composition in %s:    %d\n",nummod2str[i],
                    criteres.compo[i]);
        }
    }

if (criteres.bloc > 1)
    {
    for (i = 0; i != criteres.bloc; i++)
        {
        printf("\nBOX %d\n",i+1);
        printf("Min length:                    %d\n",criteres.longbloc[i].min);
        if (criteres.longbloc[i].max == 0)
            printf("Max length:                    MAX\n");
        else
            printf("Max length:                    %d\n",
                    criteres.longbloc[i].max);
        printf("Max number of subst.:          %d\n",
                criteres.maxerrblocs[i]);
        if (i != criteres.bloc-1)
            {
            printf("Min spacer length:             %d\n",criteres.saut[i].min);
            printf("Max spacer length:             %d\n",criteres.saut[i].max);
            if(criteres.delta != 0)
                printf("Delta     :                    %d\n",criteres.delta[i]);
            else
                printf("Delta     :                    NON\n");
            }

        if (criteres.flag_compobloc[i])
            {
            for (j = 0; j != nbSymbMod; j++)
                    if (criteres.compobloc[i][j] != -1)
                        printf("Max composition in %s:         %d\n",
                                nummod2str[j],criteres.compobloc[i][j]);
            }

        if (criteres.palindrom[i]!=-1)
            printf("Palindrom of box:              %d\n", criteres.palindrom[i]+1);
        }
    }


fprintf(stderr,"\n                   ------ CHECK THESE PARAMETERS! ------\n");

criteres.maxinter   = 0;
criteres.maxsaut    = 0;
nbfiles             = 1;

for(i=0; i<criteres.bloc-1; i++)
    {
    if((j=criteres.saut[i].max-criteres.saut[i].min+1) > criteres.maxsaut)
        criteres.maxsaut    = j;
    if(criteres.delta[i] == 0)
        j   = 1;
    else
        {
        if((j=(criteres.saut[i].max-criteres.saut[i].min+1-criteres.delta[i]*2))
                > criteres.maxinter)
            criteres.maxinter   = j;
        }
    nbfiles *= j;
    }

f = (FILE **) malloc (nbfiles * sizeof(FILE *) );
nbmod   = (int *) calloc(nbfiles, sizeof(int));
if(f==NULL || nbmod==NULL)
    fatalError("main: cannot allocate 'f/nbmod'");

if(!initFiles(f,argv[3],&criteres))
    fatalError("main: cannot initialize output files\n");

initTabSauts(&criteres);

/* Initialisax   de maxlongmod  pour recalcul de vraie longueur max           */
maxlongmod  = (LongSeq *) calloc(nbfiles , sizeof(LongSeq));
/* Idem maxlongbloc                                                           */
maxlongbloc = (LongSeq **) malloc(nbfiles * sizeof(LongSeq *));
if(maxlongmod==NULL || maxlongbloc==NULL)
    fatalError("main: cannot allocate 'maxlongmod/maxlongbloc'");
for(i=0; i!=nbfiles; i++)
    {
    maxlongbloc[i]  = (LongSeq *) calloc (criteres.bloc , sizeof(LongSeq));
    if(maxlongbloc[i]==NULL)
        fatalError("main: cannot allocate 'maxlongbloc[i]'");
    }


/******************************************************************************/
/******************************************************************************/
/* Fonction Principale                                                        */
PrintCpuTime(1);
doSpell(&criteres,nbtxt,arbre_suffixe);
user_time=PrintCpuTime(0);
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/* Insertion de la ligne d'information en tete des fichiers de sortie         */
deltatab    = (LongSeq *) malloc(criteres.bloc * sizeof(LongSeq));
if(deltatab==NULL)
    fatalError("main: cannot allocate 'deltatab', header of output will be wrong");

for(i=0,nbmodtot=0;i!=nbfiles;i++)
    {
    fprintf(f[i],"Nb models: %d\nUser time : %.2f sec.\n", nbmod[i], user_time);
    nbmodtot+=nbmod[i];

/* Insertion de la ligne de parametres dans le fichier de sortie              */
    rewind(f[i]);
    fprintf(f[i],"%%%%%% %d %d/%d %ld %d %d %d",criteres.bloc,
            criteres.quorum, nbtxt, criteres.nbsymb, criteres.longueur.min,
            maxlongmod[i], criteres.maxerr);

    file2Delta(i, deltatab, &criteres);

    for(j=0; j!=criteres.bloc; j++)
        {
        fprintf(f[i]," %d %d %d",criteres.longbloc[j].min, maxlongbloc[i][j],
            criteres.maxerrblocs[j]);

        if(j!=criteres.bloc-1)
            fprintf(f[i]," %d %d", criteres.saut[j].min+deltatab[j],
                    criteres.saut[j].min+deltatab[j]+2*criteres.delta[j]);
        }

/* Ecriture du nom du fichier alphabet utilise et de l'alphabet des sequences */
    fprintf(f[i], " %s %s", argv[1], alphaseq);

    fclose(f[i]);
    }
printf("\nNb models: %d\nUser time : %.2f sec.\n", nbmodtot, user_time);


for(i=0;i!=nbtxt;i++)
    free(text[i]);
free(text);

/* Liberations                                                                */
free(f);
free(nbmod);
Free_Arbre(arbre_suffixe);
return(0); 
}

/******************************************************************************/
/* PrintCpuTime                                                               */
/******************************************************************************/
static float PrintCpuTime(char initIt)
{
float         ust;
struct tms    tms;
static float  dust;

times(&tms);

ust =  (float) tms.tms_utime;

if (initIt)
    {
    dust = ust;
    return 0.0;
    }
else
    {
    ust -= dust;
    return  ust / sysconf(_SC_CLK_TCK);
    }
}

