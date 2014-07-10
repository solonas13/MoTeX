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

#include <pile_occ.h>


/* FONCTIONS PRIVEES                                                          */
#if OCC
int recAfficheOcc(FILE *f,Noeud * n, LongSeq l, P_Criteres cr, int codesaut);
#endif


extern  char    **text;
extern  int     nbSymbSeq;
extern  int     carseq2num[127];


/******************************************************************************/
/* creePileOcc                                                                */
/******************************************************************************/
P_PileOcc creePileOcc(void)
{
P_PileOcc p;

p=(P_PileOcc)malloc(sizeof(PileOcc));
if(p==NULL)
    fatalError("creePileOcc: cannot allocate 'p'\n");

p->occ=(P_occ)calloc(GRAIN, sizeof(Occ));
if(p->occ==NULL)
    fatalError("creePileOcc: cannot allocate 'p->occ'\n");

p->carte=(unsigned int *)malloc(GRAIN_SIZMOD*sizeof(unsigned int));
if(p->carte==NULL)
    fatalError("creePileOcc: cannot allocate 'p->carte'\n");

p->size=GRAIN;
p->size_carte=GRAIN_SIZMOD;
p->pos_carte=0;
p->pos=0;

ajouteDummy(p);

return(p);
}

/* MODE D'EMPLOI DES DUMMYS EN MULTIBLOC DELTA                                */
/* Si le parametre principal vaut VRAI c'est un dummy de separation entre
 * recursions. */
/* Sinon c'est un dummy de separation entre sauts. */
/* Dans les deux cas le dernier parametre indique le code d'intervalle
 * concerne par les occurrences qui suivent.                                  */

/******************************************************************************/
/* ajouteDummy                                                                */
/******************************************************************************/
void ajouteDummy(P_PileOcc p)
{
#if DEBUG_PILE
unsigned int * t=p->carte;
#endif

if(p->pos_carte>=p->size_carte)
    {
    p->size_carte+=GRAIN_SIZMOD;
#if DEBUG_PILE
    printf("J'etends carte a %d\n",p->size_carte);
#endif

    p->carte=(unsigned int *)realloc(p->carte,p->size_carte
            *sizeof(unsigned int));
    if(p->carte==NULL)
        fatalError("pile_occ.c: ajouteDummy: cannot reallocate 'p->carte'");
#if DEBUG_PILE
    if(t!=p->carte)
        printf("CHANGEMENT D'emplacement memoire de la carte\n");
#endif
    }

p->carte[p->pos_carte]=p->pos;
p->pos_carte++;

ajouteOcc2Pile(p, NULL, -1, -1, -1, -1, -1, -1);
}

/******************************************************************************/
/* getPrecDummy                                                               */
/******************************************************************************/
LongSeq getPrecDummy(P_PileOcc p)
{
/* if(p->pos_carte==0)  POSSIBLE DANGER*/ 
return(p->carte[p->pos_carte-1]);
}

/******************************************************************************/
/* ajouteInitOcc2Pile                                                         */
/******************************************************************************/
void ajouteInitOcc2Pile(P_PileOcc p, Noeud *x)
{
P_occ ptr;

ptr=p->occ+p->pos;
ptr->x=x;
ptr->num=carseq2num[(int) FINAL];    /* Symbole quelconque, non lu */
ptr->lon=1;
ptr->xerr=0;
ptr->blocerr=0;
ptr->saut=0;
ptr->codesaut=0;

p->pos++;
}

/******************************************************************************/
/* ajouteOcc2Pile                                                             */
/******************************************************************************/
void ajouteOcc2Pile(P_PileOcc p, Noeud *x, int num, LongSeq lon, LongSeq err,
        LongSeq blocerr, LongSeq saut, int codesaut)
{
P_occ ptr;

#if DEBUG_PILE
printf("j'ecris une nouvelle occurrence en %d\n", p->pos);
#endif

if(p->pos>=p->size)
    {
#if DEBUG_PILE
    printf("JE RESIZE (ajoute)\n");
#endif
    p->size+=GRAIN;
    ptr = p->occ;
    p->occ=(P_occ)realloc(p->occ, (p->size)*sizeof(Occ));
    if(p->occ==NULL)
        fatalError("ajouteOcc2Pile: cannot reallocate 'p->occ'\n");
#if DEBUG_PILE
    if(p->occ!=ptr)
        printf("changement d'emplacement memoire de pileocc\n");
#endif
    }
ptr             = p->occ+p->pos;
ptr->x          = x;
ptr->num        = num;
ptr->lon        = lon;
ptr->xerr       = err;
ptr->blocerr    = blocerr;
ptr->saut       = saut;
ptr->codesaut   = codesaut;

p->pos++;
}

/******************************************************************************/
/* copieLastOcc                                                               */
/******************************************************************************/
int copieLastOcc(P_PileOcc dest, P_PileOcc source)
{
int last_dummy, nbocc;

if(source->pos==0)
    return 0;

if(source->pos_carte==0)
    last_dummy=-1;
else
    last_dummy=source->carte[source->pos_carte-1];

nbocc   = source->pos-last_dummy;
if(nbocc==1)
    return 1;

if(((dest->pos)+(nbocc-1))>=(dest->size))
    {
#if DEBUG_PILE
    printf("JE RESIZE (copie)  posdest %d sizedest %d  possource %d\n",
            dest->pos, dest->size, source->pos);
#endif
    dest->size=(int)ceil(((double)(dest->pos+nbocc-1))/((double)GRAIN))*GRAIN;
#if DEBUG_PILE
    printf("New size %d\n",dest->size);
#endif
    dest->occ=(P_occ)realloc(dest->occ, dest->size*sizeof(Occ));
    if(dest->occ==NULL)
        fatalError("ajouteOcc2Pile: cannot reallocate 'dest->occ'\n");
    }
memcpy(dest->occ+dest->pos, source->occ+last_dummy+1, (nbocc-1)*sizeof(Occ));
dest->pos+=nbocc-1;
return nbocc-1;
}

/******************************************************************************/
/* transferePile2Pile                                                         */
/******************************************************************************/
void transferePile2Pile(P_PileOcc dest, P_PileOcc source)
{
#if DEBUG_PILE
P_occ t=dest->occ;
#endif

if(source->pos==0)
    return;

if(((dest->pos)+(source->pos))>=(dest->size))
    {
#if DEBUG_PILE
    printf("JE RESIZE (transfere)  posdest %d sizedest %d  possource %d\n",
            dest->pos, dest->size, source->pos);
#endif
    dest->size=(int)ceil(((double)(dest->pos+source->pos))/((double)GRAIN))*GRAIN;
#if DEBUG_PILE
    printf("New size %d\n",dest->size);
#endif
    dest->occ=(P_occ)realloc(dest->occ, dest->size*sizeof(Occ));
    if(dest->occ==NULL)
        fatalError("ajouteOcc2Pile: cannot reallocate 'dest->occ'\n");
#if DEBUG_PILE
    if(t!=dest->occ)
        printf("chgmnt d'empl memoire de pocc (transfere)\n");
#endif
    }
memcpy(dest->occ+dest->pos, source->occ, source->pos*sizeof(Occ));
dest->pos+=source->pos;

source->pos=0;
}

/******************************************************************************/
/* depileRec                                                                  */
/******************************************************************************/
void depileRec(P_PileOcc p)
{
if(p->pos_carte==0)
    p->pos=0;
else
    {
    p->pos_carte--;
    p->pos=p->carte[p->pos_carte];
    }
}

/******************************************************************************/
/* libereOcc                                                                  */
/******************************************************************************/
void videPile(P_PileOcc p)
{
p->pos          = 0;
p->pos_carte    = 0;
ajouteDummy(p);
}

/******************************************************************************/
/* liberePileOcc                                                              */
/******************************************************************************/
void liberePileOcc(P_PileOcc p)
{
free(p->occ);
free(p->carte);
free(p);
}

#if DEBUG_BASE
/******************************************************************************/
/* affichePileOcc                                                             */
/******************************************************************************/
void affichePileOcc(P_PileOcc p)
{
int i;

for(i=p->pos-1; i>=0; i--)
    if(p->occ[i].x==NULL)
        printf("====== DUMMY =======\n");
    else
        printf("num %p  branche %c lon %d err %d\n",p->occ[i].x,
                lettres[p->occ[i].num], p->occ[i].lon, p->occ[i].xerr);
/*         printf("num %d  branche %c lon %d err %d\n",p->occ[i].x->numero, */
/*                 lettres[p->occ[i].num], p->occ[i].lon, p->occ[i].xerr); */
printf("--------------------------------------\n");
}
#endif

#if OCC
/******************************************************************************/
/* afficheLastOcc                                                             */
/******************************************************************************/
/* Lance l'affichage de tous les motifs associes aux occurrences trouvees     */
/******************************************************************************/
void afficheLastOcc(FILE *f, P_PileOcc pocc, LongSeq l, P_Criteres cr)
{
int i=pocc->pos-1,nbocc=0;
P_occ ptr;

ptr=pocc->occ+i;
while((i>0) && (ptr->x!=NULL))
    {
    nbocc+=recAfficheOcc(f,ptr->x->fils[ptr->num], l-ptr->lon+ptr->saut, cr,
            ptr->codesaut);
    ptr--;
    i--;
    }
fprintf(f,"%d\n",nbocc);
}

/******************************************************************************/
/* recAfficheOcc                                                              */
/******************************************************************************/
/* Parcourt l'arbre recursivement pour atteindre les feuilles et affiche      */
/******************************************************************************/
int recAfficheOcc(FILE *f,Noeud * n, LongSeq l, P_Criteres cr, int codesaut)
{
int i,nbocc=0;

/* #if DEBUG_BASE */
/* printf("J'entre dans recAffiche avec long %d\n",l); */
/* #endif */

/* Si on a atteint une feuille                                                */
if (n->debut & LEAF_BIT)
    {
    nbocc   += Print_Positions(f, (Feuille *) n, l, cr, codesaut);
    }
else
    for(i=0; i!=nbSymbSeq+1; i++)
        if(n->fils[i])
        {
/* #if DEBUG_BASE */
/*             printf("recAffiche: je passe par %c\n",lettres[i]); */
/* #endif */
            nbocc   += recAfficheOcc(f,n->fils[i], l+n->fin-n->debut, cr,
                        codesaut);
        }
return nbocc;
}

/******************************************************************************/
/* affOcc                                                                     */
/******************************************************************************/
int afficheOcc(FILE *f, P_occ o, LongSeq longmod, P_Criteres cr)
{
if(o->x != NULL)
    return(recAfficheOcc(f,o->x->fils[o->num], longmod-o->lon+o->saut, cr,
            o->codesaut));
return 0;
}
#endif

#if DEBUG_BASE
/******************************************************************************/
/* afficheOldOcc                                                              */
/******************************************************************************/
/* Affiche les occurrences d'un niveau -n dans la pile                        */
/******************************************************************************/
void afficheOldOcc(P_PileOcc p, LongSeq l)
{
int pos=p->pos-1;
P_occ   tmpocc=p->occ+pos;

printf("=======HAUT=PILE=========\n");

while(pos >=0 && tmpocc != NULL && tmpocc->lon >= 0)
    {
    printf("*** num %d  lon %d  saut %d  codesaut %d  lmod %d\n",tmpocc->num,
            tmpocc->lon,tmpocc->saut,tmpocc->codesaut, l);
    afficheOcc(stdout,tmpocc,l,0);
    tmpocc--;
    pos--;
    }
printf("=*=*=*=*=DUMMY=*=*=*=*=*=\n");
}
#endif
