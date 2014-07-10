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

/* (C) by "coward" <coward -AT- ii.uib.no> from 1996-1999 */
/*   modified and extended by Laurent Marsan <lama@prism.uvsq.fr> 1999-2004.  */

/******************************************************************************/
/* SIGSTAT - Version avec grep+                                               */
/* (Shuffling par Eivind Coward)                                              */
/******************************************************************************/
#include <sigstat.h>

char            *alphaseq;       /* Alphabet des sequences */
extern  char    *nummod2str[127];

/******************************************************************************/
/* Prototypes prives                                                          */
/******************************************************************************/
void    afficheStats(FILE *outfile, Mot *modeles, int nbmodeles, int nbseq,
            int maxlongmod);

void    calculeStats( int nbmodeles, int nbtests, int nbseq,
                        float **res, Mot *modeles,   char    flag);


/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/
int main(int argc, char **argv)
{
Mot             *modeles;
int             nbseq, nbmodeles, i, nbtests, maxsizeseq, j, maxlongmod=0;
unsigned int    nbocc, nboccex;
float           **resultats;
FILE            *fasta, *res, *outfile, *f;
/* Shufflet's variables                                                       */
int             *nver,  **count, **vdeg, *first, *last, *count1,
                *vdeg1, nklets, nk1lets, *seqlen, order, nbseqalloc;
int             *lastedge;  /* array to indicate last edge from each vertex*/
char            **seqstart=NULL; /* to store the beginning of sequences */
char            **seq, line[BUF], *alphaseq, **origseq;
/* int             k; */

Arbre           arbre;
Criteres        cr;

#if OCC
float           **resocc;
#endif



initCriteres(&cr);

if(argc!=6)
    {
    printf("Usage: ~ Fic.Fasta Fic.res Fic.sortie nb_shufflings order\n\n");
    return 1;
    }


/** LECTURE DES PARAMETRES ****************************************************/
fasta   = openFile(argv[1], "r");

res     = openFile(argv[2], "r");

outfile = openFile(argv[3], "w");

nbtests = atoi(argv[4]);

order   = atoi(argv[5]);
/******************************************************************************/

printCpuTime(NULL);


/* Allocations necessaires a shufflet                                         */
count   = (int **) malloc(sizeof(int*)*GRAINSEQ);
seqlen  = (int *) malloc(sizeof(int)*GRAINSEQ);
vdeg    = (int **) malloc(sizeof(int*)*GRAINSEQ);      
first   = (int *) malloc(sizeof(int)*GRAINSEQ);
last    = (int *) malloc(sizeof(int)*GRAINSEQ);
nver    = (int *) malloc(sizeof(int)*GRAINSEQ);
nbseqalloc  = GRAINSEQ;
assert(count != NULL);
assert(vdeg != NULL);
assert(seqlen != NULL && first != NULL && last != NULL);
assert(nver != NULL);

if(order == 1)
    {
    origseq = (char **) malloc(GRAINSEQ*sizeof(char *));
    assert(origseq != NULL);
    }
else
	{
	seqstart = (char **) malloc(sizeof(char *)*GRAINSEQ);
	assert(seqstart != NULL);
	}

/******************************************************************************/


/* Lecture des criteres                                                       */
fgets(line, BUF, res);
if(!chargeCriteres(&cr, line))
    {
    fprintf(stderr, "File '%s' is corrupted.\n",argv[2]);
    return(1);
    }

if(valideCriteres(&cr) == FAUX)
        return 1;

/******************************************************************************/
/* Chargement alphabet sequences                                              */
initAlphabet();
alphaseq   = chargeAlphaSeq(NULL, 0, cr.alphaseq);

/******************************************************************************/
/* Chargement alphabet sequences et modeles                                   */
if(!(f=fopen(cr.ficalph,"r")))
    fatalError("main: cannot open alphabet file\n");

chargeAlphaMod(f);

fclose(f);
/******************************************************************************/


/**LECTURE DU FICHIER FASTA ***************************************************/
fprintf(stderr,"** Reading composition of the sequences to shuffle **\n");
nbseq = readseq(order,fasta,&nver,&count,&vdeg,&first,&last,&seqstart,&seqlen,
&maxsizeseq, nbseqalloc, cr.alphaseq, &origseq);
fclose(fasta);

if(nbseq<=0)
    fatalError("no sequences in file\n");


/**LECTURE DES RESULTATS DU PROGRAMME******************************************/
/* Lecture des modeles                                                        */
fprintf(stderr,"** Reading extracted models **\n");
modeles = (Mot *) lectureFichierRes(res, nbseq, &nbmodeles, &maxlongmod);
fclose(res);

if ( modeles == NULL )
    {
    fprintf(stderr, "STOP: File '%s' is empty.\n", argv[2]);
    return 0;
    }


/******************************************************************************/
/* ALLOCATIONS                                                                */
/* Allocation des tableaux de stockage des statistiques                       */
resultats   = (float **)malloc(nbtests*sizeof(float*));
assert(resultats!=NULL);
#if OCC
resocc      = (float **)malloc(nbtests*sizeof(float*));
assert(resocc!=NULL);
#endif

for ( i=0; i<nbtests; i++ )
    {
    resultats[i]    = (float*)  calloc(nbmodeles,sizeof(float));
    assert(resultats[i]!=NULL);
#if OCC
    resocc[i]       = (float*)  calloc(nbmodeles,sizeof(float));
    assert(resocc[i]!=NULL);
#endif
    }

/* Allocation des structures pour shuffling                                   */
nklets  = pow(strlen(alphaseq)-1,order);
nk1lets = nklets/(strlen(alphaseq)-1);
count1  = (int *) malloc(sizeof(int)*nklets);
vdeg1   = (int *) malloc(nk1lets*sizeof(int));
lastedge= (int *) malloc(sizeof(int)*nk1lets);

seq     = (char **) malloc(nbseq*sizeof(char*));
assert(count1!=NULL && vdeg1!=NULL && lastedge!=NULL);

for (i=0; i<nbseq; i++)
    {
    seq[i]  = (char *) malloc((seqlen[i]+2+order)*sizeof(char));  
    assert(seq[i]!=NULL);
    }




/**BOUCLE PRINCIPALE **********************************************************/
fprintf(stderr,"\n**** Starting processus ****\n");
barre(nbtests);
for( i=0; i<nbtests; i++ )
    {
/* Pour ordre 1, on recopie les sequences originales pour monoshuffling       */
    if(order==1)
        {
        for(j=0;j!=nbseq;j++)
            strncpy(seq[j],origseq[j],seqlen[j]);
        }

/*     fprintf(stderr,"\n**** Set %d/%d ****\n",i+1,nbtests); */
    generateseq(order,nver,count,vdeg,first,last,seqstart,nbseq,seqlen,
        count1, vdeg1, nklets, nk1lets, lastedge, seq);

/* Affiche les sequences shufflees                                            */
/*     for(j=0; j!=nbseq;j++) */
/*         { */
/*         printf("> %d\n",j); */
/*         k=0; */
/*         while(seq[j][k]!='\0') */
/*             { */
/*             printf("%c",seq[j][k]); */
/*             k++; */
/*             if(k%50 ==0) */
/*                 printf("\n"); */
/*             } */
/*         printf("\n");; */
/*         } */


/* Creation de l'arbre suffixe                                                */
    creeArbreSuffixeFromArray(&arbre, seq, nbseq, maxLongMod(cr), alphaseq);
/*     fprintf(stderr,"** Searching for occurrences of the %d models **\n", */
/*             nbmodeles); */

/* Gestion de la barre                                                        */
/*     if(nbmodeles>=100) */
/*         { */
/*         step    = nbmodeles/100; */
/*         barre(100); */
/*         } */
/*     else */
/*         { */
/*         step    = nbmodeles; */
/*         barre(nbmodeles); */
/*         } */


    for(j=0; j!=nbmodeles; j++)
        {
/*         fprintf(stderr,"Je cherche le mot %s\n",modeles[j].mot);  */
/*         getchar(); */
        chercheMot(arbre,&cr, modeles[j].codes, &nbocc, &nboccex, NULL);
/*         printf("J'en ai trouve %d et %d\n",nbocc, nboccex); */

/*         if(nbmodeles>=100) */
/*             { */
/*             if(j%step==0) */
/*                 barre(0); */
/*             } */
/*         else */
/*             barre(0); */

        resultats[i][j]             = (float) nbocc;
        modeles[j].moyenne_shuffle  += ((float) nbocc)/nbtests;
/*         printf("Sommeres -> %f\n",modeles[j].moyenne_shuffle); */

#if OCC
        resocc[i][j]                  = (float) nboccex;
        modeles[j].moyenne_shuffle_occ+= ((float) nboccex)/nbtests;
#endif
        }

/* Liberation de l'arbre suffixe                                              */
    libereArbreSuffixeFromArray(arbre);

    barre(0);
    }

/******************************************************************************/
/* Liberations                                                                */
for(i=0; i<nbseq; i++)
    free(seq[i]);
free(seq);

free(count);
free(vdeg);
free(first);
free(last);
free(nver);

if (order>1)
	free(seqstart);

/******************************************************************************/
/* CALCULS STATISTIQUES                                                       */

/* for( i=0; i<nbmodeles; i++ ) */
/*     printf("Le modele %s a une moyenne %f\n",modeles[i].mot,
       modeles[i].moyenne_shuffle); */

calculeStats(nbmodeles, nbtests, nbseq, resultats, modeles, 0);

#if OCC
calculeStats(nbmodeles, nbtests, cr.nbsymb, resocc, modeles,1);
#endif

/* Affichage des resultats                                                    */
afficheStats(outfile, modeles, nbmodeles, nbseq, maxlongmod);
printCpuTime(outfile);
fclose(outfile);

/* for(j=0;j!=4;j++) */
/*     { */
/*     fprintf(outfile, "%s\n",modeles[j].mot); */
/*     for(i=0; i!=nbtests; i++) */
/*         fprintf(outfile, "%f ",resultats[i][j]); */
/*     fprintf(outfile, "\n\n"); */
/*     } */
/******************************************************************************/
/* Desallocations des structures stats                                        */
free(seqlen);
for(i=0;i<nbtests;i++)
    {
    free(resultats[i]);
#if OCC
    free(resocc[i]);
#endif
    }
free(resultats);

#if OCC
free(resocc);
#endif

return 0;
}
    
    
int comparModeles(Mot *a, Mot *b)
{
if(a->zscore == b->zscore)
    return 0;
if(a->zscore > b->zscore)
    return -1;
return 1;
}

#if  OCC
int comparModelesOcc(Mot *a, Mot *b)
{
if(a->zscore_occ == b->zscore_occ)
    return 0;
if(a->zscore_occ > b->zscore_occ)
    return -1;
return 1;
}
#endif

/******************************************************************************/
/* espace                                                                     */
/******************************************************************************/
void    espace(FILE *f, int nb)
{
static char space[31]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};

fprintf(f,"%s",space+30-nb);
}

    
/******************************************************************************/
/* afficheStats                                                               */
/******************************************************************************/
void afficheStats(FILE *outfile, Mot *modeles, int nbmodeles, int nbseq,
        int maxlongmod)
{
int     i;
Mot     m;

fprintf(outfile,"STATISTICS ON THE NUMBER OF SEQUENCES HAVING AT LEAST ONE OCCURRENCE\n");
fprintf(outfile,"Model");
espace(outfile,maxlongmod-2);
fprintf(outfile," %%right\t #right\t %%shfl.\t #shfl.\tSigma\tChi2\tZ-score\n");
fprintf(outfile,"=================================================================================\n");

qsort(modeles, nbmodeles, sizeof(Mot),
        (int (*)(const void *, const void *)) comparModeles);

for( i=0; i<nbmodeles; i++ )
    {
    m   = modeles[i];

    fprintf(outfile,"%s", m.mot);
    espace(outfile, maxlongmod+4-strlen(m.mot));
    fprintf(outfile,"%3.2f%%\t%5d\t%3.2f%%\t%5.2f\t%4.2f\t%3.2f\t%3.2f\n",
        m.quorum_reel, m.nbseq_vrai, m.moyenne_shuffle*100/nbseq,
        m.moyenne_shuffle, m.sigma, m.khi2, m.zscore);
    }

#if OCC
fprintf(outfile,"\nSTATISTICS ON THE TOTAL NUMBER OF OCCURRENCES\n");
fprintf(outfile,"Model");
espace(outfile,maxlongmod-2);
fprintf(outfile," #right\t #shfl.\tSigma\tChi2\tZ-score\n");
fprintf(outfile,"=================================================================================\n");

qsort(modeles, nbmodeles, sizeof(Mot),
        (int (*)(const void *, const void *))comparModelesOcc);

for( i=0; i<nbmodeles; i++ )
    {
    m   = modeles[i];
    fprintf(outfile,"%s",m.mot);
    espace(outfile, maxlongmod+4-strlen(m.mot));
    fprintf(outfile, "%5d\t%5.2f\t%4.2f\t%3.2f\t%3.2f\n",
            m.nboccex_vrai, m.moyenne_shuffle_occ, m.sigma_occ, m.khi2_occ,
            m.zscore_occ);
    }
#endif

}



/******************************************************************************/
/* calculeStats                                                                */
/******************************************************************************/
void   calculeStats(    int     nbmodeles,  int     nbtests,       int  nbseq,
                        float   ** resultats,  Mot *modeles,       char flag)
{
float   Pv, Pf, Av, Af, nP, nA, tmp, tmp2, sigma;
int i,j;

for ( i = 0; i < nbmodeles; i++ )
    {
/*     printf("Le modele %s a %d %d occs\n", modeles[i].mot,modeles[i].nbseq, */
/*             modeles[i].nboccex); */

#if OCC
    if ( flag == 1 )
        {
        Pv  = (float)   modeles[i].nboccex_vrai;
        Pf  = modeles[i].moyenne_shuffle_occ;
        }
    else
#endif
    {
    Pv  = (float)   modeles[i].nbseq_vrai;
    Pf  = modeles[i].moyenne_shuffle; 
    }

    Av  =   (float) nbseq-Pv;
    Af  =   (float) nbseq-Pf;
    nP  =   Pv+Pf;
    nA  =   Av+Af;


/* CALCUL DU KHI2                                                             */
#if OCC
    if(flag == 1)
        modeles[i].khi2_occ =    2*CARRE(Pv*Af-Pf*Av)/(nP*nA*nbseq);
    else
#endif
        modeles[i].khi2 =    2*CARRE(Pv*Af-Pf*Av)/(nP*nA*nbseq);
/*     printf("Khi2: %f\n",modeles[i].khi2); */


/* CALCUL DU Z-SCORE                                                          */
    tmp =   0;
    for ( j = 0; j < nbtests; j++ )
        {
        tmp2    =   CARRE(resultats[j][i]-Pf)/(nbtests-1);
        tmp +=  tmp2;
        }
    sigma    =   (float) sqrt((double) tmp);
/* printf("Sigma: %f\n", sigma); */
    if ( sigma == 0.000000 )
        {
#if  OCC
        if(flag == 1)
            {
            modeles[i].sigma_occ    =   sigma;
            modeles[i].zscore_occ   =   UINT_MAX;
            }
        else
#endif
            {
            modeles[i].sigma    =   sigma;
            modeles[i].zscore   =   UINT_MAX;
            }
        }
    else
#if  OCC
        if(flag == 1)
            {
            modeles[i].sigma_occ    =   sigma;
            modeles[i].zscore_occ   =   (Pv-Pf)/sigma;
/*     printf("Z-score: %f\n",modeles[i].zscore_occ); */
            }
        else
#endif
            {
            modeles[i].sigma    =   sigma;
            modeles[i].zscore   =   (Pv-Pf)/sigma;
/*     printf("Z-score: %f\n",modeles[i].zscore); */
            }
    }
}
