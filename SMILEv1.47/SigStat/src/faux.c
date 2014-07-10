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

/******************************************************************************/
/* FAUX    - Version P motifs avec grep+                                      */
/******************************************************************************/
#include <faux.h>

extern  char    *nummod2str[127];

/******************************************************************************/
/* Prototypes prives                                                          */
/******************************************************************************/
void   afficheKhi2(Mot *modeles, int nbmodeles,  int nbseq_vrai,
                int nbseq_faux, long int nbsymbV, long int nbsymbF, FILE *out,
                int maxlongmod);


/******************************************************************************/
/* MAIN                                                                       */
/******************************************************************************/
int main(int argc, char **argv)
{
Mot             *modeles;
int             nbseq_vrai, nbseq_faux, nbmodeles, step, j, maxlongmod=0;
FILE            *res, *outfile, *f;
unsigned int    nbocc, nboccex;
long int        nbsymbF=0;
char            line[BUF];
Arbre           arbre;
Criteres        cr;

#if  OCC
char            *ptr;
#endif



initCriteres(&cr);

if(argc!=4)
    {
    printf("Usage: ~ Fic.Fasta.Faux Fic.res Fic.sortie\n\n");
    return 1;
    }

printCpuTime(NULL);

/** LECTURE DES PARAMETRES ****************************************************/
res     = openFile(argv[2], "r");

outfile = openFile(argv[3], "w");
/******************************************************************************/



/**LECTURE DES RESULTATS DU PROGRAMME******************************************/
/* Lecture des criteres                                                       */
fgets(line, BUF, res);
if(!chargeCriteres(&cr, line))
    {
    fprintf(stderr, "File '%s' corrupted.\n",argv[2]);
    return(1);
    }

if(valideCriteres(&cr) == FAUX)
        return 1;



/******************************************************************************/
/* Creation de l'arbre suffixe                                                */
initAlphabet();
creeArbreSuffixeFromFile(&arbre, argv[1], maxLongMod(cr), NULL);
nbseq_faux = arbre.nbtxt;

/******************************************************************************/
/* Chargement alphabet des modeles                                            */
if(!(f=fopen(cr.ficalph,"r")))
    fatalError("main: cannot open alphabet file\n");

chargeAlphaMod(f);

fclose(f);
/******************************************************************************/


/******************************************************************************/
/* Lecture des modeles                                                        */
nbseq_vrai   = cr.nbtotseq;
modeles = (Mot *) lectureFichierRes(res, nbseq_vrai, &nbmodeles, &maxlongmod);
/* printf("nbmodeles %d\n",nbmodeles); */
fclose(res);

if ( modeles == NULL )
    {
    fprintf(stderr, "Warning: File '%s' is empty.\n", argv[2]);
    exit(1);
    }




/**BOUCLE PRINCIPALE **********************************************************/
fprintf(stderr,"** Searching for occurrences of the %d models **\n",
    nbmodeles);

if(nbmodeles>=100)
    {
    step    = nbmodeles/100;
    barre(100);
    }
else
    {
    step    = nbmodeles;
    barre(nbmodeles);
    }

for(j=0; j!=nbmodeles; j++)
    {
/*         printf("Je cherche le mot %s\n",modeles[j].mot);  */
    chercheMot(arbre,&cr, modeles[j].codes, &nbocc, &nboccex, NULL);
/*         printf("J'en ai trouve %d et %d\n",nbocc, nboccex); */
    
    if(nbmodeles>=100)
        {
        if(j%step==0)
            barre(0);
        }
    else
        barre(0);
    
    modeles[j].nbseq_faux   = nbocc;
    
#if OCC
    modeles[j].nboccex_faux = nboccex;
#endif
    }




/******************************************************************************/
/* CALCULS STATISTIQUES                                                       */
/******************************************************************************/

fprintf(outfile, "Original sequences (%d seq) against '%s' (%d seq)\n", 
        nbseq_vrai, argv[1], nbseq_faux);

#if OCC
nbsymbF = 0;
for(j=0; j!=nbseq_faux; j++)
    {
    ptr = arbre.text[j];
    while(*ptr != FINAL)
        {
/*         printf("%c\n",*ptr); */
        ptr++;
        nbsymbF++;
        }
    }
/* printf("%d\n",nbsymbF); */
#endif

afficheKhi2(modeles, nbmodeles, nbseq_vrai, nbseq_faux, cr.nbsymb, nbsymbF,
        outfile, maxlongmod);

/* Liberation de l'arbre suffixe                                              */
libereArbreSuffixeFromFile(arbre);

printCpuTime(outfile);
fclose(outfile);

free(modeles);

return 0;
}

int comparModeles(Mot *a, Mot *b)
{
if(a->khi2 == b->khi2)
    return 0;
if(a->khi2 > b->khi2)
    return -1;
return 1;
}

/******************************************************************************/
/* espace                                                                     */
/******************************************************************************/
void    espace(FILE *f, int nb)
{
static char space[31]={' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ','\0'};

fprintf(f,"%s",space+30-nb);
}

/******************************************************************************/
/* afficheKhi2                                                                */
/******************************************************************************/
void   afficheKhi2(  Mot *modeles, int nbmodeles,  int nbseq_vrai,
                int nbseq_faux, long int nbsymbV, long int nbsymbF, FILE *out,
                int maxlongmod)
{
Mot m;
float   Pv, Pf, Av, Af, nP, nA, Fv=0, Ff=0, nbseqtot;
int     i;
#if  OCC
float   Pv_ex, Pf_ex, Av_ex, Af_ex, nP_ex, nA_ex, Fv_ex=0, Ff_ex=0;
long int    nbsymbtot;
#endif

for(i=0;i!=nbmodeles;i++)
    {
    m   = modeles[i];
    Pv  = (float) m.nbseq_vrai;
    Fv  = (float) Pv / nbseq_vrai;
    Pf  = (float) m.nbseq_faux;
    Ff  = (float) Pf / nbseq_faux;
    Av  = (float) nbseq_vrai-Pv;
    Af  = (float) nbseq_faux-Pf;
    nP  = (float) Pv+Pf;
    nA  = (float) Av+Af;
    nbseqtot    = nbseq_vrai + nbseq_faux;

#if  DEBUG
    fprintf(out,"%s\n",modeles[i].mot);
    fprintf(out,"Pv %f\tFv %f\tPf %f\tFf %f\n",Pv,Fv,Pf,Ff);
    fprintf(out,"Av %f\tAf %f\tnP %f\tnA %f\n",Av,Af,nP,nA);
    fprintf(out,"Nbseq vrai %d\tNbseq faux %d\n",nbseq_vrai,nbseq_faux);
#endif
    
    modeles[i].khi2    = nbseqtot * CARRE(Pv*Af-Pf*Av) / (nP*nA*nbseq_vrai*nbseq_faux);
    modeles[i].sign    = Fv>Ff?'+':Fv==Ff?'=':'-';
    
#if  OCC
    Pv_ex = (float) m.nboccex_vrai;
    Fv_ex = (float) Pv_ex / nbsymbV;
    Pf_ex = (float) m.nboccex_faux;
    Ff_ex = (float) Pf_ex / nbsymbF;
    Av_ex = (float) nbsymbV-Pv_ex;
    Af_ex = (float) nbsymbF-Pf_ex;
    nP_ex = (float) Pv_ex+Pf_ex;
    nA_ex = (float) Av_ex+Af_ex;
    nbsymbtot   = nbsymbV + nbsymbF;
#if  DEBUG
    fprintf(out,"Pv_ex %f\tFv_ex %f\tPf_ex %f\tFf_ex %f\n",Pv_ex,Fv_ex,Pf_ex,
            Ff_ex);
    fprintf(out,"Av_ex %f\tAf_ex %f\tnP_ex %f\tnA_ex %f\n",Av_ex,Af_ex,nP_ex,
            nA_ex);
    fprintf(out,"Nbsymb vrai %d\tNbsymb faux %d\n",nbsymbV,nbsymbF);
#endif
    modeles[i].khi2_occ    = CARRE(Pv_ex*Af_ex-Pf_ex*Av_ex)
                    / (nP_ex*nA_ex*nbsymbV*nbsymbF) * nbsymbtot;
    modeles[i].sign_occ    = Fv_ex>Ff_ex?'+':Fv_ex==Ff_ex?'=':'-';
#endif
    }

qsort(modeles, nbmodeles, sizeof(Mot),
        (int (*)(const void *, const void *)) comparModeles);


fprintf(out, "Model");
espace(out, maxlongmod-2);
fprintf(out,"#right\t#wrong");
#if  OCC
fprintf(out, "\t#rightT\t#wrongT");
#endif
fprintf(out, "\t\tChi^2");
#if  OCC
fprintf(out, "\t\tChi^2 T");
#endif
fprintf(out, "\n");
fprintf(out, "================================================================================\n");

for(i=0; i!=nbmodeles; i++)
    {
    m   = modeles[i];
    fprintf(out, "%s", m.mot);
    espace(out, maxlongmod+4-strlen(m.mot));
    fprintf(out, "%5d\t%5d", m.nbseq_vrai, m.nbseq_faux);

#if  OCC
    fprintf(out, "\t%5d\t%5d",m.nboccex_vrai, m.nboccex_faux);
#endif

    fprintf(out, "\t\t%3.2f %c", m.khi2,m.sign);
#if OCC
    fprintf(out, "\t\t%3.2f %c",m.khi2_occ,m.sign_occ);
#endif

    fprintf(out,"\n");
    }
}
