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

#include <io.h>

extern  int     numSAUT, numJOKER, nbSymbMod;

/******************************************************************************/
/* lectureFichierRes                                                          */
/******************************************************************************/
Mot* lectureFichierRes(FILE *res, int nbseq, int *nbmodeles, int *maxlongmod)
{
char    car, line[BUF], ret, modele[BUF], codes[BUF], *ptr;
signed char tmp;
int     i, nbligneslues, j, l;
Mot     *modeles;

/** Recherche du nb de modeles (en fin de fichier) **/
fseek(res, 0, SEEK_END);

i=0; /* Je cherche le 2eme ':' en partant de la fin (premier = CPU) */
do
    {
    fseek(res, -2, SEEK_CUR);
    fscanf(res,"%c",&car);
    if(car == ':')
        i++;
    }
while(i!=2);

fscanf(res,"%d",nbmodeles);
/* printf("Nb modeles = %d\n",*nbmodeles); */

if ( *nbmodeles <= 0 )
    return NULL;

/* Positionnement premier modele                                              */
fseek(res, 320, SEEK_SET);

modeles=(Mot *)calloc(*nbmodeles,sizeof(Mot));
if(modeles==NULL)
    fatalError("lectureFichierRes: cannot allocate 'modeles'");


/* LECTURE DES MODELES DANS LE FICHIER RESULTAT */
i=0;
fgets(line, BUF, res);
/* printf("J'ai lu la ligne : %s\n",line); */
do
    {
#if !OCC
    ret = sscanf(line, "%s %s %d", modele, codes, &(modeles[i].nbseq_vrai));
#else
    ret = sscanf(line, "%s %s %d %d\n", modele, codes, &(modeles[i].nbseq_vrai),
            &(modeles[i].nboccex_vrai));
#endif

/* Stockage du modele alphabetiques                                           */
    l   = strlen(modele);
    if(l>*maxlongmod)
        *maxlongmod = l;
    if(!(modeles[i].mot=(char*)malloc((l+1)*sizeof(char))))
        fatalError("lectureFichierRes: cannot allocate 'modeles[i]'");
    strcpy(modeles[i].mot, modele);

/* Lecture et stockage du modele numerique                                    */
    if(!(modeles[i].codes=(char*)malloc((strlen(codes)+1)*sizeof(char))))
        fatalError("lectureFichierRes: cannot allocate 'modeles[i].codes'");

    j=0;
    ptr = codes;

    while(*ptr != '\0')
        {
        if(*ptr == JOKERinterne)
            modeles[i].codes[j] = numJOKER;
        else if (*ptr == SAUTinterne)
            modeles[i].codes[j] = numSAUT;
        else
            {
            tmp = *ptr - SHIFTALPHA;
            if(tmp < 0 || tmp >= nbSymbMod)
                {
                fprintf(stderr,">> Error: Model '%s' corrupted (unknown symbol at position %d)\n",
                    modele, j+1);
                exit(1);
                }
            modeles[i].codes[j] = tmp;
            }
        
        j++;
        ptr++;
        }
    modeles[i].codes[j]   = -1;


    modeles[i].quorum_reel = (float) (modeles[i].nbseq_vrai)*100.0/nbseq;
/*  printf("pour le modele %s j'ai lu %d => %f\n",modeles[i].mot,modeles[i].nbseq_vrai,modeles[i].quorum_reel); */

    if (ret == 4)
        fgets(line, BUF, res);
    else
        {
        nbligneslues    = 0;
        do
            {
            fgets(line, BUF, res);
            nbligneslues++;
            }
        while(!strncmp(line, "Seq", 3));
/* printf("J'ai lu la ligne : %s\n",line); */

#if  OCC
        if(nbligneslues==1)
            {
            fprintf(stderr,
                "Error: output file contains no occurrences number,\n");
            fprintf(stderr,
                "    altough statistics on total number of occurrences have been requested.\n");
            return (NULL);
            }
        else
            sscanf(line, "%d\n", &(modeles[i].nboccex_vrai));
#endif
        if(nbligneslues != 1)
            fgets(line, BUF, res);
        }

/* printf("J'ai lu : %s %d et %d\n",modeles[i].mot,modeles[i].nbseq_vrai,modeles[i].nboccex_vrai); */
/*     printf("J'ai lu la ligne: '%s'\n",line); */
    i++;
    }
while( strncmp(line, "Nb models", 9) );

return(modeles);
}

/******************************************************************************/
/* openFile                                                                   */
/******************************************************************************/
FILE    *openFile(char *nom, char *mode)
{
FILE *f;

f   = fopen(nom,mode);

if(f == NULL)
    {
    fprintf(stderr,"Impossible to open file '%s'.\n",nom);
    exit(1);
    }

return f;
}

/******************************************************************************/
/* PrintCpuTime                                                               */
/******************************************************************************/
void printCpuTime(FILE *f)
{
float         ust;
struct tms    tms;
static float  dust;

times(&tms);

ust =  (float) tms.tms_utime;

if (f==NULL)
    {
    dust = ust;
    }
else
    {
    ust -= dust;
    printf("User time : %.2f sec.\n", ust / sysconf(_SC_CLK_TCK));
    fprintf(f,"User time : %.2f sec.\n", ust / sysconf(_SC_CLK_TCK));
    }
}
