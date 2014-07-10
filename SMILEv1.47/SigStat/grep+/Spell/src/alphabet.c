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

#include <alphabet.h>



/******************************************************************************/
/* VARIABLES GLOBALES                                                         */
/******************************************************************************/
Flag    TabSymb[MAXSYMBMOD][MAXSYMBMOD];      /* Table d'equivalences */
/* char    *alphMod[MAXSYMBMOD]={0};      Alphabet des modeles */
Symbole alphSeq[MAXSYMBMOD]={0};       /* Alphabet des sequences */
int     carseq2num[MAXSYMBMOD];        /* Conversion caractere => indice ds alphabet */
char    *nummod2str[MAXSYMBMOD];       /* Conversion num modeles => symboles */
int     nbSymbMod   = 0,        /* Nb de symboles de l'alphabet de modeles*/
        nbSymbSeq   = 0,        /* idem sequences */
        numSAUT     = -1,       /* code utilise pour SAUT ds nummod2str */
        numJOKER    = -1;       /* code utilise pour JOKER ds nummod2str */

enum    {DNA, PROTEINS, UNKNOWN} type;      /* Type de la sequence lue */


/******************************************************************************/
/* initAlphabet:                                                              */
/* initialise les variables globales de la classe.                            */
/******************************************************************************/
void    initAlphabet(void)
{
int     i, j;

for(i=0; i!=MAXSYMBMOD; i++)
    {
    carseq2num[i]   = -1;

    for(j=0; j!=MAXSYMBMOD; j++)
        TabSymb[i][j]   = 0;
    }

alphSeq[0]  = '\0';
}


/******************************************************************************/
/* chargeAlphaSeq:                                                            */
/* chargement de l'alphabet du texte passe en parametre                       */
/******************************************************************************/
Symbole *   chargeAlphaSeq(Symbole **seq, NbSeq  nbseq, char *alphaseq)
{
int     i;
Symbole *j;
char tmp[MAXSYMBMOD]={0}, s;

if(alphaseq!=NULL)
    {
    strcpy(alphSeq, alphaseq);
    nbSymbSeq   = strlen(alphSeq)-1;
    for(i=0; i!=nbSymbSeq; i++)
        carseq2num[(int) *(alphSeq+i)]   = i;
    carseq2num[(int) FINAL] = nbSymbSeq;

    return alphSeq;
    }


/* PARCOURS DU TEXTE pour recherche de l'alphabet utilise.                    */
for(i=0; i!=nbseq; i++)
    {
    for(j=seq[i]; (*j)!=FINAL; j++)
        {
        s   = *j;
        if(s<32 || s>=MAXSYMBMOD)
            {
            fprintf(stderr, ">> Error: Seq %d Pos %d control character %d ('%c') forbidden\n",
                    i, (int)(j-seq[i]), s, s);
            exit(1);
            }
        else if(s == JOKER || s == SAUT)
            {
            fprintf(stderr, ">> Error: Seq %d Pos %d character %d ('%c') forbidden in sequences\n",
                    i, (int)(j-seq[i]), s, s); 
            exit(1);
            }
        else if(!isalnum((int)s))
            fprintf(stderr, "> Warning: Seq %d Pos %d non-alphanumeric character '%c'\n",
                    i, (int)(j-seq[i]), s); 
        else  if(islower((int) s)) /* Mise en majuscules de la sequence */
            {
            s   = (Symbole) toupper((int) s);
            *j  = s;
            }
        
        tmp[(int) s]    = 1;
        }
    }

/* Construction de la chaine a fournir pour construire l'arbre                */
for(i=32, nbSymbSeq  = 0; i!=MAXSYMBMOD; i++)
    {
    if(tmp[i]==1)
        {
        alphSeq[nbSymbSeq]  = i;
        carseq2num[i]       = nbSymbSeq;
        nbSymbSeq++;
        }
    }
carseq2num[(int) FINAL] = nbSymbSeq;
alphSeq[nbSymbSeq]      = FINAL;
alphSeq[nbSymbSeq+1]    = '\0';

printf("** Text alphabet: %s (%d symbols + terminator) **\n",
        alphSeq, nbSymbSeq);
    
return alphSeq;
}



/******************************************************************************/
/* chargeAlphaMod:                                                            */
/* lecture du fichier alphabet et construction de la matrice d'equivalence.   */
/******************************************************************************/
void   chargeAlphaMod(FILE *f)
{
int     i, k;
Symbole *j;
char s, line[512];


/* LECTURE DU FICHIER ALPHABET                                                */
fgets(line, 512, f);

/* Determination du type d'alphabet                                           */
if(strstr(line, "Nucleotide"))
    {
    type    = DNA;
    fprintf(stderr, "** Models alphabet: Nucleotides **\n");
    }
else if(strstr(line, "Protein"))
    {
    type    = PROTEINS;
    fprintf(stderr, "** Models alphabet: Amino acids **\n");
    }
else
    {
    type    = UNKNOWN;
    fprintf(stderr, "** Models alphabet: Unknown type **\n");
    }

/* Lecture des lignes de l'alphabet des modeles                               */
nbSymbMod   = 0;
while(fgets(line, 512, f))
        {
        j   = (Symbole *) line;

        if(*j == '\n')
            continue;

/*         alphMod[nbSymbMod]  = (char *) malloc((strlen(line)+1)*sizeof(char)); */
/*         strcpy(alphMod[nbSymbMod], line); */
        while(*j != '\0' && *j != '\n')
            {
            s   = *j;
            if(s<=32 || s>=MAXSYMBMOD)
                {
                fprintf(stderr, ">> Error: controle character '%c' forbidden in alphabet file\n", s); 
                exit(1);
                }
            else if(s == JOKER)
                {
                if(j!=(Symbole *)line || (*(j+1)!='\0' && *(j+1)!='\n'))
                    {
                    fprintf(stderr, ">> Error: JOKER character '%c' must be alone\n",
                            JOKER); 
                    exit(1);
                    }

                if(numJOKER != -1)
                    fatalError("JOKER defined 2 times in alphabet file\n");

                numJOKER    = nbSymbMod;
                for(k=0; k!=nbSymbSeq; k++)
                    TabSymb[numJOKER][k]   = 1;
                }
            else if(s == FINAL || s == SAUT)
                {
                fprintf(stderr, ">> Error: character %d ('%c') forbidden in alphabet file\n",
                        s, s); 
                exit(1);
                }

            else if(!isalnum((int)s))
                {
                fprintf(stderr, "Warning: non-alphanumeric charactere '%c' in alphabet file\n", s); 
                TabSymb[nbSymbMod][carseq2num[(int)s]]   = 1; 
                }
            else
                {
                if(islower((int) s))
                  *j = s = (Symbole) toupper((int) s);

                TabSymb[nbSymbMod][carseq2num[(int)s]]   = 1; 
                }

            j++;
            }

        if(*j == '\n')
            *j  = '\0';
        
        if(!(nummod2str[nbSymbMod]   = (char *)
                    malloc((strlen(line)+4)*sizeof(char))))
            fatalError("chargeAlphabet: cannot allocate 'nummod2str[i]'\n");
        strcpy(nummod2str[nbSymbMod], line);


        
        nbSymbMod++;
        }

/* Ajout des symboles speciaux dans nummod2str                                */
numSAUT = nbSymbMod;
if(!(nummod2str[numSAUT]   = (char *)
            malloc(2*sizeof(char))))
    fatalError("chargeAlphabet: cannot allocate 'nummod2str[i]'\n");
nummod2str[nbSymbMod][0]    = SAUT;
nummod2str[nbSymbMod][1]    = '\0';


/* Affiche l'alphabet des modeles                                             */
/* for(i=0; i!=nbSymbMod; i++) */
/*     { */
/*     printf("SymbMod %d\t%s\n",i,nummod2str[i]); */
/*     } */

/* Info sur l'alphabet du texte                                               */
for(i=0; i!=nbSymbSeq; i++)
        {
        s   = 0;
        for(k=0; k!=nbSymbMod; k++)
            s   |= (char) TabSymb[k][i];

        if(!s)
            fprintf(stderr,"> Warning: text symbol '%c' isn't recognized by any model's symbol in alphabet file.\n",
                    alphSeq[i]);
        }
/* transAlphMod(); */
}




/******************************************************************************/
/* strshfl - teste si deux chaines sont le shuffling l'une de l'autre         */
/******************************************************************************/
Flag    strshfl(char *  a, char * b)
{
char    * p;

if(strlen(a) != strlen(b))
    return FAUX;

p   = a;
while(*p!='\0')
    {
    if(!strchr(b, *p))
            return FAUX;
    p++;
    }
return VRAI;
}

/******************************************************************************/
/* transAlphMod                                                               */
/******************************************************************************/
/* void    transAlphMod(void) */
/* { */
/* int  i,j; */
/* char tmp[512], joker[2]; */
/*  */
/* joker[0]    = JOKER; */
/* joker[1]    = '\0'; */
/*  */
/* for(i=0; i!=nbSymbMod; i++) */
/*     { */
/*     if(type == DNA) */
/*         { */
/*         if( strshfl(nummod2str[i],"ARN")) */
/*             { */
/*             fprintf(stderr, "Symbole ARN ->> A\n"); */
/*             sprintf(nummod2str[i],"A"); */
/*             } */
/*         else if( strshfl(nummod2str[i],"GRN")) */
/*             { */
/*             fprintf(stderr, "Symbole GRN ->> G\n"); */
/*             sprintf(nummod2str[i],"G"); */
/*             } */
/*         else if( strshfl(nummod2str[i],"CYN")) */
/*             { */
/*             fprintf(stderr, "Symbole CYN ->> C\n"); */
/*             sprintf(nummod2str[i],"C"); */
/*             } */
/*         else if( strshfl(nummod2str[i],"TYN")) */
/*             { */
/*             fprintf(stderr, "Symbole TYN ->> T\n"); */
/*             sprintf(nummod2str[i],"T"); */
/*             } */
/*         else if(strshfl(nummod2str[i],"AG")) */
/*             { */
/*             fprintf(stderr, "Symbole AG ->> R\n"); */
/*             sprintf(nummod2str[i],"R"); */
/*             } */
/*         else if (strshfl(nummod2str[i],"CT")) */
/*             { */
/*             fprintf(stderr, "Symbole CT ->> Y\n"); */
/*             sprintf(nummod2str[i],"Y"); */
/*             } */
/*         else if(nummod2str[i][1]=='\0') */
/*             { */
/*             if(!strcmp(nummod2str[i][0],joker)) */
/*                 { */
/*                 fprintf(stderr, "Symbole '%c' ->> N\n", JOKER); */
/*                 sprintf(nummod2str[i],"N"); */
/*                 } */
/*             } */
/*         else */
/*             { */
/*             strcpy(tmp, nummod2str[i]); */
/*             sprintf(nummod2str[i],"[%s]",tmp); */
/*             sprintf(nummod2str[i],"%s",tmp); */
/*             } */
/*         } */
/*     else if(type == PROTEINS) */
/*         { */
/*         if(!strcmp(nummod2str[i],joker)) */
/*             { */
/*             fprintf(stderr, "Symbole '%c' ->> X\n", JOKER); */
/*             sprintf(nummod2str[i],"X"); */
/*             } */
/*         else */
/*             { */
/*             strcpy(tmp, nummod2str[i]); */
/*             sprintf(nummod2str[i],"[%s]",tmp); */
/*             sprintf(nummod2str[i],"%s",tmp); */
/*             } */
/*         } */
/*     else */
/*         { */
/*         if(strcmp(nummod2str[i],joker)) */
/*             { */
/*             strcpy(tmp, nummod2str[i]); */
/*             sprintf(nummod2str[i],"[%s]",tmp); */
/*             sprintf(nummod2str[i],"%s",tmp); */
/*             } */
/*         } */
/*     } */
/* for(i=0; i!=nbSymbMod; i++) */
/*     { */
/*     for(j=0;j<=nbSymbSeq;j++) */
/*         if(TabSymb[i][j]) */
/*                 printf("Je fais matcher %s avec %c\n",alphMod[i],alphSeq[j]); */
/*     } */
/* for(j=0;j<=nbSymbSeq;j++) */
/*     if(TabSymb[numSAUT][j]) */
/*         printf("Je fais matcher %s avec %c\n",alphMod[numSAUT],alphSeq[j]); */
/* for(j=0;j<=nbSymbSeq;j++) */
/*     if(TabSymb[numJOKER][j]) */
/*         printf("Je fais matcher %s avec %c\n",alphMod[numJOKER],alphSeq[j]); */
/*  */
/* } */
