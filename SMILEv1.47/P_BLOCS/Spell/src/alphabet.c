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

void chercheComp(char *s1, char *s2);


/******************************************************************************/
/* VARIABLES GLOBALES                                                         */
/******************************************************************************/
Flag    TabSymb[MAXSYMBMOD][MAXSYMBMOD];      /* Table d'equivalences */
char    *alphMod[MAXSYMBMOD];    /* Alphabet original des modeles */
Symbole alphSeq[MAXSYMBMOD]={0}; /* Alphabet des sequences */
int     carseq2num[MAXSYMBMOD];  /* Conversion caractere => indice ds alphabet*/
char    *nummod2str[MAXSYMBMOD]; /* Conversion num modeles => symboles */
int     comp[MAXSYMBMOD];       /* Conversion complément pour palindrome */
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
/* chargeAlphabet:                                                            */
/* lecture du fichier alphabet et construction de la matrice d'equivalence.   */
/* PREND le fichier alphabet, le texte a traiter.                             */
/* RENVOIE un ptr sur l'alphabet des sequences a traiter.                     */
/******************************************************************************/
Symbole *   chargeAlphabet(FILE *f, Symbole **seq, NbSeq  nbseq)
{
int     i, k;
Symbole *j;
char tmp[MAXSYMBMOD]={0}, s, line[512];




/* PARCOURS DU TEXTE pour recherche de l'alphabet utilise.                    */
for(i=0; i!=nbseq; i++)
    {
    for(j=seq[i]; (*j)!=FINAL; j++)
        {
        s   = *j;
        if(s<32 || s>=MAXSYMBMOD)
            {
            fprintf(stderr, ">> Error: Seq %d Pos %d Control character %d ('%c')\n",
                    i, (int)(j-seq[i]), s, s);
            exit(1);
            }
        else if(s == JOKER || s == SAUT)
            {
            fprintf(stderr, ">> Error: Seq %d Pos %d Forbidden character %d ('%c') in the sequences\n",
                    i, (int)(j-seq[i]), s, s); 
            exit(1);
            }
        else if(!isalnum((int)s))
            fprintf(stderr, "> Warning: Seq %d Pos %d Non alphanumeric character '%c'\n",
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

printf("\n** Text alphabet: %s (%d symbols + terminator) **\n",
        alphSeq, nbSymbSeq);



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
    fprintf(stderr, "** Models alphabet: Amino-acids **\n");
    }
else
    {
    type    = UNKNOWN;
    fprintf(stderr, "** Models alphabet: unknown **\n");
    }

/* Lecture des lignes de l'alphabet des modeles                               */
nbSymbMod   = 0;
while(fgets(line, 512, f))
        {
        j   = (Symbole *) line;

        if(*j == '\n')
            continue;

        alphMod[nbSymbMod]  = (char *) malloc((strlen(line)+1)*sizeof(char));
        if(alphMod[nbSymbMod]==NULL)
            fatalError("alphabet.c: chargeAlphabet: cannot allocate 'alphMod[i]'\n");
        strcpy(alphMod[nbSymbMod], line);
        while(*j != '\0' && *j != '\n')
            {
            s   = *j;
            if(s<=32 || s>=MAXSYMBMOD)
                {
                fprintf(stderr, ">> Error: Control character '%c' in the alphabet file\n", s); 
                exit(1);
                }
            else if(s == JOKER)
                {
                if(j!=(Symbole *)line || (*(j+1)!='\0' && *(j+1)!='\n'))
                    {
                    fprintf(stderr, ">> Error: JOKER character '%c' not alone\n",
                            JOKER); 
                    exit(1);
                    }

                if(numJOKER != -1)
                    fatalError("JOKER defined 2 times in the alphabet file\n");

                numJOKER    = nbSymbMod;
                for(k=0; k!=nbSymbSeq; k++)
                    TabSymb[numJOKER][k]   = 1;
                }
            else if(s == FINAL || s == SAUT)
                {
                fprintf(stderr, ">> Error: Forbidden character %d ('%c') in the alphabet file\n",
                        s, s); 
                exit(1);
                }

            else if(!isalnum((int)s))
                {
                fprintf(stderr, "Warning: Non alphanumeric character '%c' in the alphabet file\n", s); 
                TabSymb[nbSymbMod][carseq2num[(int)s]]   = 1; 
                }
            else
                {
                if(islower((int) s))
                   *j = s = (Symbole) toupper((int) s);

                TabSymb[nbSymbMod][carseq2num[(int)s]]   = 1; 
/*                 printf("Je fais matcher Mod %s %d et symbole %c %d\n",line, nbSymbMod,s,s); */
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
fprintf(stderr, "Models alphabet's size: %d\n",nbSymbMod);

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
            fprintf(stderr,"> Warning: text symbol '%c' is not recognized by any models symbol.\n",
                    alphSeq[i]);
        }

return alphSeq;
}


/******************************************************************************/
/* estSymbMod                                                                 */
/******************************************************************************/
int    str2nummod(char *str)
{
int i;

for(i=0; i!=nbSymbMod; i++)
    if(!strcmp(nummod2str[i], str))
        return i;

return -1;
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
void    transAlphMod(Flag pal)
{
int  i,j;
char tmp[512];

for(i=0; i!=nbSymbMod; i++)
    {
    if(type == DNA)
        {
        if(strshfl(nummod2str[i],"AR") || strshfl(nummod2str[i],"AN") || strshfl(nummod2str[i],"ARN"))
            {
            fprintf(stderr, "Symbol %s ->> A\n", nummod2str[i]);
            sprintf(nummod2str[i],"A");
            }
        else if( strshfl(nummod2str[i],"CY") || strshfl(nummod2str[i],"CN") || strshfl(nummod2str[i],"CYN"))
            {
            fprintf(stderr, "Symbol %s ->> C\n",nummod2str[i]);
            sprintf(nummod2str[i],"C");
            }
        else if( strshfl(nummod2str[i],"GR") || strshfl(nummod2str[i],"GN") || strshfl(nummod2str[i],"GRN"))
            {
            fprintf(stderr, "Symbol %s ->> G\n", nummod2str[i]);
            sprintf(nummod2str[i],"G");
            }
        else if( strshfl(nummod2str[i],"TY") || strshfl(nummod2str[i],"TN") || strshfl(nummod2str[i],"TYN"))
            {
            fprintf(stderr, "Symbol %s ->> T\n", nummod2str[i]);
            sprintf(nummod2str[i],"T");
            }
        else if(strshfl(nummod2str[i],"AG") || strshfl(nummod2str[i],"AGR") || strshfl(nummod2str[i],"AGN") || strshfl(nummod2str[i],"AGRN")) 
            {
            fprintf(stderr, "Symbol %s ->> R\n",nummod2str[i]);
            sprintf(nummod2str[i],"R");
            }
        else if (strshfl(nummod2str[i],"CT") || strshfl(nummod2str[i],"CTY") || strshfl(nummod2str[i],"CTN") || strshfl(nummod2str[i],"CTYN"))
            {
            fprintf(stderr, "Symbol %s ->> Y\n", nummod2str[i]);
            sprintf(nummod2str[i],"Y");
            }
        else if(nummod2str[i][1]=='\0')
            {
            if(nummod2str[i][0]==JOKER)
                {
                fprintf(stderr, "Symbol '%c' ->> N\n", JOKER);
                sprintf(nummod2str[i],"N");
                }
            }
        else
            {
            strcpy(tmp, nummod2str[i]);
            sprintf(nummod2str[i],"[%s]",tmp);
            }
        }
    else if(type == PROTEINS)
        {
        if(nummod2str[i][1]=='\0')
            {
            if(nummod2str[i][0]==JOKER)
                {
                fprintf(stderr, "Symbol '%c' ->> X\n", JOKER);
                sprintf(nummod2str[i],"X");
                }
            }
        else
            {
            strcpy(tmp, nummod2str[i]);
            sprintf(nummod2str[i],"[%s]",tmp);
            }
        }
    else
        {
        if(nummod2str[i][1]=='\0')
            {
            if(nummod2str[i][0]==JOKER)
                {
                strcpy(tmp, nummod2str[i]);
                sprintf(nummod2str[i],"[%s]",tmp);
                }
            }
        }

/* Verification de collision de symboles                                      */
    for(j=0; j!=i; j++)
        if(!strcmp(nummod2str[j],nummod2str[i]))
            {
            fprintf(stderr,">> Error: possible confusion between symbols %s and %s of the alphabet.\nModify the alphabet to avoid conflict.\n",
                    alphMod[i],alphMod[j]);
            exit(1);
            }
    }

/* Positionnement des complémentaires                                         */
if (pal)
    {
    if(type != DNA)
        {
        fprintf(stderr,">> Error: palindroms can only be used with a nucleotide alphabet\n");
        exit(1);
        }

    for(i=0; i!=nbSymbMod; i++)
        comp[i] = -1;

        
    chercheComp("A", "T");
    chercheComp("C", "G");
    chercheComp("R", "Y");
/* ...and so on?                                                              */

    
/* Complementaire du joker                                                    */
    if(numJOKER!=-1)
        comp[numJOKER] = numJOKER;


/* Verification que tous les symboles de l'alphabet ont bien un palindrome    */
    for(i=0; i!=nbSymbMod; i++)
        {
        if(comp[i] == -1)
            {
            fprintf(stderr,">> Error: some symbols of the models alphabet misses their complemtentary symbol, cannot use the palindromic option\n");
            exit(1);
            }
        }
    }
}
    

void chercheComp(char *s1, char *s2)
{
int i,j;

/* Recherche de s1 et s2 dans alph modeles                                    */
    for(i=0; i!=nbSymbMod; i++)
            if (!strcmp(nummod2str[i],s1))
                    break;

    for(j=0; j!=nbSymbMod; j++)
            if (!strcmp(nummod2str[j],s2))
                    break;

/* s1 et s2 sont presents dans l'alphabet => ils sont complementaires         */
    if(i!=nbSymbMod && j!=nbSymbMod)
        {
        comp[i] = j;
        comp[j] = i;
        }
/* Si un seul d'entre eux est present => probleme                             */
    else if((i==nbSymbMod && j!=nbSymbMod) || (i!=nbSymbMod && j==nbSymbMod))
        {
        fprintf(stderr,">> Error: the models alphabet misses some symbols to be used with the palindromic option ('%s' and '%s' must appear together)\n", s1,s2);
        exit(1);
        }
/* Et si aucun... pas probleme                                                */
}
