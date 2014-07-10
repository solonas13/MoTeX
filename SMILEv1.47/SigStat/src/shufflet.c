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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
/* #include <barre.h> */
#include <shufflet.h>


/* global variables */

char wdir[];    /* directory to write to */
int debugflag;  /* 1: write debug file */
FILE *dlog;     /* debug log */
char    alpha[128];

int seqno, sizealpha;
char seqname[NAMELEN];




void usage(void)
{
fprintf(stderr,"Usage: shufflet [OPTIONS] NSEQ ORDER <INFILE >OUTFILE\n");
exit(1);
}


/******************************************************************************/
/* readseq                                                                    */
/******************************************************************************/
int readseq(int k, FILE * infile, int **nver, int ***count, int ***vdeg,
		int **first,int **last, char ***seqstart, int **seqlen, int *maxsizeseq,
        int nbseqalloc, char *alphaseq, char ***origseq)
{
                       /* k= ORDER: length of k-mers to conserve */
char *seq;             /* string to hold sequence */
char errstr[80];       /* string to hold error and warning messages */
int nklets;            /* number of k-let combinations */
int nk1lets;           /* number of (k-1)-let combinations */
int j, maxseqlenalloc;
int letter[128];

*maxsizeseq = 0;
strcpy(alpha, alphaseq);
strcpy(wdir,".");

if (k <= 0)
	usage();
if (k > MAXORDER)
	{
	fprintf(stderr,"Maximum order: %i\n",MAXORDER);
	exit(2);
	}

/* Construction du convertisseur letter                                       */
for(sizealpha=0; sizealpha!=128; sizealpha++)
    letter[sizealpha]   = -1;


for(sizealpha=0; *(alpha+sizealpha)!=FINAL; sizealpha++)
    letter[(int) *(alpha+sizealpha)] = sizealpha;

nklets = (int)pow((double)sizealpha,(double)k);

nk1lets = nklets/sizealpha;
if (nklets > MAXNKLETS)
	Error(2,"Max size of k-let table exceeded (choose smaller k)\n");


seq = (char *) malloc(GRAINLENSEQ*sizeof(char));
assert(seq!=NULL);
maxseqlenalloc  = GRAINLENSEQ;


/* main loop */
seqno=0;
while (((*seqlen)[seqno] =readnextseq(infile,seqname,&seq, &maxseqlenalloc,
                errstr)) >=0)
	/* continue until End Of File (-1) or error (<-1) */
	{
	if ( (*seqlen)[seqno] > *maxsizeseq )
		*maxsizeseq	= (*seqlen)[seqno];

	if ((*seqlen)[seqno] < k) /* too short sequence */
		{     /* print warning and continue with next sequence */
		if ((*seqlen)[seqno] == 0)
			Warning("Empty sequence");
		else
			Warning("Too short sequence");
		continue;  
		}

/* Stockage des sequences originales pour ordre 1                             */
    if(k==1)
        {
        (*origseq)[seqno]   = (char *) malloc(((*seqlen)[seqno]+2)*sizeof(char));
        assert((*origseq)[seqno]!=NULL);
        strncpy((*origseq)[seqno],seq,(*seqlen)[seqno]);
        }

	/* extract alphabet from sequence and convert to sequence of indices */
	indexseq(seq,(*seqlen)[seqno], letter);
	
	(*count)[seqno] = (int *) malloc(sizeof(int)*nklets);
	(*vdeg)[seqno] = (int *) malloc(sizeof(int)*nk1lets);
    assert((*count)[seqno]!=NULL && (*vdeg)[seqno]!=NULL);

	if (k>1)
		{
		(*seqstart)[seqno] = (char *) malloc(sizeof(char)*(k-1));
		assert((*seqstart)[seqno] != NULL);
		for(j=0;j<k-1;j++)
			(*seqstart)[seqno][j] = seq[j];
		}
	
	/* count the k-lets in the indexed sequence */
	if (debugflag)
		fprintf(dlog,"kletcount\n"); 
	kletcount(seq,(*seqlen)[seqno],k,sizealpha,(*count)[seqno]); 
	
	if (debugflag)
		fprintf(dlog,"edgecount\n");
	
	(*first)[seqno] = ind2hash(seq,k-1,sizealpha);  /* first vertex */
	/* last vertex */
	(*last)[seqno] = ind2hash(seq+(*seqlen)[seqno]-k+1,k-1,sizealpha);
	/* (k-1)-lets = #edges in/out */
	(*nver)[seqno] = edgecount(k,sizealpha,(*last)[seqno],(*count)[seqno],
                               (*vdeg)[seqno]);

	if (debugflag)
		{
		kletoutput(dlog,k,sizealpha,(*count)[seqno]); 
		kletoutput(dlog,k-1,sizealpha,(*vdeg)[seqno]);
		} 
	
	seqno++;

/* 	printf(" a traiter: seqno %d\n",seqno); */
    if(seqno == nbseqalloc)
        {
        nbseqalloc  += GRAINSEQ;
        *count   = (int **) realloc(*count, sizeof(int*)*nbseqalloc);
        *seqlen  = (int *) realloc(*seqlen, sizeof(int)*nbseqalloc);
        *vdeg    = (int **) realloc(*vdeg, sizeof(int*)*nbseqalloc);      
        *first   = (int *) realloc(*first, sizeof(int)*nbseqalloc);
        *last    = (int *) realloc(*last, sizeof(int)*nbseqalloc);
        *nver    = (int *) realloc(*nver, sizeof(int)*nbseqalloc);
        assert(*count != NULL);
        assert(*vdeg != NULL);
        assert(*seqlen != NULL);
        assert(*first != NULL);
        assert(*last != NULL);
        assert(*nver != NULL);

        if(k==1)
            {
            *origseq    = (char **) realloc(*origseq, sizeof(char *)*nbseqalloc);
            assert(*origseq != NULL);
            }
        else
            {
            *seqstart = (char **) realloc(*seqstart, sizeof(char *)*nbseqalloc);
            assert(*seqstart != NULL);
            }
        }
	}
	

return(seqno);
}


/******************************************************************************/
/* GenerateSeq                                                                */
/******************************************************************************/
void generateseq(int k,int *nver,int **count, int **vdeg, int *first,
		int *last, char **seqstart, int nbseq, int *seqlen,
		int *count1, int *vdeg1, int nklets, int nk1lets,
		int *lastedge, char **seq)
{
char seedfilename[NAMELEN]; /* file path for random generating seed */
int i,j;



/* initialize random number generator by reading seed from file*/
sprintf(seedfilename,"%s/%s",wdir,SEEDFILENAME);
readseed(seedfilename); 

/* produce shuffled sequences */
/* fprintf(stderr,"** Generating %d shuffled sequences **\n",nbseq); */

/* if(nbseq>=100) */
/*     { */
/*     step    = nbseq/100; */
/*     barre(100); */
/*     } */
/* else */
/*     barre(nbseq); */

if(k==1)
    {
    for(i=0; i!=nbseq; i++)
        {
/*         printf("Je shuffle (%d)...\n%s\nen...\n",seqlen[i], seq[i]); */
        monoshuffle(seqlen[i],seq[i]);
/*         printf("%s\n\n",seq[i]); */
        seq[i][seqlen[i]]   = FINAL;
        seq[i][seqlen[i]+1] = '\0';

/*         if(nbseq>=100) */
/*             { */
/*             if(i%step==0) */
/*                 barre(0); */
/*             } */
/*         else */
/*             barre(0); */
        }
    }
else
    for (i=0; i!=nbseq; i++)
        {
        memcpy(seq[i],seqstart[i],k-1);
    /* 		{ */
    /* 		for (j=0; j<k-1; j++) */
    /* 			seq[i][j] = seqstart[i][j]; */
    /* 		} */
        
        
        /* make working copies of count and vdeg */
        for (j = 0; j < nklets; j++)
            count1[j] = count[i][j];
        for (j = 0; j < nk1lets; j++)
            vdeg1[j] = vdeg[i][j];
        


        shuffle(sizealpha,k,nver[i],count1,vdeg1,first[i],last[i], lastedge,seq[i]);


        seq[i][seqlen[i]]   = FINAL;
        seq[i][seqlen[i]+1] = '\0';
        
/*         printf("SEQUENCE %d\n",i); */
        if (!kletverify(seq[i],seqlen[i],k,sizealpha,count[i],count1))
            Error(99,"k-let count mismatch in shuffled seq!!!");

        seqoutput(seq[i],seqlen[i],LINE);

/*         if(nbseq>=100) */
/*             { */
/*             if(i%step==0) */
/*                 barre(0); */
/*             } */
/*         else */
/*             barre(0); */
        }

writeseed(seedfilename);  /* store seed for random number generator */
}



/******************************************************/
/* error: print error message and quit                */
/* input: message: message text                       */

void Error(int code, char message[])
{
  fprintf(stderr,"\nError: sequence #%i \"%s\"",seqno,seqname);
/*   if (seqlen > 0) */
/*     fprintf(stderr," (length %i, alphabet size %i)",seqlen,m); */
  fprintf(stderr,":\n");
  fprintf(stderr,"%s\n",message);
/*   exit(code); */
}

 
/******************************************************/
/* warning: print warning message                     */
/* input: message: message text                       */
 
void Warning(char message[])
{
  fprintf(stderr,"\nWarning: sequence #%i \"%s\"",seqno,seqname);
/*   if (seqlen > 0) */
/*     fprintf(stderr," (length %i, alphabet size %i)",seqlen,m); */
  fprintf(stderr,":\n");
  fprintf(stderr,"%s\n",message);
}

/******************************************************/
/* generatename: generate name for shuffled sequence  */
/* input: seqname: name of original sequence          */
/*              n: number of shuffling                */
/* output: shfseqname: name of shuffled sequence      */

void generatename(char seqname[], int n, int nseq, char shfseqname[])
{
if (seqname[0] == '\0' || strcmp(seqname,NONAME) == 0)
	{
	if (nseq > 1)
		sprintf(shfseqname,"SHF%i",n);
	else
		shfseqname[0] = '\0'; /* 1 shuffling of unnamed sequence: no name */
	}
else
	sprintf(shfseqname,"%s_SHF%i",seqname,n);
}


