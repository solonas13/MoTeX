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


#ifndef SHUFFLET_H
#define SHUFFLET_H

#include <stdlib.h>
#include <assert.h>

#define FINAL   '$'
#define MAXORDER 10      /* max k-let order (length) */
#define MAXALPHA 27      /* max alphabet size */
#define MAXNKLETS 8000   /* max number of k-let combinations */
#define NAMELEN 64       /* max length of i.e. filenames and sequence names */
#define LINE 60          /* line length for sequence output */
#define NONAME "(unnamed)" /* assigned to unnamed sequences */
#define DLOGNAME "dlog"  /* name of debug log file */
#define ERRFILENAME "errout" /* file name for redirecting stderr */
#define SEEDFILENAME "seed"  /* name of seed file */

#define GRAINSEQ	500
#define GRAINLENSEQ 2000   /* max sequence length */


#define randomint(N) ((int)((N)*unif01())+1)  /* random integer 1..N */

extern	char wdir[NAMELEN];        /* directory to write to */
extern	int debugflag;             /* 1: write debug file */
extern	FILE *dlog;                /* debug log */
extern  char    alpha[128];      /* Alphabet des sequences */
        /* Alphabet des sequences */
/* shufflet.c */
void generatename(char seqname[], int n, int nseq, char shfseqname[]);
int  readseq(int k, FILE * infile, int **nver, int ***count,
	int ***vdeg, int **first, int **last, char ***seqstart, int **seqlen,
	int *maxsizeseq, int maxseqlenalloc, char *alphaseq, char ***origseq);
void generateseq(int k,int *nver,int **count, int **vdeg, int *first,
	int *last, char **seqstart, int nbseq, int *seqlen,
	int *count1, int *vdeg1, int nklets, int nk1lets, int *lastedge,
	char **seq);

void Error(int code, char message[]);
void Warning(char message[]);

/* euler.c */
void indexseq(signed char seq[], int seqlen, int letter[128]);
void kletcount(signed char seq[], int seqlen, int k, int m, int count[]);
int kletverify(signed char seq[], int seqlen, int k, int m, int count0[],
        int count1[]);
int edgecount(int k, int m, int last, int count[], int vdeg[]);
void kletoutput(FILE *fp, int k, int m, int count[]);
char *hash2str(int hash, int k, int m, char klet[]);
int ind2hash(char klet[], int k, int m);
void shuffle(int m, int k, int nver, int count[], int vdeg[], int first, 
	     int last, int lastedge[], char seq[]);
void monoshuffle(int seqlen, char seq[]);
void arborescence(int m, int nk1lets, int nver, int count[], int vdeg[], 
		  int first, int root, int branch[]);
void randomtrail(int m, int k, int count[], int vdeg[], int first,
		int lastedge[], char seq[]);


/* seqio.c */
int readnextseq(FILE *fp, char seqname[21], char **seq, int * maxseqlenalloc,
		char errmsg[80]);
void seqoutput(char seq[], int seqlen, int linelen);

/* random.c */
double unif01(void);
unsigned int readseed(char filename[]);
int writeseed(char filename[]);

#endif
