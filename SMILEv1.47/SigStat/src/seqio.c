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


/* seqio.c: Reading and writing sequences on fasta or plain ascii format */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <shufflet.h>

/**********************************************************/
/* readnextseq: read next sequence from file              */
/* Fasta format, or plain ascii (single sequence)         */
/* Skip first line if it starts with >                    */
/* Continue until *, /, EOF, or another >                 */
/* Convert all letters to uppercase, ignore all other     */
/* characters.                                            */
/* input: fp: file pointer                                */
/* output: seqname: sequence name (first 20 characters)   */
/*         seq: sequence string                           */
/*         errmsg: string containing error message        */
/* return value: sequence length if sequence found (>=0)  */
/*               -1: End Of File, <-1: error occurred     */

int readnextseq(FILE *fp, char seqname[21], char **seq,
        int *maxseqlenalloc, char errmsg[80])
{
  int c, pos = 0;

  (*seq)[0] = '\0';
  sprintf(seqname,NONAME);

  c = fgetc(fp); /* first character determine further processing */
/*   printf("first character: %c\n",c); */
  if (c == EOF)  /* already at End Of File */
    return -1;
  else if (c == '>')
    { /* fasta format */
      c = fgetc(fp);
      ungetc(c,fp);
      if(c!='\n' && c!=EOF)
          fscanf(fp,"%20s",seqname); /* a name must follow the > */
/* 	  printf("nameseq: %s\n",seqname); */
      while ((c = fgetc(fp)) != '\n' && c != EOF); /* skip rest of line */
    }
  else /* no initial > : plain ascii format */
    ungetc(c,fp);

  /* main loop for reading sequence */
  while ((c = fgetc(fp)) != '>' && c != '/' && c != '*' && c != EOF)
    if (isalpha(c)) /* ignore nonalphabetic characters */
      {
	if (pos == *maxseqlenalloc)
	  {
        *maxseqlenalloc += GRAINLENSEQ;
        *seq = (char *) realloc(*seq, sizeof(char)*(*maxseqlenalloc));
        assert(*seq != NULL);
	  }

	(*seq)[pos++] = toupper(c);  /* all letters converted to upper case */
    }

  if (c == '>')
    ungetc(c,fp);  /* put > back for reading next sequence */
  (*seq)[pos] = '\0'; /* nice end of string */
  if (pos == 0 && c == EOF)
	{
/* 	printf("Fini...\n"); */
    return -1; /* still empty sequence at End Of File */
	}
  else
    return pos;
}


/******************************************************/
/* seqoutput: write formatted sequence to file        */
/* NB! seq is a sequence of indices, not letters!     */
/* input: fp: output file pointer                     */
/*       seq: sequence of indices                     */
/*    seqlen: sequence length                         */
/*   linelen: number of characters per line           */

void seqoutput(char seq[], int seqlen, int linelen)
{
  int i;
  
  for (i = 0; i < seqlen; i++)
      seq[i]	=	alpha[(int)seq[i]];
}



