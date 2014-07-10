/*
 *  Copyright (c) Atelier de BioInformatique
 *  
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2 of the License, or (at your option) any later version.
 *  
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA 02111-1307, USA
 *  
 *  For questions, suggestions, bug-reports, enhancement-requests etc.
 *  I may be contacted at: Alain.Viari@inrialpes.fr
 */

/* #ifdef THINK_C */
#include <ctype.h>
/* #else */
/* #include <sys/types.h> */
/* #endif */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Gtypes.h"
#include "libsysk.h"
#include "libfasta.h"

#define CHECK   0
#define DEBUG_FASTA	0

#define READ_NEXT Faux
#define PUSH_BACK Vrai

#define SERIAL    Vrai
#define INDEXED   Faux

#ifdef THINK_C
#define LINE_FEED '\r'
#else
#define LINE_FEED '\n'
#endif

/* -------------------------------------------- */
/* @static: lecture bufferisee			*/
/* -------------------------------------------- */
static char * sNextIOBuffer(FILE *streamin, Bool retain, Bool serial)
{
/* 	Int32  lenbuf; */
	char   *buf, *end;
	
	static char sBuffer[BUFSIZ];	/* in <stdio.h>	  */
	static Bool sRetained = Faux;
	
	buf = (((retain || sRetained) && serial) 
		? sBuffer 
		: fgets(sBuffer, sizeof(sBuffer), streamin));	 	

	if (buf) {
	   end = buf + strlen(buf) - 1;
	   if (*end == LINE_FEED) *end = '\000';
	}

	sRetained = retain;
	
	return buf;
}

/* -------------------------------------------- */
/* compte le nombre de caracteres alpha	dans	*/
/* un buffer					*/
/* -------------------------------------------- */
Int32 CountAlpha(char *buf)
{
	Int32 count;
	
	for (count = 0 ; *buf ; buf++)
	    if (isalpha((int)*buf))
		count++;
	
	return count;
}


/* -------------------------------------------- */
/* copy only alpha chars from s2 to s1		*/
/* -------------------------------------------- */
char * StrcpyAlpha(char *s1, char *s2)
{
	for( ; *s2 ; s2++)
	    if (isalpha((int)*s2))
	    	*s1++ = *s2;

	*s1 = '\000';

	return s1;
}

/* -------------------------------------------- */
/* skip to next space in buffer			*/
/* -------------------------------------------- */
char * NextSpace(char *buffer)
{
	for (; *buffer ; buffer++)
	   if (isspace((int)*buffer))
	   	return buffer;
	
	return NULL;
}

/* -------------------------------------------- */
/* returns sequence name (FASTA)		*/
/* -------------------------------------------- */
char *GetFastaName(char *buffer)
{
	static char name[FASTA_NAMLEN];

	buffer[FASTA_NAMLEN] = '\000';
	
	if (sscanf(buffer + 1, "%s", name) != 1)
	    strcpy(name, "<no Name>");

	return name;
}

/* -------------------------------------------- */
/* returns sequence comment (FASTA)		*/
/* -------------------------------------------- */
char *GetFastaComment(char *buffer)
{
	char   *space;
	static char comment[FASTA_COMLEN];

	buffer[FASTA_COMLEN] = '\000';
	
	space = NextSpace(buffer);
	
	strcpy(comment, (space ? space + 1 : "<no comment>"));

	return comment;
}

/* -------------------------------------------- */
/* liberation d'une sequence			*/
/* -------------------------------------------- */
FastaSequencePtr FreeFastaSequence(FastaSequencePtr seq)
{
	if (seq) {
	    if (seq->seq)  FREE(seq->seq);
	    FREE(seq);
	}

	return NULL;
}
	
/* -------------------------------------------- */
/* allocation d'une sequence			*/
/* -------------------------------------------- */
FastaSequencePtr NewFastaSequence(void)
{
	FastaSequencePtr seq;
	
	if (! (seq = NEW(FastaSequence)))
	    return NULL;
	   
	seq->length   = 0;

	if (! (seq->seq = NEWN(char,  BUFSIZ)))
	    return FreeFastaSequence(seq);

	seq->bufsize = BUFSIZ;

	*(seq->name)    = '\000';
	*(seq->comment) = '\000';

	seq->ok = Vrai;
	
	return seq;
}

/* -------------------------------------------- */
/* lecture/redimensionnement d'une sequence au	*/
/* format Fasta	Lecture en serie		*/
/* returns : Faux -> last sequence		*/
/*	     Vrai -> more to read		*/
/*           <but> you must check seq->ok !	*/
/* -------------------------------------------- */
Bool ReadFastaSequence(FILE *streamin, FastaSequencePtr seq)
{
	Int32	readlen, buflen;
	char 	*buffer, *tbuf;

	seq->ok = Faux;				/* assume error		*/

	buflen = seq->length = 0; 
	
	seq->offset = ftell(streamin);

	buffer = sNextIOBuffer(streamin, READ_NEXT, SERIAL);

	if (! (buffer && (*buffer == '>'))) 	/* sync error		*/
	    return Faux;			/* last sequence	*/
	
	if (seq->offset)
	    seq->offset -= (strlen(buffer) + 1);

	strcpy(seq->name,    GetFastaName(buffer));
	
	strcpy(seq->comment, GetFastaComment(buffer));
	
	while ((buffer = sNextIOBuffer(streamin, READ_NEXT, SERIAL))) {

	    if (*buffer == '>') {
	    	(void) sNextIOBuffer(streamin, PUSH_BACK, SERIAL); /* push it back */
		break;
	    }

#if CHECK	    
	    readlen = CountAlpha(buffer);
#else
	    readlen = strlen(buffer);
#endif
	
	    buflen +=  readlen;
   
	    if (buflen >= seq->bufsize) {
	    
	    	if (! (tbuf = REALLOC(char, seq->seq, 2 * buflen + 1)))
	    	   return Vrai;			/* but seq->ok is Faux	*/

		seq->seq = tbuf;
		
		seq->bufsize = 2 * buflen + 1;
		
	    }		
#if CHECK
	    StrcpyAlpha(seq->seq + seq->length, buffer);
#else
	    memcpy(seq->seq + seq->length, buffer, readlen);
#endif
	
	    seq->length = buflen;
	
	}

	seq->seq[seq->length] = '\000';

	return (seq->ok = Vrai);
}

/* -------------------------------------------- */
/* lecture/redimensionnement d'une sequence au	*/
/* format Fasta	Lecture indexee			*/
/* returns : Faux -> last sequence		*/
/*	     Vrai -> more to read		*/
/*           <but> you must check seq->ok !	*/
/* -------------------------------------------- */
Bool GetFastaSequence(FILE *streamin, FastaSequencePtr seq)
{
	Int32	readlen, buflen;
	char 	*buffer, *tbuf;

	seq->ok = Faux;				/* assume error		*/

	buflen = seq->length = 0; 
	
	fseek(streamin, seq->offset, SEEK_SET);

	buffer = sNextIOBuffer(streamin, READ_NEXT, INDEXED);

	if (! (buffer && (*buffer == '>'))) 	/* sync error		*/
	    return Faux;			/* last sequence	*/
	
	if (seq->offset)
	    seq->offset -= (strlen(buffer) + 1);

	strcpy(seq->name,    GetFastaName(buffer));
	
	strcpy(seq->comment, GetFastaComment(buffer));
	
	while ((buffer = sNextIOBuffer(streamin, READ_NEXT, INDEXED))) {

	    if (*buffer == '>')
		break;

#if CHECK	    
	    readlen = CountAlpha(buffer);
#else
	    readlen = strlen(buffer);
#endif
	
	    buflen +=  readlen;
   
	    if (buflen >= seq->bufsize) {
	    
	    	if (! (tbuf = REALLOC(char, seq->seq, 2 * buflen + 1)))
	    	   return Vrai;			/* but seq->ok is Faux	*/

		seq->seq = tbuf;
		
		seq->bufsize = 2 * buflen + 1;
		
	    }		
#if CHECK
	    StrcpyAlpha(seq->seq + seq->length, buffer);
#else
	    memcpy(seq->seq + seq->length, buffer, readlen);
#endif
	
	    seq->length = buflen;
	
	}

	seq->seq[seq->length] = '\000';

	return (seq->ok = Vrai);
}

/* -------------------------------------------- */
/* ecriture d'une sequence au format Fasta	*/
/* -------------------------------------------- */
void WriteFastaSequence(FILE *streamou, FastaSequencePtr seq, Int32 char_per_line)
{
	Int32 i, nlines, rest;
	char *buf, *end, tempo;

	fputc('>', streamou);
	fputs((*(seq->name)    ? seq->name    : "<no name>")   , streamou);
	fputc(' ', streamou);
	fputs((*(seq->comment) ? seq->comment : "<no comment>"), streamou);
	fputc(LINE_FEED, streamou);

	nlines = seq->length / char_per_line;

	buf = seq->seq;

	for (i = 0 ; i < nlines ; i++) {
	    end = buf + char_per_line;
	    tempo = *end;
	    *end = '\000';
	    fputs(buf, streamou);
	    fputc(LINE_FEED , streamou);
	    *end = tempo;
	    buf += char_per_line;
	}

	if ((rest = (seq->length % char_per_line))) {
	   end = buf + rest;
	   tempo = *end;
	   *end = '\000';
	   fputs(buf, streamou);
	   fputc(LINE_FEED , streamou);
	   *end = tempo;
	}
}

