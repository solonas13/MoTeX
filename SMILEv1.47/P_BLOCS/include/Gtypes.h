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

#define _H_Gtypes

#ifndef NULL
#include <stdio.h>			/* is the official NULL	here ?	*/
#endif

/* ==================================================== */
/* constantes						*/
/* ==================================================== */

#ifndef PROTO
#define PROTO	1			/* prototypes flag		*/
#endif

#ifdef THINK_C
#define Vrai	true			/* TC boolean values		*/
#define Faux	false			/*				*/
#else
#define Vrai	0x1			/* bool values	= TRUE		*/
#define Faux	0x0			/*              = FALSE		*/
#endif

#define Nil	NULL		 	/* nil pointer			*/

#define	kBigInt16	0x7fff	    	/* plus grand 16 bits signe	*/
#define	kBigInt32	0x7fffffff  	/* plus grand 32 bits signe	*/
#define	kBigUInt16	0xffff	    	/* plus grand 16 bits ~signe	*/
#define	kBigUInt32	0xffffffff  	/* plus grand 32 bits ~signe	*/

#define kBitsPerLong	32		/*  long = 32 bits		*/
#define kMaxShftLong	31		/*  BitsPerLong - 1 max shift	*/
#define kLog2BitLong	5		/*  =log2(BitsPerLong)		*/

#ifdef THINK_C
/* ==================================================== */
/*  Types (for Macintosh ThinK C)			*/
/* ==================================================== */

typedef long		Long;		/* plus grand mot signe		*/
typedef unsigned long	ULong;		/* plus grand mot signe		*/
typedef long		Int32;		/* Int32  = 32 bits signe	*/
typedef unsigned long	UInt32;		/* UInt32 = 32 bits ~signe	*/
typedef short		Int16;		/* Int16  = 16 bits signe	*/
typedef unsigned short	UInt16;		/* UInt32 = 16 bits ~signe	*/
typedef char		Int8;		/* Int8   = 8 bits signe	*/
typedef unsigned char	UInt8;		/* UInt8  = 8 bits ~signe	*/

typedef Boolean		Bool;		/* booleen  			*/

#else
/* ==================================================== */
/*  Types (for Sun & Iris)				*/
/* ==================================================== */

typedef long		Long;		/* plus grand mot signe		*/
typedef unsigned long	ULong;		/* plus grand mot signe		*/
typedef int	        Int32;		/* Int32  = 32 bits signe	*/
typedef unsigned int    UInt32;		/* UInt32 = 32 bits ~signe	*/
typedef short		Int16;		/* Int16  = 16 bits signe	*/
typedef unsigned short	UInt16;		/* UInt32 = 16 bits ~signe	*/
typedef char		Int8;		/* Int8   = 8 bits signe	*/
typedef unsigned char	UInt8;		/* UInt8  = 8 bits ~signe	*/

typedef int		Bool;		/* booleen  (int for ANSI)	*/
		
typedef void 		*Ptr;		/* pointeur			*/
#endif

/* ==================================================== */
/*  special macro for prototypes			*/
/* ==================================================== */

#if PROTO
#define		P(s) 	s
#else
#define 	P(s) 	()
#endif
