/**
    libFLASM
    Copyright (C) 2016 Lorraine A. K. Ayad, Solon P. Pissis and Ahmad Retha

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <stdlib.h>
#include "libflasm.h"

/**
 * This is the libFLASM edit distance function.
 *
 * @param t The text (haystack) to search in
 * @param n The length of t
 * @param x The pattern which has factors that may be present in t
 * @param m The length of x
 * @param factor_length The length of a factor (needle)
 * @param max_error The maximum distance between the factor and a position in t to report
 * @return The discovered positions are returned in a set that can be iterated over
 */
ResultTupleSet flasm_ed ( unsigned char * t, unsigned int n, unsigned char * x, unsigned int m, unsigned int factor_length, unsigned int max_error )
{
	unsigned char * h;
	h = ( unsigned char * ) calloc ( factor_length + 1, sizeof ( unsigned char ) );

	int k = -max_error;
	unsigned int error = 0;
	unsigned int pos_t = 0;
	unsigned int pos_x = 0;

	ResultTupleSet results;

	CharString haystack = (CharString) t;
	CharString needle;
	Finder<CharString> finder( haystack );
	Pattern<CharString, Myers<> > pattern;

	unsigned int i;
	for ( i = 0; i < m - factor_length + 1; i++ )
	{
		memcpy( &h[0], &x[i], factor_length );
		h[factor_length] = '\0';

		needle = (CharString) h;

		setNeedle( pattern, needle );

		while ( find( finder, pattern, k ) )
		{
			pos_t = (unsigned int) endPosition( finder ) - 1;

			pos_x = i + factor_length - 1;

			error = (unsigned int) abs( getScore( pattern ) );

			ResultTuple match = {pos_t, pos_x, error};

			results.insert( match );
		}
		
		clear( finder );
		goBegin ( finder );

	}

	free ( h );

	return results;
}

/**
 * This is the libFLASM Hamming distance function.
 *
 * @param t The text (haystack) to search in
 * @param n The length of t
 * @param x The pattern which has factors that may be present in t
 * @param m The length of x
 * @param factor_length The length of a factor (needle)
 * @param max_error The maximum distance between the factor and a position in t to report
 * @return The discovered positions are returned in a set that can be iterated over
 */
ResultTupleSet flasm_hd ( unsigned char * t, unsigned int n, unsigned char * x, unsigned int m, unsigned int factor_length, unsigned int max_error )
{
	ResultTupleSet results;

	unsigned int i, j, k, err;

    	_Limit lim;
    	lim = init_limit ( factor_length, lim );
    	WORD * ones;
    	if ( ( ones = ( WORD * ) calloc ( lim.words , sizeof ( WORD ) ) ) == NULL )
    	{
        	fprintf( stderr, " Error: ow could not be allocated!\n");
        	return results;
    	}
    
    	//initialise 2 line matrix
    	WORD ** M0;
    	WORD ** M1;
    	if ( ( M0 = ( WORD ** ) calloc ( ( n + 1 ) , sizeof ( WORD * ) ) ) == NULL )
    	{
        	fprintf( stderr, " Error: M0 could not be allocated!\n");
        	return results;
    	}
    	if ( ( M1 = ( WORD ** ) calloc ( ( n + 1 ) , sizeof ( WORD * ) ) ) == NULL )
    	{
        	fprintf( stderr, " Error: M1 could not be allocated!\n");
        	return results;
    	}
    	for ( j = 0; j < n + 1; j ++ )
    	{
        	if ( ( M0[j] = ( WORD * ) calloc ( lim.words , sizeof ( WORD ) ) ) == NULL )
        	{
            		fprintf( stderr, " Error: M0J could not be allocated!\n");
            		return results;
        	}
        	if ( ( M1[j] = ( WORD * ) calloc ( lim.words , sizeof ( WORD ) ) ) == NULL )
        	{
            	fprintf( stderr, " Error: M1J could not be allocated!\n");
            	return results;
        }
    }
    
    //loop through sequences
    for ( i = 1; i < m + 1; i++ ) //loop through m
    {
    	for ( j = 0; j < n + 1; j++ ) //loop through t
        {
        	//make ones
        	if ( j == 0  && i <= factor_length ) 
		{
                	ones = shift_words ( ones, lim.words );
                	ones[lim.words - 1] = ones[lim.words - 1] + 1;
            	}
            
        	switch ( i % 2 ) 
		{
                	case 0:
                
                    	if ( j == 0 )
                    	{
                        	//fill up the first column with ones up to length h
                        	memcpy ( M1[j], ones, lim.words * sizeof ( WORD ) );
                    	}
                    	else
                    	{
                        	//copy values from diagonal cell into current cell
                        	for ( k = 0; k < lim.words; k++ ) 
				{
                            		M1[j][k] = M0[j - 1][k];
                        	}

                        	//shift things along one and clear left most bit
                        	M1[j] = shiftc_words ( M1[j], lim );

                        	//set last bit on right to hamming distance of the characters
                        	M1[j][lim.words - 1] = M1[j][lim.words - 1] | delta ( x[i - 1], t[j - 1] );
                    	}

                    	if ( j >= factor_length && i >= factor_length )
                    	{                
                       		err = popcount_words ( M1[j], lim.words );

                        	if ( err <= max_error )
                        	{
		     			ResultTuple match = { j - 1, i - 1, err };
                            		results.insert ( match );
                        	}
                        
                    	}
                    	break;
                    
                	case 1:
                    
                    	if ( j == 0 )
                    	{
                        	//fill up the first column with ones up to length h
                        	memcpy ( M0[j], ones, lim.words * sizeof ( WORD ) );
                    	}
                    	else
                    	{
                        	//copy values from diagonal cell into current cell
                        	for ( k = 0; k < lim.words; k++ ) 
				{
                            		M0[j][k] = M1[j - 1][k];
                        	}

                        	//shift things along one and clear left most bit
                        	M0[j] = shiftc_words ( M0[j], lim );

                        	//set last bit on right to hamming distance of the characters
                        	M0[j][lim.words - 1] = M0[j][lim.words - 1] | delta ( x[i - 1], t[j - 1] );
                    	}

                    	if ( j >= factor_length && i >= factor_length )
                    	{                
                        	err = popcount_words ( M0[j], lim.words );

                        	if ( err <= max_error )
                        	{
      					ResultTuple match = { j - 1, i - 1, err };
                            		results.insert ( match );
				}
                    	}
                    	break;
            	}
            }
	}
    
    	for ( j = 0; j < n + 1; j ++ ) 
	{
        	free ( M0[j] );
        	free ( M1[j] );
    	}
    	free ( M0 );
    	free ( M1 );
    	free ( ones );
    	
	return results;
}

inline _Limit init_limit ( unsigned int h, struct _Limit lim )
{
    double WSd = (double) WORD_SIZE;
    unsigned int WSi = (unsigned int) WORD_SIZE;
    
    lim.h = h;
    
    WORD yWord = ULONG_MAX; //111111
    unsigned int shift_places = ( WSi - ( h % WSi ) ) % WSi;
    yWord = yWord >> shift_places; //if shift_places is 2: 001111
    lim.yWord = yWord;
    
    lim.words = (unsigned int) ceil ( (double) lim.h / WSd );
    lim.yIndex = 0;
    
    return lim;
}

/**
 * If you supply a WORD array it returns the sum of the popcount on them
 *  
 * See <a href="http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html">Faster popcount</a>
 * 
 * @param words WORD Array
 * @param length Number of elements in the WORD array
 * @return Sum of popcounts on a WORD array
 */
inline unsigned int popcount_words ( WORD * words, int length )
{
    unsigned int count = 0;
    unsigned int i;
    for ( i = 0; i < length; i++ ) {
        count += __builtin_popcountl ( words[i] );
    }
    return count;
}

/**
 * Shifts bits in an array of WORDs one position to the left
 * 
 * @param words
 * @param length Number of elements in the WORD array
 * @return 
 */
inline WORD * shift_words ( WORD * words, int length )
{
    WORD mask = (WORD) 1 << ( (int) WORD_SIZE - 1 );
    WORD carried_bit = 0;
    WORD temp;
    int i;
    for ( i = length - 1; i > -1; i-- ) {
        temp = words[i];
        words[i] = (WORD) ( ( words[i] << 1 ) | carried_bit );
        carried_bit = (WORD) ( ( temp & mask ) != 0 );
    }
    return words;
}

/**
 * Shifts bits left one position then truncates left-most bits on most
 * significant WORD of array using Limit.yWord mask
 * 
 * @param words WORDs array
 * @param lim An initialised Limit structure
 * @return 
 */
inline WORD * shiftc_words ( WORD * words, struct _Limit lim )
{
    words = shift_words ( words, lim.words );
    words[lim.yIndex] = words[lim.yIndex] & lim.yWord;
    return words;
}

/**
 * Returns the hamming distance (1 or 0) for two characters a and b
 * 
 * @param a
 * @param b
 * @return Hamming distance between characters given
 */
inline WORD delta ( char a, char b )
{
    return (WORD)( a != b );
}
