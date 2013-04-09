/**
    MoTeX: an HPC word-based tool for MoTif eXtraction
    Copyright (C) 2012 Solon P. Pissis. 

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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#define LINE_SIZE 500
#define MAIN_ALLOC_SIZE              100

int main ( int argc, char **argv )
{

        char **         seqs    = NULL;         // the sequences in memory
	char line[LINE_SIZE];
        unsigned int    num_seqs = 0;           // the total number of sequences
        unsigned int    total_length;           // the total length of the sequences
        unsigned int    i, j;
	FILE *          in1_fd;                  // the input file descriptor
        char *          input1_filename;         // the input file
	FILE *          in2_fd;                  // the input file descriptor
        char *          input2_filename;         // the input file

	if ( argc != 3 )
	{
		fprintf ( stderr, " Error: wrong number of arguments!\n" );
		fprintf ( stderr, " %s <input motifs> <MoTeX's outputt>\n", argv[0] );
		return ( 1 );
	}
        
        input1_filename = argv[1];	          // input
        input2_filename = argv[2];		  // input

	if ( ! ( in1_fd = fopen ( input1_filename, "r") ) )
	{
		fprintf ( stderr, "Cannot open file %s\n", input1_filename );
		return ( 1 );
	}

	if ( ! ( in2_fd = fopen ( input2_filename, "r") ) )
	{
		fprintf ( stderr, "Cannot open file %s\n", input2_filename );
		return ( 1 );
	}

	while ( fgets ( line, LINE_SIZE, in2_fd ) && line[0] != '\n' ) {};	

	unsigned int max_alloc_seq = 0;
	unsigned int cur_alloc_seq = 0;

	while ( fgets ( line, LINE_SIZE, in2_fd ) && line[0] != '\n' )
	{
		
		char * pch;
		pch = strtok ( line, " " );
		//fprintf( stderr, "%s %d", pch, strlen(pch)); getchar();
		if ( cur_alloc_seq >= max_alloc_seq )
           	{
             		seqs   = ( char ** ) realloc ( seqs,   ( max_alloc_seq + MAIN_ALLOC_SIZE ) * sizeof ( char * ) ); 
             		max_alloc_seq += MAIN_ALLOC_SIZE;
           	}

          	seqs[ cur_alloc_seq ]   = strdup ( pch );
          	++ cur_alloc_seq;
	}

	seqs   = ( char ** ) realloc ( seqs, ( cur_alloc_seq + 1 ) * sizeof ( char * ) );
        seqs[ cur_alloc_seq ]   = NULL;
	num_seqs = cur_alloc_seq;

	if ( fclose ( in2_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	unsigned int succ = 0;
	unsigned int total = 0;
	while ( ! feof ( in1_fd ) )
	{	
		while ( fgets ( line, LINE_SIZE, in1_fd ) && line[0] != '\n' )
		{
			
			char * pch;
			pch = strtok ( line, "\n" );
			//fprintf( stderr, "%s %d", pch, strlen(pch)); getchar();
			for ( i = 0; i < cur_alloc_seq; ++ i )
			{
				if ( strcmp ( pch, seqs[i] ) == 0 )
				{
					succ++;
					break;
				}
			}
			total++;
		}
	}

	if ( fclose ( in1_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	fprintf( stderr, " Accuracy: %d/%d\n", succ, total);

	for ( i = 0; i < cur_alloc_seq; ++ i )
        	free ( ( void * ) seqs[i] );
        free ( seqs );
}
