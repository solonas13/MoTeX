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

#define ALLOC_SIZE              512
#define LINE_SIZE               512
#define DEL                     '$' 

#define DNA                     "ACGTN"                         //DNA alphabet
#define PROT                    "ARNDCQEGHILKMFPSTWYV"          //Proteins alphabet
#define USR                     "0123456789"                    //User-defined alphabet, e.g. integer alphabet

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

int main ( int argc, char **argv )
{

        char*           t       = NULL;         // the sequences concatenated for broadcast and checking 
        char **         seqs    = NULL;         // the sequences in memory
        unsigned int    num_seqs = 0;           // the total number of sequences
        unsigned int    total_length;           // the total length of the sequences
	char            * alphabet;
        unsigned int    i, j;
	FILE *          in_fd;                    // the input file descriptor
        char *          input_filename;           // the input file
	FILE *          boxes_in_fd;              // the input file descriptor
        char *          boxes_in_filename;     // the input boxes file
	FILE *          out1_fd;                  // the input file descriptor
        char *          output1_filename;         // the input file
	FILE *          out2_fd;                  // the input file descriptor
        char *          output2_filename;         // the input file

	if ( argc != 9 )
	{
		fprintf ( stderr, " Error: wrong number of arguments!\n" );
		fprintf ( stderr, " %s <input> <k> <e_max> <boxes> <q> <# motifs> <sequences FASTA output> <implanted motifs output>\n", argv[0] );
		return ( 1 );
	}
        
	alphabet = DNA;
        input_filename = argv[1];	          // input
        unsigned int k = atoi ( argv[2] );        // motif length         
	unsigned int eMAX = atoi ( argv[3] );     // (unsigned int)(rand())%( eMAX + 1 );        // number of errors
	boxes_in_filename = argv[4];     	  // boxes
	unsigned int q = atoi ( argv[5] );	  // proportion of sequences into which we will insert the motif
        unsigned int num_mot = atoi ( argv[6] );  // number of motifs to be inserted         
        output1_filename = argv[7];		  // output
        output2_filename = argv[8];		  // output

	/* read the (multi)FASTA file with the sequences to t -- record the number num_seq of sequences */
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, "Cannot open file %s\n", input_filename );
		return ( 1 );
	}

	unsigned int max_alloc, len;
	max_alloc = len = 0;
	char c;
	c = fgetc( in_fd );

	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' );

		num_seqs++;
		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( c == '\n' ) continue;

			c = toupper( c );

			if ( len >= max_alloc )
			{
				t = ( char * ) realloc ( t,   ( max_alloc + ALLOC_SIZE ) * sizeof ( char ) );
				max_alloc += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
				t[ len++ ] = c;
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character!\n", input_filename );
				return ( 1 );
			}

		}
		t[ len++ ] = DEL;

	} while( c != EOF );

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	t[ len - 1 ] = '\0';
	total_length = len - 1;

	seqs   = ( char ** ) calloc ( ( num_seqs ) , sizeof ( char * ) );
        for ( i = 0, j = 0; i < num_seqs; i++ )
        {
                unsigned int max_alloc, len;
                max_alloc = len = 0;
                char * seq = NULL;
                char c;

                while ( ( c = t[ j++ ] ) && c != DEL && c != '\0' )
                {
                        if ( len >= max_alloc )
                        {
                                seq = ( char * ) realloc ( seq,   ( max_alloc + ALLOC_SIZE ) * sizeof ( char ) );
                                max_alloc += ALLOC_SIZE;
                        }
                        seq[ len++ ] = c;
                }

                seq[ len++ ] = '\0';
                seqs[i]   = strdup ( seq );

                free ( seq );
        }
        free ( t );

	/* read the boxes file in memory */
	if ( ! ( boxes_in_fd = fopen ( boxes_in_filename, "r") ) ) 
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", boxes_in_filename );
		return ( 1 );
	}

	unsigned int nb_gaps = 0;

	if ( fscanf ( boxes_in_fd, "%d", &nb_gaps ) != 1 )
	{
		fprintf ( stderr, " Error: file %s is not in the correct format!\n", boxes_in_filename );
		return ( 1 );
	}

	unsigned int nb_boxes = nb_gaps + 1;

	unsigned int total_mot_length = k;
	unsigned int * bgaps_min = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );
	unsigned int * bgaps_max = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );
	unsigned int * blens = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );
	unsigned int * berrs = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );

	for ( i = 0; i < nb_gaps; i ++ )
	{
		if ( fscanf ( boxes_in_fd, "%d", &bgaps_min[i] ) != 1 )
		{
			fprintf ( stderr, " Error: file %s is not in the correct format!\n", boxes_in_filename );
			return ( 1 );
		}
		total_mot_length += bgaps_min[i];
		if ( fscanf ( boxes_in_fd, "%d", &bgaps_max[i] ) != 1 )
		{
			fprintf ( stderr, " Error: file %s is not in the correct format!\n", boxes_in_filename );
			return ( 1 );
		}
		if ( fscanf ( boxes_in_fd, "%d", &blens[i] ) != 1 )
		{
			fprintf ( stderr, " Error: file %s is not in the correct format!\n", boxes_in_filename );
			return ( 1 );
		}
		total_mot_length += blens[i];
		if ( fscanf ( boxes_in_fd, "%d", &berrs[i] ) != 1 )
		{
			fprintf ( stderr, " Error: file %s is not in the correct format!\n", boxes_in_filename );
			return ( 1 );
		}
	}

	if ( fclose ( boxes_in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if ( ! ( out2_fd = fopen ( output2_filename, "w") ) )
	{
		fprintf ( stderr, "Cannot open file %s\n", output2_filename );
		return ( 1 );
	}

	for ( i = 0; i < num_mot; i ++ )
	{
		//Create the motif
                char * motif = NULL;
                motif = calloc ( ( total_mot_length + 1 ) , sizeof( char ) );
                char * out_motif = NULL;
                out_motif = calloc ( ( total_mot_length + 1 ) , sizeof( char ) );

		unsigned int b;
		unsigned int index = 0;
		unsigned int ell;
		unsigned int gap;

		for ( j = 0; j < total_mot_length; j++ )
		{
			unsigned int DNAchar = ( unsigned int ) ( rand() ) %( 4 );        // number of errors
			if ( DNAchar == 0 ) motif[j] = 'A';
			if ( DNAchar == 1 ) motif[j] = 'C';
			if ( DNAchar == 2 ) motif[j] = 'G';
			if ( DNAchar == 3 ) motif[j] = 'T';
		}
                motif[ j ] = '\0';

		unsigned int offset = 0;
		for ( b = 0; b < nb_boxes; b ++ )
		{
			/* Extract the structured motif from the sequence */
			if ( b == 0 )
			{
				ell = k;
				for ( j = 0; j < ell; j++, index++ )
					out_motif[index] = motif[j];
				out_motif [ index ] = '_';
				index ++;
			}
			else
			{
				offset += bgaps_min[b - 1] + ell;
				ell = blens[b - 1];
				for ( j = 0; j < ell; j++, index++ )
					out_motif[ index ] = motif[offset + j];
				out_motif [ index ] = '_';
				index ++;
			}
		}
		index--;
                out_motif[ index ] = '\0';
		//fprintf ( stderr, "%s\n", motif );	
		//fprintf ( stderr, "%s\n", out_motif );
		fprintf ( out2_fd, "%s\n", out_motif );	

		//Compute the number of sequences where we will implant motifs
		unsigned int alt_seqs = ceil ( ( double ) ( q * num_seqs ) / ( double ) 100 );
 		unsigned int * alt_seqs_arr;
                alt_seqs_arr = ( unsigned int * ) calloc ( ( num_seqs ) , sizeof ( unsigned int ) );
		unsigned int succ = 0;

		while ( succ < alt_seqs )
		{
			unsigned int pos = ( unsigned int ) ( rand() ) % ( num_seqs + 1 );
			if ( alt_seqs_arr[ pos ] == 0 )
			{
				alt_seqs_arr[ pos ] = 1;
				succ ++;
			}
		}

		//Implant the structured motifs by also inserting errors in them
		for ( j = 0; j < num_seqs; j ++ )
		{
			if ( alt_seqs_arr[ j ] == 1 )
			{
				//fprintf ( stderr, "%s\n", motif ); 
                		char * motif_cpy = ( char * ) calloc ( ( total_mot_length + 1 ) , sizeof ( char ) );
				strcpy ( motif_cpy, motif );

				unsigned int offset = 0;
				for ( b = 0; b < nb_boxes; b ++ )
				{
					/* Extract the structured motif from the sequence */
					if ( b == 0 )
					{
						ell = k;
						unsigned int e = ( unsigned int ) ( rand() ) % ( eMAX + 1 );
						unsigned int jj;
						for ( jj = 0; jj < e; jj++ )
						{
							unsigned int pos = ( unsigned int ) ( rand() ) % ( ell );
							unsigned int DNAchar = ( unsigned int ) ( rand() ) %( 4 );        
							if ( DNAchar == 0 ) motif_cpy[pos] = 'A';
							if ( DNAchar == 1 ) motif_cpy[pos] = 'C';
							if ( DNAchar == 2 ) motif_cpy[pos] = 'G';
							if ( DNAchar == 3 ) motif_cpy[pos] = 'T';
						}
					}
					else
					{
						offset += bgaps_min[b - 1] + ell;
						ell = blens[b - 1];
						unsigned int e = ( unsigned int ) ( rand() ) % ( berrs[b - 1] + 1 );
						unsigned int jj;
						for ( jj = 0; jj < e; jj++ )
						{
							unsigned int pos = ( unsigned int ) ( rand() ) % ( ell );
							unsigned int DNAchar = ( unsigned int ) ( rand() ) %( 4 );        
							if ( DNAchar == 0 ) motif_cpy[pos + offset] = 'A';
							if ( DNAchar == 1 ) motif_cpy[pos + offset] = 'C';
							if ( DNAchar == 2 ) motif_cpy[pos + offset] = 'G';
							if ( DNAchar == 3 ) motif_cpy[pos + offset] = 'T';
						}
					}
				}
				//fprintf ( stderr, "%s\n", motif ); 
				//fprintf ( stderr, "%s\n", motif_cpy ); 
				//fprintf ( stderr, "%s\n", seqs[j] ); 

				char * new_seq = ( char * ) calloc ( strlen( motif_cpy ) + strlen ( seqs[ j ] ) + 1, sizeof ( char ) );
				sprintf ( new_seq, "%s%s", seqs[j], motif_cpy );
				new_seq [ strlen( motif_cpy ) + strlen ( seqs[ j ] ) ] = '\0'; 
                		free ( seqs[j] );
                		seqs[j] = strdup ( new_seq );
				//fprintf ( stderr, "%s\n", seqs[j] ); 
				//getchar();
                		free ( motif_cpy );
                		free ( new_seq );
			}
		}

                free ( motif );
                free ( out_motif );
                free ( alt_seqs_arr );
	}

	if ( fclose ( out2_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	if ( ! ( out1_fd = fopen ( output1_filename, "w") ) )
	{
		fprintf ( stderr, "Cannot open file %s\n", output1_filename );
		return ( 1 );
	}

        for ( i = 0; i < num_seqs; i ++ )
	{
		fprintf ( out1_fd, ">%d\n", i );	
		fprintf ( out1_fd, "%s\n", seqs[i] );	
	}
	fprintf ( stderr, "%d structured motifs were implanted; each in %d distinct sequences.\n", num_mot,  ( int ) ceil ( ( double ) ( q * num_seqs ) / ( double ) 100 ) );	

	if ( fclose ( out1_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

        /* Deallocate & Finalize */
        for ( i = 0; i < num_seqs; i ++ )
                free ( ( void * ) seqs[i] );
        free ( seqs );
	free ( bgaps_min );
	free ( bgaps_max );
	free ( blens );
	free ( berrs );

}
