/**
    MoTeX: A word-based HPC tool for MoTif eXtraction
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

#ifdef _USE_MPI
#include <mpi.h>
#endif

#ifdef _USE_OMP
#include <omp.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "motexdefs.h"

int main ( int argc, char **argv )
{
	struct TSwitch  sw;

	unsigned int 	l;			// the length for motifs
	unsigned int 	e;			// the max number of errors
	unsigned int 	d;			// the distance
	char *          alphabet;		// the alphabet
	char * 		input_filename;		// the input file
	char * 		output_filename;	// the output file
	char * 		background_in_filename;	// the background file
	char * 		un_in_filename;		// the input file of the unmatched motifs
	char * 		un_out_filename;	// the output file of the unmatched motifs
	char * 		smile_out_filename;	// the output file of the reported motifs a la SMILE
	char * 		un_smile_out_filename;	// the output file of the unmatched motifs a la SMILE
	char * 		boxes_in_filename;	// the input file for the boxes
	FILE *		in_fd;			// the input file descriptor
	FILE *		un_in_fd;		// the input file descriptor for the unmatched motifs
	FILE *		boxes_in_fd;		// the input file descriptor for the boxes

        char* 		t	= NULL;		// the sequences concatenated for broadcast and checking 
        char const ** 	seqs    = NULL; 	// the sequences in memory
	unsigned int 	num_seqs = 0;		// the total number of sequences
	unsigned int 	total_length;		// the total length of the sequences
	unsigned int	total_mot_length;   	// the min total length of the motif

	unsigned int 	nb_boxes;		// the number of boxes
	unsigned int 	nb_gaps;		// the number of gaps
	unsigned int 	nb_structs;		// the number of different motifs based on the structure of the distance intervals
   	unsigned int  * bgaps_min;		//   the min gaps
   	unsigned int  * bgaps_max;		//   the max gaps
   	unsigned int  * blens;			//   the lengths of boxes
   	unsigned int  * berrs;			//   the errors	in boxes
	unsigned int ** S;			// the array that stores the motifs structs

	unsigned int 	i, j;

	unsigned int ** g_all_occur;		// global all-occurrences vector
	unsigned int ** g_occur;		// global occurrences vector

	double start, end;

	int P = 1;				// the number of processors 

	#ifdef _USE_MPI
	int rank;				// the rank of the processor
	unsigned int long_seq;			// few and long sequences or not
	unsigned int * l_all_occur;
	unsigned int * l_occur;
	unsigned int * sendv;
	unsigned int * recv;
	int * nb_elements;
	int * disp;
	#endif

	#ifdef _USE_OMP
	unsigned int threads;
	#endif

	/* decodes the arguments */
   	i = decode_switches ( argc, argv, &sw );

   	if ( i < 15 ) 
    	{
      		usage ();
      		return ( 1 );
    	}
   	else 
    	{
		if ( sw . q > 100 )
       		{
         		fprintf ( stderr, " Error: quorum argument q (%%) should be less or equal to 100!\n" );
         		return ( 1 );
       		}

		if ( sw . Q > 100 )
       		{
         		fprintf ( stderr, " Error: max-quorum argument Q (%%) should be less or equal to 100!\n" );
         		return ( 1 );
       		}

		if ( sw . l > sizeof( unsigned int ) * CHAR_BIT - 1 )
                {
                        fprintf( stderr, " Error: the fixed-length of motifs k must be less or equal than %d!\n", ( unsigned int ) ( sizeof( unsigned int ) * CHAR_BIT - 1 ) );
                        return ( 1 );
                }
      		else 	l = sw . l;

      		e   	= sw . e;	

		if ( e >= l )
       		{
         		fprintf ( stderr, " Error: error threshold e must be less than the fixed-length of motifs k!\n" );
         		return ( 1 );
       		}

      		if ( sw . d == 0 )       d = sw . d;
      		else if ( sw . d == 1 )  d = sw . d;
      		else
       		{
         		fprintf ( stderr, " Error: distance argument d should be `0' for Hamming distance or `1' for edit distance!\n" );
         		return ( 1 );
       		}

      		if      ( ! strcmp ( "DNA", sw . alphabet ) )	alphabet = DNA;
      		else if ( ! strcmp ( "PROT", sw . alphabet ) ) 	alphabet = PROT;
      		else if ( ! strcmp ( "USR", sw . alphabet ) )  	alphabet = USR;
      		else
       		{
         		fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!\n" );
         		return ( 1 );
       		}

		//sw . alphabet = alphabet;
		input_filename      	= sw . input_filename;
		output_filename     	= sw . output_filename;
    	}

	background_in_filename  = sw . background_in_filename;
	un_in_filename          = sw . un_in_filename;
	un_out_filename         = sw . un_out_filename;
	smile_out_filename      = sw . smile_out_filename;
	un_smile_out_filename   = sw . un_smile_out_filename;
	boxes_in_filename       = sw . boxes_in_filename;
	
	/* If a background file is NOT given as input */
	if ( background_in_filename == NULL )
       	{
		if ( un_out_filename != NULL )
		{
        		fprintf ( stderr, " Error: `-u' option must be used with `-b' option!\n" );
        		return ( 1 );
		}

		if ( un_smile_out_filename != NULL )
		{
        		fprintf ( stderr, " Error: `-U' option must be used with `-b' option!\n" );
        		return ( 1 );
		}
       	}

	/* If a background file is given as input */
	if ( background_in_filename != NULL ) 
       	{
		if ( un_out_filename == NULL )
		{
			fprintf ( stderr, " Error: `-b' option must be used with `-u' option!\n" );
			return ( 1 );
		}

		/* If SMILE-like output is used, then it must be used for the unmatched motifs as well */
		if ( un_smile_out_filename != NULL  && smile_out_filename == NULL )
		{
			fprintf ( stderr, " Error: `-U' option must be used with `-S' option!\n" );
			return ( 1 );
		}
		
		if ( un_in_filename != NULL )
		{
        		fprintf ( stderr, " Error: `-I' option cannot be used with `-b' option!\n" );
        		return ( 1 );
		}	
       	}

	if ( boxes_in_filename != NULL )
       	{
		if ( un_in_filename != NULL )
		{
			fprintf ( stderr, " Error: `-I' option cannot be used with `-s' option!\n" );
			return ( 1 );
		}
	}

	#ifdef _USE_OMP
	/* set the num of threads to be used */
	threads = sw . t;
   	omp_set_num_threads( threads );
	#endif

	#ifdef _USE_MPI
        MPI_Init( NULL, NULL );
	MPI_Comm_size( MPI_COMM_WORLD, &P );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	long_seq = sw . L;

	if ( long_seq && un_in_filename != NULL )
       	{
        	fprintf ( stderr, " Error: `-I' option cannot be used with `-L' option!\n" );
        	return ( 1 );
       	}

	/* The master processor reads and broadcasts the input data as concatenated text t */
	if ( rank == 0 )
	{
	#endif
		/* read the (multi)FASTA file with the sequences to t -- record the number num_seq of sequences */
        	if ( ! ( in_fd = fopen ( input_filename, "r") ) ) 
        	{
                	fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
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

			unsigned int seq_len = 0;

        		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
        		{
				if( seq_len == 0 && c == '\n' ) 
				{
                			fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
					c = fgetc( in_fd );
					break;
				}
				if( c == '\n' ) continue;

				c = toupper( c );

                		if ( len >= max_alloc )
                		{
                        		t = ( char * ) realloc ( t,   ( max_alloc + ALLOC_SIZE ) * sizeof ( char ) ); 
                       	 		max_alloc += ALLOC_SIZE;
                		}

				if( strchr ( alphabet, c ) )
				{
                			t[ len++ ] = c;
					seq_len++;
				}
				else
        			{
                			fprintf ( stderr, " Error: input file %s contains an unexpected character!\n", input_filename );
                			return ( 1 ); 
        			}

        		}

			if( seq_len != 0 )
			{
                		if ( len >= max_alloc )
                		{
                        		t = ( char * ) realloc ( t,   ( max_alloc + ALLOC_SIZE ) * sizeof ( char ) ); 
                       	 		max_alloc += ALLOC_SIZE;
                		}
                		t[ len++ ] = DEL;
				num_seqs++;
			}

		} while( c != EOF );

		if ( fclose ( in_fd ) )
        	{
                	fprintf( stderr, " Error: file close error!\n");
                	return ( 1 );
        	}

        	t[ len - 1 ] = '\0';
        	total_length = len - 1;

		/* read the structure of the boxes */
		if ( boxes_in_filename != NULL )
		{
			/* read the boxes file in memory */
			if ( ! ( boxes_in_fd = fopen ( boxes_in_filename, "r") ) ) 
			{
				fprintf ( stderr, " Error: Cannot open file %s!\n", boxes_in_filename );
				return ( 1 ); 
			}

			if ( fscanf ( boxes_in_fd, "%d", &nb_gaps ) != 1 )
			{
				fprintf ( stderr, " Error: file %s is not in the correct format!\n", boxes_in_filename );
				return ( 1 ); 
			}

			nb_boxes = nb_gaps + 1;

			total_mot_length = l;
			bgaps_min = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );
			bgaps_max = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );
			blens = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );
			berrs = ( unsigned int * ) calloc ( nb_boxes , sizeof( unsigned int ) );

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
		}
		else
		{
			nb_gaps = 0;
			nb_boxes = 1;
			total_mot_length = l;
		}


	#ifdef _USE_MPI
		/* broadcast the data */
		MPI_Bcast ( &total_length, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send the total length
		MPI_Bcast ( &nb_gaps, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send the total number of boxes
    		MPI_Bcast ( &t[0], total_length + 1, MPI_CHAR, 0, MPI_COMM_WORLD); 	//send the actual (concatenated) text
		MPI_Bcast ( &num_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send num_seqs
		MPI_Bcast ( &total_mot_length, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send num_seqs
		MPI_Bcast ( &bgaps_min[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//send the boxes structure
		MPI_Bcast ( &bgaps_max[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//send the boxes structure
		MPI_Bcast ( &blens[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//
		MPI_Bcast ( &berrs[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//
	}
	else
	{
		/* receive the data */
		MPI_Bcast ( &total_length, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive the total length of data
		MPI_Bcast ( &nb_gaps, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive the total number of boxes

		/* allocate space for the data */
		t = ( char * ) calloc ( total_length + 1 , sizeof( char ) );
		if ( !t )
        	{
        		fprintf( stderr, " Error: the text could not be allocated!\n" );
                	return ( 1 );
        	}
		
		nb_boxes = nb_gaps + 1;	
		bgaps_min = ( unsigned int * ) calloc ( nb_gaps , sizeof( unsigned int ) );
		bgaps_max = ( unsigned int * ) calloc ( nb_gaps , sizeof( unsigned int ) );
		blens = ( unsigned int * ) calloc ( nb_gaps , sizeof( unsigned int ) );
		berrs = ( unsigned int * ) calloc ( nb_gaps , sizeof( unsigned int ) );
		
		/* receive the data into t */
    		MPI_Bcast ( &t[0], total_length + 1, MPI_CHAR, 0, MPI_COMM_WORLD ); 	//receive the actual text
		MPI_Bcast ( &num_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive num_seqs
		MPI_Bcast ( &total_mot_length, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send num_seqs
		MPI_Bcast ( &bgaps_min[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive the boxes structure
		MPI_Bcast ( &bgaps_max[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive the boxes structure
		MPI_Bcast ( &blens[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//
		MPI_Bcast ( &berrs[0], nb_gaps, MPI_INT, 0, MPI_COMM_WORLD ); 		//
	}
	#endif

	/* Each processor splits the concatenated text into the respective sequences */
        total_length = 0;
	seqs   = ( char const ** ) calloc ( ( num_seqs ) , sizeof ( char * ) ); 
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
			total_length++;	
           	}

                seq[ len++ ] = '\0';
           	seqs[i]   = strdup ( seq ); 
		
		free ( seq );
	}
	free ( t ); 

	sw . total_length = total_length;

	if ( nb_gaps )
	{
		unsigned int ** A;
		A = ( unsigned int ** ) malloc ( ( nb_gaps ) * sizeof ( unsigned int * ) );
		unsigned int * y = malloc( ( nb_gaps ) * sizeof ( unsigned int * ) );

		nb_structs = 1;	
		for ( i = 0; i < nb_gaps; i++ )
		{
			unsigned int nb_intervals = bgaps_max[i] - bgaps_min[i] + 1;
			nb_structs *= nb_intervals;
			y[i] = nb_intervals;
			A[i] = ( unsigned int * ) calloc ( ( nb_intervals ) , sizeof ( unsigned int ) );
			for ( j = 0; j < nb_intervals; j ++ )
				A[i][j] = bgaps_min[i] + j;
		}

		sw . nb_gaps = nb_gaps;
		sw . nb_boxes = nb_gaps + 1;
		sw . bgaps_min = bgaps_min;
		sw . bgaps_max = bgaps_max;
		sw . blens = blens;
		sw . berrs = berrs;
		sw . nb_structs = nb_structs;

		S = malloc( sizeof( unsigned int * ) * nb_structs );
		for ( i = 0 ; i < nb_structs; i++ ) 
			S[i] = calloc( nb_gaps, sizeof( unsigned int ) );

		unsigned int xi = 0;
		unsigned int * tmp_prod = calloc( nb_gaps, sizeof( unsigned int ) );

		recurse( nb_gaps, 0, y, A, &xi, tmp_prod, S );
	
		sw . S = S;

		free( tmp_prod );
		free( y );
		for ( i = 0; i < nb_gaps; i++ )
			free ( A[i] );
		free ( A );
	}
	else
	{
		sw . nb_boxes = nb_boxes;
		nb_structs = 1;
		sw . nb_structs = nb_structs;
	}


	/* The algorithm for motif extraction */
	if ( un_in_filename == NULL )	//This is, in the normal case, when we search for motifs in a set of sequences (MultiFASTA file)
	{
		start = gettime();

		#if defined _USE_CPU || defined _USE_OMP
		g_all_occur   = ( unsigned int ** ) calloc ( ( num_seqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_all_occur )
		{
			fprintf( stderr, " Error: the global vector could not be allocated!\n");
			return ( 1 );
		}

		g_occur   = ( unsigned int ** ) calloc ( ( num_seqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_occur )
		{
			fprintf( stderr, " Error: the global vector could not be allocated!\n");
			return ( 1 );
		}

		#ifdef _USE_OMP
		#pragma omp parallel for private ( i, j )
		#endif
		for ( i = 0; i < num_seqs; i++ )
		{
			unsigned int m = strlen ( seqs[i] );

			/* allocate space for vectors lv and gv */
			g_all_occur[i] = ( unsigned int * ) calloc ( m * nb_structs, sizeof( unsigned int ) );
			if ( ! g_all_occur[i] )
			{
				fprintf( stderr, " Error: the global all-occurrences vector could not be allocated!\n");
				exit ( 1 );
			}

			g_occur[i] = ( unsigned int * ) calloc ( m * nb_structs, sizeof( unsigned int ) );
			if ( ! g_occur[i] )
			{
				fprintf( stderr, " Error: the global occurrences vector could not be allocated!\n");
				exit ( 1 );
			}

			for ( j = 0; j < num_seqs; j++ )
			{
				unsigned int n = strlen ( seqs[j] );

				/* check if the length of the sequence satisfies the restrictions set by the algorithm */
				if ( total_mot_length > n || total_mot_length > m )
				{
					fprintf ( stderr, " Omitting short sequence in file %s!\n", input_filename );
					continue;
				}

				if ( nb_gaps == 0 )
				{

					if ( d == 0 )
					{
						if ( ! motifs_extraction_hd ( seqs[i], m, seqs[j], n, sw, g_occur[i], g_all_occur[i] ) )
						{
							fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
							exit ( 1 );
						}
					}	
					else
					{
						if ( ! motifs_extraction_ed ( seqs[i], m, seqs[j], n, sw, g_occur[i], g_all_occur[i] ) )
						{
							fprintf( stderr, " Error: motifs_extraction_ed() failed!\n");                        
							exit ( 1 );
						}
					}
				}
				else
				{
					if ( d == 0 )
					{	
						if ( ! structured_motifs_extraction_hd ( seqs[i], m, seqs[j], n, sw, g_occur[i], g_all_occur[i] ) )
						{
							fprintf( stderr, " Error: structured_motifs_extraction_hd() failed!\n");                        
							exit ( 1 );
						}
					}	
					else
					{
						if ( ! structured_motifs_extraction_ed ( seqs[i], m, seqs[j], n, sw, g_occur[i], g_all_occur[i] ) )
						{
							fprintf( stderr, " Error: structured_motifs_extraction_ed() failed!\n");                        
							exit ( 1 );
						}
					}
				}
			}
		}
		#endif

		#ifdef _USE_MPI
		g_all_occur   = ( unsigned int ** ) calloc ( ( num_seqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_all_occur )
		{
			fprintf( stderr, " Error: the global vector could not be allocated!\n");
			return ( 1 );
		}

		g_occur   = ( unsigned int ** ) calloc ( ( num_seqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_occur )
		{
			fprintf( stderr, " Error: the global vector could not be allocated!\n");
			return ( 1 );
		}

		MPI_Barrier( MPI_COMM_WORLD );

		if( rank == 0 )	start = MPI_Wtime();
	 
		if ( long_seq )
		{
			/* The algorithm for motif extraction */
			for ( i = 0; i < num_seqs; i++ )
			{
				unsigned int m = strlen ( seqs[i] );

				/* allocate space for vectors lv and gv */
				g_all_occur[i] = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
				if ( ! g_all_occur[i] )
				{
					fprintf( stderr, " Error: the global all-occurrences vector could not be allocated!\n");
					exit ( 1 );
				}

				g_occur[i] = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
				if ( ! g_occur[i] )
				{
					fprintf( stderr, " Error: the global occurrences vector could not be allocated!\n");
					exit ( 1 );
				}

				int first, last, lovcnt;
			
				vec_allocation ( rank, P, m, &first, &last, &lovcnt );

				l_all_occur = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
				if ( ! l_all_occur )
				{
					fprintf( stderr, " Error: the local all-occurrences vector could not be allocated!\n");
					return ( 1 );
				}

				l_occur = ( unsigned int * ) calloc ( lovcnt , sizeof( unsigned int ) );
				if ( ! l_occur )
				{
					fprintf( stderr, " Error: the local occurrences vector could not be allocated!\n");
					return ( 1 );
				}

				nb_elements = ( int * ) calloc ( P , sizeof( int ) );
				if ( ! nb_elements )
				{
					fprintf( stderr, " Error: the nb_elements vector could not be allocated!\n");
					return ( 1 );
				}
				MPI_Gather ( &lovcnt, 1, MPI_UNSIGNED, nb_elements, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD ); 

				disp = ( int * ) calloc ( P , sizeof( int ) );
				if ( ! disp )
				{
					fprintf( stderr, " Error: the nb_elements vector could not be allocated!\n");
					return ( 1 );
				}
				MPI_Gather ( &first, 1, MPI_UNSIGNED, disp, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD ); 

				for ( j = 0; j < num_seqs; j++ )
				{
					unsigned int n = strlen ( seqs[j] );

					/* check if the length of the sequence satisfies the restrictions set by the algorithm */
					if ( total_mot_length > n || total_mot_length > m )
					{
						fprintf ( stderr, " Omitting short sequence in file %s!\n", input_filename );
						continue;
					}

					sendv = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
					if ( ! sendv )
					{
						fprintf( stderr, " Error: the local vector could not be allocated!\n");
						return ( 1 );
					}

					recv = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
					if ( ! recv )
					{
						fprintf( stderr, " Error: the local vector could not be allocated!\n");
						return ( 1 );
					}

					if ( d == 0 )
					{
						if ( ! motifs_extraction_opasm_hd ( seqs[i], m, seqs[j], n, sw, l_all_occur, sendv, rank, P ) )
						{
							fprintf( stderr, " Error: motifs_extraction_opasm_hd() failed!\n");                        
							return ( 1 );
						}
					}	
					else
					{
						if ( ! motifs_extraction_opasm_ed ( seqs[i], m, seqs[j], n, sw, l_all_occur, sendv, rank, P ) )
						{
							fprintf( stderr, " Error: motifs_extraction_opasm_ed() failed!\n");                        
							return ( 1 );
						}
					}
				
					MPI_Allreduce( sendv, recv, m, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD ); 

					int ii, jj; 
					for ( ii = first, jj = 0; ii <= last; ii++, jj++ )	if ( recv[ii] > 0 ) l_occur[jj] = l_occur [jj] + 1;

					free ( sendv );
					free ( recv );

				}

				/* Collective communication */
				MPI_Reduce ( l_all_occur, g_all_occur[i], m, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Gatherv ( l_occur, lovcnt, MPI_UNSIGNED, g_occur[i], nb_elements, disp, MPI_UNSIGNED, 0, MPI_COMM_WORLD ); 
				free ( l_all_occur );
				free ( l_occur );
				free ( nb_elements );
				free ( disp );
			}
		}
		else
		{
			for ( i = 0; i < num_seqs; i++ )
			{
				unsigned int m = strlen ( seqs[i] );

				/* allocate space for vectors */
				g_all_occur[i] = ( unsigned int * ) calloc ( m * nb_structs, sizeof( unsigned int ) );
				if ( ! g_all_occur[i] )
				{
					fprintf( stderr, " Error: the global all-occurrences vector could not be allocated!\n");
					return ( 1 );
				}

				g_occur[i] = ( unsigned int * ) calloc ( m * nb_structs, sizeof( unsigned int ) );
				if ( ! g_occur[i] )
				{
					fprintf( stderr, " Error: the global occurrences vector could not be allocated!\n");
					return ( 1 );
				}

				l_all_occur = ( unsigned int * ) calloc ( m * nb_structs, sizeof( unsigned int ) );
				if ( ! l_all_occur )
				{
					fprintf( stderr, " Error: the local all-occurrences vector could not be allocated!\n");
					return ( 1 );
				}

				l_occur = ( unsigned int * ) calloc ( m * nb_structs, sizeof( unsigned int ) );
				if ( ! l_occur )
				{
					fprintf( stderr, " Error: the local occurrences vector could not be allocated!\n");
					return ( 1 );
				}

				int first, last, lonum_seqs;
				vec_allocation ( rank, P, num_seqs, &first, &last, &lonum_seqs );

				for ( j = first; j <= last; j++ )
				{
					unsigned int n = strlen ( seqs[j] );

					/* check if the length of the sequence satisfies the restrictions set by the algorithm */
					if ( total_mot_length > n || total_mot_length > m )
					{
						fprintf ( stderr, " Omitting short sequence in file %s!\n", input_filename );
						continue;
					}

					if ( nb_gaps == 0 )
					{
						if ( d == 0 )
						{
							if ( ! motifs_extraction_hd ( seqs[i], m, seqs[j], n, sw, l_occur, l_all_occur ) )
							{
								fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
								return ( 1 );
							}
						}	
						else
						{
							if ( ! motifs_extraction_ed ( seqs[i], m, seqs[j], n, sw, l_occur, l_all_occur ) )
							{
								fprintf( stderr, " Error: motifs_extraction_ed() failed!\n");                        
								exit ( 1 );
							}
						}
					}
					else
					{
						if ( d == 0 )
						{
							if ( ! structured_motifs_extraction_hd ( seqs[i], m, seqs[j], n, sw, l_occur, l_all_occur ) )
							{
								fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
								return ( 1 );
							}
						}	
						else
						{
							if ( ! structured_motifs_extraction_ed ( seqs[i], m, seqs[j], n, sw, l_occur, l_all_occur ) )
							{
								fprintf( stderr, " Error: motifs_extraction_ed() failed!\n");                        
								exit ( 1 );
							}
						}

					}	
				}

				/* Collective communication */
				MPI_Reduce ( l_all_occur, g_all_occur[i], m * nb_structs, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Reduce ( l_occur, g_occur[i], m * nb_structs, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
				free ( l_all_occur );
				free ( l_occur );
			}
		}

		MPI_Barrier( MPI_COMM_WORLD );
		#endif

		#ifdef _USE_MPI
		/* write the motifs */
		if( rank == 0 )
		{
			end = MPI_Wtime();
		#endif

		#if defined _USE_CPU || defined _USE_OMP
		end = gettime();
		#endif

		#ifdef _USE_OMP
		P = sw . t;
		#endif
			if ( nb_gaps == 0 )
			{
				if ( background_in_filename == NULL )
				{
					write_motifs ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
					if ( smile_out_filename != NULL )
						write_motifs_smile ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start );
				}
				else
				{
					write_motifs_back ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
					if ( un_smile_out_filename != NULL )
						write_motifs_back_smile ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
				}

			}
			else
			{
				if ( background_in_filename == NULL )
				{
					write_structured_motifs ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
					if ( smile_out_filename != NULL )
						write_structured_motifs_smile ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start );
				}
				else
				{
					write_structured_motifs_back ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
					if ( un_smile_out_filename != NULL )
						write_structuted_motifs_back_smile ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
				}
			}
				
		#ifdef _USE_MPI
		}
		#endif

		/* Deallocate */
		for ( i = 0; i < num_seqs; i ++ ) 
		{
			free ( ( void * ) g_all_occur[i] );
			free ( ( void * ) g_occur[i] );
		}
		free ( g_all_occur );	
		free ( g_occur );	

	}
	else //This is in the case when we search for foreground motifs, which were previously not exactly matched with background motifs, in a set of sequences (MultiFASTA file)
	{
		/* open the file with the unmatched motifs for reading */
        	if ( ! ( un_in_fd = fopen ( un_in_filename, "r") ) ) 
        	{
                	fprintf ( stderr, " Error: Cannot open file %s!\n", un_in_filename );
                	return ( 1 ); 
        	}

                struct 		Tdata * fdata = NULL;
		unsigned int 	num_fseqs;
		unsigned int 	max_alloc;
		char 		line[LINE_SIZE];
		char const	** fseqs      = NULL;

		/* read the motifs in memory */
		max_alloc = num_fseqs = 0;

		/* This is to avoid the header */
		while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' );
 
		while ( fgets ( line, LINE_SIZE, un_in_fd ) && line[0] != '\n' )
		{
			if ( num_fseqs >= max_alloc )
                	{
                		fdata   = ( struct Tdata * ) realloc ( fdata,   ( max_alloc + ALLOC_SIZE ) * sizeof ( struct Tdata ) );
             			fseqs   = ( char const ** ) realloc ( fseqs,   ( max_alloc + ALLOC_SIZE ) * sizeof ( char * ) ); 
                        	max_alloc += ALLOC_SIZE;
               		}
			
			char * pch;
			pch = strtok ( line, " " );
    			fseqs[ num_fseqs  ]   = strdup ( pch );

			pch = strtok ( NULL, " " );
    			fdata[ num_fseqs ] . u = atoi( pch );

			pch = strtok ( NULL, " " );
    			fdata[ num_fseqs ] . n = atoi( pch );

			pch = strtok ( NULL, " " );
    			fdata[ num_fseqs ] . r = atof( pch );

    			pch = strtok ( NULL, " ");
    			fdata[ num_fseqs ] . v = atoi( pch );

			num_fseqs++;
		}

		if ( fclose ( un_in_fd ) ) 
		{
      			fprintf( stderr, " Error: file close error!\n");				      
			return ( 1 );
		}

        	start = gettime();
		
		/* find whether the above read sequences are motifs in the set of sequences (MultiFASTA file)*/
		#if defined _USE_CPU || defined _USE_OMP
		g_all_occur   = ( unsigned int ** ) calloc ( ( num_fseqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_all_occur )
        	{
        		fprintf( stderr, " Error: the global vector could not be allocated!\n");
               		return ( 1 );
        	}

		g_occur   = ( unsigned int ** ) calloc ( ( num_fseqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_occur )
        	{
        		fprintf( stderr, " Error: the global vector could not be allocated!\n");
               		return ( 1 );
        	}

		#ifdef _USE_OMP
		#pragma omp parallel for private ( i, j )
		#endif
		for ( i = 0; i < num_fseqs; i++ )
		{
			unsigned int m = strlen ( fseqs[i] );
			l = m;

			/* allocate space for vectors lv and gv */
			g_all_occur[i] = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
			if ( ! g_all_occur[i] )
        		{
        			fprintf( stderr, " Error: the global all-occurrences vector could not be allocated!\n");
                		exit ( 1 );
        		}

			g_occur[i] = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
			if ( ! g_occur[i] )
        		{
        			fprintf( stderr, " Error: the global occurrences vector could not be allocated!\n");
                		exit ( 1 );
        		}

			for ( j = 0; j < num_seqs; j++ )
			{
				unsigned int n = strlen ( seqs[j] );

				/* check if the length of the sequence satisfies the restrictions set by the algorithm */
				if ( total_mot_length > n )
				{
					fprintf ( stderr, " Omitting short sequence in file %s!\n", input_filename );
					continue;
				}

				if ( d == 0 )
				{
					if ( ! motifs_extraction_hd ( fseqs[i], m, seqs[j], n, sw, g_occur[i], g_all_occur[i] ) )
        				{
              					fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
              					exit ( 1 );
        				}
        			}	
				else
				{
					if ( ! motifs_extraction_ed ( fseqs[i], m, seqs[j], n, sw, g_occur[i], g_all_occur[i] ) )
        				{
              					fprintf( stderr, " Error: motifs_extraction_ed() failed!\n");                        
              					exit ( 1 );
        				}
        			}
			}
		}
		#endif

		#ifdef _USE_MPI
		g_all_occur   = ( unsigned int ** ) calloc ( ( num_fseqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_all_occur )
        	{
        		fprintf( stderr, " Error: the global vector could not be allocated!\n");
            		return ( 1 );
        	}

		g_occur   = ( unsigned int ** ) calloc ( ( num_fseqs ) , sizeof ( unsigned int * ) ); 
		if ( ! g_occur )
        	{
        		fprintf( stderr, " Error: the global vector could not be allocated!\n");
         		return ( 1 );
        	}

		MPI_Barrier( MPI_COMM_WORLD );

		if( rank == 0 )	start = MPI_Wtime();

		for ( i = 0; i < num_fseqs; i++ )
		{
			unsigned int m = strlen ( fseqs[i] );
			l = m;

			/* allocate space for vectors */
			g_all_occur[i] = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
			if ( ! g_all_occur[i] )
        		{
        			fprintf( stderr, " Error: the global all-occurrences vector could not be allocated!\n");
                		return ( 1 );
        		}

			g_occur[i] = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
			if ( ! g_occur[i] )
        		{
        			fprintf( stderr, " Error: the global occurrences vector could not be allocated!\n");
                		return ( 1 );
        		}

			l_all_occur = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
			if ( ! l_all_occur )
        		{
        			fprintf( stderr, " Error: the local all-occurrences vector could not be allocated!\n");
                		return ( 1 );
        		}

			l_occur = ( unsigned int * ) calloc ( m , sizeof( unsigned int ) );
			if ( ! l_occur )
        		{
        			fprintf( stderr, " Error: the local occurrences vector could not be allocated!\n");
                		return ( 1 );
        		}

			int first, last, lonum_seqs;
			vec_allocation ( rank, P, num_fseqs, &first, &last, &lonum_seqs );

			for ( j = first; j <= last; j++ )
			{
				unsigned int n = strlen ( seqs[j] );

				/* check if the length of the sequence satisfies the restrictions set by the algorithm */
				if ( total_mot_length > n )
				{
					fprintf ( stderr, " Omitting short sequence in file %s!\n", input_filename );
					continue;
				}

				if ( d == 0 )
				{
					if ( ! motifs_extraction_hd ( fseqs[i], m, seqs[j], n, sw, l_occur, l_all_occur ) )
        				{
              					fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
              					return ( 1 );
        				}
        			}	
				else
				{
					if ( ! motifs_extraction_ed ( fseqs[i], m, seqs[j], n, sw, l_occur, l_all_occur ) )
        				{
              					fprintf( stderr, " Error: motifs_extraction_ed() failed!\n");                        
              					exit ( 1 );
        				}
        			}
			}

			/* Collective communication */
			MPI_Reduce ( l_all_occur, g_all_occur[i], m, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
			MPI_Reduce ( l_occur, g_occur[i], m, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
			free ( l_all_occur );
			free ( l_occur );
		}

		MPI_Barrier( MPI_COMM_WORLD );
		#endif
 
		/* write the output */
		#ifdef _USE_MPI
		if( rank == 0 )
		{
	 		end = MPI_Wtime();
		#endif

		#if defined _USE_CPU || defined _USE_OMP
        	end = gettime();
		#endif

		#ifdef _USE_OMP
		P = sw . t;
		#endif

                write_motifs_fore ( sw, num_fseqs, fseqs, g_occur, g_all_occur, end - start, P, num_seqs, fdata );
		#ifdef _USE_MPI
		}
		#endif   	

		/* Deallocate */
		for ( i = 0; i < num_fseqs; i ++ ) 
		{
			free ( ( void * ) fseqs[i] );
			free ( ( void * ) g_all_occur[i] );
			free ( ( void * ) g_occur[i] );
		}
		free ( fseqs );	
		free ( g_all_occur );	
		free ( g_occur );	
		free ( fdata );
        	fdata = NULL;
	}

	/* Deallocate & Finalize */
	for ( i = 0; i < num_seqs; i ++ ) 
		free ( ( void * ) seqs[i] );
	free ( seqs );	
	free ( sw . input_filename );
   	free ( sw . output_filename );
   	free ( sw . alphabet );
	if ( background_in_filename )    free ( sw . background_in_filename );
	if ( un_in_filename )  free ( sw . un_in_filename );
	if ( un_out_filename ) free ( sw . un_out_filename );
	if ( smile_out_filename )    free ( sw . smile_out_filename );
	if ( un_smile_out_filename )    free ( sw . un_smile_out_filename );
	if ( nb_gaps )
	{ 
		free ( sw . boxes_in_filename );
		free ( sw . bgaps_min );
		free ( sw . bgaps_max );
		free ( sw . blens );
		free ( sw . berrs );
		for ( i = 0 ; i < nb_structs; i++ ) 
			free ( S[i] );
		free ( S );
	}

	#ifdef _USE_MPI
	MPI_Finalize();
	#endif

	return ( 0 );
}
