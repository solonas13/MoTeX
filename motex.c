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

	unsigned int 	l;			// the length of motifs
	unsigned int 	q;			// the quorum
	unsigned int 	e;			// the max number of errors
	unsigned int 	d;			// the distance
	char *          alphabet;		// the alphabet
	char * 		input_filename;		// the input file
	char * 		output_filename;	// the output file
	char * 		background_filename;	// the background file
	FILE *		in_fd;			// the input file descriptor

        char* 		t	= NULL;		// the sequences concatenated for broadcast and checking 
        char const ** 	seqs    = NULL; 	// the sequences in memory
	unsigned int 	num_seqs = 0;		// the total number of sequences
	unsigned int 	total_length;		// the total length of the sequences

	unsigned int 	n, m;			
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
      		q   	= sw . q;	
      		l   	= sw . l;	
      		e   	= sw . e;	

      		if ( sw . d == 0 )       d = sw . d;
      		else if ( sw . d == 1 )  d = sw . d;
      		else
       		{
         		fprintf ( stderr, "Error: distance argument d should be `0' for Hamming distance or `1' for edit distance!!!\n" );
         		return ( 1 );
       		}

      		if      ( ! strcmp ( "DNA", sw . alphabet ) )	alphabet = DNA;
      		else if ( ! strcmp ( "PROT", sw . alphabet ) ) 	alphabet = PROT;
      		else if ( ! strcmp ( "USR", sw . alphabet ) )  	alphabet = USR;
      		else
       		{
         		fprintf ( stderr, "Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences or `USR' for sequences over a user-defined alphabet!!!\n" );
         		return ( 1 );
       		}

		input_filename      = sw . input_filename;
		output_filename     = sw . output_filename;
    	}

	background_filename = sw . background_filename;
        
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

	/* The master processor reads and broadcasts the input data as concatenated text t */
	if ( rank == 0 )
	{
	#endif
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

        	fclose ( in_fd );
        	t[ len - 1 ] = '\0';
        	total_length = len - 1;

		/* l should be less than sizeof( type ) * CHAR_BIT */
		if ( l > sizeof( unsigned int ) * CHAR_BIT - 1 )
                {
                        fprintf( stderr, " Error: the fixed-length of motifs must be less or equal than %d!\n", ( unsigned int ) ( sizeof( unsigned int ) * CHAR_BIT - 1 ) );
                        return ( 1 );
                }

	#ifdef _USE_MPI
		/* broadcast the data */
		MPI_Bcast ( &total_length, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send the total length
    		MPI_Bcast ( &t[0], total_length + 1, MPI_CHAR, 0, MPI_COMM_WORLD); 	//send the actual (concatenated) text
		MPI_Bcast ( &num_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//send num_seqs
	}
	else
	{
		/* receive the data */
		MPI_Bcast ( &total_length, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive the total length of data

		/* allocate space for the data */
		t = ( char * ) calloc ( total_length + 1 , sizeof( char ) );
		if ( !t )
        	{
        		fprintf( stderr, " Error: the text could not be allocated!\n" );
                	return ( 1 );
        	}
		
		/* receive the data into t */
    		MPI_Bcast ( &t[0], total_length + 1, MPI_CHAR, 0, MPI_COMM_WORLD ); 	//receive the actual text
		MPI_Bcast ( &num_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD ); 		//receive num_seqs
	}
	#endif
	
	/* Each processor splits the concatenated text into the respective sequences */
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
           	}

                seq[ len++ ] = '\0';
           	seqs[i]   = strdup ( seq ); 
		
		free ( seq );
	}
	free ( t );

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

        start = gettime();

	/* The algorithm for motif extraction */
	for ( i = 0; i < num_seqs; i++ )
	{
		m = strlen ( seqs[i] );

		/* allocate space for vectors lv and gv */
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

		#ifdef _USE_OMP
		#pragma omp parallel for private ( j )
		#endif
		for ( j = 0; j < num_seqs; j++ )
		{
			n = strlen ( seqs[j] );

			/* check if the length of the sequence satisfies the restrictions set by the algorithm */
			if ( l > n || l > m )
			{
        			fprintf( stderr, " Error: the fixed-length of motifs must be less or equal to the length of the sequences!\n");
                		exit ( 1 );
        		}

			if ( d == 0 )
			{
				if ( ! motifs_extraction_hd ( seqs[i], m, seqs[j], n, l, e, g_occur[i], g_all_occur[i] ) )
        			{
              				fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
              				exit ( 1 );
        			}
        		}	
			else
			{
				if ( ! motifs_extraction_ed ( seqs[i], m, seqs[j], n, l, e, g_occur[i], g_all_occur[i] ) )
        			{
              				fprintf( stderr, " Error: motifs_extraction_ed() failed!\n");                        
              				exit ( 1 );
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
			m = strlen ( seqs[i] );

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
				n = strlen ( seqs[j] );

				/* check if the length of the sequence satisfies the restrictions set by the algorithm */
				if ( l > n || l > m )
				{
        				fprintf( stderr, " Error: the fixed-length of motifs must be less or equal to the length of the sequences!\n");
                			exit ( 1 );
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
					if ( ! motifs_extraction_opasm_hd ( seqs[i], m, seqs[j], n, l, e, l_all_occur, sendv, rank, P ) )
        				{
              					fprintf( stderr, " Error: motifs_extraction_opasm_hd() failed!\n");                        
              					return ( 1 );
        				}
        			}	
				else
				{
					if ( ! motifs_extraction_opasm_ed ( seqs[i], m, seqs[j], n, l, e, l_all_occur, sendv, rank, P ) )
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
			m = strlen ( seqs[i] );

			/* allocate space for vectors */
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
			vec_allocation ( rank, P, num_seqs, &first, &last, &lonum_seqs );

			for ( j = first; j <= last; j++ )
			{
				n = strlen ( seqs[j] );

				/* check if the length of the sequence satisfies the restrictions set by the algorithm */
				if ( l > n || l > m )
				{
        				fprintf( stderr, " Error: the fixed-length of motifs must be less or equal to the length of the sequences!\n");
                			exit ( 1 );
        			}

				if ( d == 0 )
				{
					if ( ! motifs_extraction_hd ( seqs[i], m, seqs[j], n, l, e, l_occur, l_all_occur ) )
        				{
              					fprintf( stderr, " Error: motifs_extraction_hd() failed!\n");                        
              					exit ( 1 );
        				}
        			}	
				else
				{
					if ( ! motifs_extraction_ed ( seqs[i], m, seqs[j], n, l, e, l_occur, l_all_occur ) )
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

		if ( background_filename == NULL )
                	write_motifs ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
		else
                	write_motifs_back ( sw, num_seqs, seqs, g_occur, g_all_occur, end - start, P );
	#ifdef _USE_MPI
	}
	#endif

	/* Deallocate & Finalize */
	for ( i = 0; i < num_seqs; i ++ ) 
	{
		free ( ( void * ) g_all_occur[i] );
		free ( ( void * ) g_occur[i] );
		free ( ( void * ) seqs[i] );
	}
	free ( g_all_occur );	
	free ( g_occur );	
	free ( seqs );	

	free ( sw . input_filename );
   	free ( sw . output_filename );
   	free ( sw . alphabet );

	#ifdef _USE_MPI
	 MPI_Finalize();
	#endif

	return ( 0 );
}
