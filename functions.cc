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

#include <time.h>
#include <sys/time.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <datrie/trie.h>
#include <libflasm.h>
#include "motexdefs.h"

using namespace std;
using namespace libflasm;

#if 0
#ifdef _USE_MPI
/*
returns the ith coordinate of element x at the diagonal step
*/
inline unsigned int i_coord ( unsigned int step, unsigned int n, unsigned int m, unsigned int x )
{
	if ( step <= n )		return step - x;
	else				return n - x;
}

/*
returns the jth coordinate of element x at the diagonal step
*/
inline unsigned int j_coord ( unsigned int step, unsigned int n, unsigned int m, unsigned int x )
{
        if ( step <= n )                return x;
        else		                return x + step - n;
}

/*
returns the coordinate of the upper left cell at the diagonal step
*/
inline unsigned int upper_left ( unsigned int step, unsigned int n, unsigned int j )
{
	if ( step <= n )		return j - 1;
	else if ( step == n + 1 )	return j;
	else  /*if ( step > n + 1 )*/	return j + 1;
}

/*
returns the coordinate of the upper cell at the diagonal step
*/
inline unsigned int upper ( unsigned int step, unsigned int n, unsigned int j )
{
	if ( step <= n )                return j;
        else /*if ( step >= n + 1 )*/       return j + 1;
}

/*
returns the coordinate of the left cell at the diagonal step
*/
inline unsigned int left ( unsigned int step, unsigned int n, unsigned int j )
{
        if ( step <= n )                return j - 1;
        else /*if ( step >= n + 1 )*/       return j;

}
/*
returns the number of elements at the diagonal step 
*/
inline unsigned int elements( unsigned int step, unsigned int n, unsigned int m )
{
	if ( m <= n )
	{
		if ( step <= m )	return ( step + 1 );
        	else if ( step < n)	return ( m + 1 );
        	else	                return ( m + 1 - step + n );
	}
	else
	{
		if ( step <= n )	return ( step + 1 );
        	else if ( step < m)	return ( n + 1 );
        	else	                return ( m + 1 - step + n );
	}
}

/*
allocates the amount of cells for each processor at the diagonal step
*/
inline void allocation ( int rank, int P, unsigned int n, unsigned int m, unsigned int step, int *first, int* last )
{
	unsigned int d = elements ( step, n, m ); 
	unsigned int modulo = d % P;

	if( rank < modulo )
	{
		( *first ) = rank * ( ceil ( (double) ( d ) / P ) );
		( *last ) = ( *first ) + ceil ( (double) ( d ) / P ) - 1;
	}
	else
	{
		( *first ) = rank * ( floor ( (double) ( d ) / P ) ) + modulo;
		( *last ) = ( *first ) + floor ( (double) ( d ) / P ) - 1;
	}
}

/*
allocates the workload for each processor given the total number of jobs m
*/
inline void vec_allocation ( int rank, int P, unsigned int m, int *first, int* last, int* count )
{
	unsigned int d = m; 
	unsigned int modulo = d % P;

	if( rank < modulo )
	{
		( *first ) = rank * ( ceil ( (double) ( d ) / P ) );
		( *last ) = ( *first ) + ceil ( (double) ( d ) / P ) - 1;
	}
	else
	{
		( *first ) = rank * ( floor ( (double) ( d ) / P ) ) + modulo;
		( *last ) = ( *first ) + floor ( (double) ( d ) / P ) - 1;
	}

	( *count ) = ( *last ) - ( *first ) + 1;
}

/*
handles the p2p boundary cells swaps between the processors
*/
inline void communication ( int P, int rank, unsigned int step, unsigned int n, unsigned int m, WORD1 * D, int first, int last, int first_n, int last_n)
{
	MPI_Status status;
	
        unsigned int a = elements( step, n, m );
	unsigned int b = elements( step + 1, n, m );
	unsigned int sLj = first;
	unsigned int sRj = last;
	unsigned int rLj = sLj - 1;
        unsigned int rRj = sRj + 1;
	unsigned int prev_rank_smaller_b = ( rank - 1 < b ) ? 1 : 0;
	unsigned int next_rank_smaller_b = ( rank + 1 < b ) ? 1 : 0;
	unsigned int prev_rank_smaller_a = ( rank - 1 < a ) ? 1 : 0;
	unsigned int next_rank_smaller_a = ( rank + 1 < a ) ? 1 : 0;
	unsigned int not_last_rank = ( rank != P - 1 ) ? 1 : 0;
	unsigned int not_first_rank = ( rank != 0 ) ? 1 : 0;

	if ( rank < a ) //send
	{
		if ( step < n )
		{
			if ( ! ( first_n <= last + 1 && last + 1 <= last_n ) && not_last_rank && next_rank_smaller_b ) //to send right
			{
				MPI_Send( &D[sRj] , 1 , MPI_UNSIGNED, rank + 1, 0 , MPI_COMM_WORLD );
		       		//printf ( "\nP[%d] sends D[%d] to P[%d] in step %d", rank , sRj , rank + 1 , step );
			}
			if ( ! (first_n <= first && first <= last_n ) && not_first_rank && prev_rank_smaller_b )	//to send left
			{
				MPI_Send( &D[sLj] , 1 , MPI_UNSIGNED , rank - 1 , 0 , MPI_COMM_WORLD);
		               	//printf ( "\nP[%d] sends D[%d] to P[%d] in step %d", rank , sLj , rank - 1 , step );
			}
		}
		else    // step: n..n+m-1
		{
			if ( ! (first_n <= first - 1 && first - 1 <= last_n ) && not_first_rank && prev_rank_smaller_b )//to send left	
			{
				MPI_Send( &D[sLj] , 1 , MPI_UNSIGNED , rank - 1 , 0 , MPI_COMM_WORLD);
		               	//printf ( "\nP[%d] sends D[%d] to P[%d] in step %d", rank , sLj , rank - 1 , step );
			}	
			if ( ! ( first_n <= last && last <= last_n ) && not_last_rank && next_rank_smaller_b ) //to send right
			{
				MPI_Send( &D[sRj] , 1 , MPI_UNSIGNED, rank + 1, 0 , MPI_COMM_WORLD );
		         	//printf ( "\nP[%d] sends D[%d] to P[%d] in step %d", rank , sRj , rank + 1 , step );
			}
		}
        }

	if ( rank < b )
	{
		if ( step < n  )
		{
			if ( ! (first <= first_n - 1 && first_n - 1 <= last) && not_first_rank && prev_rank_smaller_a ) //to receive from left
			{
		           	MPI_Recv( &D[rLj] , 1 , MPI_UNSIGNED , rank - 1 , 0 , MPI_COMM_WORLD, &status);
		               	//printf("\nP[%d] receives D[%d] from P[%d] in step %d" , rank , rLj , rank - 1 , step );
			}
			if ( ! (first <= last_n && last_n <= last) && not_last_rank && next_rank_smaller_a ) //receive from right
			{
				MPI_Recv( &D[rRj], 1, MPI_UNSIGNED, rank + 1, 0, MPI_COMM_WORLD, &status);
				//printf ( "\nP[%d] receives D[%d] from P[%d] in step %d", rank ,  rRj , rank + 1 , step );
			}
		}	
		else    // step: n..n+m-1
		{
			if ( ! (first <= last_n + 1 && last_n + 1 <= last) && not_last_rank && next_rank_smaller_a ) //receive from right
			{
				MPI_Recv( &D[rRj], 1, MPI_UNSIGNED, rank + 1, 0, MPI_COMM_WORLD, &status);
				//printf ( "\nP[%d] receives D[%d] from P[%d] in step %d", rank ,  rRj , rank + 1 , step );
			}
			if ( ! (first <= first_n && first_n <= last) && not_first_rank && prev_rank_smaller_a ) //to receive from left
			{
				MPI_Recv( &D[rLj] , 1 , MPI_UNSIGNED , rank - 1 , 0 , MPI_COMM_WORLD, &status);
		           	//printf("\nP[%d] receives D[%d] from P[%d] in step %d" , rank , rLj , rank - 1 , step );
			}
		}
	}

}

/*
the dynamic programming algorithm under the hamming distance model for one (or a few) very large sequence(s) for MPI
*/
unsigned int motifs_extraction_opasm_hd ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v, int rank, int P )
{
	unsigned int step;
	int first, last;
	int first_n;
	int last_n;
	int x;

	WORD1 y; 
	WORD1 *D0, *D1, *D2; 	//matrix B in anti-diagonal vectors
	unsigned int i;
	unsigned int j;
	unsigned int ul;

  	/* Memory Allocation */
	if ( ( D0 = ( WORD1 * ) calloc ( ( m + 1 ) , sizeof( WORD1 ) ) ) == NULL )
	{
		fprintf( stderr, " Error: D0 could not be allocated!\n");
		return ( 0 );
	}

	if ( ( D1 = ( WORD1 * ) calloc ( ( m + 1 ) , sizeof( WORD1 ) ) ) == NULL )
	{
		fprintf( stderr, " Error: D1 could not be allocated!\n");
		return ( 0 );
	}

	if ( ( D2 = ( WORD1 * ) calloc ( ( m + 1 ) , sizeof( WORD1 ) ) ) == NULL )
	{
		fprintf( stderr, " Error: D2 could not be allocated!\n");
		return ( 0 );
	}

	y = ( ( WORD1 ) ( 1 ) << ( sw . l - 1 ) ) - 1;

        step = 0;
	allocation ( rank, P , n , m, step , &first , &last);		//Processor Allocation
	allocation ( rank, P , n , m, step + 1 , &first_n , &last_n);   //Processor Allocation

	while ( step < n + m + 1 ) 
	{  
		#if 1
		for( x = ( first ); x <= ( last ) ; x++ )			
		{
			i = i_coord ( step, n, m, x );
			j = j_coord ( step, n, m, x );
			if( j == 0 ) continue;
			ul = upper_left ( step, n, x );

			switch ( step % 3 )
			{
				case 0 :

				if ( i == 0 )
					D0[x] = ( ( WORD1 ) 2 << ( min ( j , sw . l ) - 1 ) ) - 1;
				else
					D0[x] = shiftc ( D1[ul], y ) | delta1 ( t[i - 1], p[j - 1] );

				if ( popcount ( D0[x] ) <= sw . e && j >= sw . l )	
				{
					v[j - 1] = v[j - 1] + 1; 
					u[j - 1] = u[j - 1] + 1; 
				}

				break;
				
				case 1 :
				
				if( i == 0 )
					D1[x] = ( ( WORD1 ) 2 << ( min ( j , sw . l ) - 1 ) ) - 1;
				else
					D1[x] = shiftc ( D2[ul], y ) | delta1 ( t[i - 1], p[j - 1] );

				if ( popcount ( D1[x] ) <= sw . e && j >= sw . l )	
				{
					v[j - 1] = v[j - 1] + 1; 
					u[j - 1] = u[j - 1] + 1; 
				}

				break;
			
				case 2 :
 
				if( i == 0 )
					D2[x] = ( ( WORD1 ) 2 << ( min ( j , sw . l ) - 1 ) ) - 1;
				else
					D2[x] = shiftc ( D0[ul], y ) | delta1 ( t[i - 1], p[j - 1] );

				if ( popcount ( D2[x] ) <= sw . e && j >= sw . l )	
				{
					v[j - 1] = v[j - 1] + 1; 
					u[j - 1] = u[j - 1] + 1; 
				}

				break;
				
				default:
				fprintf ( stderr, " Error: this should never happen!\n");
				break;		
			} 
			
		}
		
		if ( step < n + m )
		{
			switch ( step % 3 )
			{
	      			case 0 :
	           		communication ( P, rank, step, n, m, D0, first, last, first_n, last_n );
				break;
		
	      			case 1 :
		    		communication ( P, rank, step, n, m, D1, first, last, first_n, last_n );
				break;	
	
	      			case 2 :
		         	communication ( P, rank, step, n, m, D2, first, last, first_n, last_n );
				break;
		
				default:
				fprintf ( stderr, " Error: this should never happen!\n");
				break;	
			}	
		}
		#endif
		step++;
		first = first_n; last = last_n;
		allocation ( rank, P , n , m, step + 1 , &first_n , &last_n);   //Processor Allocation
		
	}

	free ( D0 );
	free ( D1 );
	free ( D2 );

	return  ( 1 );
}

/*
the dynamic programming algorithm under the hamming distance model for one (or a few) very large sequence(s) for MPI
*/
unsigned int motifs_extraction_opasm_ed ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v, int rank, int P )
{
	unsigned int step;
	int first, last;
	int first_n;
	int last_n;
	int x;

	WORD1 y; 
	WORD1 *D0, *D1, *D2; 	//matrix B in anti-diagonal vectors
	
	unsigned int i;
	unsigned int j;
	unsigned int ul;
	unsigned int up;
	unsigned int lft;

  	// Memory Allocation 
	if ( ( D0 = ( WORD1 * ) calloc ( ( m + 1 ) , sizeof( WORD1 ) ) ) == NULL )
	{
		fprintf( stderr, " Error: D0 could not be allocated!\n");
		return ( 0 );
	}

	if ( ( D1 = ( WORD1 * ) calloc ( ( m + 1 ) , sizeof( WORD1 ) ) ) == NULL )
	{
		fprintf( stderr, " Error: D1 could not be allocated!\n");
		return ( 0 );
	}

	if ( ( D2 = ( WORD1 * ) calloc ( ( m + 1 ) , sizeof( WORD1 ) ) ) == NULL )
	{
		fprintf( stderr, " Error: D2 could not be allocated!\n");
		return ( 0 );
	}

	y = ( ( WORD1 ) ( 1 ) << ( sw . l - 1 ) ) - 1;
        step = 0;
	allocation ( rank, P , n , m, step , &first , &last);		//Processor Allocation
	allocation ( rank, P , n , m, step + 1 , &first_n , &last_n);   //Processor Allocation

	while ( step < n + m + 1 ) 
	{  
		for( x = ( first ); x <= ( last ) ; x++ )			
		{
			i = i_coord ( step, n, m, x );
			j = j_coord ( step, n, m, x );
			if( j == 0 ) continue;
			ul = upper_left ( step, n, x );
			up = upper ( step, n, x );
			lft = left ( step, n, x );

			switch ( step % 3 )
			{
				case 0 :

				if ( i == 0 )
					D0[x] = ( ( WORD1 ) 2 << ( min ( j , sw . l ) - 1 ) ) - 1;
				else if ( j <= sw . l )
					D0[x] = bitminmax ( shift ( D2[up] ) | 1, shift( D2[lft] ) | 1, shift( D1[ul] ) | delta1 ( t[i - 1], p[j - 1] ) );
				else
					D0[x] = bitminmax ( shiftc( D2[lft] , y ) | 1, shift ( D2[up] ) | 1, shiftc ( D1[ul], y ) | delta1 ( t[i - 1], p[j - 1] ) );
				if ( popcount ( D0[x] ) <= sw . e && j >= sw . l )	
				{
					v[j - 1] = v[j - 1] + 1; 
					u[j - 1] = u[j - 1] + 1; 
				}
				break;
				
				case 1 :
				
				if( i == 0 )
					D1[x] = ( ( WORD1 ) 2 << ( min ( j , sw . l ) - 1 ) ) - 1;
				else if ( j <= sw . l )
					D1[x] = bitminmax ( shift ( D0[up] ) | 1, shift( D0[lft] ) | 1, shift( D2[ul] ) | delta1 ( t[i - 1], p[j - 1] ) );
				else
					D1[x] = bitminmax ( shiftc( D0[lft] , y ) | 1, shift ( D0[up] ) | 1, shiftc ( D2[ul], y ) | delta1 ( t[i - 1], p[j - 1] ) );
				if ( popcount ( D1[x] ) <= sw . e && j >= sw . l )	
				{
					v[j - 1] = v[j - 1] + 1; 
					u[j - 1] = u[j - 1] + 1; 
				}
				break;
			
				case 2 :
 
				if( i == 0 )
					D2[x] = ( ( WORD1 ) 2 << ( min ( j , sw . l ) - 1 ) ) - 1;
				else if ( j <= sw . l )
					D2[x] = bitminmax ( shift ( D1[up] ) | 1, shift( D1[lft] ) | 1, shift( D0[ul] ) | delta1 ( t[i - 1], p[j - 1] ) );
				else
					D2[x] = bitminmax ( shiftc( D1[lft] , y ) | 1, shift ( D1[up] ) | 1, shiftc ( D0[ul], y ) | delta1 ( t[i - 1], p[j - 1] ) );
				if ( popcount ( D2[x] ) <= sw . e && j >= sw . l  )	
				{
					v[j - 1] = v[j - 1] + 1; 
					u[j - 1] = u[j - 1] + 1; 
				}
				break;
				
				default:
				fprintf ( stderr, " Error: this should never happen!\n");
				break;		
			} 
		}
		
		if ( step < n + m )
		{
			switch ( step % 3 )
			{
	      			case 0 :
	           		communication ( P, rank, step, n, m, D0, first, last, first_n, last_n );
				break;
		
	      			case 1 :
		    		communication ( P, rank, step, n, m, D1, first, last, first_n, last_n );
				break;	
	
	      			case 2 :
		         	communication ( P, rank, step, n, m, D2, first, last, first_n, last_n );
				break;
		
				default:
				fprintf ( stderr, " Error: this should never happen!\n");
				break;	
			}	
		}

		step++;
		first = first_n; last = last_n;
		allocation ( rank, P , n , m, step + 1 , &first_n , &last_n);   //Processor Allocation
		
	}

	free ( D0 );
	free ( D1 );
	free ( D2 );

	return  ( 1 );
}
#endif
#endif

/*
the dynamic programming algorithm under the hamming distance model
*/
unsigned int motifs_extraction_hd ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v )
{
	unsigned int * occ;

  	/* Memory Allocation */
	if ( ( occ = ( unsigned int* ) calloc ( ( m ) , sizeof( unsigned int ) ) ) == NULL )
	{
		fprintf( stderr, " Error: occ could not be allocated!\n");
		return ( 0 );
	}


	ResultTupleSet pattern_matching_hd = flasm_hd( ( unsigned char * ) t, n, ( unsigned char * ) p, m, sw . l, sw . e, true);
	
	ResultTupleSetIterator it;
	
	for ( it = pattern_matching_hd.begin(); it != pattern_matching_hd.end(); it++ )
	{
		ResultTuple result = *it;

		if ( result . error <= sw . e )	
		{
			if ( occ[result . pos_x] == 0 )
			{
				u[result . pos_x] = u[result . pos_x] + 1;
				occ[result . pos_x] = 1;
			} 
			v[result . pos_x] = v[result . pos_x] + 1; 
		}
	}
	
	free ( occ );
	return  ( 1 );
}

/*
the dynamic programming algorithm under the hamming distance model for structured motifs extraction
*/
unsigned int structured_motifs_extraction_hd ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v )
{
	unsigned int *** B; 		//matrices B
	unsigned int ** occ;
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int b;
	unsigned int g;

	/* 3d dynamic memory allocation for matrices B*/
   	if( ( B = ( unsigned int *** ) malloc ( ( sw . nb_boxes ) * sizeof( unsigned int ** ) ) ) == NULL )
    	{
      		fprintf( stderr, " Error: B could not be allocated!\n" );
      		return 0;
    	} 
	
   	for ( i = 0; i < sw . nb_boxes; i ++ )
    	{
      		if( ( B[i] = ( unsigned int ** ) malloc ( ( n + 1 ) * sizeof( unsigned int * ) ) ) == NULL )
       		{
         		fprintf( stderr, " Error: B could not be allocated!\n" );
         		return 0;
       		}

		for ( j = 0; j < n + 1; j++ )
		{
			if( ( B[i][j] = ( unsigned int * ) malloc ( ( m + 1 ) * sizeof( unsigned int ) ) ) == NULL )
       			{
         			fprintf( stderr, " Error: B could not be allocated!\n" );
         			return 0;
       			}
			
			for( k = 0; k< m + 1; k++ )
				B[i][j][k] = m ;
		}
    	}

	if ( ( occ = ( unsigned int ** ) malloc ( sw . nb_structs * sizeof ( unsigned int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: R could not be allocated!\n");
                return ( 0 );
        }
	for ( j = 0; j < sw . nb_structs; ++ j )
        {	
		if ( ( occ[j] = ( unsigned int * ) calloc ( ( m + 1 ) , sizeof ( unsigned int ) ) ) == NULL )
        	{
                	fprintf( stderr, " Error: R could not be allocated!\n");
                	return ( 0 );
        	}
	}

	/* Computing the DP matrices */
	for ( b = 0; b < sw . nb_boxes; b++ ) 
	{
		unsigned int ell;
		unsigned int err;

		if ( b != 0 )
		{	
			ell = sw . blens[b - 1];
			err = sw . berrs[b - 1];
		}
		else
		{
			ell = sw . l;
			err = sw . e;
		}

		ResultTupleSet pattern_matching_hd = flasm_hd( ( unsigned char * ) t, n, ( unsigned char * ) p, m, ell, err, true );
	
		ResultTupleSetIterator it;
	
		for ( it = pattern_matching_hd.begin(); it != pattern_matching_hd.end(); it++ )
		{
			ResultTuple result = *it;

			B[b][result.pos_t+1][result.pos_x+1] = result . error;
		}
	}

	/* Merging the single-motif occurrences into structured-motif occurrences */
	
	for ( i = 0; i < n + 1; i++ ) 
	{
		for( j = 0; j < m + 1 ; j++ )			
		{
			if ( B[0][i][j] <= sw . e )
			{
				for ( k = 0; k < sw . nb_structs; k ++ )
				{
					unsigned int l;
					for ( l = 0; l < sw . nb_structs; l ++ )
					{
						unsigned int array_offset = k * m;
						unsigned int dist_offset_j = j;
						unsigned int dist_offset_i = i;
						unsigned int r = 1;
						for ( g = 0; g < sw . nb_gaps; g++ )
						{
							dist_offset_j += ( sw . S[k][g] + sw . blens[g] );
							dist_offset_i += ( sw . S[l][g] + sw . blens[g] );
							
							if ( dist_offset_i < n + 1 && dist_offset_j < m + 1 )
							{
								if ( B[g + 1][dist_offset_i][dist_offset_j] <= sw . berrs[g] )
								{
									r++;
									if ( r == sw . nb_boxes )
									{
										if ( occ[k][j - 1] == 0 )
										{
											u[ array_offset + j - 1] = u[ array_offset + j - 1] + 1;
											occ[k][j- 1] = 1;
										} 
										v[ array_offset + j - 1] = v[ array_offset + j - 1] + 1;
									}
								}	
								else
									break;
							}
						}
					}
				}
			}
		}
	}

	for ( i = 0;  i < sw . nb_boxes; i++ )
	{
		for ( j = 0; j < n + 1; j++ )
			free( B[i][j] );
      		free ( B[i] );
	}
   	free ( B );

	for ( j = 0; j < sw . nb_structs; ++ j )
		free( occ[j] );
	free ( occ );

	return  ( 1 );
}

/*
A recursive function that creates the cartesian product of multiple (sets) arrays
*/
void recurse( unsigned int x, unsigned int xi, unsigned int * y, unsigned int ** A, unsigned int * newx, unsigned int * temp_inner, unsigned int ** B ) 
{
	unsigned int i;	
	if ( xi < x ) 
	{
    		for ( i = 0; i < y[xi]; i++ ) 
		{
      			temp_inner[xi] = A[xi][i];
      			recurse( x, xi + 1, y, A, newx, temp_inner, B );
    		}
  	} 
	else 
	{
    		for ( i = 0; i < x; i++ ) 
		{
      			B[*newx][i] = temp_inner[i];
    		}
    		( *newx )++;
  	}
}

/*
the dynamic programming algorithm under the edit distance model
*/
unsigned int motifs_extraction_ed ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v )
{

	unsigned int * occ;	

  	// Memory Allocation 
	if ( ( occ = ( unsigned int* ) calloc ( ( m ) , sizeof( unsigned int ) ) ) == NULL )
	{
		fprintf( stderr, " Error: occ could not be allocated!\n");
		return ( 0 );
	}

	ResultTupleSet pattern_matching_ed = flasm_ed( ( unsigned char * ) t, n, ( unsigned char * ) p, m, sw . l, sw . e, true );
	
	ResultTupleSetIterator it;
	
	for ( it = pattern_matching_ed.begin(); it != pattern_matching_ed.end(); it++ )
	{
		ResultTuple result = *it;

		if ( result . error <= sw . e )	
		{
			if ( occ[result . pos_x ] == 0 )
			{
				u[result . pos_x ] = u[result . pos_x ] + 1;
				occ[result . pos_x ] = 1;
			} 
			v[result . pos_x ] = v[result . pos_x ] + 1; 
	
		}
	}
	
	free ( occ );
	return  ( 1 );
}

/*
the dynamic programming algorithm under the edit distance model for structured motifs extraction
*/
unsigned int structured_motifs_extraction_ed ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v )
{
	unsigned int *** B; 		//matrices B
	unsigned int ** occ;
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int b;
	unsigned int g;

	/* 3d dynamic memory allocation for matrices B*/
		if( ( B = ( unsigned int *** ) malloc ( ( sw . nb_boxes ) * sizeof( unsigned int ** ) ) ) == NULL )
    	{
      		fprintf( stderr, " Error: B could not be allocated!\n" );
      		return 0;
    	} 
	
   	for ( i = 0; i < sw . nb_boxes; i ++ )
    	{
      		if( ( B[i] = ( unsigned int ** ) malloc ( ( n + 1 ) * sizeof( unsigned int * ) ) ) == NULL )
       		{
         		fprintf( stderr, " Error: B could not be allocated!\n" );
         		return 0;
       		}

		for ( j = 0; j < n + 1; j++ )
		{
			if( ( B[i][j] = ( unsigned int * ) malloc ( ( m + 1 ) * sizeof( unsigned int ) ) ) == NULL )
       			{
         			fprintf( stderr, " Error: B could not be allocated!\n" );
         			return 0;
       			}
			
			for( k = 0; k< m + 1; k++ )
				B[i][j][k] = m;
		}
    	}

	if ( ( occ = ( unsigned int ** ) malloc ( sw . nb_structs * sizeof ( unsigned int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: R could not be allocated!\n");
                return ( 0 );
        }
	for ( j = 0; j < sw . nb_structs; ++ j )
        {	
		if ( ( occ[j] = ( unsigned int * ) calloc ( m , sizeof ( unsigned int ) ) ) == NULL )
        	{
                	fprintf( stderr, " Error: R could not be allocated!\n");
                	return ( 0 );
        	}
	}


	/* Computing the DP matrices */
	for ( b = 0; b < sw . nb_boxes; b++ ) 
	{
		unsigned int ell;
		unsigned int err;

		if ( b != 0 )
		{	
			ell = sw . blens[b - 1];
			err = sw . berrs[b - 1];
		}
		else
		{
			ell = sw . l;
			err = sw . e;
		}

		ResultTupleSet pattern_matching_ed = flasm_ed( ( unsigned char * ) t, n, ( unsigned char * ) p, m, ell, err, true );
	
		ResultTupleSetIterator it;
	
		for ( it = pattern_matching_ed.begin(); it != pattern_matching_ed.end(); it++ )
		{
			ResultTuple result = *it;

			B[b][result.pos_t+1][result.pos_x+1] = result . error;

		}
	}

	/* Merging the single-motif occurrences into structured-motif occurrences */
	for ( i = 0; i < n + 1; i++ ) 
	{
		for( j = 0; j < m + 1 ; j++ )			
		{
			if ( B[0][i][j] <= sw . e )
			{
				for ( k = 0; k < sw . nb_structs; k ++ )
				{
					unsigned int l;
					for ( l = 0; l < sw . nb_structs; l ++ )
					{
						unsigned int array_offset = k * m;
						unsigned int dist_offset_j = j;
						unsigned int dist_offset_i = i;
						unsigned int r = 1;
						for ( g = 0; g < sw . nb_gaps; g++ )
						{
							dist_offset_j += ( sw . S[k][g] + sw . blens[g] );
							dist_offset_i += ( sw . S[l][g] + sw . blens[g] );
							
							if ( dist_offset_i < n + 1 && dist_offset_j < m + 1 )
							{
								if ( B[g + 1][dist_offset_i][dist_offset_j] <= sw . berrs[g] )
								{
									r++;
									if ( r == sw . nb_boxes )
									{
										if ( occ[k][j - 1] == 0 )
										{
											u[ array_offset + j - 1] = u[ array_offset + j - 1] + 1;
											occ[k][j- 1] = 1;
										} 
										v[ array_offset + j - 1] = v[ array_offset + j - 1] + 1;
									}
								}	
								else
									break;
							}
						}
					}
				}
			}
		}
	}

	for ( i = 0;  i < sw . nb_boxes; i++ )
	{
		for ( j = 0; j < n + 1; j++ )
			free( B[i][j] );
      		free ( B[i] );
	}
   	free ( B );

	for ( j = 0; j < sw . nb_structs; ++ j )
		free( occ[j] );
	free ( occ );

	return  ( 1 );
}

#if 0
/*
given integers a, b, c this function returns one of the integers a, b, c
with the property that it has the least number of 1's (bits set on). If there is 
a draw then it returns the maximum of the two when viewed as decimal integers
*/
inline WORD1 bitminmax ( WORD1 a, WORD1 b, WORD1 c )
{
        unsigned int x , y , z; 
	WORD1 minimum, maximum;

        x = popcount ( a );
        y = popcount ( b );
        z = popcount ( c );

        minimum = min( x , y );

        if( z < minimum )
        {
                return c;
        }
        else if ( z == minimum && ( x != y ) )
        {
                if( x < y )
                        maximum = max ( c , a );
                else
                        maximum = max ( c , b );

                return maximum;
        }
        else if ( ( z == x ) && ( x == y ) )
        {
                maximum = max ( b , c );
                maximum = max ( a , maximum );
                return maximum;
        }
        else if ( z > minimum &&  ( x != y ) )
        {
                if( x < y )
                        return a;
                else
                        return b;

        }
        else /*if ( z > minimum && ( x == y ) )*/
        {
                maximum = max ( a , b );
                return maximum;
        }

}


/*
moves the bits one position to the left and enters zeros from the right
*/
inline WORD1 shift ( WORD1 a ) 
{
	return ( ( WORD1 ) a << 1 );
}

/*
shifts and truncates the leftmost bit
*/
inline WORD1 shiftc( WORD1 a, WORD1 x ) 
{
	return shift ( a & x );
}

/*
returns the hamming distance for two characters a and b 
*/
inline unsigned int delta1( char a, char b )
{
	if	( a == b )			return 0;
	else					return 1;
}

/*
returns the number of 1's in an integer x
*/
inline unsigned int popcount ( WORD1 x )
{
	//__builtin_popcount: 32-bit
	//__builtin_popcountll: 128-bit
        //http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html
	return __builtin_popcountl( x );
}
#endif

/*
write the header of the output
*/
unsigned int write_motex_header ( FILE * out_fd, struct TSwitch sw, unsigned int num_seqs, double exectime, int P )
{
	time_t               t;
   	time ( &t );
	unsigned int         i;

	if ( sw . nb_boxes == 1 )
	{
		fprintf ( out_fd, "####################################\n" );
		fprintf ( out_fd, "# Program: MoTeX\n" );
		fprintf ( out_fd, "# Rundate: %s", ctime ( &t ) );
		fprintf ( out_fd, "# Input file: %s\n", sw . input_filename );
		fprintf ( out_fd, "# Output file: %s\n", sw . output_filename );
		fprintf ( out_fd, "# For %d input sequences\n", num_seqs );
		fprintf ( out_fd, "#     d = %d\n", sw . d );
		fprintf ( out_fd, "#     q = %d\n", sw . q );
		if ( sw . Q != 100 )
		fprintf ( out_fd, "#     Q = %d\n", sw . Q );
		if ( sw . n != 1 )
		fprintf ( out_fd, "#     n = %d\n", sw . n );
		if ( sw . N != 10000 )
		fprintf ( out_fd, "#     N = %d\n", sw . N );
		fprintf ( out_fd, "#     l = %d\n", sw . l ) ;
		fprintf ( out_fd, "#     e = %d\n", sw . e );
		fprintf ( out_fd, "# Run on %d proc(s) in %lf secs\n", P, exectime );
		fprintf ( out_fd, "####################################\n\n" );
	}
	else
	{
		fprintf ( out_fd, "####################################\n" );
		fprintf ( out_fd, "# Program: MoTeX\n" );
		fprintf ( out_fd, "# Rundate: %s", ctime ( &t ) );
		fprintf ( out_fd, "# Input file: %s\n", sw . input_filename );
		fprintf ( out_fd, "# Output file: %s\n", sw . output_filename );
		fprintf ( out_fd, "# For %d input sequences\n", num_seqs );
		fprintf ( out_fd, "#     d = %d\n", sw . d );
		fprintf ( out_fd, "#     q = %d\n", sw . q );
		if ( sw . Q != 100 )
		fprintf ( out_fd, "#     Q = %d\n", sw . Q );
		if ( sw . n != 1 )
		fprintf ( out_fd, "#     n = %d\n", sw . n );
		if ( sw . N != 10000 )
		fprintf ( out_fd, "#     N = %d\n", sw . N );
		fprintf ( out_fd, "# Number of boxes: %d\n", sw . nb_boxes );
		fprintf ( out_fd, "#     l0 = %d\n", sw . l ) ;
		fprintf ( out_fd, "#     e0 = %d\n", sw . e );
		for ( i = 0; i < sw . nb_gaps; i ++ )	
		{
			fprintf ( out_fd, "#     g%d = [%d,%d]\n", i + 1, sw . bgaps_min[i], sw . bgaps_max[i] ) ;
			fprintf ( out_fd, "#     l%d = %d\n", i + 1, sw . blens[i] ) ;
			fprintf ( out_fd, "#     e%d = %d\n", i + 1, sw . berrs[i] );
		}
		fprintf ( out_fd, "# Boxes input file: %s\n", sw . boxes_in_filename );
		fprintf ( out_fd, "# Run on %d proc(s) in %lf secs\n", P, exectime );
		fprintf ( out_fd, "####################################\n\n" );
	}

	return ( 1 );
}

/*
write the output
*/
unsigned int write_motifs ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P )
{
	FILE * 		out_fd;				// file with the motifs
	unsigned int i, j, k;
	unsigned int valid = 0;

	AlphaMap *	alphabet = NULL;
        _Trie *      	trie = NULL;

	if ( ( out_fd = fopen( sw . output_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

	/* Write the header */
	write_motex_header ( out_fd, sw, num_seqs, exectime, P );

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );


	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
                        if (  (  ( ( double ) u[i][j] / num_seqs ) * 100 ) >= sw . q && v[i][j] >= sw . n )
                        	if (  (  ( ( double ) u[i][j] / num_seqs ) * 100 ) <= sw . Q && v[i][j] <= sw . N )
				{
					AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( sw . l + 1 ) , sizeof( AlphaChar ) );
					for ( k = 0; k < sw. l; k ++ )
						ACmotif[k] = ( AlphaChar )seqs[i][k + j - ( sw . l ) + 1];
					ACmotif [ k ] = 0;
					if ( trie_retrieve ( trie, ACmotif, NULL ) != TRUE )
					{
						trie_store ( trie, ACmotif, 0 );

						char      * motif   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );
						memcpy ( motif, &seqs[i][j - ( sw . l ) + 1 ], sw . l );
						fprintf ( out_fd, "%s %d %d %lf %d\n", motif, u[i][j], num_seqs, (  ( ( double ) u[i][j] / num_seqs ) ), v[i][j] );
						valid ++;
						free ( motif );
					}
					free ( ACmotif );
				}
		}
	}

	fprintf ( out_fd, "\n#A total of %d valid motifs were extracted.\n\n", valid );
                
	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;

	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output for structured motifs
*/
unsigned int write_structured_motifs ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P )
{
	FILE * 		out_fd;				// file with the motifs
	unsigned int i, j, k, l;
	unsigned int valid = 0;

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

	if ( ( out_fd = fopen( sw . output_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

	/* Write the header */
	write_motex_header ( out_fd, sw, num_seqs, exectime, P );

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	/* Calculate the actual length of the structured motif --- excluding the gaps */
	unsigned int total_length_sm = sw . l;
	unsigned int b;
	for ( b = 0; b < sw . nb_gaps; b ++ )	
		total_length_sm += sw . blens[b];

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
                	for ( k = 0; k < sw . nb_structs; k++ )
			{
				if (  (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100 ) >= sw . q && v[i][k * m + j] >= sw . n )
					if (  (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100 ) <= sw . Q && v[i][k * m + j] <= sw . N )
					{
						AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( AlphaChar ) );  //AAAA_ATAT_TATTT
						char      * motif   = ( char * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( char ) );
						unsigned int ell;
						unsigned int offset;
						unsigned int index = 0;
						unsigned int jj = j - ( sw . l  ) + 1;

						for ( b = 0; b < sw . nb_boxes; b ++ ) 
						{
							/* Extract the structured motif from the sequences */
							if ( b == 0 )	
							{
								offset = 0;
								ell = sw . l;
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + l];
									motif[index]   = seqs[i][jj + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}	
							else
							{
								offset += sw . S[k][b - 1] + ell;
								ell = sw . blens[b - 1];
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + offset + l];
									motif[index]   = seqs[i][jj + offset + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}
						}

						index--;
						ACmotif [ index ] = 0;
						motif [ index ] = 0;

						/* Insert the extracted structured motif into the trie */
						if ( trie_retrieve ( trie, ACmotif, NULL ) != TRUE )
						{
							trie_store ( trie, ACmotif, 0 );
							fprintf ( out_fd, "%s %d %d %lf %d\n", motif, u[i][k * m + j], num_seqs, (  ( ( double ) u[i][k * m + j] / num_seqs ) ), v[i][k * m + j] );
							valid ++;
						}
						free ( ACmotif );
						free ( motif );
					}
                        }
		}
	}

	fprintf ( out_fd, "\n#A total of %d valid structured motifs were extracted.\n\n", valid );
                
	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;

	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the header of the output a la SMILE
*/
unsigned int write_smile_header ( FILE * out_fd, struct TSwitch sw, unsigned int num_seqs, char * alphabet_str )
{
	unsigned int i;
	char line1[LINE_SIZE];

	if ( sw . nb_boxes == 1 )
	{
		// %%% 1 120/1200 1200000 8 8 1 alphabet ACGT$
		sprintf ( line1, "%%%%%% 1 %d/%d %d %d %d %d alphabet %s$", 
					( int ) ceil ( ( double ) ( ( double )sw . q * num_seqs ) / ( double ) 100 ), 
					num_seqs, 
					sw . total_length, 
					sw . l, 
					sw . l, 
					sw . e, 
					alphabet_str );
	}
   	else
	{
		//nb of boxes, absolute quorum/number of sequences, total number of symbols in the input sequences, total min length, 
		//total max length, total substitutions, description given for each box in the input parameters, and the alphabet)
		//%%% 2 128/1062 196736 12 12 2 6 6 1 16 18 6 6 1 alphabet ACGT$
		unsigned int min_len = sw . l;
		unsigned int max_len = sw . l;
		unsigned int total_errs = sw . e;
		unsigned int total_length_sm = sw . l;

		unsigned int b;
		for ( b = 0; b < sw . nb_gaps; b ++ )
		{
			total_length_sm += sw . blens[b];
			min_len += sw . blens[b];
			max_len += sw . blens[b];
			total_errs += sw . berrs[b];
		} 		

		line1[0] = 0;
		sprintf ( line1, "%%%%%% %d %d/%d %d %d %d %d %d %d %d ", 
				 sw . nb_boxes, 
				 ( int ) ceil ( ( double ) ( ( double )sw . q * num_seqs ) / ( double ) 100 ), 
				 num_seqs, 
				 sw . total_length, 
				 min_len, 
				 max_len, 
				 total_errs, 
				 sw . l, 
				 sw . l, 
				 sw . e );

		char line2[LINE_SIZE];
		line2[0] = 0;
		for ( b = 0; b < sw . nb_gaps; b ++ )
		{
			char buffer[LINE_SIZE];
			buffer[0] = 0;
			sprintf( buffer, "%d %d %d %d %d ", sw . bgaps_min[b], sw . bgaps_max[b], sw . blens[b], sw . blens[b], sw . berrs[b] );
			strncat( line2, buffer, strlen ( buffer ) );
		}

		char line3[LINE_SIZE];
		line3[0] = 0;
		sprintf ( line3, "alphabet %s$", alphabet_str);

		strcat ( line1, line2 );
		strcat ( line1, line3 );
	}

     	fprintf ( out_fd, "%s", line1 );

	if      ( ! strcmp ( "DNA", sw . alphabet ) )
	{
        	for( i = strlen ( line1 ); i < 79; i++ ) 
     			fprintf ( out_fd, " " );

		// two empty lines
		fprintf ( out_fd, "\n" );
        	for( i = 0; i != 2; i++ ) 
    		{ 
    			fprintf( out_fd, "                                        " ); 
    			fprintf( out_fd, "                                       \n" ); 
    		}
	}
 
        if ( ! strcmp ( "PROT", sw . alphabet ) )
	{
        	for( i = strlen ( line1 ); i < 159; i++ ) 
     			fprintf ( out_fd, " " );

		// one empty line
		fprintf ( out_fd, "\n" );
        	for( i = 0; i != 1; i++ ) 
    		{ 
    			fprintf( out_fd, "                                        " ); 
    			fprintf( out_fd, "                                       \n" ); 
    		}
	}
	// one line with `='
	fprintf( out_fd, "========================================" ); 
	fprintf( out_fd, "=======================================\n" );
	return ( 1 );
}

/*
write the output a la SMILE
*/
unsigned int write_motifs_smile ( struct TSwitch sw, unsigned int num_seqs, char const ** seqs, unsigned int ** u, unsigned int ** v, double exectime )
{
	FILE * 		out_fd;				// file with the motifs
	unsigned int i, j, k;
	unsigned int valid = 0;
        char * alphabet_str = NULL;

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

	if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet_str = ( char * ) DNA;
        else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet_str = ( char * ) PROT;
        else if ( ! strcmp ( "USR", sw . alphabet ) )   alphabet_str = ( char * ) USR;

	if ( ( out_fd = fopen( sw . smile_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

	write_smile_header ( out_fd, sw, num_seqs, alphabet_str );

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
                        if (  (  ( ( double ) u[i][j] / num_seqs ) * 100 ) >= sw . q && v[i][j] >= sw . n )
                        	if (  (  ( ( double ) u[i][j] / num_seqs ) * 100 ) <= sw . Q && v[i][j] <= sw . N )
				{
					AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( sw . l + 1 ) , sizeof( AlphaChar ) );
					for ( k = 0; k < sw. l; k ++ )
						ACmotif[k] = ( AlphaChar )seqs[i][k + j - ( sw . l ) + 1];
					ACmotif [ k ] = 0;
					if ( trie_retrieve ( trie, ACmotif, NULL ) != TRUE )
					{
						trie_store ( trie, ACmotif, 0 );

						char      * motif   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );
						memcpy ( motif, &seqs[i][j - ( sw . l ) + 1 ], sw . l );

						char      * id   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );

						unsigned int ii, jj;
						for ( ii = 0; ii < sw . l; ii++ )
						{
							for ( jj = 0; jj < strlen ( alphabet_str ); jj++ )
								if ( motif[ii] == alphabet_str[jj] ) break;

							char * s   = ( char * ) calloc ( ( sw . l ) , sizeof( char ) );
							sprintf( s, "%c", ( char ) ( '0' + jj ) );
							id[ii] = s[0] ;
							free ( s );
						}
						id[sw . l] = '\0';

						fprintf ( out_fd, "%s %s %d\t%d\n", motif, id, u[i][j], v[i][j] );

						valid ++;
						free ( motif );
						free ( id );
					}
					free ( ACmotif );
				}
		}
	}

	fprintf ( out_fd, "Nb models: %d\n", valid );
	fprintf ( out_fd, "User time: %lf sec.\n", exectime );
                
	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;

	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output a la SMILE
*/
unsigned int write_structured_motifs_smile ( struct TSwitch sw, unsigned int num_seqs, char const ** seqs, unsigned int ** u, unsigned int ** v, double exectime )
{
	FILE * 		out_fd;				// file with the motifs
	unsigned int i, j, k, l;
	unsigned int valid = 0;

        char * alphabet_str = NULL;

	if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet_str = ( char * ) DNA;
        else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet_str = ( char * ) PROT;
        else if ( ! strcmp ( "USR", sw . alphabet ) )   alphabet_str = ( char * ) USR;

	if ( ( out_fd = fopen( sw . smile_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

	write_smile_header ( out_fd, sw, num_seqs, alphabet_str );

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	unsigned int b;
	unsigned int total_length_sm = sw . l;
	for ( b = 0; b < sw . nb_gaps; b ++ )
		total_length_sm += sw . blens[b];

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
			for ( k = 0; k < sw . nb_structs; k ++ )
			{
				if (  (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100 ) >= sw . q && v[i][k * m + j] >= sw . n )
					if (  (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100 ) <= sw . Q && v[i][k * m + j] <= sw . N )
					{
						AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( AlphaChar ) );  //AAAA_ATAT_TATTT
						char      * motif   = ( char * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( char ) );
						unsigned int ell;
						unsigned int offset;
						unsigned int index = 0;
						unsigned int jj = j - ( sw . l  ) + 1;

						for ( b = 0; b < sw . nb_boxes; b ++ ) 
						{
							/* Extract the structured motif from the sequences */
							if ( b == 0 )	
							{
								offset = 0;
								ell = sw . l;
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + l];
									motif[index]   = seqs[i][jj + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}	
							else
							{
								offset += sw . S[k][b - 1] + ell;
								ell = sw . blens[b - 1];
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + offset + l];
									motif[index]   = seqs[i][jj + offset + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}
						}

						index--;
						ACmotif [ index ] = 0;
						motif [ index ] = 0;

						/* Insert the extracted structured motif into the trie */
						if ( trie_retrieve ( trie, ACmotif, NULL ) != TRUE )
						{
							trie_store ( trie, ACmotif, 0 );

							char      * id   = ( char * ) calloc ( ( total_length_sm + sw . nb_boxes ) , sizeof( char ) );

							unsigned int ii, jj;
							for ( ii = 0; ii < total_length_sm + sw . nb_gaps; ii++ )
							{
								for ( jj = 0; jj < strlen ( alphabet_str ); jj++ )
									if ( motif[ii] == alphabet_str[jj] ) break;

								char * s   = ( char * ) calloc ( strlen ( alphabet_str ) , sizeof( char ) );
								if ( jj < strlen ( alphabet_str ) )	
								{
									sprintf( s, "%c", ( char ) ( '0' + jj ) );
									id[ii] = s[0] ;
								}
								else
									id[ii] = '-';
									 
								free ( s );
							}
							id[total_length_sm + sw . nb_gaps] = '\0';
							                   
							fprintf ( out_fd, "%s %s %d\t%d\n", motif, id, u[i][k * m + j], v[i][k * m + j] );

							free ( id );
							valid ++;
						}
						free ( ACmotif );
						free ( motif );
					}
                        }
		}
	}

	fprintf ( out_fd, "Nb models: %d\n", valid );
	fprintf ( out_fd, "User time: %lf sec.\n", exectime );
                
	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;

	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output considering a background file as input
*/
unsigned int write_motifs_back ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P )
{
	FILE * 		out_fd;				// file with the motifs
	FILE * 		un_out_fd;			// file with the unmatched fg motifs
	FILE * 		in_fd;				// file with the bg motifs
	unsigned int i, j, k;
	char line[LINE_SIZE];

	unsigned int pvalid = 0;
	unsigned int uvalid = 0;
	unsigned int num_bseqs;
	unsigned int max_alloc;

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

        struct Tdata * bdata = NULL;
	
	#ifdef _USE_MPFR
	mpfr_t mpfr_bin;
	mpfr_init2(mpfr_bin, ACC);

	mpfr_t * mpfr_factLUT = ( mpfr_t * ) malloc( ( num_seqs + 1 ) * sizeof( mpfr_t ) );
	mpfr_fillTable( mpfr_factLUT, ( unsigned long int ) num_seqs );
	#else
	long double bin_cdf = 0.;
	long double *log_factLUT = calloc( (num_seqs + 1 ), sizeof(long double) );
	fillTable ( log_factLUT, num_seqs );
	#endif
	char background_size_string[100];
	char background_quorum_size_string[100];

	int background_size = 0, background_quorum_size = 0;

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	if ( ( out_fd = fopen( sw . output_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

        write_motex_header ( out_fd, sw, num_seqs, exectime, P );

	if ( ( un_out_fd = fopen( sw . un_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

        write_motex_header ( un_out_fd, sw, num_seqs, exectime, P );

	if ( ( in_fd = fopen( sw . background_in_filename, "r") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

	/* Here we need to read the num of seqs and the quorum of the bf file */
	unsigned int line_num = 0;
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		line_num++;
		if ( line_num == 6 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			strcpy ( background_size_string, pch );
		}
		if ( line_num == 8 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, "\n" );
			strcpy ( background_quorum_size_string, pch );
		}
	}
	
	background_size = atoi( background_size_string );
	
	background_quorum_size = atoi( background_quorum_size_string );
	
	assert( background_size > 0 && background_quorum_size > 0 );
	
	max_alloc = num_bseqs = 0;
 
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		if ( num_bseqs >= max_alloc )
                {
                	bdata = ( struct Tdata * ) realloc ( bdata,   ( max_alloc + ALLOC_SIZE ) * sizeof ( struct Tdata ) );
                        max_alloc += ALLOC_SIZE;
               	}
                AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( sw . l + 1 ) , sizeof( AlphaChar ) );
		char * pch;
		pch = strtok ( line, " " );
		for ( k = 0; k < sw. l; k ++ )
			ACmotif[k] = ( AlphaChar ) pch[k];

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . u = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . n = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . r = atof( pch );

    		pch = strtok ( NULL, " ");
    		bdata[ num_bseqs ] . v = atoi( pch );

		TrieData data = num_bseqs;
		trie_store ( trie, ACmotif, data );

                free ( ACmotif );
		num_bseqs++;
	}
        
	if ( fclose ( in_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
                        if (  ( ( double ) u[i][j] / num_seqs ) * 100  >= sw . q && v[i][j] >= sw . n )
                        	if (  ( ( double ) u[i][j] / num_seqs ) * 100  <= sw . Q && v[i][j] <= sw . N )
				{
					AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( sw . l + 1 ) , sizeof( AlphaChar ) );
					char      * motif   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );
					memcpy ( motif, &seqs[i][j - ( sw . l ) + 1 ], sw . l );
					for ( k = 0; k < sw. l; k ++ )
						ACmotif[k] = ( AlphaChar )seqs[i][k + j - ( sw . l ) + 1];
					ACmotif [ k ] = 0;
					TrieData data;
					if ( trie_retrieve ( trie, ACmotif, &data ) != TRUE ) 	//it does not exist as bg motif; add it as a fg motif
					{
						data = -1;
						trie_store ( trie, ACmotif, data );
										
						#ifdef _USE_MPFR
						mpfr_binomial_cdf_less_than(mpfr_bin, (unsigned long int)u[i][j], (unsigned long int)num_seqs, (long double)background_quorum_size/background_size, mpfr_factLUT);
						
						mpfr_ui_sub(mpfr_bin, 1, mpfr_bin, GMP_RNDU);
						
						fprintf ( out_fd, "%s %d %d %lf %d %d %d ", 
							  	motif, 
							  	u[i][j],
							  	num_seqs,
							  	(double) u[i][j]/num_seqs,
							  	v[i][j],
							  	0,
							  	0
							  	);
						
						mpfr_out_str( out_fd, 10, 10, mpfr_bin, GMP_RNDU);
						
						fprintf( out_fd, "\n");
						
						#else
						bin_cdf = binomial_cdf_less_than (u[i][j], num_seqs, ( long double ) background_quorum_size/background_size, log_factLUT );
							
						fprintf ( out_fd, "%s %d %d %lf %d %d %d %Le\n", 
								motif, 
								u[i][j], 
								num_seqs, 
								(  ( ( double ) u[i][j] / num_seqs ) ), 
								v[i][j], 
								0, 
								0, 
								1. - bin_cdf  );
							
						#endif
							
						fprintf ( un_out_fd, "%s %d %d %lf %d\n", 	//write it out as an unmatched fg motif
								motif, 
								u[i][j], 
								num_seqs, 
								(  ( ( double ) u[i][j] / num_seqs ) ), 
								v[i][j] );
						
						uvalid ++;
					}
					else 							//it already exists 
					{
						if ( data != -1 )				//as a bg motif; add it as a fg motif matching this bg motif
						{
						    
							#ifdef _USE_MPFR

							mpfr_binomial_cdf_less_than(mpfr_bin, (unsigned long int)u[i][j], (unsigned long int) num_seqs, ( long double ) bdata[ data ] . r , mpfr_factLUT);
							    
							mpfr_ui_sub(mpfr_bin, 1, mpfr_bin, GMP_RNDU);
							   
							fprintf ( out_fd, "%s %d %d %lf %d %lf %d ", 
									motif, 
								      	u[i][j],
								      	num_seqs,
								      	(double) u[i][j]/num_seqs,
								      	v[i][j],						      
								      	bdata[data] . r,
								      	bdata[data] . v
								      	);
							    
							mpfr_out_str( out_fd, 10, 10, mpfr_bin, GMP_RNDU);
							fprintf(out_fd, "\n");
							    
							#else

							bin_cdf = binomial_cdf_less_than ( u[i][j], num_seqs, bdata[ data ].u, log_factLUT);

							fprintf ( out_fd, "%s %d %d %lf %d %lf %d %Le\n", 
									motif, 
									u[i][j], 
									num_seqs,
									(  ( ( double ) u[i][j] / num_seqs ) ), 
									v[i][j], 
									bdata[ data ]  . r, 
									bdata [ data ] . v,
									1. - bin_cdf
							    );

							#endif

							data = -1;
							trie_store ( trie, ACmotif, data );
							pvalid ++;
						}
							//as a fg motif; do nothing as it is a duplicate	
					}
					free ( motif );
					free ( ACmotif );
				}
		}
	}

	fprintf ( out_fd, "\n#A total of %d background motifs were read.\n", num_bseqs );
	fprintf ( out_fd, "#A total of %d valid motifs were extracted:\n", pvalid + uvalid );
	fprintf ( out_fd, "#  %d of them were matched with background motifs\n", pvalid );
	fprintf ( out_fd, "#  %d of them were not matched with any background motif.\n\n", uvalid );
	
	fprintf ( un_out_fd, "\n#A total of %d background motifs were read.\n", num_bseqs );
	fprintf ( un_out_fd, "#A total of %d valid motifs were extracted:\n", pvalid + uvalid );
	fprintf ( un_out_fd, "#  %d of them were not matched with any background motif.\n\n", uvalid );
	

	#ifdef _USE_MPFR
	
	#else
	free(log_factLUT);
	#endif
	
	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;
	free ( bdata );
	bdata = NULL;

	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	if ( fclose ( un_out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output for structured motifs considering a background file as input
*/
unsigned int write_structured_motifs_back ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P )
{
	FILE * 		out_fd;				// file with the motifs
	FILE * 		un_out_fd;			// file with the unmatched fg motifs
	FILE * 		in_fd;				// file with the bg motifs
	unsigned int i, j, k, l;
	char line[LINE_SIZE];

	unsigned int pvalid = 0;
	unsigned int uvalid = 0;
	unsigned int num_bseqs;
	unsigned int max_alloc;

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

        struct Tdata * bdata = NULL;
	
	#ifdef _USE_MPFR
	mpfr_t mpfr_bin;
	mpfr_init2(mpfr_bin, ACC);

	mpfr_t * mpfr_factLUT = ( mpfr_t * ) malloc( ( num_seqs + 1 ) * sizeof( mpfr_t ) );
	mpfr_fillTable( mpfr_factLUT, ( unsigned long int ) num_seqs );
	#else
	long double bin_cdf = 0.;
	long double *log_factLUT = calloc( (num_seqs + 1 ), sizeof(long double) );
	fillTable ( log_factLUT, num_seqs );
	#endif
	char background_size_string[100];
	char background_quorum_size_string[100];

	int background_size = 0, background_quorum_size = 0;

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	if ( ( out_fd = fopen( sw . output_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

        write_motex_header ( out_fd, sw, num_seqs, exectime, P );

	if ( ( un_out_fd = fopen( sw . un_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

        write_motex_header ( un_out_fd, sw, num_seqs, exectime, P );

	if ( ( in_fd = fopen( sw . background_in_filename, "r") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

	/* Here we need to read the num of seqs and the quorum of the bf file */
	unsigned int line_num = 0;
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		line_num++;
		if ( line_num == 6 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			strcpy ( background_size_string, pch );
		}
		if ( line_num == 8 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, "\n" );
			strcpy ( background_quorum_size_string, pch );
		}
	}
	
	background_size = atoi( background_size_string );
	
	background_quorum_size = atoi( background_quorum_size_string );
	
	assert( background_size > 0 && background_quorum_size > 0 );

	unsigned int b;
	unsigned int total_length_sm = sw . l;
	for ( b = 0; b < sw . nb_gaps; b ++ )
		total_length_sm += sw . blens[b];
	
	max_alloc = num_bseqs = 0;

	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		if ( num_bseqs >= max_alloc )
                {
                	bdata = ( struct Tdata * ) realloc ( bdata,   ( max_alloc + ALLOC_SIZE ) * sizeof ( struct Tdata ) );
                        max_alloc += ALLOC_SIZE;
               	}
		AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( AlphaChar ) );  //AAAA_ATAT_TATTT
		char * pch;
		pch = strtok ( line, " " );
		for ( k = 0; k < total_length_sm + sw . nb_gaps; k ++ )
			ACmotif[k] = ( AlphaChar ) pch[k];

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . u = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . n = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . r = atof( pch );

    		pch = strtok ( NULL, " ");
    		bdata[ num_bseqs ] . v = atoi( pch );

		TrieData data = num_bseqs;
		trie_store ( trie, ACmotif, data );

                free ( ACmotif );
		num_bseqs++;
	}
 
	if ( fclose ( in_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );
                for ( j = sw . l - 1; j < m; j ++ )
		{
			for ( k = 0; k < sw . nb_structs; k ++ )
			{
				if (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100  >= sw . q && v[i][k * m + j] >= sw . n )
					if (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100  <= sw . Q && v[i][k * m + j] <= sw . N )
					{
						AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( AlphaChar ) );  //AAAA_ATAT_TATTT
						char      * motif   = ( char * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( char ) );
						unsigned int ell;
						unsigned int offset;
						unsigned int index = 0;
						unsigned int jj = j - ( sw . l  ) + 1;

						for ( b = 0; b < sw . nb_boxes; b ++ ) 
						{
							/* Extract the structured motif from the sequences */
							if ( b == 0 )	
							{
								offset = 0;
								ell = sw . l;
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + l];
									motif[index]   = seqs[i][jj + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}	
							else
							{
								offset += sw . S[k][b - 1] + ell;
								ell = sw . blens[b - 1];
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + offset + l];
									motif[index]   = seqs[i][jj + offset + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}
						}

						index--;
						ACmotif [ index ] = 0;
						motif [ index ] = 0;

						TrieData data;

						if ( trie_retrieve ( trie, ACmotif, &data ) != TRUE ) 	//it does not exist as bg motif; add it as a fg motif
						{
							data = -1;
							trie_store ( trie, ACmotif, data );

							#ifdef _USE_MPFR
							mpfr_binomial_cdf_less_than( 	mpfr_bin, 
											(unsigned long int) u[i][k * m + j], 
											(unsigned long int) num_seqs, 
											(long double) background_quorum_size / background_size, 
											mpfr_factLUT);
							
							mpfr_ui_sub(mpfr_bin, 1, mpfr_bin, GMP_RNDU);

							fprintf ( out_fd, "%s %d %d %lf %d %d %d ", 
									motif, 
									u[i][k * m + j],
									num_seqs,
									(double) u[i][k * m + j]/num_seqs,
									v[i][k * m + j],
									0,
									0);
							
							mpfr_out_str( out_fd, 10, 10, mpfr_bin, GMP_RNDU);

							fprintf( out_fd, "\n");
							
							#else
							bin_cdf = binomial_cdf_less_than ( u[i][k * m + j], num_seqs, ( long double ) background_quorum_size/background_size, log_factLUT );
								
							fprintf ( out_fd, "%s %d %d %lf %d %d %d %Le\n", 
									motif, 
									u[i][k * m + j], 
									num_seqs, 
									(  ( ( double ) u[i][k * m + j] / num_seqs ) ), 
									v[i][k * m + j], 
									0, 
									0, 
									1. - bin_cdf  );
							#endif
								
							fprintf ( un_out_fd, "%s %d %d %lf %d\n", 	//write it out as an unmatched fg motif
									motif, 
									u[i][k * m + j], 
									num_seqs, 
									(  ( ( double ) u[i][k * m + j] / num_seqs ) ), 
									v[i][k * m + j] );
							uvalid ++;
						}
						else 							//it already exists 
						{
							if ( data != -1 )				//as a bg motif; add it as a fg motif matching this bg motif
							{
								#ifdef _USE_MPFR
							
								mpfr_binomial_cdf_less_than(mpfr_bin, (unsigned long int) u[i][k * m + j], (unsigned long int) num_seqs, ( long double ) bdata[ data ] . r , mpfr_factLUT);
							    
								mpfr_ui_sub(mpfr_bin, 1, mpfr_bin, GMP_RNDU);

								fprintf ( out_fd, "%s %d %d %lf %d %lf %d ", 
										motif, 
										u[i][k * m + j],
										num_seqs,
										(double) u[i][k * m + j]/num_seqs,
										v[i][k * m + j],						      
										bdata[data] . r,
										bdata[data] . v
										);
								    
								mpfr_out_str( out_fd, 10, 10, mpfr_bin, GMP_RNDU);
								fprintf(out_fd, "\n");
								    
								#else

								bin_cdf = binomial_cdf_less_than ( u[i][k * m + j], num_seqs, bdata[ data ].u, log_factLUT);

								fprintf ( out_fd, "%s %d %d %lf %d %lf %d %Le\n", 
										motif, 
										u[i][k * m + j], 
										num_seqs,
										(  ( ( double ) u[i][k * m + j] / num_seqs ) ), 
										v[i][k * m + j], 
										bdata[ data ]  . r, 
										bdata [ data ] . v,
										1. - bin_cdf
								    );

								#endif

								data = -1;
								trie_store ( trie, ACmotif, data );
								pvalid ++;
							}
								//as a fg motif; do nothing as it is a duplicate	
						}
						free ( motif );
						free ( ACmotif );					
					}
			}
		}
	}

	fprintf ( out_fd, "\n#A total of %d background motifs were read.\n", num_bseqs );
	fprintf ( out_fd, "#A total of %d valid motifs were extracted:\n", pvalid + uvalid );
	fprintf ( out_fd, "#  %d of them were matched with background motifs\n", pvalid );
	fprintf ( out_fd, "#  %d of them were not matched with any background motif.\n\n", uvalid );

	fprintf ( un_out_fd, "\n#A total of %d background motifs were read.\n", num_bseqs );
	fprintf ( un_out_fd, "#A total of %d valid motifs were extracted:\n", pvalid + uvalid );
	fprintf ( un_out_fd, "#  %d of them were not matched with any background motif.\n\n", uvalid );
	

	#ifdef _USE_MPFR
	
	#else
	free(log_factLUT);
	#endif
	
	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;
	free ( bdata );
	bdata = NULL;

	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	if ( fclose ( un_out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output a la SMILE considering a motex background file as input
*/
unsigned int write_motifs_back_smile ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P )
{
	FILE * 		out_fd;				// file with the motifs
	FILE * 		un_out_fd;			// file with the motifs
	FILE * 		in_fd;				// file with the bg motifs
	unsigned int i, j, l, k;
	char line[LINE_SIZE];
        char * alphabet_str = NULL;

	unsigned int pvalid = 0;
	unsigned int uvalid = 0;
	unsigned int num_bseqs;
	unsigned int max_alloc;

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

        struct Tdata * bdata = NULL;
	
	if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet_str = ( char * ) DNA;
        else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet_str = ( char * ) PROT;
        else if ( ! strcmp ( "USR", sw . alphabet ) )   alphabet_str = ( char * ) USR;

	if ( ( un_out_fd = fopen( sw . un_smile_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file %s!\n", sw . un_smile_out_filename );
		return  ( 0 );
	}

	write_smile_header ( un_out_fd, sw, num_seqs, alphabet_str );

	if ( ( out_fd = fopen( sw . smile_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file %s!\n", sw . smile_out_filename );
		return  ( 0 );
	}

	write_smile_header ( out_fd, sw, num_seqs, alphabet_str );

	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	if ( ( in_fd = fopen( sw . background_in_filename, "r") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file %s!\n", sw . background_in_filename );
		return  ( 0 );
	}

	unsigned int line_num = 0;
	char background_size_string[100];
	char background_quorum_size_string[100];

	/* Here we need to read the num of seqs and the quorum of the bf file */
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		line_num++;
		if ( line_num == 6 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			strcpy ( background_size_string, pch );
		}
		if ( line_num == 8 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, "\n" );
			strcpy ( background_quorum_size_string, pch );
		}
	}

	max_alloc = num_bseqs = 0;
 
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		if ( num_bseqs >= max_alloc )
                {
                	bdata = ( struct Tdata * ) realloc ( bdata,   ( max_alloc + ALLOC_SIZE ) * sizeof ( struct Tdata ) );
                        max_alloc += ALLOC_SIZE;
               	}
                AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( sw . l + 1 ) , sizeof( AlphaChar ) );
		char * pch;
		pch = strtok ( line, " " );
		for ( k = 0; k < sw. l; k ++ )
			ACmotif[k] = ( AlphaChar ) pch[k];

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . u = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . n = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . r = atof( pch );

    		pch = strtok ( NULL, " ");
    		bdata[ num_bseqs ] . v = atoi( pch );

		TrieData data = num_bseqs;
		trie_store ( trie, ACmotif, data );

                free ( ACmotif );
		num_bseqs++;
	}

	if ( fclose ( in_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
			if (  ( ( double ) u[i][j] / num_seqs ) * 100  >= sw . q && v[i][j] >= sw . n )
				if (  ( ( double ) u[i][j] / num_seqs ) * 100  <= sw . Q && v[i][j] <= sw . N )
				{
					AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( sw . l + 1 ) , sizeof( AlphaChar ) );
					char      * motif   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );
					memcpy ( motif, &seqs[i][j - ( sw . l ) + 1 ], sw . l );

					for ( k = 0; k < sw. l; k ++ )
						ACmotif[k] = ( AlphaChar )seqs[i][k + j - ( sw . l ) + 1];
					ACmotif [ k ] = 0;

					TrieData data;

					if ( trie_retrieve ( trie, ACmotif, &data ) != TRUE ) 	// ACmotif does not exist as bg motif --- output it to the unmatched ones
					{
						data = -1;

						trie_store ( trie, ACmotif, data );

						char      * id   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );

						unsigned int ii, jj;
						for ( ii = 0; ii < sw . l; ii++ )
						{
							for ( jj = 0; jj < strlen ( alphabet_str ); jj++ )
								if ( motif[ii] == alphabet_str[jj] ) break;

							char * s   = ( char * ) calloc ( ( sw . l ) , sizeof( char ) );
							sprintf( s, "%c", ( char ) ( '0' + jj ) );
							id[ii] = s[0] ;
							free ( s );
						}
						id[sw . l] = '\0';
						fprintf ( un_out_fd, "%s %s %d\t%d\n", motif, id, u[i][j], v[i][j] );
						fprintf ( out_fd,    "%s %s %d\t%d\n", motif, id, u[i][j], v[i][j] );

						free ( id );
						uvalid ++;
					}
					else 							// it already exists 
					{
						if ( data != -1 )				// as a bg motif; add it as a fg motif matching this bg motif
						{
							char      * id   = ( char * ) calloc ( ( sw . l + 1 ) , sizeof( char ) );

							unsigned int ii, jj;
							for ( ii = 0; ii < sw . l; ii++ )
							{
								for ( jj = 0; jj < strlen ( alphabet_str ); jj++ )
									if ( motif[ii] == alphabet_str[jj] ) break;

								char * s   = ( char * ) calloc ( ( sw . l ) , sizeof( char ) );
								sprintf( s, "%c", ( char ) ( '0' + jj ) );
								id[ii] = s[0] ;
								free ( s );
							}
							id[sw . l] = '\0';
							
							fprintf ( out_fd, "%s %s %d\t%d\n", motif, id, u[i][j], v[i][j] );
							free ( id );
							data = -1;
							trie_store ( trie, ACmotif, data );
							pvalid ++;
						}
					}
					free ( motif );
					free ( ACmotif );
				}
		}
	}

	fprintf ( out_fd, "Nb models: %d\n", pvalid + uvalid );
	fprintf ( out_fd, "User time: %lf sec.\n", exectime );
	fprintf ( un_out_fd, "Nb models: %d\n", uvalid );
	fprintf ( un_out_fd, "User time: %lf sec.\n", exectime );

	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;
	free ( bdata );
	bdata = NULL;

	if ( fclose ( un_out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output for structured motifs a la SMILE considering a motex background file as input
*/
unsigned int write_structuted_motifs_back_smile ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P )
{
	FILE * 		out_fd;				// file with the motifs
	FILE * 		un_out_fd;			// file with the motifs
	FILE * 		in_fd;				// file with the bg motifs
	unsigned int i, j, k, l;
	char line[LINE_SIZE];
        char * alphabet_str = NULL;

	unsigned int pvalid = 0;
	unsigned int uvalid = 0;
	unsigned int num_bseqs;
	unsigned int max_alloc;

	AlphaMap *	alphabet = NULL;
        _Trie  *      	trie = NULL;

        struct Tdata * bdata = NULL;
	
	if      ( ! strcmp ( "DNA", sw . alphabet ) )   alphabet_str = ( char * ) DNA;
        else if ( ! strcmp ( "PROT", sw . alphabet ) )  alphabet_str = ( char * ) PROT;
        else if ( ! strcmp ( "USR", sw . alphabet ) )   alphabet_str = ( char * ) USR;

	if ( ( un_out_fd = fopen( sw . un_smile_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file %s!\n", sw . un_smile_out_filename );
		return  ( 0 );
	}

	write_smile_header ( un_out_fd, sw, num_seqs, alphabet_str );

	if ( ( out_fd = fopen( sw . smile_out_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file %s!\n", sw . smile_out_filename );
		return  ( 0 );
	}

	write_smile_header ( out_fd, sw, num_seqs, alphabet_str );


	/* Create an empty alphabet */
	alphabet = alpha_map_new();

	/* Define the alphabet's range */
        alpha_map_add_range ( alphabet, 0, 127 );

	/* Create an empty trie based on the alphabet */
        trie = trie_new ( alphabet );

	if ( ( in_fd = fopen( sw . background_in_filename, "r") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file %s!\n", sw . background_in_filename );
		return  ( 0 );
	}

	unsigned int line_num = 0;
	char background_size_string[100];
	char background_quorum_size_string[100];

	/* Here we need to read the num of seqs and the quorum of the bf file */
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		line_num++;
		if ( line_num == 6 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			strcpy ( background_size_string, pch );
		}
		if ( line_num == 8 )
		{
				
			char * pch;
			pch = strtok ( line, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, " " );
			pch = strtok ( NULL, "\n" );
			strcpy ( background_quorum_size_string, pch );
		}
	}

	unsigned int b;
	unsigned int total_length_sm = sw . l;
	for ( b = 0; b < sw . nb_gaps; b ++ )
		total_length_sm += sw . blens[b];
	
	max_alloc = num_bseqs = 0;
 
	while ( fgets ( line, LINE_SIZE, in_fd ) && line[0] != '\n' )
	{
		if ( num_bseqs >= max_alloc )
                {
                	bdata = ( struct Tdata * ) realloc ( bdata,   ( max_alloc + ALLOC_SIZE ) * sizeof ( struct Tdata ) );
                        max_alloc += ALLOC_SIZE;
               	}
		AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( AlphaChar ) );  //AAAA_ATAT_TATTT
		char * pch;
		pch = strtok ( line, " " );
		for ( k = 0; k < total_length_sm + sw . nb_gaps; k ++ )
			ACmotif[k] = ( AlphaChar ) pch[k];

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . u = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . n = atoi( pch );

    		pch = strtok ( NULL, " " );
    		bdata[ num_bseqs ] . r = atof( pch );

    		pch = strtok ( NULL, " ");
    		bdata[ num_bseqs ] . v = atoi( pch );

		TrieData data = num_bseqs;
		trie_store ( trie, ACmotif, data );

                free ( ACmotif );
		num_bseqs++;
	}

	if ( fclose ( in_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}

	for ( i = 0; i < num_seqs; i ++ ) 
	{
		unsigned int m = strlen ( seqs[i] );

                for ( j = sw . l - 1; j < m; j ++ )
		{
			for ( k = 0; k < sw . nb_structs; k ++ )
			{
				if (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100  >= sw . q && v[i][k * m + j] >= sw . n )
					if (  ( ( double ) u[i][k * m + j] / num_seqs ) * 100  <= sw . Q && v[i][k * m + j] <= sw . N )
					{
						AlphaChar * ACmotif = ( AlphaChar * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( AlphaChar ) );  //AAAA_ATAT_TATTT
						char      * motif   = ( char * ) calloc ( ( total_length_sm + sw . nb_gaps + 1 ) , sizeof( char ) );
						unsigned int ell;
						unsigned int offset;
						unsigned int index = 0;
						unsigned int jj = j - ( sw . l  ) + 1;

						for ( b = 0; b < sw . nb_boxes; b ++ ) 
						{
							/* Extract the structured motif from the sequences */
							if ( b == 0 )	
							{
								offset = 0;
								ell = sw . l;
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + l];
									motif[index]   = seqs[i][jj + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}	
							else
							{
								offset += sw . S[k][b - 1] + ell;
								ell = sw . blens[b - 1];
								for ( l = 0; l < ell; l++, index++ )
								{
									ACmotif[index] = ( AlphaChar ) seqs[i][jj + offset + l];
									motif[index]   = seqs[i][jj + offset + l];
								}
								ACmotif [ index ] = '_'; 
								motif[index] = '_';
								index ++;
							}
						}

						index--;
						ACmotif [ index ] = 0;
						motif [ index ] = 0;

						TrieData data;

						if ( trie_retrieve ( trie, ACmotif, &data ) != TRUE ) 	// ACmotif does not exist as bg motif --- output it to the unmatched ones
						{
							data = -1;
							trie_store ( trie, ACmotif, data );
							char      * id   = ( char * ) calloc ( ( total_length_sm + sw . nb_boxes ) , sizeof( char ) );

							unsigned int ii, jj;
							for ( ii = 0; ii < total_length_sm + sw . nb_gaps; ii++ )
							{
								for ( jj = 0; jj < strlen ( alphabet_str ); jj++ )
									if ( motif[ii] == alphabet_str[jj] ) break;

								char * s   = ( char * ) calloc ( strlen ( alphabet_str ) , sizeof( char ) );
								if ( jj < strlen ( alphabet_str ) )	
								{
									sprintf( s, "%c", ( char ) ( '0' + jj ) );
									id[ii] = s[0] ;
								}
								else
									id[ii] = '-';
									 
								free ( s );
							}
							id[total_length_sm + sw . nb_gaps] = '\0';
							fprintf ( un_out_fd, "%s %s %d\t%d\n", motif, id, u[i][k * m + j], v[i][k * m + j] );
							fprintf ( out_fd,    "%s %s %d\t%d\n", motif, id, u[i][k * m + j], v[i][k * m + j] );
							free ( id );
							uvalid ++;
						}
						else 							// it already exists 
						{
							if ( data != -1 )				// as a bg motif; add it as a fg motif matching this bg motif
							{
								char      * id   = ( char * ) calloc ( ( total_length_sm + sw . nb_boxes ) , sizeof( char ) );

								unsigned int ii, jj;
								for ( ii = 0; ii < total_length_sm + sw . nb_gaps; ii++ )
								{
									for ( jj = 0; jj < strlen ( alphabet_str ); jj++ )
										if ( motif[ii] == alphabet_str[jj] ) break;

									char * s   = ( char * ) calloc ( strlen ( alphabet_str ) , sizeof( char ) );
									if ( jj < strlen ( alphabet_str ) )	
									{
										sprintf( s, "%c", ( char ) ( '0' + jj ) );
										id[ii] = s[0] ;
									}
									else
										id[ii] = '-';
										 
									free ( s );
								}
								id[total_length_sm + sw . nb_gaps] = '\0';
								fprintf ( out_fd, "%s %s %d\t%d\n", motif, id, u[i][k * m + j], v[i][k * m + j] );
								free ( id );
								data = -1;
								trie_store ( trie, ACmotif, data );
								pvalid ++;
							}
						}
						free ( motif );
						free ( ACmotif );
					}
			}
		}
	}

	fprintf ( out_fd, "Nb models: %d\n", pvalid + uvalid );
	fprintf ( out_fd, "User time: %lf sec.\n", exectime );
	fprintf ( un_out_fd, "Nb models: %d\n", uvalid );
	fprintf ( un_out_fd, "User time: %lf sec.\n", exectime );

	trie_free ( trie );    
        trie = NULL;
        alpha_map_free ( alphabet );
        alphabet = NULL;
	free ( bdata );
	bdata = NULL;

	if ( fclose ( un_out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}

/*
write the output considering a foreground file of previously unmatched motifs as input
*/
unsigned int write_motifs_fore ( struct TSwitch sw, unsigned int num_fseqs, char const ** fseqs, unsigned int ** u, unsigned int ** v, double exectime, int P, unsigned int num_seqs, struct Tdata * fdata )
{
	FILE * 		out_fd;				// file with the motifs
	unsigned int i, j, k;
	unsigned int valid = 0;

#ifdef _USE_MPFR
	mpfr_t mpfr_bin;
	mpfr_init2(mpfr_bin, ACC);

	mpfr_t * mpfr_factLUT = ( mpfr_t * ) malloc( (num_seqs +1) * sizeof(mpfr_t));
	mpfr_fillTable(mpfr_factLUT, (unsigned long int)num_seqs);
#else
	long double *log_factLUT = calloc( (num_seqs + 1), sizeof(long double) );
	fillTable ( log_factLUT, num_seqs );
#endif
	if ( ( out_fd = fopen( sw . output_filename, "w") ) == NULL) 
	{	 
		fprintf( stderr, " Error: cannot open file!\n");
		return  ( 0 );
	}

        write_motex_header ( out_fd, sw, num_seqs, exectime, P );


	for ( i = 0; i < num_fseqs; i ++ ) 
	{
		unsigned int m = strlen ( fseqs[i] );

                for ( j = m - 1; j < m; j ++ )
		{
                        if (  (  ( ( double ) u[i][j] / num_seqs ) * 100 ) >= sw . q && v[i][j] >= sw . n )
                        	if (  (  ( ( double ) u[i][j] / num_seqs ) * 100 ) <= sw . Q && v[i][j] <= sw . N )
				{
					char      * motif   = ( char * ) calloc ( ( m + 1 ) , sizeof( char ) );
					memcpy ( motif, &fseqs[i][j - ( m ) + 1 ], m );
					
	#ifdef _USE_MPFR
					mpfr_binomial_cdf_less_than(mpfr_bin, (unsigned long int) fdata[i] . u, (unsigned long int) fdata[i] . n, (long double) u[i][j] / num_seqs, mpfr_factLUT);
					
					mpfr_ui_sub(mpfr_bin, 1, mpfr_bin, GMP_RNDU);
					
					fprintf ( out_fd, "%s %d %d %lf %d %lf %d ", 
						  motif, 
						  fdata[i] . u,
						  fdata[i] . n,
						  fdata[i] . r,
						  fdata[i] . v,
						  (double) u[i][j]/ num_seqs,
						  v[i][j]
						  );
					mpfr_out_str( out_fd, 10, 10, mpfr_bin, GMP_RNDU);
					fprintf( out_fd, "\n");
	#else
					
					long double bin_cdf = binomial_cdf_less_than ( fdata[i] . u, fdata[i] . n, ( long double ) u[i][j] / num_seqs, log_factLUT );

					fprintf ( out_fd, "%s %d %d %lf %d %lf %d %Le\n", 
						  motif, 
						  fdata[i] . u,
						  fdata[i] . n,
						  fdata[i] . r,
						  fdata[i] . v,
						  (double) u[i][j] /num_seqs,
						  v[i][j],
						  1. - bin_cdf
						  );
	#endif
					valid ++;
					free ( motif );
				}
		}
	}

	fprintf ( out_fd, "\n#A total of %d valid motifs were extracted.\n\n", valid );
                
	if ( fclose ( out_fd ) ) 
	{
      		fprintf( stderr, " Error: file close error!\n");				      
		return ( 0 );
	}
	return ( 1 );
}


static struct option long_options[] =
 {
   { "alphabet",           	required_argument, NULL, 'a' },
   { "background-in-file",    	required_argument, NULL, 'b' },
   { "input-file",         	required_argument, NULL, 'i' },
   { "output-file",        	required_argument, NULL, 'o' },
   { "quorum",         	   	required_argument, NULL, 'q' },
   { "max-quorum",         	required_argument, NULL, 'Q' },
   { "motifs-length",      	required_argument, NULL, 'k' },
   { "distance",   	   	required_argument, NULL, 'd' },
   { "errors",   	   	required_argument, NULL, 'e' },
   { "num-of-occ", 	        required_argument, NULL, 'n' },
   { "max-num-of-occ", 	        required_argument, NULL, 'N' },
   { "num-of-threads",     	required_argument, NULL, 't' },
//   { "long-sequences", 	  	required_argument, NULL, 'L' },
   { "boxes-in-file",		required_argument, NULL, 's' },
   { "unmatched-in-file",	required_argument, NULL, 'I' },
   { "unmatched-out-file",	required_argument, NULL, 'u' },
   { "SMILE-out-file",		required_argument, NULL, 'S' },
   { "unmatched-SMILE-out-file",required_argument, NULL, 'U' },
   { "help",               	no_argument,       NULL, 'h' },
   { NULL,                 	0,                 NULL, 0   }
 };

/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet			= NULL;
   sw -> background_in_filename		= NULL;
   sw -> input_filename 		= NULL;
   sw -> output_filename		= NULL;
   sw -> un_in_filename 		= NULL;
   sw -> un_out_filename		= NULL;
   sw -> smile_out_filename		= NULL;
   sw -> un_smile_out_filename		= NULL;
   sw -> boxes_in_filename		= NULL;
   sw -> q 				= 0;
   sw -> Q 				= 100;
   sw -> l        			= 0;
   sw -> d       			= 0;
   sw -> e       			= 0;
   sw -> n       			= 1;
   sw -> N       			= 10000;
   sw -> t       			= 4;
//   sw -> L       			= 0;

   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:b:i:o:q:Q:k:d:e:n:N:s:t:I:S:u:U:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
	   args ++;
           break;
         
         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
	   args ++;
           break;

         case 'b':
           sw -> background_in_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> background_in_filename, optarg );
           break;

         case 'I':
           sw -> un_in_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> un_in_filename, optarg );
           break;

         case 'u':
           sw -> un_out_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> un_out_filename, optarg );
           break;

         case 'S':
           sw -> smile_out_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> smile_out_filename, optarg );
           break;

         case 'U':
           sw -> un_smile_out_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> un_smile_out_filename, optarg );
           break;

         case 's':
           sw -> boxes_in_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> boxes_in_filename, optarg );
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
	   args ++;
           break;

         case 'q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> q = val;
	   args ++;
           break;

         case 'Q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> Q = val;
           break;

         case 'k':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
	   args ++;
           break;

         case 'n':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> n = val;
           break;

         case 'N':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> N = val;
           break;

         case 'd':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> d = val;
	   args ++;
           break;

         case 't':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> t = val;
           break;

	#if 0
         case 'L':
           val = strtol ( optarg, &ep, 10 );
           if ( ! ep || ep == optarg )
            {
              return ( 0 );
            }
           sw -> L = val;
           break;
	#endif

         case 'e':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> e = val;
	   args ++;
           break;

         case 'h':
           return ( 0 );
       }
    }
   
   if ( args < 7 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }

/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, "\nm    m       mmmmmmm        m    m\n" );
   fprintf ( stdout, "##  ##  mmm     #     mmm    #  #\n" );
   fprintf ( stdout, "# ## # #\" \"#    #    #\"  #    ##\n" );
   fprintf ( stdout, "# \"\" # #   #    #    #\"\"\"\"   m\"\"m\n" );
   fprintf ( stdout, "#    # \"#m#\"    #    \"#mm\"  m\"  \"m\n\n" );

   fprintf ( stdout, " Usage: motexCPU|motexOMP|motexMPI <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'\n"
                     "                                      for protein  sequences. You may use `USR'\n"
                     "                                      for user-defined alphabet; edit the file\n" 
                     "                                      motexdefs.h accordingly.\n" ); 
   fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     MoTeX output filename.\n" );
   fprintf ( stdout, "  -d, --distance            <int>     The  distance  used  for extracting  the\n"
                     "                                      motifs. It can be  either 0 (for Hamming\n" 
	             "                                      distance) or 1 (for edit distance).\n"); 
   fprintf ( stdout, "  -k, --motifs-length       <int>     The length for motifs.\n");
   fprintf ( stdout, "  -e, --errors              <int>     Limit the  max number  of errors to this\n"
                     "                                      value.\n" );
   fprintf ( stdout, "  -q, --quorum              <int>     The quorum is the minimum percentage (%%)\n"
                     "                                      of sequences in which a motif must occur.\n\n" );
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, "  -Q, --max-quorum          <int>     The maximum percentage (%%) of sequences\n"
                     "                                      in which a motif can occur (default: 100).\n" );
   fprintf ( stdout, "  -n, --num-of-occ          <int>     The minimum  number of  occurrences of a\n"
                     "                                      reported  motif in any  of the sequences\n"
                     "                                      (default: 1).\n" );
   fprintf ( stdout, "  -N, --max-num-of-occ      <int>     The maximum  number of  occurrences of a\n"
                     "                                      reported  motif in any  of the sequences\n"
                     "                                      (default: 10000).\n" );
   fprintf ( stdout, "  -s, --structured-motifs   <str>     Input filename  for the structure of the\n"
                     "                                      boxes in the case of structured motifs.\n" );
   fprintf ( stdout, "  -S, --SMILE-out-file      <str>     SMILE-like output filename to be used by\n"
                     "                                      SMILE.\n" );
   fprintf ( stdout, "  -b, --background-in-file  <str>     MoTeX background filename for statistical\n"
                     "                                      evaluation passed as input.\n" );
   fprintf ( stdout, "  -t, --threads             <int>     Number of threads to be used by the OMP\n"
                     "                                      version (default: 4).\n" );
   #if 0
   fprintf ( stdout, "  -L, --long-sequences      <int>     If the number of input sequences is less\n"
                     "                                      than  the number of  processors  used by\n" 
                     "                                      the MPI version, this should be set to 1\n"
                     "                                      (default: 0); useful  for a few (or one)\n"
                     "                                      very long sequence(s), e.g. a chromosome.\n" );
   #endif	
   fprintf ( stdout, "  -u, --un-out-file         <str>     Output filename for foreground motifs not\n"
                     "                                      matched exactly with any background motif\n"
                     "                                      in the file passed with the `-b' option.\n" );
   fprintf ( stdout, "  -I, --un-in-file          <str>     Input filename of the aforementioned file\n"
                     "                                      with the unmatched  motifs. These  motifs\n"
                     "                                      will be approximately searched  as motifs\n"
                     "                                      in the file passed with the `-i' option.\n" );
   fprintf ( stdout, "  -U, --SMILE-un-out-file   <str>     SMILE-like output filename for foreground\n"
                     "                                      motifs  not  matched  exactly  with  any\n"
                     "                                      background motif in the file passed with\n"
                     "                                      the `-b' option.\n" );
 }

/*
Timing function
*/
double gettime( void )
{
    struct timeval ttime;
    gettimeofday(&ttime , 0);
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};
