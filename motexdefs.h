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

#define ALLOC_SIZE 		512
#define LINE_SIZE 		512
#define DEL 			'$' 

#define DNA			"ACGTN"				//DNA alphabet
#define PROT			"ARNDCQEGHILKMFPSTWYV"		//Proteins alphabet
#define USR			"0123456789"			//User-defined alphabet, e.g. integer alphabet

#define max(a,b) ((a) > (b)) ? (a) : (b)
#define min(a,b) ((a) < (b)) ? (a) : (b)

struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
   char               * background_filename;
   char               * un_in_filename;  			
   char               * un_out_filename;  			
   char               * boxes_in_filename;  			
   char               * smile_out_filename;  			
   char               * un_smile_out_filename;  			
   unsigned int         q;
   unsigned int         Q;
   unsigned int         l;
   unsigned int         d;
   unsigned int         e;
   unsigned int		n;
   unsigned int		N;
   unsigned int		t;
   unsigned int		L;
   unsigned int         total_length;
   unsigned int         nb_gaps;
   unsigned int         nb_boxes;
   unsigned int         nb_structs;
   unsigned int       * bgaps_min;
   unsigned int       * bgaps_max;
   unsigned int       * blens;
   unsigned int       * berrs;
   unsigned int	     ** S;
 };

struct Tdata
 {
   unsigned int	        u;
   unsigned int	        n;
   double               r;
   unsigned int         v;
 };

int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
void usage ( void );

void recurse( unsigned int x, unsigned int xi, unsigned int * y, unsigned int ** A, unsigned int * newx, unsigned int * temp_inner, unsigned int ** B ); 
inline unsigned int bitminmax ( unsigned int a, unsigned int b, unsigned int c );
inline unsigned int shift ( unsigned int a );
inline unsigned int shiftc ( unsigned int a, unsigned int x );
inline unsigned int delta ( char a, char b );
inline unsigned int popcount ( unsigned int x );

#ifdef _USE_MPFR
#define ACC 200
#include <mpfr.h>
void mpfr_fillTable(mpfr_t * a, unsigned long int n);
void mpfr_binomial_cdf_less_than( mpfr_t t, unsigned long int x, unsigned long int N, long double p, mpfr_t *LUT );
#else
void fillTable( long double * a, int n );
long double binomial_cdf_less_than( int x, int N, long double p, long double * LUT );
#endif

#ifdef _USE_MPI
inline unsigned int i_coord ( unsigned int step, unsigned int n, unsigned int m, unsigned int x );
inline unsigned int j_coord ( unsigned int step, unsigned int n, unsigned int m, unsigned int x );
inline unsigned int upper_left ( unsigned int step, unsigned int n, unsigned int j);
inline unsigned int upper ( unsigned int step, unsigned int n, unsigned int j );
inline unsigned int left ( unsigned int step, unsigned int n, unsigned int j );
inline unsigned int elements( unsigned int step, unsigned int n, unsigned int m );
inline void allocation ( int rank, int P, unsigned int n, unsigned int m, unsigned int step, int* first, int* last );
inline void vec_allocation ( int rank, int P, unsigned int m, int *first, int* last, int* count );
inline void communication ( int P, int rank, unsigned int step, unsigned int n, unsigned int m, unsigned int* D, int first, int last, int first_n, int last_n );
unsigned int motifs_extraction_opasm_hd ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int*  v, int rank, int P );
unsigned int motifs_extraction_opasm_ed ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v, int rank, int P );
#endif

unsigned int motifs_extraction_hd ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v );
unsigned int motifs_extraction_ed ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v );
unsigned int structured_motifs_extraction_hd ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v );
unsigned int structured_motifs_extraction_ed ( const char * p, unsigned int m, const char * t, unsigned int n, struct TSwitch sw, unsigned int * u, unsigned int * v );

unsigned int write_motex_header ( FILE * out_fd, struct TSwitch sw, unsigned int num_seqs, double exectime, int P );
unsigned int write_smile_header ( FILE * out_fd, struct TSwitch sw, unsigned int num_seqs, char * alphabet_str );
unsigned int write_motifs ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P );
unsigned int write_structured_motifs ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P );
unsigned int write_motifs_smile ( struct TSwitch sw, unsigned int num_seqs, char const ** seqs, unsigned int ** u, unsigned int ** v, double exectime );
unsigned int write_structured_motifs_smile ( struct TSwitch sw, unsigned int num_seqs, char const ** seqs, unsigned int ** u, unsigned int ** v, double exectime );
unsigned int write_motifs_back ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P );
unsigned int write_structuted_motifs_back_smile ( struct TSwitch sw, unsigned int num_seqs, char const   ** seqs, unsigned int ** u, unsigned int ** v, double exectime, int P );
unsigned int write_motifs_fore ( struct TSwitch sw, unsigned int num_fseqs, char const ** fseqs, unsigned int ** u, unsigned int ** v, double exectime, int P, unsigned int num_seqs, struct Tdata * fdata );
double gettime( void );

