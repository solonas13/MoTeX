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
#include "motexdefs.h"

#ifdef _USE_MPFR

void mpfr_fillTable( mpfr_t * a, unsigned long int n )
{

	unsigned long int i = 0;
    	mpfr_init2( a[0], ACC );
  	mpfr_set_ui( a[0], 1, GMP_RNDU );
  
 	for( i = 1; i <= n; ++i )
    	{
     	 	mpfr_init2( a[i], ACC );
      		mpfr_mul_ui( a[i], a[i-1], i, GMP_RNDU );   
    	}
}

void mpfr_binomial_cdf_less_than( mpfr_t t, unsigned long int x, unsigned long int N, long double p, mpfr_t *LUT )
{
	unsigned long int i = 0;
  	mpfr_t mp, oneminusmp, pr, tmp1, tmp2;
  
  	mpfr_init2(mp, ACC);
  	mpfr_init2(oneminusmp, ACC);  
  	mpfr_init2(pr, ACC);
  	mpfr_init2(tmp1, ACC);
  	mpfr_init2(tmp2, ACC);

	mpfr_set_ld(t, 0.0L, GMP_RNDU);
	  
	if( p >= 1.0 )	return;
	  
	mpfr_set_ld(mp, p, GMP_RNDU);
	mpfr_ui_sub(oneminusmp, 1, mp, GMP_RNDU);
	  
	for( i = 0; i < x; ++i )
	{
		mpfr_pow_ui(tmp1, mp, i, GMP_RNDU);
	      	mpfr_pow_ui(tmp2, oneminusmp, (N-i), GMP_RNDU);
	      	mpfr_set(pr, LUT[N], GMP_RNDU);
	      	mpfr_div(pr, pr, LUT[i], GMP_RNDU);
	      	mpfr_div(pr, pr, LUT[N-i], GMP_RNDU);
	      	mpfr_mul(pr, pr, tmp1, GMP_RNDU);
	      	mpfr_mul(pr, pr, tmp2, GMP_RNDU);
	      	mpfr_add(t, t, pr, GMP_RNDU);
	}

	mpfr_clear(mp);
	mpfr_clear(oneminusmp);  
	mpfr_clear(pr);
	mpfr_clear(tmp1);
	mpfr_clear(tmp2);
}

#else

void fillTable ( long double * a, int n )
{ 
	int i = 0;
  	a[0] = 1;
  
  	for( i = 1; i <= n; ++i )	a[i] = a[i-1] * i;
}


long double binomial_cdf_less_than ( int x, int N, long double p, long double *LUT )
{
	int i = 0;
  	long double prob = 0.;
  	long double pr = 0.;

  	if ( p >= 1.0 )		return 0.;
  
  	for ( i = 0; i < x; ++i )
    	{
		pr = (LUT[N]/LUT[i])/LUT[N-i] * pow(p, i) * pow(1. - p, N-i); 
	  	prob += pr;
      		if (x > N )	break;
    	}

  	return ( prob );
}

#endif
