/***************************************************************************
                          common.cpp  -  description
                             -------------------
    begin                : Thu Apr 11 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


// This file contains common functions, that may be used throughout the library

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "commonPK.h"
#include "common.h"
#include "externs.h"
#include "externsPK.h"
#include "structs.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_internal_loop.h"
//#include "s_specific_functions.h"
#include "constantsPK.h"
//#include "common.h"     // July 17 - removed - unnecessary
//#include "constants.h"  // July 17 - removed - unnecessary
//#include "structsPK.h"  // July 17 - removed - unnecessary
//#include "externs.h"    // July 17 - removed - unnecessary
//#include "externsPK.h"  // July 17 - removed - unnecessary

#include "Stack.h"

#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_internal_loop.h"
#include "s_specific_functions.h"
#include "common.h"     // July 17 - removed - unnecessary
//#include "constants.h"  // July 17 - removed - unnecessary
//#include "structsPK.h"  // July 17 - removed - unnecessary
//#include "externs.h"    // July 17 - removed - unnecessary
//#include "externsPK.h"  // July 17 - removed - unnecessary
//#include "computeEnergy.h"

#include "params.h"
#include "paramsPK.h"

PARAMTYPE LEstacked_pair_energy (int i, int j, int *sequence)
{
//	printf("INSIDE LE_STACK = %f\n", s_stacked_pair::get_energy (i, j, sequence));
	return s_stacked_pair::get_energy (i, j, sequence);
}


PARAMTYPE LEhairpin_loop_energy (int i, int j, int* sequence, char *csequence)
// PRE:  seq is the sequence of the loop; important for tetraloops
//       I assume i-j can pair
// POST: Help function to compute_hairpin_matrix
//       Return the energy of a hairpin
{
	//PARAMTYPE testvar = s_hairpin_loop::get_energy(i, j, sequence, csequence);
	//printf("******* %f", testvar);
	//return testvar;
//	printf("Inside LE_HAIR = %f\n", s_hairpin_loop::get_energy(i, j, sequence, csequence));
	return s_hairpin_loop::get_energy(i, j, sequence, csequence);
}


PARAMTYPE LEinternal_loop_energy (int i, int j, int ip, int jp, int *sequence)
// PRE:  The energy matrix was calculated
// POST: Read from the read pool, write to nodes;
//       Store the node and return the energy.
{
//	printf("Inside LE_internal = %f\n", s_internal_loop::get_energy(i,j,ip,jp,sequence));
	return s_internal_loop::get_energy(i,j,ip,jp,sequence);
}

/*
void substr (char *source, int begin, int end, char *dest)
// PRE:  begin and end are smaller than strlen(source)
// POST: Put in dest what is in source between position begin and position end
{
	// in src/common/common.cpp
//    common::substr(source, begin, end, dest);
	substr(source, begin, end, dest);

}
int penalty_by_size (int size, char type)
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop
{
	// in src/common/common.cpp
	return penalty_by_size(size, type);

}

void giveup (char *string1, char *string2)
// to add: variable nb of parameters, as in scanf, printf
{
	// in src/common/common.cpp
	giveup(string1, string2);
}

void empty_string (char * str)
// PRE:  str is a string
// POST: Put '\0' at all positions
{
	// in src/common/common.cpp
	empty_string(str);
}
*/


// Cristina: changed this to match the simfold prototype

PARAMTYPE dangling_energy (int *sequence, char* structure, int i1, int i2, int i3, int i4)
//      (   )...(   )
//      i1..i2..i3..i4
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling energy between i2 and i3
{
	if (no_pk_dangling_ends == 0)
		return s_dangling_energy(sequence, structure, i1, i2, i3, i4);
	return 0;
}


PARAMTYPE dangling_energy_left (int *sequence, char* structure, int i1, int i2, int i3, int i4)
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i1 and i3
{
	if (no_pk_dangling_ends == 0)
		return s_dangling_energy_left(sequence, structure, i1, i2, i3, i4);
	return 0;
}


PARAMTYPE dangling_energy_right (int *sequence, char* structure, int i1, int i2, int i3, int i4)
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2
{
	if (no_pk_dangling_ends == 0)
		return s_dangling_energy_right(sequence, structure, i1, i2, i3, i4);
	return 0;
}


// COAXIAL STACKING

int LEcoax_stack_energy_flush_a (int i, int j, int ip, int jp, int *sequence)
// Pre: ip is i+1 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
// Post: Returns coaxial stacking energy from Walter et al 1999, or Turner parameters
//       if the given combination is not calculated in Walter et al. 1999, or 0 if
//       the returned value is INF (base pair is not canonical)
// Energy returned in 10cal/mol
{
    int coaxial = coaxstack_f_a [sequence[i]]
                  	   		[sequence[j]]
                    	 	[sequence[ip]]
                     		[sequence[jp]];

	if (DEBUG2)
		printf("LEcoax_stack_energy_flush_a: coaxstack[%d][%d]{%d][%d] = %d\n", i, j, ip, jp, coaxial);
	                     		
	if (coaxial >= INF)
	{
		coaxial = stack [sequence[i]]
	                     [sequence[j]]
    	                 [sequence[ip]]
        	             [sequence[jp]];

		if (DEBUG2)
			printf("LEcoax_stack_energy_flush_a: stack[%d][%d]{%d][%d] = %d\n", i, j, ip, jp, coaxial);

		if (coaxial >= INF)
			return 0;
	}
	
	return coaxial;
}

int LEcoax_stack_energy_flush_b (int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int other_j, int other_jp, int *sequence, int ignore_dangles)
// Pre: ip is i+1 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
//		flag = COAX_MULTI, COAX_PSEUDO, or COAX_OTHER
//		dangle_i corresponds to the base pair i.j
//		dangle_ip corresponds to the base pair ip.jp
// Post: Returns coaxial stacking energy from Turner lab / Mathews parameters (coaxial.dat).
//		 Currently, this is simply the Turner parameters for a regular stacked pair, or 0 if
//       the returned value is INF (base pair is not canonical)
// Energy returned in 10cal/mol
{
    int coaxial = coaxstack_f_b [sequence[i]]
                  	   		[sequence[j]]
                    	 	[sequence[ip]]
                     		[sequence[jp]];

	int dangle_en1 = 0;
	int dangle_en2 = 0;

	if (flag == COAX_MULTI)  // between rightmost child and outer base pair of multiloop (here, ip > jp, so dangle_ip = ip - 1 and use dangle_bot, dangle 3')
	{
	// outer base pair and leftmost child stack
	// dangle_ip   dangle_i
	//         |    |
	//       ( .    . (    )  ) 
	//       jp       j    i  ip

	// outer base pair and rightmost child stack
	//      dangle_ip  dangle_i
	//             |    |
	//  (   (    ) .    . )
	//  i   ip   jp      j

		if (ignore_dangles == 0)
		{
			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
															  [sequence[other_j]]
															  [sequence[j-1]]);	                          
				}
			}

			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
											[sequence[ip]]
											[sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// take dangle_en2 as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_f_b, multi: dangle_en1 = %d, dangle_en2 = %d\n", dangle_en1, dangle_en2);

		if (dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return 0;
	}
	else if (flag == COAX_OTHER)  //e.g two children in a multiloop (here, ip < jp here so dangle_ip is jp+1, so use dangle_top, dangle 5')
	{
	// dangle_i    dangle_ip
	//  |                |
	//  . (    )  (    ) .
	//    j    i  ip   jp

		if (ignore_dangles == 0)
		{
			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)  // the next base pair to the left of j shares the dangling end
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
					                          		[sequence[other_j]]
					                          		[sequence[j-1]]);
					else
			    		dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
			                         						[sequence[other_j]]
			                          						[sequence[j-1]]);
				}
			}
			
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
		
				if (other_jp > 0)  // the next base pair to the right of jp shares the dangling end
				{
					if (simple_dangling_ends)
					{
						// return the current dangle_en2
					}		
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_f_b, other: dangle_en1 = %d, dangle_en2 = %d\n", dangle_en1, dangle_en2);

		if (dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return 0;
	}
	else  //  if (flag == COAX_PSEUDO)
	{
	// dangle_ip          dangle_i
	//       |              |
	//     ( .    (    )    . ) 
	//     jp     i    ip     j

		if (ignore_dangles == 0)
		{
			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// leave as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_f_b, pseudo: dangle_en1 = %d, dangle_en2 = %d\n", dangle_en1, dangle_en2);

		if (dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return 0;
	}

	if (DEBUG2)
		printf("LEcoax_stack_energy_flush_b: coaxstack[%d][%d]{%d][%d] = %d\n", i, j, ip, jp, coaxial);
	                     		
	if (coaxial >= INF)
	{
		return 0;
	}
	
	return coaxial;
}

int LEcoax_stack_energy_mismatch(int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int *sequence, char* structure, int other_j, int other_jp, int& stack_dangle, int ignore_dangles)
// Pre: ip is i+2 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
//		flag = COAX_MULTI, COAX_PSEUDO, or COAX_OTHER
//		dangle_i corresponds to the base pair i.j
//		dangle_ip corresponds to the base pair ip.jp
// Post: Returns coaxial stacking mismatch energy from Turner lab / Mathews parameters.
//		 stack_dangle is set to 0 or 1 to indicate whether the i or ip dangle, respectively, is involved in the stacking
//
// Note: "structure" does not have to be the true representation (i.e. use different paranetheses for
// different bands. (In particular, do not use < or > since these have a special interpretation in
// the dangling energy function.)
//
// Energy returned in 10cal/mol
{
	int coaxial_i = 0;
	int coaxial_ip = 0;
	if (dangle_i > 0)  // involve dangle_i in coaxial stack only if it's a free base
    		coaxial_i = ( coaxstack_m1 [sequence[i]] [sequence[j]] [sequence[i+1]] [sequence[dangle_i]] +
    				  coaxstack_m2 [sequence[i+1]] [sequence[dangle_i]] [sequence[ip]] [sequence[jp]] );
	if (dangle_ip > 0)  // involve dangle_ip in coaxial stack only if it's a free base
		coaxial_ip = ( coaxstack_m1 [sequence[jp]] [sequence[ip]] [sequence[dangle_ip]] [sequence[i+1]] +
    				   coaxstack_m2 [sequence[i]] [sequence[j]] [sequence[i+1]] [sequence[dangle_ip]] );
	
	if (DEBUG2)
	{
		printf("LEcoax_stack_energy_mismatch: i.j, dangle_i %d(%d).%d(%d) dangle %d(%d); ip.jp, dangle_ip %d(%d).%d(%d), dangle %d(%d); i+1 = %d(%d)\n", 
					i, sequence[i], j, sequence[j], dangle_i, sequence[dangle_i], ip, sequence[ip], jp, sequence[jp], dangle_ip, sequence[dangle_ip], i+1, sequence[i+1]);
		printf("  coaxial_i = ");
		if (dangle_i > 0)
			printf("%d + %d = ", coaxstack_m1 [sequence[i]] [sequence[j]][sequence[i+1]] [sequence[dangle_i]],
			                        coaxstack_m2 [sequence[i+1]] [sequence[dangle_i]] [sequence[ip]] [sequence[jp]]);
		printf("%d; ", coaxial_i);
		printf("coaxial_ip = ");
		if (dangle_ip > 0)
			printf("%d + %d = ", coaxstack_m1 [sequence[jp]] [sequence[ip]] [sequence[dangle_ip]] [sequence[i+1]],
			                        coaxstack_m2 [sequence[i]] [sequence[j]] [sequence[i+1]] [sequence[dangle_ip]]);
		printf("%d \n", coaxial_ip);					
	}
	
	int coaxial = 0;

	if (coaxial_i <= coaxial_ip)
	{
		coaxial = coaxial_i;
		stack_dangle = 0;  // i is part of the continuous backbone stack
	}
	else
    {
		coaxial = coaxial_ip;
		stack_dangle = 1;  // ip is part of the continuous backbone stack
    }    				  

	int dangle_en0 = 0;
	int dangle_en1 = 0;
	int dangle_en2 = 0;

	if (flag == COAX_MULTI)  // between rightmost child and outer base pair of multiloop (here, ip > jp, so dangle_ip = ip - 1 and use dangle_bot, dangle 3')
	{
	// outer base pair and leftmost child stack
	// dangle_ip   dangle_i
	//         |    |
	//       ( .    . (    ) . ) 
	//       jp       j    i   ip

	// outer base pair and rightmost child stack
	//      dangle_ip  dangle_i
	//             |    |
	//  ( . (    ) .    . ) 
	//  i   ip   jp      j

		if (ignore_dangles == 0)
		{
			// dangling energy of base between i and ip
			if (ip < jp)  // outer base pair and leftmost child stacking
				dangle_en0 = s_dangling_energy_left(sequence, structure, i, j, ip, jp);
			else  // outer base pair and rightmost child stacking
				dangle_en0 = s_dangling_energy_right(sequence, structure, ip, j, i, jp);

			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
															  [sequence[other_j]]
															  [sequence[j-1]]);	                          
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
											[sequence[ip]]
											[sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// take dangle_en2 as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_mismatch, multi: dangle_en0 = %d, dangle_en1 = %d, dangle_en2 = %d\n", dangle_en0, dangle_en1, dangle_en2);

		if (dangle_en0 + dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return 0;
	}
	else if (flag == COAX_OTHER)  //e.g two children in a multiloop (here, ip < jp here so dangle_ip is jp+1, so use dangle_top, dangle 5')
	{
	// dangle_i    dangle_ip
	//  |                |
	//  . (    ) . (    ) .
	//    j    i   ip   jp

		if (ignore_dangles == 0)
		{
			// dangling energy of base between i and ip
			dangle_en0 = s_dangling_energy(sequence, structure, j, i, ip, jp);

			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)  // the next base pair to the left of j shares the dangling end
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
					                          		[sequence[other_j]]
					                          		[sequence[j-1]]);
					else
			    		dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
			                         						[sequence[other_j]]
			                          						[sequence[j-1]]);
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
		
				if (other_jp > 0)  // the next base pair to the right of jp shares the dangling end
				{
					if (simple_dangling_ends)
					{
						// return the current dangle_en2
					}		
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_mismatch, other: dangle_en0 = %d, dangle_en1 = %d, dangle_en2 = %d\n", dangle_en0, dangle_en1, dangle_en2);

		if (dangle_en0 + dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return 0;
	}
	else  //  if (flag == COAX_PSEUDO)
	{
	// ignore danling, calculate just coaxial energy
	// dangle_ip          dangle_i
	//       |              |
	//     ( .    (    )    . ) 
	//     jp     i    ip     j

		if (ignore_dangles == 0)
		{
			// dangling energy of base between i and ip
			dangle_en0 = s_dangling_energy(sequence, structure, j, i, ip, jp);
			// NOTE: even though this is the dangling energy between two stems of a pseudoknot,
			//       simply using the schema shown below for the parameters of s_dangling_energy
			//       yields the correct result, since the dangling bases are now on the inside of 
			//		 base pairs i1.i2, i3.i4 in the implementation of s_dangling_energy
			// (    ( ... )    ) 
			// i4   i2    i3   i1


			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// leave as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_mismatch, pseudo: dangle_en0 = %d, dangle_en1 = %d, dangle_en2 = %d\n", dangle_en0, dangle_en1, dangle_en2);

		if (dangle_en0 + dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return 0;
	}

	if (coaxial >= INF)
	{
		return 0;
	}
	
	return coaxial;

}


//// FOR PARAMETER TUNING //////

void count_LEcoax_stack_energy_flush_b (int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int other_j, int other_jp, int *sequence, double* counter, int ignore_dangles)
// Pre: ip is i+1 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
//		flag = COAX_MULTI, COAX_PSEUDO, or COAX_OTHER
//		dangle_i corresponds to the base pair i.j
//		dangle_ip corresponds to the base pair ip.jp
// Post: Returns coaxial stacking energy from Turner lab / Mathews parameters (coaxial.dat).
//		 Currently, this is simply the Turner parameters for a regular stacked pair, or 0 if
//       the returned value is INF (base pair is not canonical)

// NOTE: all the flags work, not just COAX_PSEUDO
{
    int coaxial = coaxstack_f_b [sequence[i]]
                  	   		[sequence[j]]
                    	 	[sequence[ip]]
                    		[sequence[jp]];

	int num_params_DP_and_pkfree = get_num_params_PK_DP();

	char type[100];
	sprintf(type, "coaxstack_f_b[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[ip], sequence[jp]);
	int index_coax = structure_type_index_PK_CC(type);

	int dangle_en1 = 0;
	int dangle_en2 = 0;

	if (flag == COAX_MULTI)  // between rightmost child and outer base pair of multiloop (here, ip > jp, so dangle_ip = ip - 1 and use dangle_bot, dangle 3')
	{
	// outer base pair and leftmost child stack
	// dangle_ip   dangle_i
	//         |    |
	//       ( .    . (    )  ) 
	//       jp       j    i  ip

	// outer base pair and rightmost child stack
	//      dangle_ip  dangle_i
	//             |    |
	//  (   (    ) .    . )
	//  i   ip   jp      j

		if (ignore_dangles == 0)
		{
			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
															  [sequence[other_j]]
															  [sequence[j-1]]);	                          
				}
			}

			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
											[sequence[ip]]
											[sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// take dangle_en2 as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_f_b, multi: dangle_en1 = %d, dangle_en2 = %d\n", dangle_en1, dangle_en2);

		if (dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return;
	}
	else if (flag == COAX_OTHER)  //e.g two children in a multiloop (here, ip < jp here so dangle_ip is jp+1, so use dangle_top, dangle 5')
	{
	// dangle_i    dangle_ip
	//  |                |
	//  . (    )  (    ) .
	//    j    i  ip   jp

		if (ignore_dangles == 0)
		{
			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)  // the next base pair to the left of j shares the dangling end
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
					                          		[sequence[other_j]]
					                          		[sequence[j-1]]);
					else
			    		dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
			                         						[sequence[other_j]]
			                          						[sequence[j-1]]);
				}
			}
			
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
		
				if (other_jp > 0)  // the next base pair to the right of jp shares the dangling end
				{
					if (simple_dangling_ends)
					{
						// return the current dangle_en2
					}		
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_f_b, other: dangle_en1 = %d, dangle_en2 = %d\n", dangle_en1, dangle_en2);

		if (dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return;
	}
	else  //  if (flag == COAX_PSEUDO)
	{
	// dangle_ip          dangle_i
	//       |              |
	//     ( .    (    )    . ) 
	//     jp     i    ip     j

		if (ignore_dangles == 0)
		{
			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// leave as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_f_b, pseudo: dangle_en1 = %d, dangle_en2 = %d\n", dangle_en1, dangle_en2);

		if (dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return;
	}

	if (DEBUG2)
		printf("LEcoax_stack_energy_flush_b: coaxstack[%d][%d]{%d][%d] = %d\n", i, j, ip, jp, coaxial);
	                     		
	if (coaxial >= INF)
	{
		return;
	}

	counter[num_params_DP_and_pkfree + index_coax - 1]++;	
	return;  // only case where we return coaxial (ie. use coaxial stacking)
}

void count_LEcoax_stack_energy_mismatch(int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int *sequence, char* structure, int other_j, int other_jp, int& stack_dangle, double* counter, int ignore_dangles)
// Pre: ip is i+2 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
//		flag = COAX_MULTI, COAX_PSEUDO, or COAX_OTHER
//		dangle_i corresponds to the base pair i.j
//		dangle_ip corresponds to the base pair ip.jp
// Post: Returns coaxial stacking mismatch energy from Turner lab / Mathews parameters.
//		 stack_dangle is set to 0 or 1 to indicate whether the i or ip dangle, respectively, is involved in the stacking
//
// NOTE: all the flags work, not just COAX_PSEUDO
// Note: "structure" does not have to be the true representation (i.e. use different paranetheses for
// different bands. (In particular, do not use < or > since these have a special interpretation in
// the dangling energy function.)
{
	char type[100];
	int index_coax1 = -1;
	int index_coax2 = -1;
	int index_coax_i1 = -1;
	int index_coax_ip1 = -1;
	int index_coax_i2 = -1;
	int index_coax_ip2 = -1;
	int num_params_DP_and_pkfree = get_num_params_PK_DP();

	int coaxial_i = 0;
	int coaxial_ip = 0;
	if (dangle_i > 0)  // involve dangle_i in coaxial stack only if it's a free base
	{
    		coaxial_i = ( coaxstack_m1 [sequence[i]] [sequence[j]] [sequence[i+1]] [sequence[dangle_i]] +
    				  coaxstack_m2 [sequence[i+1]] [sequence[dangle_i]] [sequence[ip]] [sequence[jp]] );
		sprintf(type, "coaxstack_m1[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[i+1], sequence[dangle_i]);
		index_coax_i1 = structure_type_index_PK_CC(type);
		sprintf(type, "coaxstack_m2[%d][%d][%d][%d]", sequence[i+1], sequence[dangle_i], sequence[ip], sequence[jp]);
		index_coax_i2 = structure_type_index_PK_CC(type);
	}
	if (dangle_ip > 0)  // involve dangle_ip in coaxial stack only if it's a free base
	{
		coaxial_ip = ( coaxstack_m1 [sequence[jp]] [sequence[ip]] [sequence[dangle_ip]] [sequence[i+1]] +
    				   coaxstack_m2 [sequence[i]] [sequence[j]] [sequence[i+1]] [sequence[dangle_ip]] );
		sprintf(type, "coaxstack_m1[%d][%d][%d][%d]", sequence[jp], sequence[ip], sequence[dangle_ip], sequence[i+1]);
		index_coax_ip1 = structure_type_index_PK_CC(type);
		sprintf(type, "coaxstack_m2[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[i+1], sequence[dangle_ip]);
		index_coax_ip2 = structure_type_index_PK_CC(type);
	}
	
	if (DEBUG2)
	{
		printf("LEcoax_stack_energy_mismatch: i.j, dangle_i %d(%d).%d(%d) dangle %d(%d); ip.jp, dangle_ip %d(%d).%d(%d), dangle %d(%d); i+1 = %d(%d)\n", 
					i, sequence[i], j, sequence[j], dangle_i, sequence[dangle_i], ip, sequence[ip], jp, sequence[jp], dangle_ip, sequence[dangle_ip], i+1, sequence[i+1]);
		printf("  coaxial_i = ");
		if (dangle_i > 0)
			printf("%d + %d = ", coaxstack_m1 [sequence[i]] [sequence[j]][sequence[i+1]] [sequence[dangle_i]],
			                        coaxstack_m2 [sequence[i+1]] [sequence[dangle_i]] [sequence[ip]] [sequence[jp]]);
		printf("%d; ", coaxial_i);
		printf("coaxial_ip = ");
		if (dangle_ip > 0)
			printf("%d + %d = ", coaxstack_m1 [sequence[jp]] [sequence[ip]] [sequence[dangle_ip]] [sequence[i+1]],
			                        coaxstack_m2 [sequence[i]] [sequence[j]] [sequence[i+1]] [sequence[dangle_ip]]);
		printf("%d \n", coaxial_ip);					
	}
	
	int coaxial = 0;

	if (coaxial_i <= coaxial_ip)
	{
		coaxial = coaxial_i;
		stack_dangle = 0;  // i is part of the continuous backbone stack
		index_coax1 = index_coax_i1;
		index_coax2 = index_coax_i2;
	}
	else
	{
		coaxial = coaxial_ip;
		stack_dangle = 1;  // ip is part of the continuous backbone stack
		index_coax1 = index_coax_ip1;
		index_coax2 = index_coax_ip2;
	}    				  

	int dangle_en0 = 0;
	int dangle_en1 = 0;
	int dangle_en2 = 0;

	if (flag == COAX_MULTI)  // between rightmost child and outer base pair of multiloop (here, ip > jp, so dangle_ip = ip - 1 and use dangle_bot, dangle 3')
	{
	// outer base pair and leftmost child stack
	// dangle_ip   dangle_i
	//         |    |
	//       ( .    . (    ) . ) 
	//       jp       j    i   ip

	// outer base pair and rightmost child stack
	//      dangle_ip  dangle_i
	//             |    |
	//  ( . (    ) .    . ) 
	//  i   ip   jp      j

		if (ignore_dangles == 0)
		{
			// dangling energy of base between i and ip
			if (ip < jp)  // outer base pair and leftmost child stacking
				dangle_en0 = s_dangling_energy_left(sequence, structure, i, j, ip, jp);
			else  // outer base pair and rightmost child stacking
				dangle_en0 = s_dangling_energy_right(sequence, structure, ip, j, i, jp);

			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
															  [sequence[other_j]]
															  [sequence[j-1]]);	                          
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
											[sequence[ip]]
											[sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// take dangle_en2 as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_mismatch, multi: dangle_en0 = %d, dangle_en1 = %d, dangle_en2 = %d\n", dangle_en0, dangle_en1, dangle_en2);

		if (dangle_en0 + dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return;
	}
	else if (flag == COAX_OTHER)  //e.g two children in a multiloop (here, ip < jp here so dangle_ip is jp+1, so use dangle_top, dangle 5')
	{
	// dangle_i    dangle_ip
	//  |                |
	//  . (    ) . (    ) .
	//    j    i   ip   jp

		if (ignore_dangles == 0)
		{
			// dangling energy of base between i and ip
			dangle_en0 = s_dangling_energy(sequence, structure, j, i, ip, jp);

			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)  // the next base pair to the left of j shares the dangling end
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
					                          		[sequence[other_j]]
					                          		[sequence[j-1]]);
					else
			    		dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
			                         						[sequence[other_j]]
			                          						[sequence[j-1]]);
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
		
				if (other_jp > 0)  // the next base pair to the right of jp shares the dangling end
				{
					if (simple_dangling_ends)
					{
						// return the current dangle_en2
					}		
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_mismatch, other: dangle_en0 = %d, dangle_en1 = %d, dangle_en2 = %d\n", dangle_en0, dangle_en1, dangle_en2);

		if (dangle_en0 + dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return;
	}
	else  //  if (flag == COAX_PSEUDO)
	{
	// ignore danling, calculate just coaxial energy
	// dangle_ip          dangle_i
	//       |              |
	//     ( .    (    )    . ) 
	//     jp     i    ip     j

		if (ignore_dangles == 0)
		{
			// dangling energy of base between i and ip
			dangle_en0 = s_dangling_energy(sequence, structure, j, i, ip, jp);
			// NOTE: even though this is the dangling energy between two stems of a pseudoknot,
			//       simply using the schema shown below for the parameters of s_dangling_energy
			//       yields the correct result, since the dangling bases are now on the inside of 
			//		 base pairs i1.i2, i3.i4 in the implementation of s_dangling_energy
			// (    ( ... )    ) 
			// i4   i2    i3   i1


			// dangling energy of unpaired base outside of i.j
			if (dangle_i > 0)
			{
				dangle_en1 = MIN (0, dangle_bot[sequence[i]]
	    	        						  [sequence[j]]
    	    		        				  [sequence[dangle_i]]);		
				if (other_j > 0)
				{
					if (simple_dangling_ends)
						dangle_en1 = MIN (0, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
					else
						dangle_en1 = MIN (dangle_en1, dangle_top[sequence[j-2]]
													  [sequence[other_j]]
													  [sequence[j-1]]);
				}
			}
			// dangling energy of unpaired base outside of ip.jp
			if (dangle_ip > 0)
			{
				dangle_en2 = MIN (0, dangle_top[sequence[jp]]
	    	        						  [sequence[ip]]
    	    		        				  [sequence[dangle_ip]]);
				if (other_jp > 0)
				{
					if (simple_dangling_ends)
					{
						// leave as is
					}
					else
						dangle_en2 = MIN (dangle_en2, dangle_bot[sequence[other_jp]]
															  [sequence[jp+2]]
															  [sequence[jp+1]]);
				}
			}
		}

		if (DEBUG2)
			printf("LEcoax_stack_energy_mismatch, pseudo: dangle_en0 = %d, dangle_en1 = %d, dangle_en2 = %d\n", dangle_en0, dangle_en1, dangle_en2);

		if (dangle_en0 + dangle_en1 + dangle_en2 <= coaxial)  // if dangling bases yield lower energy, take dangling, include no coaxial stacking
			return;
	}

	if (coaxial >= INF)
	{
		return;
	}

	if (index_coax1 == -1)
		printf("ERROR: commonPK.cpp:: count_LEcoax_stack_energy_mismatch - index_coax1 == -1\n");
	if (index_coax2 == -1)
		printf("ERROR: commonPK.cpp:: count_LEcoax_stack_energy_mismatch - index_coax2 == -1\n");

	counter[num_params_DP_and_pkfree + index_coax1 - 1]++;	
	counter[num_params_DP_and_pkfree + index_coax2 - 1]++;	
	return;  // only case where we return coaxial (ie. use coaxial stacking)

}

///// END FOR PARAMETER TUNING ///////


// similar to simfold/src/s_specific_functions.cpp
int dangling_energy_res (int *sequence, int i1, int i2, int i3, int i4, int d_12, int d_34)
//      (   )...(   )
//      i1..i2..i3..i4
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
//		 d_12 = dangling end corresponding to base pair i1.i2
//		 d_34 = dangling end corresponding to base pair 13.i4
//       similar to dangling_energy() if d_12 or d_13 are 1, they will be consider as restricted 
//			and not included in the dangling calculation
// POST: return dangling energy between i2 and i3
{
    int energy;
    int d_top, d_bot;
    d_top = 0;
    d_bot = 0;


	if (d_12 == 1)
		d_top = 0;
	else
		d_top = MIN (0, dangle_top[sequence[i2]]
                          [sequence[i1]]
                          [sequence[i2+1]]);

	if (d_34 == 1)
		d_bot = 0;
	else
		d_bot = MIN (0, dangle_bot[sequence[i4]]
                          [sequence[i3]]
                          [sequence[i3-1]]);

	if (DEBUG2)
		printf("dangling_energy_res: d_top(%d) = %d, d_bot(%d) = %d\n", i2+1, d_top, i3-1, d_bot);

	if (i2+1 == i3-1)     // see which is smaller
    {
        if (simple_dangling_ends)
            energy = d_top;
        else
            energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i2+1 < i3-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}

// similar to simfold/src/s_specific_functions.cpp
int dangling_energy_left_res (int *sequence, int i1, int i2, int i3, int i4, int d_12, int d_34)
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
//		 d_12 = dangling end corresponding to base pair i1.i2
//		 d_34 = dangling end corresponding to base pair 13.i4
//       similar to dangling_energy_left() if d_12 or d_13 are 1, they will be consider as restricted 
//		   and not included in the dangling calculation
// POST: return dangling energy between i1 and i3
{
    int energy;
    int d_top, d_bot;
    d_top = 0;
    d_bot = 0;

	if (d_12 == 1)
		d_top = 0;
	else
	    d_top = MIN (0, dangle_top[sequence[i1]]
                          [sequence[i2]]
                          [sequence[i1+1]]);

	if (d_34 == 1)
		d_bot = 0;
	else
	    d_bot = MIN (0, dangle_bot[sequence[i4]]
                          [sequence[i3]]
                          [sequence[i3-1]]);

	if (DEBUG2)
		printf("dangling_energy_left_res: d_top(%d) = %d, d_bot(%d) = %d\n", i1+1, d_top, i3-1, d_bot);

    if (i1+1 == i3-1)     // see which is smaller
    {
        if (simple_dangling_ends)
            energy = d_top;
        else    
            energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i1+1 < i3-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}

// similar to simfold/src/s_specific_functions.cpp
int dangling_energy_right_res (int *sequence, int i1, int i2, int i3, int i4, int d_12, int d_34)
//      (    (    )...)
//      i1   i3  i4  i2
//		 d_12 = dangling end corresponding to base pair i1.i2
//		 d_34 = dangling end corresponding to base pair 13.i4
//       similar to dangling_energy_right() if d_12 or d_13 are 1, they will be consider as restricted 
//		   and not included in the dangling calculation
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2
{
    int energy;
    int d_top, d_bot;
    d_top = 0;

    d_bot = 0;

	if (d_34 == 1)
		d_top = 0;
	else
	    d_top = MIN (0, dangle_top[sequence[i4]]
                          [sequence[i3]]
                          [sequence[i4+1]]);

	if (d_12 == 1)
		d_bot = 0;
	else
		d_bot = MIN (0, dangle_bot[sequence[i1]]
                          [sequence[i2]]
                          [sequence[i2-1]]);

	if (DEBUG2)
		printf("dangling_energy_right_res: d_top(%d) = %d, d_bot(%d) = %d\n", i4+1, d_top, i2-1, d_bot);

	if (i4+1 == i2-1)     // see which is smaller
    {
        if (simple_dangling_ends)
            energy = d_top;
        else    
            energy = d_top < d_bot ? d_top : d_bot;
    }
    else if (i4+1 < i2-1)
    {
        energy = d_top + d_bot;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    return energy;
}


// NOT USED: only handles two types of brackets
void detect_original_PKed_pairs_limited(char *structure, int *p_table)
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair
{
        int i, j, struct_len;
        stack_ds st; //stach used for (
        stack_ds st_brack; // stack used for [
        h_init (&st);
        h_init (&st_brack);
    /******************************************
    **    ASK MIRELA
    *******************************************/
        //remove_space (structure);
        //printf("remove_space done \n");
    /*****************************************/
        struct_len = strlen (structure);
        for (i=0; i <= struct_len; i++)
          {
			if (i == 0){
			  p_table[i] = 0;
			}
            else if (structure[i-1] == '.'){
              p_table[i] = 0;
            }
            else if (structure[i-1] == '_'){
              p_table[i] = 0;
            }
            else if (structure[i-1] == '('){
              h_push (&st, i);
            }
            else if (structure[i-1] == '['){
              h_push (&st_brack, i);
            }
            else if (structure[i-1] == ')')
            {
                j = h_pop (&st);
                p_table[i] = j;
                p_table[j] = i;
            }
            else if (structure[i-1] == ']')
            {
                j = h_pop (&st_brack);
                p_table[i] = j;
                p_table[j] = i;
            }
          }


        if (st.top != 0 || st_brack.top != 0)
        {
            fprintf (stderr, "detect_orignal_PKed_pairs_limited::The given structure is not valid: %d more left parenthesis than right parentheses\n", st.top);
            exit (1);
        }
}

// handles 30 types of brackets
void detect_original_PKed_pairs_many (char *structure, short* p_table)
{
        int i, j, struct_len;
	// Assume brackets appear in the following order:
	//my @pk_left  = ('(','[','{','<','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z');
	//my @pk_right = (')',']','}','>','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z');
	char bl[] = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	char br[] = ")]}>abcdefghijklmnopqrstuvwxyz";

    stack_ds st1;  // (
    stack_ds st2;  // [
    stack_ds st3;  // {
    stack_ds st4;  // <
    stack_ds st5;  // a
    stack_ds st6;  // b
    stack_ds st7;  // c
    stack_ds st8;  // d
    stack_ds st9;  // e
    stack_ds st10;  // f
    stack_ds st11;  // g
    stack_ds st12;  // h
    stack_ds st13;  // i
    stack_ds st14;  // j
    stack_ds st15;  // k
    stack_ds st16;  // l
    stack_ds st17;  // m
    stack_ds st18;  // n
    stack_ds st19;  // o
    stack_ds st20;  // p
    stack_ds st21;  // q
    stack_ds st22;  // r
    stack_ds st23;  // s
    stack_ds st24;  // t
    stack_ds st25;  // u
    stack_ds st26;  // v
    stack_ds st27;  // w
    stack_ds st28;  // x
    stack_ds st29;  // y
    stack_ds st30;  // z

    h_init (&st1);
    h_init (&st2);
    h_init (&st3);
    h_init (&st4);
    h_init (&st5);
    h_init (&st6);
    h_init (&st7);
    h_init (&st8);
    h_init (&st9);
    h_init (&st10);
    h_init (&st11);
    h_init (&st12);
    h_init (&st13);
    h_init (&st14);
    h_init (&st15);
    h_init (&st16);
    h_init (&st17);
    h_init (&st18);
    h_init (&st19);
    h_init (&st20);
    h_init (&st21);
    h_init (&st22);
    h_init (&st23);
    h_init (&st24);
    h_init (&st25);
    h_init (&st26);
    h_init (&st27);
    h_init (&st28);
    h_init (&st29);
    h_init (&st30);

	struct_len = strlen (structure);
	for (i=0; i <= struct_len; i++)
	{
//printf("Analyzing %d, character %c\n",i,structure[i-1]);
		if (i == 0){
			p_table[i] = 0;
		}
		else if (structure[i-1] == '.'){
              		p_table[i] = 0;
            	}
           	else if (structure[i-1] == '_'){
			p_table[i] = 0;
            	}
            	else if (structure[i-1] == bl[0]){
//printf("Pushing %c ( for i = %d\n",structure[i-1],i);
              		h_push (&st1, i);
            	}
            	else if (structure[i-1] == bl[1]){
//printf("Pushing %c [ for i = %d\n",structure[i-1],i);
              		h_push (&st2, i);
            	}
            	else if (structure[i-1] == bl[2]){
              		h_push (&st3, i);
            	}
            	else if (structure[i-1] == bl[3]){
              		h_push (&st4, i);
            	}
            	else if (structure[i-1] == bl[4]){
              		h_push (&st5, i);
            	}
            	else if (structure[i-1] == bl[5]){
              		h_push (&st6, i);
            	}
            	else if (structure[i-1] == bl[6]){
              		h_push (&st7, i);
            	}
            	else if (structure[i-1] == bl[7]){
              		h_push (&st8, i);
            	}
            	else if (structure[i-1] == bl[8]){
              		h_push (&st9, i);
            	}
            	else if (structure[i-1] == bl[9]){
              		h_push (&st10, i);
            	}
            	else if (structure[i-1] == bl[10]){
              		h_push (&st11, i);
            	}
            	else if (structure[i-1] == bl[11]){
              		h_push (&st12, i);
            	}
            	else if (structure[i-1] == bl[12]){
              		h_push (&st13, i);
            	}
            	else if (structure[i-1] == bl[13]){
              		h_push (&st14, i);
            	}
            	else if (structure[i-1] == bl[14]){
              		h_push (&st15, i);
            	}
            	else if (structure[i-1] == bl[15]){
              		h_push (&st16, i);
            	}
            	else if (structure[i-1] == bl[16]){
              		h_push (&st17, i);
            	}
            	else if (structure[i-1] == bl[17]){
              		h_push (&st18, i);
            	}
            	else if (structure[i-1] == bl[18]){
              		h_push (&st19, i);
            	}
            	else if (structure[i-1] == bl[19]){
              		h_push (&st20, i);
            	}
            	else if (structure[i-1] == bl[20]){
              		h_push (&st21, i);
            	}
            	else if (structure[i-1] == bl[21]){
              		h_push (&st22, i);
            	}
            	else if (structure[i-1] == bl[22]){
              		h_push (&st23, i);
            	}
            	else if (structure[i-1] == bl[23]){
              		h_push (&st24, i);
            	}
            	else if (structure[i-1] == bl[24]){
              		h_push (&st25, i);
            	}
            	else if (structure[i-1] == bl[25]){
              		h_push (&st26, i);
            	}
            	else if (structure[i-1] == bl[26]){
              		h_push (&st27, i);
            	}
            	else if (structure[i-1] == bl[27]){
              		h_push (&st28, i);
            	}
            	else if (structure[i-1] == bl[28]){
              		h_push (&st29, i);
            	}
            	else if (structure[i-1] == bl[29]){
              		h_push (&st30, i);
            	}
            	else if (structure[i-1] == br[0]){
//printf("Popping %c )",structure[i-1]);
			j = h_pop (&st1);
//printf(" for pair %d.%d\n",i,j);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[1]){
//printf("Popping %c ]",structure[i-1]);
			j = h_pop (&st2);
//printf(" for pair %d.%d\n",i,j);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[2]){
			j = h_pop (&st3);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[3]){
			j = h_pop (&st4);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[4]){
			j = h_pop (&st5);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[5]){
			j = h_pop (&st6);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[6]){
			j = h_pop (&st7);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[7]){
			j = h_pop (&st8);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[8]){
			j = h_pop (&st9);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[9]){
			j = h_pop (&st10);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[10]){
			j = h_pop (&st11);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[11]){
			j = h_pop (&st12);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[12]){
			j = h_pop (&st13);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[13]){
			j = h_pop (&st14);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[14]){
			j = h_pop (&st15);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[15]){
			j = h_pop (&st16);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[16]){
			j = h_pop (&st17);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[17]){
			j = h_pop (&st18);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[18]){
			j = h_pop (&st19);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[19]){
			j = h_pop (&st20);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[20]){
			j = h_pop (&st21);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[21]){
			j = h_pop (&st22);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[22]){
			j = h_pop (&st23);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[23]){
			j = h_pop (&st24);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[24]){
			j = h_pop (&st25);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[25]){
			j = h_pop (&st26);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[26]){
			j = h_pop (&st27);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[27]){
			j = h_pop (&st28);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[28]){
			j = h_pop (&st29);
			p_table[i] = j;
			p_table[j] = i;
            	}
            	else if (structure[i-1] == br[29]){
			j = h_pop (&st30);
			p_table[i] = j;
			p_table[j] = i;
            	}
        }

        if (st1.top != 0 || st2.top != 0 || st3.top != 0 || st4.top != 0 || st5.top != 0 || st6.top != 0 || st7.top != 0
		 || st8.top != 0 || st9.top != 0 || st10.top != 0 || st11.top != 0 || st12.top != 0 || st13.top != 0
		 || st14.top != 0 || st15.top != 0 || st16.top != 0 || st17.top != 0 || st18.top != 0 || st19.top != 0
		 || st20.top != 0 || st21.top != 0 || st22.top != 0 || st23.top != 0 || st24.top != 0 || st25.top != 0
		 || st26.top != 0 || st27.top != 0 || st28.top != 0 || st29.top != 0 || st30.top != 0)
        {
            fprintf (stderr, "detect_original_PKed_pairs_many::The given structure is not valid: more left parenthesis than right parentheses\n");
            exit (1);
        }
}


void h_init (stack_ds *st)
// PRE:  None
// POST: Initialize the stack st
{
        st->top = 0;
}

void h_push (stack_ds *st, int el)
// PRE:  st is a valid stack
// POST: Push an element to the stack
{
	st->top = st->top +1;
    st->elem[st->top] = el;
}

int h_pop (stack_ds *st)
// PRE:  st is a valid stack, that is not empty
// POST: pop an element from the stack and return it
{
    if (st->top <= 0)
    {
        fprintf (stderr, "h_pop::The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
    int result = st->elem[st->top];
    st->top = st->top -1 ;
    return result;
}

int structure_type_index_PK(char type[])
{
	// FOR DP ENERGY MODEL
	// Assume the input type[] (left column) means the following (right column):
	// ps: 		# double Ps;		// pseudoloop initiation energy
	// psm: 	# double Psm;		// penalty for introducing pseudoknot inside a multiloop
	// psp: 	# double Psp;		// penalty for introducting pseudoknot inside a pseudoloop
	// pb:		# double Pb;		// penalty for band
	// pup: 	# double Pup;		// penalty for unpaired base in pseudoloop or band
	// pps: 	# double Pps;		// penalty for nested closed region inside a pseudoloop
	// stp: 	# double stP;		// multiplicative penalty for stacked pair that spans a band 
	// intp: 	# double intP;		// multiplicative penalty for internal loop that spans a band
	// a:		# double a;
	// b:		# double b;
	// c:		# double c;
	// a_p: 	# double a_p;		// penalty for introducing a multiloop that spans a band
	// b_p: 	# double b_p;		// penalty for multiloop base pair when the multiloop spans a band
	// c_p: 	# double c_p;		// penalty for unpaired base in multiloop that spans a band

	if (type[0] == 'p')
	{
		if (type[1] == 's')
		{
			if (type[2] == '\0')
				return 1;
			if (type[2] == 'm')
				return 2;
			if (type[2] == 'p')
				return 3;
		}
		if (type[1] == 'b')
			return 4;
		if (type[1] == 'u')
			return 5;
		if (type[1] == 'p')
			return 6;
	}

	if (type[0] == 's')
		return 7;  // stP
	if (type[0] == 'i')
		return 8;  // intP
	if (type[0] == 'a')
	{
		if (type[1] == '\0')
			return 9;  // a
		else
			return 12;  // a_p
	}
	if (type[0] == 'b')
	{
		if (type[1] == '\0')
			return 10;  // b
		else
			return 13;  // b_p
	}
	if (type[0] == 'c')
	{
		if (type[1] == '\0')
			return 11;  // c
		else
			return 14;  // c_p
	}
	
	printf("ERROR: invalid input to structure_type_index_PK: %s\n", type);
	return 0;
}

// fills an array string_params_PK_CC with readable form of the parameters
// used in the CC2006b model, in order that they appear in a list of
// parameters used for parameter tuning
void create_string_params_PK_CC()
{
	// coaxial stacking parameters are first
	// assumes parameters are in this order:
	// coaxstack_f_b
	// coaxstack_m1
	// coaxstack_m2

	int index = 0;
	int i,j,k,l;
	for (i=0; i < NUCL; i++) {
		for (j=0; j < NUCL; j++) {
			for (k=0; k < NUCL; k++) {
				for (l=0; l < NUCL; l++)
				{
					if (coaxstack_f_b[i][j][k][l] < INF)
					{
						sprintf(string_params_PK_CC[index],"coaxstack_f_b[%d][%d][%d][%d]",i,j,k,l);
						index++;
					}
				}
			}
		}
	}
	for (i=0; i < NUCL; i++) {
		for (j=0; j < NUCL; j++) {
			for (k=0; k < NUCL; k++) {
				for (l=0; l < NUCL; l++)
				{
					if (coaxstack_m1[i][j][k][l] < INF)
					{
						sprintf(string_params_PK_CC[index],"coaxstack_m1[%d][%d][%d][%d]",i,j,k,l);
						index++;
					}
				}
			}
		}
	}
	for (i=0; i < NUCL; i++) {
		for (j=0; j < NUCL; j++) {
			for (k=0; k < NUCL; k++) {
				for (l=0; l < NUCL; l++)
				{
					if (coaxstack_m2[i][j][k][l] < INF)
					{
						sprintf(string_params_PK_CC[index],"coaxstack_m2[%d][%d][%d][%d]",i,j,k,l);
						index++;
					}
				}
			}
		}
	}

	// read the CC2006b energy model parameters
	// assumes the parameters are in this order:
	// cc2006_s2_l1, cc2006_s1_l2, cc2006_s2_formula, cc2006_s1_formula

	for (i=0; i < CC2006_STEMSIZE; i++) {
		for (j=0; j < CC2006_LOOPSIZE; j++) {
			if (cc2006_s2_l1[i][j] < INF)
			{
				sprintf(string_params_PK_CC[index],"cc2006_s2_l1[%d][%d]",i,j);
				index++;
			}
		}
	}
	for (i=0; i < CC2006_STEMSIZE; i++) {
		for (j=0; j < CC2006_LOOPSIZE; j++) {
			if (cc2006_s1_l2[i][j] < INF)
			{
				sprintf(string_params_PK_CC[index],"cc2006_s1_l2[%d][%d]",i,j);
				index++;
			}
		}
	}
	for (i=1; i < 4; i++) {  // ignore the first row (l_min, which is not a parameter to be changed)
		for (j=0; j < CC2006_STEMSIZE_FORMULA; j++) {
			if (cc2006_s2_formula[i][j] < INF)
			{
				sprintf(string_params_PK_CC[index],"cc2006_s2_formula[%d][%d]",i,j);
				index++;
			}
		}
	}
	for (i=1; i < 4; i++) {  // ignore the first row (l_min, which is not a parameter to be changed)
		for (j=0; j < CC2006_STEMSIZE_FORMULA; j++) {
			if (cc2006_s1_formula[i][j] < INF)
			{
				sprintf(string_params_PK_CC[index],"cc2006_s1_formula[%d][%d]",i,j);
				index++;
			}
		}
	}
	// Mirela: added deltaG_assemble
	sprintf(string_params_PK_CC[index],"deltaG_assemble");
        index++;
}

int structure_type_index_PK_CC(char type[])
{
	int i = 0;
	int start = 0;
	int found = 0;

	// only coaxial, CC model params, not the DP model params
	int num_params = get_num_params_PK_CC2006b() - get_num_params_PK_DP();

	// FOR CC ENERGY MODEL (DP model not included)
	// Assumes coaxial stacking params come first, starting from index 1

	for (i=start; i < num_params; i++)
    {
        if (strcmp (type, string_params_PK_CC[i]) == 0)
        {
            found = 1;
            //printf ("%s found in %d steps\n", type, i-start+1);
            break;
        }
    }
    if (!found)
    {
		printf ("ERROR: commonPK.cpp:: structure_type_index_PK_CC - type not found: %s!!!\n", type);
        exit(1);
    }
	return i+1;

}

// from simfold/init.cpp
double ascii_to_doublePK (char *string)
// PRE:  string is either in float format or it is a '.'
// POST: convert in infinity (INF) if it is '.', in a float otherwise
{
    char *ptr;
    double en;
    if (strcmp(string,".") == 0)
        return INF;
    en = strtod(string, &ptr);
    return en;
}

// from simfold/init.cpp
int ascii_to_intPK (char *string)
// PRE:  string is either in float format or it is a '.'
// POST: convert in infinity (INF) if it is '.', in a float otherwise    
{
    char *ptr;
    double en;
    if (strcmp(string,".") == 0)
        return INF;
    en = strtod(string, &ptr);
    if (en < 0) en -= EPSILON;
    else en += EPSILON;
    return (int) (en*100);
}

// from simfold/init.cpp with my own additions for CC model
// converts input to int then multiplies by 100, then multiplies by multiplier
double ascii_to_CCdouble_PK (char *string, double multiplier)
// PRE:  string is either in float format or it is a '.'
// POST: convert in infinity (INF) if it is '.', in a float otherwise    
{
    char *ptr;
    double en;
    if (strcmp(string,".") == 0)
        return INF;
    en = strtod(string, &ptr);
    if (en < 0) en -= EPSILON;
    else en += EPSILON;
    return multiplier * (int) (en*100);
}

// from simfold/init.cpp
void giveupPK (char *string1, char *string2)
// to add: variable nb of parameters, as in scanf, printf
{
    char temp[100];
    sprintf (temp, "%s %s", string1, string2);
    perror (temp);
    exit(1);
}

double compute_PK_sensitivity (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    short ptable_ref[MaxN];
    short ptable_pred[MaxN];
    int len, i;
    double sens;
    int num_correct_bp;
    int num_true_bp;
    len = strlen(ref_structure);
	// change ref_structure to bpseq
	detect_original_PKed_pairs_many (ref_structure, ptable_ref);
	// change pred_structure to bpseq
	detect_original_PKed_pairs_many (pred_structure, ptable_pred);
	num_correct_bp = 0;
    num_true_bp = 0;
    // Mirela: for some reason, detect_original_PKed_pairs_many starts from 1, and has everything incremented by 1
    // So modify below accordingly. 
    //for (i=0; i < len; i++)
    for (i=1; i <= len; i++)
    {
        //if (ptable_ref[i] > -1)    // paired base
        if (ptable_ref[i] > 0)    // paired base
        {
            num_true_bp++;
            if (ptable_pred[i] == ptable_ref[i])
                num_correct_bp++;
        }
    }
    if (num_true_bp == 0)
        return 0.0;
    // Mirela: Feb 8, 2009, changed to 0
    //    return -1.0;
    sens = num_correct_bp*1.0/num_true_bp;
    return sens;
}
double compute_PK_ppv (char *ref_structure, char *pred_structure)
// returns 0 if undefined (denominator is 0)
{
    short ptable_ref[MAXSLEN];
    short ptable_pred[MAXSLEN];
    int len, i;
    double ppv;
    int num_correct_bp;
    int num_pred_bp;
    len = strlen(ref_structure);
	// change ref_structure to bpseq
    detect_original_PKed_pairs_many (ref_structure, ptable_ref);
	// change pred_structure to bpseq
    detect_original_PKed_pairs_many (pred_structure, ptable_pred);
    num_correct_bp = 0;
    num_pred_bp = 0;
    // Mirela: for some reason, detect_original_PKed_pairs_many starts from 1, and has everything incremented by 1
    // So modify below accordingly. 
    //for (i=0; i < len; i++)
    for (i=1; i <= len; i++)
    {
        //if (ptable_ref[i] > -1 && ptable_pred[i] == ptable_ref[i])    // paired base
        if (ptable_ref[i] > 0 && ptable_pred[i] == ptable_ref[i])    // paired base
            num_correct_bp++;
        //if (ptable_pred[i] > -1)    // paired base
        if (ptable_pred[i] > 0)    // paired base
            num_pred_bp++;
    }
    if (num_pred_bp == 0)
        return 0.0;
    // Mirela: Feb 8, 2009, changed to 0
    //    return -1.0;
    ppv = num_correct_bp*1.0/num_pred_bp;
    return ppv;
}



///////////// FUNCTIONS FROM SIMFOLD ///////////////////
// count_dangling_energy from params.cpp

void count_LEdangling_energy (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter)
//      (   )...(   )
//      i1..i2..i3..i4
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling energy between i2 and i3
// Mirela: Nov 23, 2003
// Feb 28, 2008. We might have a situation like this: <   >...(   ) or like this: (   )...<   >.
//  In that case, only add the parameter dangling onto the () pair, if it's at least 2 unpaired bases away
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    char type[100];
    int index_top, index_bot;
    int first_index;    // first index should be dangle_top[0][3][0]
    first_index = structure_type_index( (char *) "dangle_top[0][3][0]");
    
    d_top = 0;
    d_bot = 0;

    if (i2 != link && structure[i2] != '>')
      {
        //d_top = MIN (0, IGINF(dangle_top[sequence[i2]][sequence[i1]][sequence[i2+1]]));
        d_top = dangle_top[sequence[i2]] [sequence[i1]] [sequence[i2+1]];
        sprintf (type, "dangle_top[%d][%d][%d]",sequence[i2], sequence[i1], sequence[i2+1]);
        index_top = structure_type_index (type);
      }
    if (i3-1 != link && structure[i3] != '<')
      {
        //d_bot = MIN (0, IGINF(dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]]));
        d_bot = dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]];
        sprintf (type, "dangle_bot[%d][%d][%d]",sequence[i4], sequence[i3], sequence[i3-1]);
        index_bot = structure_type_index (type);
      }

    if (structure[i2] == '>' && structure[i3] == '(')   // pseudoknot, ignore dangling end dangling on it
    {
        if (i3 <= i2+2)     // >.( or >(   ignore completely
            energy = 0;
        else                // >...(
        {
            energy = d_bot;
            counter[index_bot]++;
        }
    }    
    else if (structure[i2] == ')' && structure[i3] == '<')   // pseudoknot, ignore dangling end dangling on it
    {
        if (i3 <= i2+2)     // ).< or )<   ignore completely
            energy = 0;
        else                // )...<
        {
            energy = d_top;
            counter[index_top]++;
        }
    }
    else if (structure[i2] == '>' && structure[i3] == '<')  // case >..<  ignore completely
    {
        energy = 0;
    }                                
    else if (i2+1 == i3-1 && i2 == link)
      {
        energy = d_bot;
        counter[index_bot]++;
      }
    else if (i2+1 == i3-1 && i3-1 == link)
      {
        energy = d_top;
        counter[index_top]++;
      }
    else if (i2+1 == i3-1)     // see which is smaller
    {
        //energy = d_top < d_bot ? d_top : d_bot;
        // NOTE: the comparison of d_top with d_bot is not right!
        // NO! This is not right if we don't know which of d_top and d_bot is smaller        

        // if we restrict the 3' dangling ends to be less than the 5' ones, then it's ok to do what follows

        // if we fix the dangling ends to the Turner parameters, we have to count this
        //if (d_top < d_bot) counter[index_top]++;
        //else counter[index_bot]++;

        if (simple_dangling_ends)
        {
            energy = d_top;
            counter[index_top]++;
        }
        else
        {
            if (d_top < d_bot)
            {
                energy = d_top;
                counter[index_top]++;
            }
            else
            {
                energy = d_bot;
                counter[index_bot]++;
            }
        }

        // if we introduce another variable as min, we need to do this
        //counter_min_dangle[index_top-first_index][index_bot-first_index]++;
    }
    else if (i2+1 < i3-1)
    {
        energy = d_top + d_bot;
        counter[index_top]++;
        counter[index_bot]++;
    }

    else // if there is no free base between the two branches, return 0
        energy = 0;
    //    return energy;
}

void count_LEdangling_energy_left (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter)
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i1 and i3
// Mirela: Nov 23, 2003
// Feb 28, 2008. We might have a situation like this:  (....<    >   ). In that case, only add the 3' dangling end.
// If it's (.<    > ), don't add any
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    char type[100];
    int index_top, index_bot;
    d_top = 0;
    d_bot = 0;
    int first_index;    // first index should be dangle_top[0][3][0]
    first_index = structure_type_index( (char *) "dangle_top[0][3][0]");

    // this will be used in multi-loops.
    // add the dangle_top, even if it is positive
    if (i1 != link)
      {
        //d_top = MIN (0, IGINF(dangle_top[sequence[i1]] [sequence[i2]] [sequence[i1+1]]));
        d_top = dangle_top[sequence[i1]] [sequence[i2]] [sequence[i1+1]];
        sprintf (type, "dangle_top[%d][%d][%d]",sequence[i1], sequence[i2], sequence[i1+1]);
        index_top = structure_type_index (type);
      }
    // in the other parts of the multi-loop, the dangles are added only if they are negative
    if (i3-1 != link && structure[i3] != '<')
      {
        //d_bot = MIN (0, IGINF(dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]]));
        d_bot = dangle_bot[sequence[i4]] [sequence[i3]] [sequence[i3-1]];
        sprintf (type, "dangle_bot[%d][%d][%d]",sequence[i4], sequence[i3], sequence[i3-1]);
        index_bot = structure_type_index (type);
      }

    if (structure[i3] == '<')   // pseudoknot inside, ignore dangling end dangling on it
    {
        if (i3 <= i1+2)     // (< or (.<, ignore completely
            energy = 0;
        else                // (....<
        {
            energy = d_top;
            counter[index_top]++;
        }
    }          
    else if (i1+1 == i3-1 && i1 == link)
      {
        energy = d_bot;
        counter[index_bot]++;
      }
    else if (i1+1 == i3-1 && i3-1 == link)
      {
        energy = d_top;
        counter[index_top]++;
      }
    else if (i1+1 == i3-1)     // see which is smaller
    {
        //energy = d_top < d_bot ? d_top : d_bot;
        // NOTE: the comparison of d_top with d_bot is not right!
        //if (d_top < d_bot) counter[index_top]++;
        //else counter[index_bot]++;
        //counter_min_dangle[index_top-first_index][index_bot-first_index]++;

        if (simple_dangling_ends)
        {
            energy = d_top;
            counter[index_top]++;
        }
        else
        {
            if (d_top < d_bot)
            {
                energy = d_top;
                counter[index_top]++;
            }
            else
            {
                energy = d_bot;
                counter[index_bot]++;
            }
        }
    }
    else if (i1+1 < i3-1)
    {
        energy = d_top + d_bot;
        counter[index_top]++;
        counter[index_bot]++;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    //    return energy;
}


void count_LEdangling_energy_right (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter)
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2
// Mirela: Nov 23, 2003
// Feb 28, 2008. We might have a situation like this:  (    <    >...)
//  In that case, only add the 5' dangling end if it's at least two unpaired bases away
{
    PARAMTYPE energy;
    PARAMTYPE d_top, d_bot;
    char type[100];
    int index_top, index_bot;
    int first_index;    // first index should be dangle_top[0][3][0]
    first_index = structure_type_index( (char *) "dangle_top[0][3][0]");
    
    d_top = 0;
    d_bot = 0;

    if (i4 != link && structure[i3] != '<')
      {
        //d_top = MIN (0, IGINF(dangle_top[sequence[i4]] [sequence[i3]] [sequence[i4+1]]));
        d_top = dangle_top[sequence[i4]] [sequence[i3]] [sequence[i4+1]];
        sprintf (type, "dangle_top[%d][%d][%d]",sequence[i4], sequence[i3], sequence[i4+1]);
        index_top = structure_type_index (type);
      }
    if (i2-1 != link)
      {
        //d_bot = MIN (0, IGINF(dangle_bot[sequence[i1]] [sequence[i2]] [sequence[i2-1]]));
        d_bot = dangle_bot[sequence[i1]] [sequence[i2]] [sequence[i2-1]];
        sprintf (type, "dangle_bot[%d][%d][%d]",sequence[i1], sequence[i2], sequence[i2-1]);
        index_bot = structure_type_index (type);
      }

    if (structure[i4] == '>')   // pseudoknot inside, ignore dangling end dangling on it
    {
        if (i2 <= i4+2)     // >.) or >)   ignore completely
            energy = 0;
        else                // >...)
        {
            energy = d_bot;
            counter[index_bot]++;
        }
    }          
    else if (i4+1 == i2-1 && i4 == link)
      {
        energy = d_bot;
        counter[index_bot]++;
      }
    else if (i4+1 == i2-1 && i2-1 == link)
      {
        energy = d_top;
        counter[index_top]++;
      }
    else if (i4+1 == i2-1)     // see which is smaller
    {
        //energy = d_top < d_bot ? d_top : d_bot;
        // NOTE: the comparison of d_top with d_bot is not right!
        //if (d_top < d_bot) counter[index_top]++;
        //else counter[index_bot]++;
        //counter_min_dangle[index_top-first_index][index_bot-first_index]++;
        
        if (simple_dangling_ends)
        {
            energy = d_top;
            counter[index_top]++;
        }
        else
        {
            if (d_top < d_bot)
            {
                energy = d_top;
                counter[index_top]++;
            }
            else
            {
                energy = d_bot;
                counter[index_bot]++;
            }
        }
    }
    else if (i4+1 < i2-1)
    {
        energy = d_top + d_bot;
        counter[index_top]++;
        counter[index_bot]++;
    }
    else // if there is no free base between the two branches, return 0
        energy = 0;
    //    return energy;
}
