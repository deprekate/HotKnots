/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots
         File: sc.cpp
         Description: 
               Interface between folding functions and LEModel library.
 
    Date        : Oct. 16, 2004
    copyright   : (C) 2004 by Jihong Ren, Baharak Rastegari  
    email       : jihong@cs.ubc.ca, baharak@cs.ubc.ca       
******************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Defines.h"
#include "score.h"
#include "stdio.h"
#include "sc.h"

extern int TRACE;

void Char_To_Intsequence(const int len, const char*const csequence, int * sequence){
	for (int i = 0 ; i < len; i++){
		switch (csequence[i]){
			case 'A':
				sequence[i] = A;
				break;
			case 'C':
				sequence[i] = C;
				break;
			case 'G':
				sequence[i] = G;
				break;
			case 'T':
				sequence[i] = T;
				break;
			case 'U':
				sequence[i] = U;
				break;
			default:		// Cristina: added default case
				printf("ERROR: sc.cpp::Char_To_Intseqeunce - wrong letter csequence[%d] = %c\n", i, csequence[i]);
				exit(1);
		}
	}
}

float score(int len, char*s0, short*p0, int TRACE, int energyModel){
	Score sc;
	char*s = (new char[len+1]);
	short*p = (new short[len+1]);
	p[0] = 0; s[0]=0;
	for(int i=1; i<=len; ++i) {
		p[len-i+1] = p0[len-i];
		//	printf("p[%d] = %t\n", len-i+1, p0[len-i]);
		switch(s0[i-1]){
			case'A': s[i]=0;  break;
			case'C': s[i]=1;  break;
			case'G': s[i]=2;  break;
			case'U': s[i]=3;  break;
			case'T': s[i]=3;  break;  // Cristina: added 'T' case and default case
			default: printf("ERROR: sc.cpp::score() - wrong RNA letter s0[%d] = %c\n", i-1, s0[i-1]); exit(1);  
		}
	}
	float score = sc.score(len, s, p, TRACE, energyModel);  // ADDED: energyModel parameter
	delete[](s);
	delete[](p);
	return score;
}


int F2(int len, char* s, int i,  int k,  int l,  int j){
	Score sc;
	int BIGINT = 999999;
	int*sequence=(new  int[len+1]);

	Char_To_Intsequence(len, s, sequence);

	int score = sc.F2(len, sequence, i, k, l, j);
	delete[](sequence);

	if (score >= BIGINT) score = -score;
	return score;
}

int sPair(int len, char* s, int i, int j){
	Score sc;
	int BIGINT = 999999;
	int*sequence=(new  int[len+1]);    
	Char_To_Intsequence(len, s, sequence);
	int score = sc.sPair(len, sequence, i, j);

	delete[](sequence);

	if (score >= BIGINT) score = -score;
	return score;
}




