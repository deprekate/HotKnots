/*****************************************************************
HotKnot: A heuristic algorithm for RNA secondary 
structure prediction including pseudoknots
File: score.cpp
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



#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "score.h"

#include "Stack.h"
#include "Loop.h"
#include "Input.h"
#include "Defines.h"
#include "Bands.h"
#include "LoopList.h"

// CHECK: why does this function use tstackh (hairpin) when it is a stacked pair?
int Score::sPair( int len, int* sequence, int i, int j){
	i = i - 1; 
	j = j - 1;

	int dangle = dangle_bot[sequence[j-1]][sequence[i+1]][sequence[i]] + 
		dangle_top[sequence[j-1]][sequence[i+1]][sequence[j]];

	// Mirela: bug??

	//if (tstackh[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]] >= INF)  // BUG!
	//if (tstackh[sequence[i+1]][sequence[j-1]][sequence[i]][sequence[j]] >= INF)
	if (dangle >= INF/2)
		return LOWINF;

	//return tstackh[sequence[i]][sequence[j]][sequence[i+1]][sequence[j-1]];  // BUG!
	//return -abs(tstackh[sequence[i+1]][sequence[j-1]][sequence[i]][sequence[j]]);
	return dangle;
}

int Score::F1( int len,  char* csequence, int* sequence, int i, int j){
	if (LEhairpin_loop_energy(i-1, j-1, sequence, csequence) >= INF)
		return LOWINF;
	return LEhairpin_loop_energy(i-1, j-1, sequence, csequence);
}

int Score::F2( int len, int* sequence, int i, int k, int l, int j){
	if ( (k== i+1) && (j== l+1)){
		//if (DEBUGH)
		//	printf("F2 stack (i.j = %d.%d) returns: %d\n", i, j, LEstacked_pair_energy (i-1, j-1, sequence));
		if (LEstacked_pair_energy (i-1, j-1, sequence) >= INF)
			return LOWINF;
		return LEstacked_pair_energy (i-1, j-1, sequence);
	} else {
		if (LEinternal_loop_energy (i-1, j-1, k-1, l-1, sequence) >= INF)
			return LOWINF;
		return LEinternal_loop_energy (i-1, j-1, k-1, l-1, sequence);
	}
}

int Score::intloop( int len, int* s, int i, int k, int l, int j){
	return F2(len, s, i, k, l, j);
}


char ccc(int i){
	switch(i){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		default:		// Cristina: added default case
			printf("ERROR: score.cpp::ccc - wrong int representation of RNA letter: %d\n", i);
			exit(1);
	};
};


bool Amazing(int i, int j){
	int x;
	if (i > j){
		x = i;
		i = j;
		j = x;
	};

	if ( (i == 2) && (j == 3))
		return 0;
	if ( (i == 0) && (j != 3))
		return 1;
	if ( (i == 1) && (j != 2))
		return 1;
	return 0;
}


long Score::score(int len, char* s, short* p, int TRACE, int energyModel){   
	ReadInput* input = new ReadInput(len, s, p);
	Stack* st = new Stack(input);
	Bands* ban = new Bands(input, st);
	Loop* tree = new Loop(0, len + 1, input, ban, st);
	int a, b;

	for (int i = 1; i <= input->Size; i++){
		if (input->BasePair(i) > 0) {
			if(st->Add(i, a, b)){
				tree->addLoop(a, b);
			}
		}
	}

	tree->countNumberOfChildren();  // Cristina: set number of children for the top loop

	float  f = tree->EnergyViaSimfold(energyModel) + tree->EnergyDanglingViaSimfold(energyModel);
	//        float  f = tree->Energy(energyModel) + tree->EnergyDangling();

	delete [] tree->DotParanthStructure;
	delete [] st->PrevInStack;
	delete input;
	delete st;
	delete ban;
	delete tree;
	return f;

}

Score::~Score(){}
Score::Score(){
	sc = 0;
}
