/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots
         File: Bands.cpp
         Description: 
             Contains functions for identifying and storing band regions.
 
    Date        : Oct. 16, 2004
    copyright   : (C) 2004 by Baharak Rastegari, Jihong Ren  
    email       : baharak@cs.ubc.ca, jihong@cs.ubc.ca        
******************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "Bands.h"

Bands::Bands(ReadInput * R, Stack * S){

	Input = R;
	St = S;
	pattern = new B_pattern[MaxN];
	memset(pattern, 0, MaxN*sizeof(B_pattern));

	int i;
	for (i=1; i <= Input->Size; i++){
		pattern[i].isLeftBorder = false;		pattern[i].prev = i-1;		pattern[i].next = i+1;
		pattern[i].OtherBorder = 0;
	}
	for (i=1; i<= Input->Size; i++) 		if (Input->BasePair(i) <= 0) {
			Update_links (pattern[i].prev, pattern[i].next);
		}
};

/*********************************************************************************
*********************************************************************************/
Bands::~Bands(){
	delete [] pattern;
};

/*********************************************************************************
Output: gives the list of band regions for the loop
*********************************************************************************/
void Bands::Output(int border1, int border2,int NumberOfBands){

	int i=border1;
	if (DEBUG)
		printf("bands for loop [%d, %d] are (Number of Bands = %d)", border1, border2,NumberOfBands);
	while (i<border2+1){
		if (DEBUG)
			printf("%d:%d | ", i, pattern[i].OtherBorder);
		i = pattern[i].next;
	}

	if (DEBUG)
		printf("\n");

}

/*********************************************************************************
Find_next_good_index: returns the left border of the next band region
*********************************************************************************/
int Bands::Find_next_good_index (int i, int border1, int border2){

	while (  (pattern[i].isLeftBorder == true) && (i < border2+1)){ 
		  i = pattern[i].next;
	}
	return i;

}

/*********************************************************************************
aux_Find_bands: returns the band regions of the pseudoknotted closed region,
Number of Bands and the pointer to the last band region.
*********************************************************************************/

// A band is the union of two regions [i,i_prime] and [j_prime,j]
// i and j_prime are considered to be left borders of the band region (isLeftBorder true)
// i_prime is considered to be the RightBorder of pattern[i]
// j is considered to be the RightBorder of pattern[j_prime]

// TODO: change OtherBorder to RightBorder

void Bands::aux_Find_bands(int border1, int border2, int* NumberOfBands, int* CurrentBandRegion){
	int i,j,i_prime,j_prime; // the band will be [i,i_prime] [j_prime,j]	int i_help, j_help;
	i = Find_next_good_index(border1, border1, border2);
	*NumberOfBands = 0;
	*CurrentBandRegion = 0;
	while ( i < border2+1){		
		*NumberOfBands = *NumberOfBands+1 ;
		j = Input->BasePair(i);
    	i_prime = i;
    	j_prime = j;
		while (true){
			i_help = pattern[i_prime].next;
			j_help = pattern[j_prime].prev;
			if (Input->BasePair(i_help) == j_help){
				i_prime = i_help;
				j_prime = j_help;
			}else
			  break;
		}
		// found the borders of the band region, now set them in pattern
		pattern[i].isLeftBorder = true;
		pattern[j_prime].isLeftBorder = true;
		pattern[i].OtherBorder = i_prime;
		pattern[j_prime].OtherBorder = j;
        if (*CurrentBandRegion < j_prime)
			 *CurrentBandRegion = j_prime;
		if (i != i_prime) { // if the band is not just a base pair
			Update_links(i, pattern[i_prime].next);  // the next possible left border of a region of a band, after i, is now after i_prime, since [i,i_prime] is now known to be a band
			Update_links(j_prime,pattern[j].next);   // the next possible left border of a region of a band, after j_prime, is now after j, since [j_prime,j] is now known to be a band
		}
		// check for the next possible basepairs in a band starting from i.next = i_prime.next (i.e. starting from the next possible base pair outside of [i,i_prime])
		i = Find_next_good_index(pattern[i].next, border1, border2); fflush(stdout);
	}
}

/*********************************************************************************
*********************************************************************************/
void Bands::Update_links(int from , int to){
	pattern[from].next = to;
	pattern[to].prev = from;
}


/*********************************************************************************
*********************************************************************************/
int Bands::Prev(int i){
	return pattern[i].prev;
};

/*********************************************************************************
*********************************************************************************/
int Bands::Next(int i){
	return pattern[i].next;
};
