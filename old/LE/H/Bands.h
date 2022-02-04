//#include "Sequence.h"
#ifndef BANDS_H
#define BANDS_H

#include "Defines.h"
#include "Input.h"
#include "LoopList.h"
#include "Stack.h"


class Bands {

public:

//FUNCTIONS	

    Bands(ReadInput * R, Stack * S);					
	~Bands();
	void Output(int border1, int border2, int NumberOfBands);
	int Find_next_good_index (int i, int border1, int border2);
	void aux_Find_bands(int border1, int border2, int* NumberOfBands, int* CurrentBandRegion);
	void Update_links(int from , int to);


	int Prev(int i);
	int Next(int i);


//ATTRIBUTES
	ReadInput * Input;
	Stack	* St;
	B_pattern * pattern;  // holds in order from left to right, the left borders of regions associated with a band
						  // (i.e. if [i,j]|_|[i',j'] is a band interleaved with [x,y]|_|[x',y'], then 
						  // pattern[i].next = x, pattern[x].next = i', and pattern[i'].next = x'
};

#endif

