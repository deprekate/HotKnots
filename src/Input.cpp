/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots
         File: Input.cpp
         Description: 
             Reading the input and construct the input-strucutres which will be used
             by the program!
 
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


#include "Input.h"
#include "LoopList.h"


/******************************************************************
AddPair: add a base pair to the incomplete-structures.
*******************************************************************/
void ReadInput::AddPair(int m, int n, char c){

//		printf("Indeed adding pair %d %d\n",m,n);

		int size = 0;
		c = toupper(c);
		CSequence[m] = c;
		switch (c){
			case 'A':
				type[m] = A;
				break;
			case 'C':
				type[m] = C;
				break;
			case 'G':
				type[m] = G;
				break;
			case 'T':
				type[m] = T;
				break;
			case 'U':
				type[m] = U;
				break;
			default:
				printf("ERROR: Input.cpp::AddPair() - wrong RNA letter (%c) at index %d\n", c, m);
				exit(1);
			};



		if (size < m)
			size = m;
		if (size < n)
			size = n;
		if ((m > 0) && (n > 0)){
			Sequence[m] = n;
			Sequence[n] = m;
		};

		if (m == 0)
			Sequence[n] = -1;
		if (n == 0)
			Sequence[m] = -1;

		if (Size < size)
			Size = size;
};

void ReadInput::clearDanglingRestriction() {
	for (int i = 0; i < MaxN; i++){
		cannot_add_dangling[i] = 0;
//		must_add_dangling[i] = 0;
	}
}

/*************************************************************************************
ReadInput: Takes as input the file names which contain RNA primary structure (fseq) and 
RNA secondary structure (fbseq) and construct the input-structures according to these files!
**************************************************************************************/
ReadInput::ReadInput(char * fbpseq){
	//INITIALIZING
	for (int i = 0; i < MaxN; i++){
		Sequence[i] = -1;
		CSequence[i] = ' ';
		Next[i] = 0;
		Prev[i] = 0;
		type[i] = 0;
		loops[i].pair = 0;
		looplists[i] = NULL;
		cannot_add_dangling[i] = 0;
		//must_add_dangling[i] = 0;  // Cristina: June 11 added
	}

	CSequence[0] = ' ';
	
	FILE * fileBpseq = fopen(fbpseq, "r+");
	int  m, n;
	
	char c[MaxN+1];
	int Size = 0;	
	int i=0;	
	while (!feof(fileBpseq)){
		fscanf(fileBpseq, "%d %c %d\n", &m, &c[i], &n);
		AddPair(m, n, c[i]);
		i++;
	}
    Size = i;
	fclose(fileBpseq);

	//REMOVING UNPAIRED BASES
	int Last = 0;
	for (int m = 1; m <= Size; m++)
		if (Sequence[m] > 0){
			loops[m].pair = BasePair(m);
			Prev[m] = Last;
			Next[Last] = m;
			Last = m;
	}
}


/*************************************************************************************
ReadInput: Takes as input the file names which contain RNA primary structure (fseq) and 
RNA secondary structure (fbseq) and construct the input-structures according to these files!
**************************************************************************************/
ReadInput::ReadInput(char * fseq, char * fbpseq){

//INITIALIZING
	for (int i = 0; i < MaxN; i++){
		Sequence[i] = -1;
		CSequence[i] = ' ';
		Next[i] = 0;
		Prev[i] = 0;
		type[i] = 0;
		loops[i].pair = 0;
		looplists[i] = NULL;
		cannot_add_dangling[i] = 0;
		//must_add_dangling[i] = 0;   // Cristina: June 11 added
	};

	CSequence[0] = ' ';


	FILE * fileSeq = fopen(fseq, "r+");
	if (fileSeq == NULL){
		printf("problem opening %s \n",fseq);	
	}
	FILE * fileBpseq = fopen(fbpseq, "r+");
	int  m, n;
	
	char c[MaxN+1];
	while(!feof(fileSeq)) {
		fscanf(fileSeq, "%s", c);
	}
	printf("Done reading sequence file \n");
	int i=0;	
	while (!feof(fileBpseq)){
		fscanf(fileBpseq, "%d %d\n", &m,  &n);		
		AddPair(m, n, c[i]);
		i++;
	};
	fclose(fileSeq);
	fclose(fileBpseq);

	//REMOVING UNPAIRED BASES
	int Last = 0;
	for (int m = 1; m <= Size; m++)
		if (Sequence[m] > 0){  // if base in sequence is paired...
			loops[m].pair = BasePair(m);  // set the closing base pair of the loop at m
			Prev[m] = Last;  // set the previous base pair
			Next[Last] = m;  // set the next base pair
			Last = m;
	};
};

/****************************************************************************************
ReadInput: Takes as input the length of the RNA strand, RNA primary structure (csequence) 
and RNA secondary structure (sequence) and construct the input-structures according to 
these information!
*****************************************************************************************/

ReadInput::ReadInput(int size, char* baseSequence, short* pairRefSequence){
	for (int i = 0; i < MaxN; i++){
		Sequence[i] = -1;
		CSequence[i] = ' ';
		looplists[i] = NULL;
		loops[i].pair = 0;
		Next[i] = 0;
		Prev[i] = 0;
		type[i] = 0;
		cannot_add_dangling[i] = 0;
		//must_add_dangling[i] = 0;   // Cristina: June 11 added
	};


	char c = toupper(baseSequence[0]);
	int isChar = 0;  // flag to detect how the char* input works

	if (DEBUG)
		printf("Size is %d\n", size);
	if (DEBUG)
		printf("Csequence in read input: ");

	if (c == 'A' || c == 'C' || c == 'G' || c == 'U' || c == 'T')
		isChar = 1;

	for (int i = 1; i <= size; i++){
		if (isChar == 1){
			c = toupper(baseSequence[i-1]);
		}else{
			switch ((int)baseSequence[i]){
				case A:
					c = 'A';
					break;
				case C:
					c = 'C';
					break;
				case G:
					c = 'G';
					break;
				case U:
					c = 'U';
					break;
				default:            // Added just one default case here.
					printf("ERROR: Input.cpp::ReadInput - wrong RNA letter (%c) at index %d\n", (int)baseSequence[i], i);
					exit(1);
				//case U:
				//default:
				//	c = 'U';
				//	break;
			};
		}

		//printf("Adding Pair %d %d\n", i, (int)pairRefSequence[i]);
		AddPair(i, (int)pairRefSequence[i], c);
	};
	Size = size;

	int Last = 0;
	for (int m = 1; m <= Size; m++){
		if (Sequence[m] > 0){
			loops[m].pair = BasePair(m);
			Prev[m] = Last;
			Next[Last] = m;
			Last = m;
		}
	}

	for (int i = 0; i <= Size; i++)
		ClosedRegions[i] = NULL;
};


/*****************************************************************************
RemoveUnpaired: Removes unpaired bases from the the secondary structure.
******************************************************************************/

void ReadInput::RemoveUnpaired(){
	int Last = 0;
	for (int m = 1; m <= Size; m++)
		if (Sequence[m] > 0){
			Prev[m] = Last;
			Next[Last] = m;
			Last = m;
	};

}


/*************************************************************************************
BasePair: Returns the base pair of an element. Returns -1 if the element is not paired.
**************************************************************************************/
int ReadInput::BasePair(int a){
	return Sequence[a];
};

/*********************************************************************************
*********************************************************************************/
ReadInput::~ReadInput(){  
};





