/*****************************************************************
		 HotKnot: A heuristic algorithm for RNA secondary 
			structure prediction including pseudoknots
		 File: Loop.cpp
		 Description: 
			 Contains functions for Constructing the tree of closedRegions-loops! 
			 Each node in the tree corresponds to a closed region in the RNA secondary structure. 
			 Each closed region is a loop (hairpin, interior, multi or pseudoknotted loops)
			 but not every loop is a closed region (interior-pseudo and multi-pseudoknotted loops).
			 Each node contains information which are essential for calculating the free energy
			 of all loops in the structure and therefore calculating the free energy
			 of a secondary structure.
			 
			 It also contains functions for calculating the free energy of the loops
			 and the free energy of the secodary structure.
			 
			 

 
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

#include <math.h>
#include <iostream>

#include "Loop.h"

#include "Defines.h"  // July 16 - added - includes common.h and commonPK.h
#include "common.h"
#include "commonPK.h"
#include "params.h"   // FOR PARAMETER TUNING
#include "paramsPK.h"
//#include "params.cpp"   // FOR PARAMETER TUNING

#include "externs.h"
#include "constants.h"
#include "string.h"

#define NUL 0

/******************************************************************
  Takes as input the right and left border of an identified closedRegion-loop (b, e),
  plus the stack S (which contains the information about the current status of the closed 
  region finding program).
  It then creates a node corresponding to this closedRegion-loop!
 *******************************************************************/
Loop::Loop(int b, int e,  ReadInput * R, Bands * B, Stack * S){
	ILoops = NUL;
	MLoops = NUL;
	St = S;
	LeftSibling   = NUL;
	RightChild    = NUL;
	Parent		= NUL;
	begin         = b;
	end           = e;
	Input			= R;
	bandpattern	= B;
	type = external;
	NumberOfChildren	= 0;
	NumberOfUnpaird	= 0;
	NumberOfUnpairedInPseudo = 0;
	NumberOfUnpairedInUnbandChild = 0;

	NumberOfUnBandChild = 0;
	NumberOfBandChild = 0; 
	nested = nothing;

	finalEnergy = 0;
	finalCoaxEnergy = 0;

	DotParanthStructure = new char[R->Size];
}

/******************************************************************
  Clears the loop structure so that it can be used again with another
  energy model. Only element cleared is number of unpaired bases in a
  pseudoloop.
 *******************************************************************/
/*
   void Loop::clearLoop(){
   NumberOfUnpairedInPseudo = 0;
   }
   */

/*********************************************************************************
 *********************************************************************************/
Loop::~Loop(){
	if (RightChild != NUL)
		delete RightChild;
}

/********************************************************************************
setDotParanthStruct:  sets the DotParanthStructure
 ********************************************************************************/
void Loop::setDotParanthStruct(char * structure) {
	for (int i = 0; i < (int) strlen(structure); i++){
		DotParanthStructure[i] = structure[i];
	}
	DotParanthStructure[strlen(structure)] = '\0';
}

/*********************************************************************************
WhereLocated: It figures out the location status of the loop (in-Band or out-Band).
span-Band loops are identified when a pseudoknotted loops is identified 
in FindInnerLoops function.
 *********************************************************************************/
void Loop::WhereLocated(){
	int previousRegion;
	//Starting from the last band region, trying to find out the location status of the loop
	while (true){
		if (( end <= bandpattern->pattern[Parent ->CurrentBandRegion].OtherBorder) && (end >= Parent->CurrentBandRegion) ) {
			if (begin >= Parent->CurrentBandRegion) {
				nested = inBand;  // completely in the band
				return;
			}
			// else it is cross-Band and will identified in FindInnerLoops function
		} else {  // moving to the previous band region
			previousRegion = bandpattern->pattern[Parent ->CurrentBandRegion].prev;
			if (( end > bandpattern->pattern[previousRegion].OtherBorder) && (end < Parent->CurrentBandRegion) ) { 
				nested = unBand;  // inside pseudoloop but not in band (i.e. between band borders)
				return;
			} else {
				Parent ->CurrentBandRegion = previousRegion;				
			}
		}
	}

}

/*********************************************************************************
FindInnerLoops: It figures out the interior-pseudoknotted and multi-pseudoknotted 
loops (span-band loops) nested in a pseudoknotted closed region.
 *********************************************************************************/
void Loop::FindInnerLoops(int border1, int border2){		
	int i=border1;
	while (i <= border2){  // check all left band borders
		int j = bandpattern->pattern[i].OtherBorder;
		if  (Input->BasePair(i) > i){  // ignore left borders of the second region associated with a band (i.e. ignore jp)
			int ip = Input->BasePair(i);
			int jp = Input->BasePair(j);
			//[i,j]|_| [jp,ip] is a band and i < i'
			int x = j;
			while (x > i){  // check if x strictly greater than i to prevent y from being set to -1
				// EITHER [y,x] |_| [xp,yp] is a candidate for an internal loop spanning the band [i,j]|_| [jp,ip]
				// OR     [y,x] |_| [xp,yp] is a candidate for a multiloop (with other closed regions in [y,x] or [xp,yp]) spanning the band [i,j]|_| [jp,ip]
				int y = St->PrevInStack[x];  // PrevInStack[x] = the left border of the closed region that x potentially belongs to, without knowing anything about what x is paired with
				// See the definition of 'the left border y of a potentially closed region from the viewpoint of x' - ENABLE USING DEBUG2
				int xp = Input->BasePair(x);
				int yp = Input->BasePair(y);
				// check that base pair x.xp spans the band [i,j]|_|[jp,ip]
				if (  (xp >= jp) && (yp >= jp) && (xp >= St->PrevInStack[yp]) ){
					// check for internal loops (i.e. no closed regions between [y,x] and [xp,yp] band)
					if ( (x == Input->Next[y]) && (yp == Input->Next[xp])){  // Next = next paired base
						// Internal loops that spans a band (interior-pseudoknotted loop) has been found and should be added to ILoops list
						if  (ip > i){
							T_IntList * IntList = new T_IntList;
							IntList->Num = y;
							IntList->Next = NUL;
							IntList->tuning_flag = 0;  // PARAMETER TUNING
							Input->looplists[y] = new LoopList(Input, y, x);
							IntList->Next = ILoops;
							ILoops = IntList;
						};
					}else{
						// Multiloop that spans a band (multi-pseudoknotted loop) has been found and should be added to MLoops list
						// Call the FindChildren function to finding out the tuples in this multi-pseudoknotted loop
						//MultiLoops->Add(y, x);
						T_IntList * IntList = new T_IntList;
						IntList->Num = y;
						IntList->Next = NUL;
						IntList->tuning_flag = 0;  // PARAMETER TUNING
						Input->looplists[y] = new LoopList(Input, y, x);
						Input->looplists[y]->FindChildren();
						IntList->Next = MLoops;
						MLoops = IntList;
					};
				};
				x = y;
			};
		}
		// set i to the next left border of all regions associated with a band 
		// (i.e. i and ip are left borders of band [i,j] |_| [ip,jp]
		i = bandpattern->pattern[i].next;
	}
}

/*********************************************************************************
PseudoNestedCheck: It checks whether the type of a closed region is pseudoknotted
and if then call appropriate function to find the location status of its children
besides calculating the number of in-Band and un-Band children. The number of
unpaired bases in a pseudoloop is also set up here (by subtracting the number of
bases in nested closed regions).
 *********************************************************************************/
void Loop::PseudoNestedCheck(){
	if (type != pseudo)
		return;
	Loop * L = RightChild;
	while (L != NUL){
		L->WhereLocated();	 // sets L->nested to inBand or unBand
		if (L->nested == inBand){
			NumberOfBandChild = NumberOfBandChild + 1;
		} else if (L->nested == unBand){
			NumberOfUnBandChild = NumberOfUnBandChild +1;
		}
		if (L->nested == unBand)
			NumberOfUnpairedInUnbandChild -= (L->end - L->begin + 1);  // number of unpaired bases in pseudoloop must not include bases 
		L = L->LeftSibling;
	};
}


/*******************************************************************************
LoopSetType: It figures out type of the loop based on the borders and the 
number/type of the children.
 ********************************************************************************/
void Loop::loopSetType(){
	if (Input->BasePair(begin) != end){
		type = pseudo;
		return;
	};
	if (!RightChild){
		type = hairpin;
		return;
	};
	// if only one child and it's not pseudoknotted...
	if ( (!RightChild->LeftSibling) && (Input->BasePair(RightChild->begin) == RightChild->end)){
		// ...and if child's base pair is stacked with this loop's pair, then stacked
		if ((RightChild->begin == begin + 1) && (RightChild->end == end - 1 )){
			type = stackloop;
			return;
		};
		// ...and if not stacked, then interior loop
		type = interior;
		return;
	};
	// if there is more than one child or the child is pseudoknotted, then multiloop
	if ( (RightChild->LeftSibling) 	||	(Input->BasePair(RightChild->begin) != RightChild->end)) {
		type = multi;
		return;
	};
	// if it's none of the above, then it's external loop
	type = external;
	return;    
}

/*******************************************************************************
LoopSetType: It assigns to each loop its type by calling setTypes() function.
 ********************************************************************************/
void Loop::setTypes(){
	loopSetType();
	Loop * L = RightChild;
	if (!L)
		return;
	while (L){
		L->setTypes();
		L = L->LeftSibling;
	};    
}

/*******************************************************************************
UnpairedBases: It calculates the number of unpaired bases in a closed region!
 ********************************************************************************/
int UnpairedBases(ReadInput* Input, int from, int to){
	int k = 0;
	for (int j = from; j < to; j++)
		if (Input->BasePair(j) <= 0)
			k++;
	return k;
}

void Loop::countNumberOfChildren(){
	int numChildren = 0;
	Loop * L1 = RightChild;
	while (L1 != NUL){
		L1 = L1->LeftSibling;
		numChildren++;
	}
	NumberOfChildren = numChildren;
}

/*******************************************************************************
addLoop: It adds an identified closedRegion/loop to the tree (of closed region)!
And calls appropriate functions for setting different attributes of the loop.
 ********************************************************************************/
void Loop::addLoop(int begin, int end){
	// T is the root of the tree
	// L1 is new loop to add
	Loop *L1 = new Loop(begin, end, Input, bandpattern, St);
	Input->ClosedRegions[begin] = L1;  // ClosedRegions[begin] = closed region/loop L1 starting at begin
	int last = end;
	int first = begin;
	// by default, set L1 to be child of T
	L1->Parent = this;

	// If L1 is the first loop to be added to the tree
	if (!RightChild){
		RightChild = L1;
		if (begin != end)
			L1->NumberOfUnpaird = UnpairedBases(Input, first, last);
		L1->loopSetType();

		L1->LoopFindBands();
		L1->PseudoNestedCheck();
		return;
	};
	Loop *L = RightChild;  // L is the rightmost child of T
	int HaveChild = (L->begin >= begin);  // is L supposed to be a child of L1 (L.begin >= L1.begin)?
	// Note: if L.begin >= L1.begin, then L.end < L1.end (i.e. L nested in L1). Otherwise, L1.end < L.begin
	// and L1 would have already been added as another child of T and a sibling of L. 

	if (HaveChild)
		L1->RightChild = L;  // set L to be a child of L1

	Loop *L2 = L;
	Input->loops[begin].num_branches = 0;
	Input->loops[begin].pseudo_num_branches = 0;

	// Calculating the number of tuples (branches) in the closed region.
	// Calculating the number of tuples by assgining the value equal to 2 to 
	// pseudoknotted tuples.
	// Calculating the number of children for the closed region while making the appropriate father-child
	// relationship between them.

	// add new node to tree: need to shift existing children (L) to become children of new node (L1)
	//                       if they are nested inside (i.e. L.begin >= L1.begin)
	while (L && (L->begin >= begin)){
		// set the base number (L->begin) of the next branch (Input->loops[begin].num_branches) in loop L1
		Input->loops[begin].bri[Input->loops[begin].num_branches++] = L->begin;
		if (L->type == pseudo)
			Input->loops[begin].pseudo_num_branches += 2;  // add 2 since have 2 base pairs in a pseudoloop
		else
			Input->loops[begin].pseudo_num_branches += 1;  // add 1 since have 1 base pair in a non-pseudoloop

		if (begin != end)
			L1->NumberOfUnpaird += UnpairedBases(Input, L->end, last);
		last = L->begin;

		// Note: if there is a 2nd, 3rd, etc child, they are added as left siblings of the 1st child (RightChild)
		L1->NumberOfChildren++;
		L2 = L;  // this is used to keep track of the last non-null L considered in the loop
		L->Parent = L1;  // udpate parent child relationship (L1 parent of L)

		L = L->LeftSibling;  // update L to be the next rightmost child of T
	};

	// Assigning the appropriate tuples (branches) to the closed region
	// The closed regions were assigned to branches in opposite order in the above while loop
	// (i.e. the rightmost child -- region with largest left border -- is assigned to branch 0).
	// This loop reverses the order so that branch 0 is the region with the smallest left border. 
	for (int index = 0; index < Input->loops[begin].num_branches; index++){
		int index2 = Input->loops[begin].num_branches - index - 1;
		if (index < index2){
			int x = Input->loops[begin].bri[index];
			Input->loops[begin].bri[index] = Input->loops[begin].bri[index2];
			Input->loops[begin].bri[index2] = x;
		};
	};

	// DONECHECK: if this is here, do you need the calls to UnpairedBases inside the while loop?
	//        i.e. could you have only one call to UnpairedBases(Input, begin, end)
	// DONE: either implementation is ok

	// this line adds the last set of unpaired bases: between begin and leftmost child's begin
	if (begin != end)
		L1->NumberOfUnpaird += UnpairedBases(Input, first, last);

	// make L1's new left sibling the last L that was considered in the while loop
	// (i.e. the left sibling from the first row, from where the new children of L1 were taken)
	if (L)
		L1->LeftSibling = L;

	// for the last valid L that was moved to become a child of L1, make L's left sibling NUL,
	// since there were no more L's moved beside it (note that the rest of the L's valid during
	// an iteration of the while loop above keep the same left siblings since they were all moved
	// to become children of L1 until this last L)
	if (L2 && HaveChild)
		L2->LeftSibling = NUL;  

	RightChild = L1;  // add L1 as rightmost child of T

	L1->loopSetType();

	L1->LoopFindBands();
	L1->PseudoNestedCheck();
	delete [] L1->DotParanthStructure;
}

/*********************************************************************************
Print: It prints out the parsing tree!
 *********************************************************************************/

void Loop::Print(int i){
	for (int j = 0; j < i; j++)
		printf(".");
	if (Parent)
		printf("[%d, %d] nested: %d\n",begin, end, nested);
	Loop * L = RightChild;
	if (L){
		do {
			L->Print(i+1);
			L = L->LeftSibling;
		} while (L != NUL);

	};
}

/*********************************************************************************
LoopFindBands: It finds the bands of a pseudoknotted region and stores them
in bandpattern.
 *********************************************************************************/
void Loop::LoopFindBands(){

	// if loop is not a pseudoknot or it begins at 0, there are no bands borders to find
	// so simply update the next potential left band border (bandpattern array) to the next base pair in preparation for next iteration
	if ( (type != pseudo) || (begin == 0)){
		if (begin > 0)
			bandpattern->Update_links(bandpattern->Prev(begin) , bandpattern->Next(end));
		return;
	};

	// if the loop is a pseudoknot, find the bands (aux_Find_Bands) and update the band borders (Update_links)
	// note: bandpattern->pattern contains the list of the next left band border of each of the two regions associated with a band
	bandpattern->aux_Find_bands(begin, end, &NumberOfBands, &CurrentBandRegion);
	bandpattern->Update_links(bandpattern->Prev(begin) , bandpattern->Next(end));


	// find internal and multiloops that span a band (interior-pseudoknotted and multi-pseudoknotted loops)
	FindInnerLoops(begin, end);
}

/*********************************************************************************
FindSpanBandLoops: Enumerate the loops that span the bands associated with closed
the closed region that this loop represents i.e. region [begin, end].
 *********************************************************************************/
void Loop::FindSpanBandLoops() {

	// band borders of this region are known

	// visit closed regions tree in postfix order, using the list of base pairs in the sequence
	// to add elements to the BL list, if they span the band

	// for each element in the BL list, mark down the base number, the base pair, prev, next, nested

	// scan BL from left to right, and for each element k, determine the type of span-band loop it is
}

/****************************************************
 *           ENERGY CALCULATION PART					*
 *													*
 *****************************************************/

/*********************************************************************************
 *********************************************************************************/
char c(int i){
	switch(i){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
	};
	return 0;
};

void Loop::printEnergyTrace(){

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)
	Loop * L = RightChild;
	while (L != NUL){
		L->printEnergyTrace();
		L = L->LeftSibling;
	};

	switch (type){
		case	stackloop:
			printf("***[%d, %d] Stacked Pair Energy = %f cal/mol\n", begin, end, finalEnergy);
			break;
		case	hairpin	:
			printf("***[%d, %d] Hairpin Energy = %f cal/mol\n", begin, end, finalEnergy);
			break;
		case	interior:
			printf("***[%d, %d] Internal Loop Energy = %f cal/mol\n", begin, end, finalEnergy);
			break;
		case	multi	:
			printf("***[%d, %d] Multiloop Energy = %f cal/mol\n", begin, end, finalEnergy);
			break;
		case	pseudo	:
			printf("***[%d, %d] Pseudoloop Energy = %f cal/mol\n", begin, end, finalEnergy);
			break;
		default:
			break;
	};

	if (RightChild != NUL && RightChild->finalCoaxEnergy != 0)
	{
		if (type == external)
			printf("***[1, %d] coaxial energy term added for all its children = %f cal/mol\n", Input->Size, RightChild->finalCoaxEnergy);
		else
			printf("\t + coaxial energy term added for all its children = %f cal/mol\n", RightChild->finalCoaxEnergy);
	}
}


/*********************************************************************************
stackEnergy: calls the appropriate function from another program (simFold) to
calculate the free energy of stacked pair. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
PARAMTYPE	Loop::stackEnergy(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = begin;

	PARAMTYPE en = LEstacked_pair_energy (i, f[i].pair, sequence);
	//	int en = s_stacked_pair::get_energy(i, f[i].pair, sequence);

	return en;
}


/*********************************************************************************
hairpinEnergy: calls the appropriate function from another program (simFold) to
calculate the free energy of a hairpin Loop. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
PARAMTYPE	Loop::hairpinEnergy(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = begin;
	char * csequence = Input->CSequence;
	PARAMTYPE en = LEhairpin_loop_energy (i, f[i].pair, sequence, csequence);

	return en;
}


/*********************************************************************************
interiorEnergy: calls the appropriate function from another program (simFold) to
calculate the free energy of an interior Loop. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
PARAMTYPE	Loop::interiorEnergy(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = begin;

	int ip, jp;
	ip = f[i].bri[0];  // the base index representing the first branch of this loop
	jp = f[f[i].bri[0]].pair;  // the base pair of ip
	Input->cannot_add_dangling[ip-1] = 1;
	Input->cannot_add_dangling[jp+1] = 1;
	PARAMTYPE en = LEinternal_loop_energy (i, f[i].pair, ip, jp, sequence);

	return en;

}


/*********************************************************************************
Energy: Calculate the free energy of RNA secondary structure by calling
getEnergy function for calculating the energy associated with the loops and
then adding some other penalties and other energy values corresponding to the
secondary structure (where some of them construct the free energy of the external loop).
The main part is calling functions from simFold program and is changed slightly to consider
pseudoknotted substructures.
 *********************************************************************************/
float Loop::Energy(int model){

	if (model == DP)
		return -10*getEnergyDP();  // returns energy in 10 x 10*cal/mol = cal/mol
	if (model == RE)
		return -10*getEnergyRE();  // returns energy in 10 x 10*cal/mol = cal/mol
	if (model == CC2006a)
		return -10*getPartialCoaxialEnergy(CC2006a) - 10*getEnergyCC2006a();  // returns energy in 10 x 10*cal/mol = cal/mol
	//		return -10*getEnergyCC2006a();
	if (model == CC2006b)
		//		return -10*getPartialCoaxialEnergy(CC2006b) - 10*getEnergyCC2006b();  // returns energy in 10 x 10*cal/mol = cal/mol
		return -10*getEnergyCC2006b();
	if (model == CC2006c)
		return -10*getPartialCoaxialEnergyAll(CC2006b) - 10*getEnergyCC2006b();  // returns energy in 10 x 10*cal/mol = cal/mol
	// *** Add new energy model code here *** //

	printf("ERROR: invalid energy model input\n");
	return 0;

	// the dangling energies and AU penalties are calculated in EnergyDangling();
}

// Energy via general Simfold function 
// (same one used in Energy(int model, double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles))
// Ignores dangling ends if no_pk_dangling_ends is set (constantsPK.cpp)
// Returns energy in cal/mol
// DP and CCb work

float Loop::EnergyViaSimfold(int model){
	double free_value = 0;
	int reset_c = 0;
	int ignore_dangles = no_pk_dangling_ends;

	//	printf("reset_c = %d, ignore_dangles = %d  model & DP = %d %d\n", reset_c, ignore_dangles, model, DP);

	if (model == DP){
		return 1000*Energy(model, NUL, NUL, free_value, reset_c, ignore_dangles);
	}
	if (model == RE)
		return -10*getEnergyRE();  // returns energy in 10 x 10*cal/mol = cal/mol
	if (model == CC2006a)
		return -10*getPartialCoaxialEnergy(CC2006a) - 10*getEnergyCC2006a();  // returns energy in 10 x 10*cal/mol = ca$
		//return -10*getEnergyCC2006a();
	if (model == CC2006b)
		return 1000*Energy(model, NUL, NUL, free_value, reset_c, ignore_dangles);
		//return -10*getPartialCoaxialEnergy(CC2006b) - 10*getEnergyCC2006b();  // returns energy in 10 x 10*cal/mol = ca$
		//return -10*getEnergyCC2006b();
	if (model == CC2006c) 
		return -10*getPartialCoaxialEnergyAll(CC2006b) - 10*getEnergyCC2006b();  // returns energy in 10 x 10*cal/mol =$
	// *** Add new energy model code here *** //

	printf("ERROR: invalid energy model input\n");
	return 0;
	// the dangling energies and AU penalties are calculated in EnergyDangling();
}


/*********************************************************************************
Energy: FOR PARAMETER TUNING Calculate the free energy of RNA secondary structure by calling
getEnergy function for calculating the energy associated with the loops and
then adding some other penalties and other energy values corresponding to the
secondary structure (where some of them construct the free energy of the external loop).
The main part is calling functions from simFold program and is changed slightly to consider
pseudoknotted substructures. Returns energy in kcal/mol. DP works.
 *********************************************************************************/
float Loop::Energy(int model, double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles){

	int num_params = 0;
	if (model == DP)
		num_params = get_num_params_PK_DP();
	else if (model == CC2006b)
		num_params = get_num_params_PK_CC2006b();

	if (num_params == 0)
		printf("WARNING: Loop.cpp::Energy(int model, double P_matrix, ...) - num_params not initialized\n");

	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	if (model == DP)
		return -getEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);  // returns energy in kcal/mol

	// RE has not been modified to work with parameter tuning
	//	if (model == RE)
	//		return 10*getEnergyRE();  // returns energy in kcal/mol

	if (model == CC2006b)
		return -getEnergyCC2006b(P_matrix, c, f, reset_c, ignore_dangles);  // returns energy in kcal/mol

	return 0;

	/* Add more energy models here */

	// the dangling energies and AU penalties are calculated in EnergyDangling();
}



int FLOOR_OVER2(int i)
{
	//	printf("FLOOR OVER 2 received %d, modulo %d\n", i, i%2);
	if ((i % 2) == 0)  // i is even
		return i/2;
	return (i-1)/2;
}

void Loop::removeDangling(int i, LoopType loop_type, pk_coax_features oneCoaxStack, int stack_dangle)
	// Pre: i is from 0 to NumberOfChildren
	//      coax_energy = allCoaxStacks[i].coax_energy
{
	if (stack_dangle != -1 && stack_dangle != 0 && stack_dangle != 1)
		printf("ERROR: stack_dangle = %d not set properly; must be 0,1, or -1\n", stack_dangle); 

	switch (loop_type)
	{
		case multi:
			if (DEBUG2)
				printf("removeDangling: multi start\n");
			if (i == 0)  // outer base pair stacks with right child
			{
				if (oneCoaxStack.loop0->end == oneCoaxStack.loop1->end + 1)  //L1,L2 flush coaxial stacking
				{
					if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
					{
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  //L2,  only change array if the dangling end was free originally
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
						if (Input->BasePair(oneCoaxStack.loop0->begin + 1) < 0)  //L1,  only change array if the dangling end was free originally
							Input->cannot_add_dangling[oneCoaxStack.loop0->begin + 1] = 1;
					}
				}
				else if ( (oneCoaxStack.loop0->end == oneCoaxStack.loop1->end + 2) && ( (Input->BasePair(oneCoaxStack.loop0->begin + 1) < 0) || (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0) ) )  // mismatch coaxial stacking
				{
					if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
					{
						// make sure the dangling end involved in stacking is not included elsewhere
						if (RES_STACK_DANGLE == 1)
						{
							if (stack_dangle == 0 && Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  //L2,  only change array if the dangling end was free originally
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
							if (stack_dangle == 1 && Input->BasePair(oneCoaxStack.loop0->begin + 1) < 0)  //L1,  only change array if the dangling end was free originally
								Input->cannot_add_dangling[oneCoaxStack.loop0->begin + 1] = 1; 
							if (stack_dangle != 0 && stack_dangle != 1)
								printf("ERROR: stack_dangle = %d not set properly for mismatch loop; must be 0 or 1\n", stack_dangle); 
						}
						else  // restrict both dangling ends
						{
							if (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  //L2,  only change array if the dangling end was free originally
							{
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
								//							printf("multi: restriction for %d\n", Input->BasePair(oneCoaxStack.loop1->end) - 1);
							}
							if (Input->BasePair(oneCoaxStack.loop0->begin + 1) < 0)  //L1,  only change array if the dangling end was free originally
							{
								Input->cannot_add_dangling[oneCoaxStack.loop0->begin + 1] = 1; 
								//							printf("multi: restriction for %d\n", oneCoaxStack.loop0->begin + 1);
							}
						}

						// make sure dangling end between the two loops (the mismatch base) is not included as a dangling end elsewhere
						if (Input->BasePair(oneCoaxStack.loop1->end + 1) < 0)
						{
							Input->cannot_add_dangling[oneCoaxStack.loop1->end + 1] = 1;
							//						printf("multi: restriction for %d\n", oneCoaxStack.loop1->end + 1);
						}
						else
						{
							printf("WARNING! removingDangling(), multi, i=0: dangling end %d is not free, yet it was seen to be in getPartialCoaxialStacking!\n", oneCoaxStack.loop1->end + 1);
						}
					}
				}
			}
			else if (i == Parent->NumberOfChildren)  // outer base pair stacks with left child
			{
				if (oneCoaxStack.loop0->begin == oneCoaxStack.loop1->begin + 1)  // flush coaxial stacking
				{
					if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
					{
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(oneCoaxStack.loop1->end - 1) < 0)  // only change array if the dangling end was free originally
							Input->cannot_add_dangling[oneCoaxStack.loop1->end - 1] = 1;
						if (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
					}
				}
				else if ( (oneCoaxStack.loop0->begin == oneCoaxStack.loop1->begin + 2) && ( (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0) || (Input->BasePair(oneCoaxStack.loop1->end - 1) < 0) ) )  // mismatch coaxial stacking (dangling ends are free bases)
				{
					if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
					{
						if (RES_STACK_DANGLE == 1)
						{
							// make sure the dangling end involved in stacking is not included elsewhere
							if (stack_dangle == 0 && Input->BasePair(oneCoaxStack.loop1->end - 1) < 0)  // only change array if the dangling end was free originally
								Input->cannot_add_dangling[oneCoaxStack.loop1->end - 1] = 1;
							if (stack_dangle == 1 && Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
							if (stack_dangle != 0 && stack_dangle != 1)
								printf("ERROR: stack_dangle = %d not set properly for mismatch loop; must be 0 or 1\n", stack_dangle); 
						}
						else  // restrict both dangling ends
						{
							if (Input->BasePair(oneCoaxStack.loop1->end - 1) < 0)  // only change array if the dangling end was free originally
							{
								Input->cannot_add_dangling[oneCoaxStack.loop1->end - 1] = 1;
								//							printf("multi: restriction for %d\n", oneCoaxStack.loop1->end - 1);
							}
							if (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
							{
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
								//							printf("multi: restriction for %d\n", Input->BasePair(oneCoaxStack.loop0->begin) + 1);
							}
						}

						// make sure dangling end between the two loops (the mismatch base) is not included as a dangling end elsewhere
						if (Input->BasePair(oneCoaxStack.loop1->begin + 1) < 0)
						{
							Input->cannot_add_dangling[oneCoaxStack.loop1->begin + 1] = 1;
							//						printf("multi: restriction for %d\n", oneCoaxStack.loop1->begin + 1);
						}
						else
						{
							printf("WARNING! removingDangling(), multi, i=NumChildren: dangling end %d is not free, yet it was seen to be in getPartialCoaxialStacking!\n", oneCoaxStack.loop1->begin + 1);
						}
					}
				}
			}
			else  // two children in a multiloop stack (not the outer base pair)
			{
				if (oneCoaxStack.loop0->begin == oneCoaxStack.loop1->end + 1)  //L1,L2  flush coaxial stacking
				{
					if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
					{
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  // only change array if the dangling end was free originally; also do array bounds check
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
						if (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
					}

				}
				else if ( (oneCoaxStack.loop0->begin == oneCoaxStack.loop1->end + 2) && ( (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0) || (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0) ) )  // mismatch coaxial stacking
				{
					if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
					{
						if (RES_STACK_DANGLE == 1)
						{
							// make sure the dangling end involved in stacking is not included elsewhere
							if (stack_dangle == 0 && Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  // only change array if the dangling end was free originally
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
							if (stack_dangle == 1 && Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
							if (stack_dangle != 0 && stack_dangle != 1)
								printf("ERROR: stack_dangle = %d not set properly for mismatch loop; must be 0 or 1\n", stack_dangle); 
						}
						else  // restrict both dangling ends
						{
							if (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  // only change array if the dangling end was free originally
							{
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
								//							printf("multi: restriction for %d\n", Input->BasePair(oneCoaxStack.loop1->end) - 1);
							}
							if (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
							{
								Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
								//							printf("multi: restriction for %d\n", Input->BasePair(oneCoaxStack.loop0->begin) + 1);
							}
						}
						// make sure dangling end between the two loops (the mismatch base) is not included as a dangling end elsewhere
						if (Input->BasePair(oneCoaxStack.loop1->end + 1) < 0)
						{
							Input->cannot_add_dangling[oneCoaxStack.loop1->end + 1] = 1;
							//						printf("multi: restriction for %d\n", oneCoaxStack.loop1->end + 1);
						}
						else
						{
							printf("WARNING! removingDangling(), multi, two children: dangling end %d is not free, yet it was seen to be in getPartialCoaxialStacking!\n", oneCoaxStack.loop1->end + 1);
						}
					}
				}
			}
			break;
		case external:
			if (DEBUG2)
				printf("removeDangling: external start\n");

			//printf("oneCoaxStack.loop0->begin = %d\n", oneCoaxStack.loop0->begin);
			//printf("oneCoaxStack.loop1->end = %d\n", oneCoaxStack.loop1->end);
			if (oneCoaxStack.loop0->begin == oneCoaxStack.loop1->end + 1)  //L1,L2  flush coaxial stacking
			{
				if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
				{
					//printf("Input->BasePair(oneCoaxStack.loop1->end) - 1 = %d", Input->BasePair(oneCoaxStack.loop1->end) - 1);
					//printf("Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) = %d", Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1));
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  // only change array if the dangling end was free originally; also do array bounds check
						Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
					if (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
				}

			}
			else if ( (oneCoaxStack.loop0->begin == oneCoaxStack.loop1->end + 2) && ( (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0) || (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0) ) )  // mismatch coaxial stacking
			{
				//			printf("here2\n");
				if (oneCoaxStack.coax_energy != 0)  // coaxial stacking favourable over dangling ends
				{
					if (RES_STACK_DANGLE == 1)
					{
						// make sure the dangling end involved in stacking is not included elsewhere
						if (stack_dangle == 0 && Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  // only change array if the dangling end was free originally
						{
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
							//						printf("external: restriction for %d\n", Input->BasePair(oneCoaxStack.loop1->end) - 1);
						} 
						if (stack_dangle == 1 && Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
						{
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
							//						printf("external: restriction for %d\n", Input->BasePair(oneCoaxStack.loop0->begin) + 1);
						}
						if (stack_dangle != 0 && stack_dangle != 1)
							printf("ERROR: stack_dangle = %d not set properly for mismatch loop; must be 0 or 1\n", stack_dangle); 
					}
					else
					{
						if (Input->BasePair(Input->BasePair(oneCoaxStack.loop1->end) - 1) < 0)  // only change array if the dangling end was free originally
						{
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop1->end) - 1] = 1;
							//						printf("external: restriction for %d\n", Input->BasePair(oneCoaxStack.loop1->end) - 1);
						}
						if (Input->BasePair(Input->BasePair(oneCoaxStack.loop0->begin) + 1) < 0)  // only change array if the dangling end was free originally
						{
							Input->cannot_add_dangling[Input->BasePair(oneCoaxStack.loop0->begin) + 1] = 1;
							//						printf("external: restriction for %d\n", Input->BasePair(oneCoaxStack.loop0->begin) + 1);
						}
					}

					// make sure dangling end between the two loops (the mismatch base) is not included as a dangling end elsewhere
					if (Input->BasePair(oneCoaxStack.loop1->end + 1) < 0)
					{
						Input->cannot_add_dangling[oneCoaxStack.loop1->end + 1] = 1;
						//					printf("external: restriction for %d\n", oneCoaxStack.loop1->end + 1);
					}
					else
					{
						printf("WARNING! removingDangling(), external: dangling end %d is not free, yet it was seen to be in getPartialCoaxialStacking!\n", oneCoaxStack.loop1->end + 1);
					}
				}
			}
			break;
		default:
			break;
	}
}


/*********************************************************************************
getPartialCoaxialEnergy: Calculate the coaxial stacking energies between children
in the loop tree (not between structural elements within those children, e.g. stems
of pseudoknots, which are included in individual energy calculations, 
e.g getPseudoEnergy()). This function allows coaxial stacking between adjacent pairs 
of stems only. If one stem stacks with another stem, it cannot stack with a third. 
(The exception is if the child is a multiloop of pseudoloop, in which case coaxial 
stacking may occur twice for the same stem, and would already be accounted for in
pseudoEnergy() or in this function on a lower row of the tree.)

MUST be called BEFORE getEnergy() since it sets restrictions on some dangling ends
which should not be included in multiloops.

Cristina: added last parameter and structure to LEcoax_stack_energy_flush_b and mismatch
 *********************************************************************************/
float Loop::getPartialCoaxialEnergy(int flag) {
	float sum = 0;
	//	float energyToAdd = 0;
	//	int * sequence = Input->type;
	//	Loop * L1;
	//	Loop * L2;

	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	// look at neighbouring children across each row of the loop tree
	// use energy minimization for children/branches in a multiloop

	Loop * L = RightChild;
	while (L != NUL){
		sum += L->getPartialCoaxialEnergy(flag);
		L = L->LeftSibling;  // go to the child of this parent

	}

	//    float sum = 0;
	float energyToAdd = 0;
	int * sequence = Input->type;
	Loop * L1;
	Loop * L2;

	int i = 0;
	int combos_start = 0;  // index into allCoaxStacks[] array
	int combos_end = 0;  // index into allCoaxStacks[] array
	int indexMinEnergy = 0;  // index into allCoaxStacks[] array which holds the lowest energy
	int nextToAdd = 0;  // ranges from 0 to NumberOfChildren; index represents one of the possible pairs
	int j = 1;  // keeps track of next spot available to add a new coaxial stacking energy possibility to allCoaxStacks[]
	int k = 2;  // number of pairs being considered

	int dangle0 = -1;  // free base associated with L1
	int dangle1 = -1;  // free base associated with L2
	pk_coax_features allCoaxStacks[MAXCOAXSTACK];

	//printf("Current Loop (Rightmost child): %d to %d\n", begin, end);

	// only perform this calculation for right children (since the left siblings are included in this calculation)
	if (Parent != NUL && this == Parent->RightChild)  
	{
		switch (Parent->type){
			case multi :
				i = 0;
				energyToAdd = 0;
				indexMinEnergy = 0;

				L1 = Parent;  // outer base pair
				L2 = this;  // stacks with rightmost child

				allCoaxStacks[i].index0 = i;
				allCoaxStacks[i].index1 = i+1;
				allCoaxStacks[i].loop0 = L1;
				allCoaxStacks[i].loop1 = L2;
				allCoaxStacks[i].stack_dangle = -1;
				allCoaxStacks[i].lastPairAdded = i;
				allCoaxStacks[i].lastArrayIndex = -1;

				dangle0 = -1;  // free base associated with L1
				dangle1 = -1;  // free base associated with L2

				if (L1->end == L2->end + 1)  // flush coaxial stacking
				{
					if (Input->BasePair(L2->end) - 1 <= 0)
					{
						dangle0 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and right child: dangling end 0 is in wrong place!\n");
					}
					else
						dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1 if base is paired
					if (L1->begin + 1 > Input->Size)
					{
						dangle1 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and right child: dangling end 1 is in wrong place!\n");
					}
					else
						dangle1 = (Input->BasePair(L1->begin + 1) > 0) ? -1 : L1->begin + 1;  // -1 if base is paired
					// NOTE: can't use L2->begin - 1 for dangle0, since L2 may be a pseudoloop and we want the base pair of the stem starting at i

					if (flag == CC2006a)
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), sequence);
					else
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->end)+2), sequence, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair & right child: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->begin - 1] = 1;
					if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->begin + 1] = 1;
					}
					*/
				}
				else if ( (L1->end == L2->end + 2) && ( (L1->begin + 1 <= Input->Size && Input->BasePair(L1->begin + 1) < 0) || (Input->BasePair(L2->end)- 1 > 0 && Input->BasePair(Input->BasePair(L2->end)- 1) < 0) ) )  // mismatch coaxial stacking
				{
					dangle0 = (Input->BasePair(Input->BasePair(L2->end)- 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1: this base can't act in coaxial stack since already paired
					dangle1 = (Input->BasePair(L1->begin + 1) > 0) ? -1 : L1->begin + 1;  // -1: this base can't act in coaxial stack since already paired
					// NOTE: can't use L2->begin - 1 for dangle0, since L2 may be a pseudoloop and we want the base pair of the stem starting at i

					allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->end)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair & right child: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->begin - 1] = 1;
					if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->begin + 1] = 1;
					}
					*/
				}
				else  // coaxial stacking not possible
				{
					if (DEBUG2)
						printf("multi loop: outer pair and right child: neither flush nor mismatch\n");

					allCoaxStacks[i].coax_energy = 0;
				}

				i++;
				L1 = this;
				L2 = LeftSibling;

				while (L2 != NUL)
				{
					allCoaxStacks[i].index0 = i;
					allCoaxStacks[i].index1 = i+1;
					allCoaxStacks[i].loop0 = L1;
					allCoaxStacks[i].loop1 = L2;
					allCoaxStacks[i].stack_dangle = -1;	
					allCoaxStacks[i].lastPairAdded = i;
					allCoaxStacks[i].lastArrayIndex = -1;

					// consider coaxial stacking between children of a multiloop (not the outer base pair of the multiloop)
					if (L1->begin == L2->end + 1)  // flush coaxial stacking
					{
						if (Input->BasePair(L2->end) - 1 <= 0)
							dangle0 = -1;
						else
							dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1 if base is paired
						if (Input->BasePair(L1->begin) + 1 > Input->Size)
							dangle1 = -1;
						else
							dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1 if base is paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						if (flag == CC2006a)
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), sequence);
						else
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), sequence, no_pk_dangling_ends);

						if (DEBUG2)
							printf("multi loop: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (L2->begin - 1 >= 0 && Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally; also do array bounds check
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->begin + 1] = 1;
						}
						*/
					}
					else if ( (L1->begin == L2->end + 2) && ( (Input->BasePair(L1->begin) + 1 <= Input->Size && Input->BasePair(Input->BasePair(L1->begin) + 1) < 0) || (Input->BasePair(L2->end) - 1 > 0 && Input->BasePair(Input->BasePair(L2->end) - 1) < 0) ) )  // mismatch coaxial stacking
					{

						dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1: this base can't act in coaxial stack since already paired
						dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1: this base can't act in coaxial stack since already paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

						if (DEBUG2)
							printf("multi loop: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->end + 1] = 1;
						}
						*/
					}
					else  // coaxial stacking not possible
					{
						if (DEBUG2)
							printf("multi loop: neither flush nor mismatch\n");

						allCoaxStacks[i].coax_energy = 0;
					}

					if (allCoaxStacks[i].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
						indexMinEnergy = i;  // find min of currently calculated stacking possibilities

					i++;

					// add coaxial stacking energy of L1 and L2
					//					energyToAdd += LEcoax_stack_energy_flush_a (L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), sequence); 

					L1 = L2;		
					L2 = L2->LeftSibling;
				}
				//printf("after children are done i = %d (should be 6)", i);

				// consider coaxial stacking between leftmost child (index i) and the outer base pair of multiloop
				// L1 is now leftmost child
				L2 = Parent;  // L2 is outer base pair of multiloop

				allCoaxStacks[i].index0 = i;
				allCoaxStacks[i].index1 = 0;
				allCoaxStacks[i].loop0 = L1;
				allCoaxStacks[i].loop1 = L2;
				allCoaxStacks[i].stack_dangle = -1;
				allCoaxStacks[i].lastPairAdded = i;
				allCoaxStacks[i].lastArrayIndex = -1;

				if (L1->begin == L2->begin + 1)  // flush coaxial stacking
				{
					if (L2->end - 1 <= 0)
					{
						dangle0 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and left child: dangling end 0 is in wrong place!\n");
					}
					else
						dangle0 = (Input->BasePair(L2->end - 1) > 0) ? -1 : L2->end - 1;  // -1 if base is paired
					if (Input->BasePair(L1->begin) + 1 > Input->Size)
					{
						dangle1 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and left child: dangling end 1 is in wrong place!\n");
					}
					else
						dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1 if base is paired
					// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip

					if (flag == CC2006a)
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), sequence);
					else
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->begin)-2), Input->BasePair(Input->BasePair(L1->begin)+2), sequence, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair & left child: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->end - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->end - 1] = 1;
					if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->end + 1] = 1;
					}
					*/
				}
				else if ( (L1->begin == L2->begin + 2) && ( (Input->BasePair(L1->begin) + 1 <= Input->Size && Input->BasePair(Input->BasePair(L1->begin) + 1) < 0) || (L2->end - 1 > 0 && Input->BasePair(L2->end - 1) < 0) ) )  // mismatch coaxial stacking (dangling ends are free bases)
				{
					dangle0 = (Input->BasePair(L2->end - 1) > 0) ? -1 : L2->end - 1;  // -1: this base can't act in coaxial stack since already paired
					dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1: this base can't act in coaxial stack since already paired
					// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip

					allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->begin)-2), Input->BasePair(Input->BasePair(L1->begin)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair and left child: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->end - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->end - 1] = 1;
					if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->end + 1] = 1;
					}
					*/
				}
				else  // coaxial stacking not possible
				{
					if (DEBUG2)
						printf("multi loop: outer pair and left child: neither flush nor mismatch\n");

					allCoaxStacks[i].coax_energy = 0;
				}

				if (allCoaxStacks[i].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
					indexMinEnergy = i;  // find min of currently calculated stacking possibilities

				//printf("multi: allCoaxStacks[%d].coax_energy = %d\n", i, allCoaxStacks[i].coax_energy);

				// up to this point, allCoaxStacks[] from 0 to NumberOfChildren has been filled;
				// next, need to consider combinations of coaxial stacking pairs;
				// note that a stem can only coaxial stack once, e.g. we can't consider both 1.2 and 2.3 stacking

				combos_start = 0;
				combos_end = Parent->NumberOfChildren+1;
				// due to outer base pair, we have NumberOfChildren+1 adjacent pairs (indexed 0 to NumberOfChildren),
				// since the outer base pair can pair with the right child or with the left child

				nextToAdd = 0;  // ranges from 0 to NumberOfChildren-1; index represents one of the possible stems
				j = 1;  // keeps track of next spot available to add a new coaxial stacking energy possibility to allCoaxStacks[]
				k = 2;  // number of pairs being considered

				// each iteration represents i pairs being considered; start from two pairs of coaxially stacked stems
				// can go up to combinations of floor((NumberOfChildren+1)/2) pairs (because a stem can't pair with two other stems)
				for (k = 2; k <= FLOOR_OVER2(Parent->NumberOfChildren+1); k++)
				{
					j = 0;

					// loop through the last set of combinations added and try to add another basic pair
					for (i = combos_start; i < combos_end; i++)
					{
						// add another pair to combination i
						nextToAdd = allCoaxStacks[i].index1 + 1;  // next stem in a pair that could be considered
						if (allCoaxStacks[i].index1 == 0)  // if we are trying to add to a combination that already uses the left child and the outer pair (stem 0)...
							nextToAdd = Parent->NumberOfChildren + 1;  // ...force the while loop below to end

						while (nextToAdd <= Parent->NumberOfChildren)  // no more to stems to add
							// the last stem possible to add is the left most child (index NumberOfChildren-1)
						{
							//printf("allCoaxStacks[%d].coax_energy = %d; adding to %d\n", nextToAdd, allCoaxStacks[nextToAdd].coax_energy, combos_end+j);

							//                             printf("i = %d, nextToAdd = %d, j = %d\n", i, nextToAdd, j);

							//printf("allCoaxStacks[%d].coax_energy = %d\n", nextToAdd, allCoaxStacks[nextToAdd].coax_energy);
							//printf("allCoaxStacks[%d].coax_energy = %d; adding to %d\n", nextToAdd, allCoaxStacks[nextToAdd].coax_energy, combos_end+j);

							if (allCoaxStacks[i].index0 == 0 && nextToAdd > Parent->NumberOfChildren - 1)
								break;  // trying to add a pair with stem #NumberOfChildren-1 and stem #0, but a pair with stem #0 already exists

							allCoaxStacks[combos_end + j].coax_energy = allCoaxStacks[i].coax_energy + allCoaxStacks[nextToAdd].coax_energy;
							allCoaxStacks[combos_end + j].index0 = allCoaxStacks[i].index0;  // the first stem ever added in allCoaxStacks[i] remains the same here
							allCoaxStacks[combos_end + j].index1 = nextToAdd + 1;  // the last stem in the added pair is the one adjacent to nextToAdd
							if (nextToAdd > Parent->NumberOfChildren - 1)
								allCoaxStacks[combos_end + j].index1 = 0;  // i.e. are pairing stem #NumberOfChildren-1 and stem #0
							allCoaxStacks[combos_end + j].stack_dangle = -1;
							allCoaxStacks[combos_end + j].lastArrayIndex = i;
							allCoaxStacks[combos_end + j].lastPairAdded = nextToAdd;

							if (allCoaxStacks[combos_end + j].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
								indexMinEnergy = combos_end + j;  // find min of currently calculated stacking possibilities

							if (DEBUG2)
							{
								printf("multi: total coaxing energy found = %d (%d + %d)\n", allCoaxStacks[combos_end + j].coax_energy, allCoaxStacks[i].coax_energy, allCoaxStacks[nextToAdd].coax_energy);
								printf("multi: current min coaxing energy = %d\n", allCoaxStacks[indexMinEnergy].coax_energy);
							}

							j++;
							nextToAdd++;  // try to add the next basic pair
						}
					}

					combos_start = combos_end;  // the first new combination that was added
					combos_end = combos_end + j;  // the last new combination that was added
				}

				// the minimum free energy coaxial stacking combination is now represented by index indexMinEnergy

				// remove the dangling base energies as appropriate for the coaxial stackings that were chosen
				i = indexMinEnergy;

				do
				{
					if (DEBUG2)
						printf("i = %d, allCoaxStacks[i].lastPairAdded = %d, allCoaxStacks[lastPairAdded].stack_dangle = %d, allCoaxStacks[i].lastArrayIndex = %d\n",
								i, allCoaxStacks[i].lastPairAdded, allCoaxStacks[allCoaxStacks[i].lastPairAdded].stack_dangle, allCoaxStacks[i].lastArrayIndex);

					removeDangling(allCoaxStacks[i].lastPairAdded, Parent->type, allCoaxStacks[allCoaxStacks[i].lastPairAdded], allCoaxStacks[allCoaxStacks[i].lastPairAdded].stack_dangle);
					i = allCoaxStacks[i].lastArrayIndex;
				} while (i != -1);  // while we haven't reached a basic pairing

				// DOES NOT APPLY TO MULTILOOP
				//				if (indexMinEnergy != -1)
				//                               {       
				energyToAdd = allCoaxStacks[indexMinEnergy].coax_energy;
				if (DEBUG2)
				{
					printf("external: final min coaxing energy = %d\n", allCoaxStacks[indexMinEnergy].coax_energy);
				}
				//                                }   
				//                                else
				//                                        energyToAdd = 0;

				sum += energyToAdd;

				break;

			case external :
				i = 0;
				energyToAdd = 0;
				indexMinEnergy = -1;

				//				pk_coax_features allCoaxStacks[MAXCOAXSTACK];

				L1 = this;
				L2 = LeftSibling;

				while (L2 != NUL)
				{
					allCoaxStacks[i].index0 = i;
					allCoaxStacks[i].index1 = i+1;
					allCoaxStacks[i].loop0 = L1;
					allCoaxStacks[i].loop1 = L2;
					allCoaxStacks[i].stack_dangle = -1;	
					allCoaxStacks[i].lastPairAdded = i;
					allCoaxStacks[i].lastArrayIndex = -1;

					// consider coaxial stacking between children of a multiloop (not the outer base pair of the multiloop)
					if (L1->begin == L2->end + 1)  // flush coaxial stacking
					{
						if (Input->BasePair(L2->end) - 1 <= 0)
							dangle0 = -1;
						else
							dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1 if base is paired
						if (Input->BasePair(L1->begin) + 1 > Input->Size)
							dangle1 = -1;
						else
							dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1 if base is paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						if (flag == CC2006a)
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), sequence);
						else
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), sequence, no_pk_dangling_ends);

						//						printf("base pair of Input->BasePair(L1->begin) + 1 = %d\n", Input->BasePair(Input->BasePair(L1->begin) + 1));

						if (DEBUG2)
							printf("external loop: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (L2->begin - 1 >= 0 && Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally; also do array bounds check
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->begin + 1] = 1;
						}
						*/
					}
					else if ( (L1->begin == L2->end + 2) && ( ( Input->BasePair(L1->begin) + 1 <= Input->Size && (Input->BasePair(Input->BasePair(L1->begin) + 1) < 0)) || (Input->BasePair(L2->end) - 1 > 0 && Input->BasePair(Input->BasePair(L2->end) - 1) < 0) ) )  // mismatch coaxial stacking
					{

						dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1: this base can't act in coaxial stack since already paired
						dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1: this base can't act in coaxial stack since already paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

						if (DEBUG2)
							printf("external loop: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->end + 1] = 1;
						}
						*/
					}
					else  // coaxial stacking not possible
					{
						if (DEBUG2)
							printf("external loop: neither flush nor mismatch\n");

						allCoaxStacks[i].coax_energy = 0;
					}

					if (indexMinEnergy == -1 && allCoaxStacks[i].coax_energy <= 0)
						indexMinEnergy = i;
					else if (allCoaxStacks[i].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
						indexMinEnergy = i;  // find min of currently calculated stacking possibilities

					if (DEBUG2)
					{
						printf("external: total coaxing energy found = %d\n", allCoaxStacks[i].coax_energy);
						printf("external: current min coaxing energy = %d\n", allCoaxStacks[indexMinEnergy].coax_energy);
					}

					i++;

					L1 = L2;		
					L2 = L2->LeftSibling;
				}

				// up to this point, allCoaxStacks[] from 0 to NumberOfChildren has been filled;
				// next, need to consider combinations of coaxial stacking pairs;
				// note that a stem can only coaxial stack once, e.g. we can't consider both 1.2 and 2.3 stacking

				if (DEBUG2)
					printf("external: starting combinations of basic pairings with NumberOfChildren = %d\n", Parent->NumberOfChildren);

				//if (Parent->type == external && LeftSibling != NUL)
				//	printf("Parent is external %d to %d, sibling exists %d to %d!\n", Parent->begin, Parent->end, LeftSibling->begin, LeftSibling->end);

				//printf("FLOOR_OVER2(Parent->NumberOfChildren+1) = %d\n", FLOOR_OVER2(Parent->NumberOfChildren+1));

				combos_start = 0;
				combos_end = Parent->NumberOfChildren-1;  // for NumberOfChildren stems, have NumberOfChildren-1 adjacent pairs

				nextToAdd = 0;  // ranges from 0 to NumberOfChildren-1; index represents one of the possible stems
				j = 1;  // keeps track of next spot available to add a new coaxial stacking energy possibility to allCoaxStacks[]
				k = 2;  // number of pairs being considered

				// each iteration represents i pairs being considered; start from two pairs of coaxially stacked stems
				// can go up to combinations of floor(NumberOfChildren/2) = ceiling((NumberOfChildren-1)/2) pairs (because a stem can't pair with two other stems)
				for (k = 2; k <= FLOOR_OVER2(Parent->NumberOfChildren); k++)
				{
					//					printf("here\n");

					j = 0;

					// loop through the last set of combinations added and try to add another basic pair
					for (i = combos_start; i < combos_end; i++)
					{
						// add another pair to combination i
						nextToAdd = allCoaxStacks[i].index1 + 1;  // next stem in a pair that could be considered (index nextToAdd and nextToAdd+1)
						// NOT NECESSARY FOR EXTERNAL LOOP
						//						if (allCoaxStacks[i].index1 == 0)
						//							nextToAdd = Parent->NumberOfChildren + 1;  // i.e. force the while loop below to end

						while (nextToAdd < Parent->NumberOfChildren - 1)  // no more pairs left to add 
							// (-1 is since nextToAdd+1 -- the second stem in the pair being added -- must be less than NumberOfChildren)
						{

							//							printf("i = %d, nextToAdd = %d, j = %d\n", i, nextToAdd, j);

							// NOT NECESSARY FOR EXTERNAL LOOP
							//							if (allCoaxStacks[i].index0 == 0 && nextToAdd >= Parent->NumberOfChildren)
							//								break;  // trying to add a pair with stem #NumberOfChildren and stem #0, but a pair with stem #0 already exists

							allCoaxStacks[combos_end + j].coax_energy = allCoaxStacks[i].coax_energy + allCoaxStacks[nextToAdd].coax_energy;
							allCoaxStacks[combos_end + j].index0 = allCoaxStacks[i].index0;  // the first stem ever added in allCoaxStacks[i] remains the same here
							allCoaxStacks[combos_end + j].index1 = nextToAdd + 1;  // the last stem is the one adjacent to nextToAdd
							// NOT NECESSARY FOR EXTERNAL LOOP
							//							if (nextToAdd + 1 > Parent->NumberOfChildren)
							//								allCoaxStacks[combos_end + j].index1 = 0;  // i.e. are pairing stem #NumberOfChildren and stem #0
							allCoaxStacks[combos_end + j].stack_dangle = -1;
							allCoaxStacks[combos_end + j].lastArrayIndex = i;
							allCoaxStacks[combos_end + j].lastPairAdded = nextToAdd;

							if (indexMinEnergy == -1 && allCoaxStacks[combos_end + j].coax_energy <= 0)
								indexMinEnergy = combos_end + j;
							else if (allCoaxStacks[combos_end + j].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
								indexMinEnergy = combos_end + j;  // find min of currently calculated stacking possibilities

							if (DEBUG2)
							{
								if (indexMinEnergy == -1)
									printf("external: current min coaxing energy = 0\n");
								else
									printf("external: current min coaxing energy = %d\n", allCoaxStacks[indexMinEnergy].coax_energy);
							}

							j++;
							nextToAdd++;  // try to add the next basic pair

						}
					}

					combos_start = combos_end;  // the first new combination that was added
					combos_end = combos_end + j;  // the last new combination that was added

					//					printf("k = %d\n", k);
				}

				// the minimum free energy coaxial stacking combination is now represented by index indexMinEnergy

				// remove the dangling base energies as appropriate for the coaxial stackings that were chosen
				i = indexMinEnergy;

				// don't do this since i could be -1
				//printf("indexMinEnergy = %d, allCoaxStacks[i].lastPairAdded = %d, allCoaxStacks[i].stack_dangle = %d, allCoaxStacks[i].lastArrayIndex = %d\n", 
				//indexMinEnergy, allCoaxStacks[i].lastPairAdded, allCoaxStacks[i].stack_dangle, allCoaxStacks[i].lastArrayIndex);

				while (i != -1)  // while we haven't reached a basic pairing, or if there was only 1 stem
				{
					if (DEBUG2)
						printf("i = %d, allCoaxStacks[i].lastPairAdded = %d, allCoaxStacks[lastPairAdded].stack_dangle = %d, allCoaxStacks[i].lastArrayIndex = %d\n", 
								i, allCoaxStacks[i].lastPairAdded, allCoaxStacks[allCoaxStacks[i].lastPairAdded].stack_dangle, allCoaxStacks[i].lastArrayIndex);

					removeDangling(allCoaxStacks[i].lastPairAdded, Parent->type, allCoaxStacks[allCoaxStacks[i].lastPairAdded], allCoaxStacks[allCoaxStacks[i].lastPairAdded].stack_dangle);
					i = allCoaxStacks[i].lastArrayIndex;
				}

				if (indexMinEnergy != -1)
				{
					energyToAdd = allCoaxStacks[indexMinEnergy].coax_energy;
					if (DEBUG2)
					{
						printf("external: final min coaxing energy = %d\n", allCoaxStacks[indexMinEnergy].coax_energy);
					}
				}
				else
					energyToAdd = 0;

				sum += energyToAdd;

				break;

			default:
				break;
		}
	}

	finalCoaxEnergy = energyToAdd;
	return sum;
}


/*********************************************************************************
getPartialCoaxialEnergyAll: Calculate the coaxial stacking energies between children
in the loop tree (not between structural elements within those children, e.g. stems
of pseudoknots, which are included in individual energy calculations, 
e.g getPseudoEnergy()). This considers coaxial stacking between ALL children (does
not taking minimum between pairs, for example). 
THIS FUNCTION NOT USED IN CURRENT MODELS.

Cristina: added last parameter and structure to LEcoax_stack_energy_flush_b and mismatch
 *********************************************************************************/
float Loop::getPartialCoaxialEnergyAll(int flag) {
	float sum = 0;
	//	float energyToAdd = 0;
	//	int * sequence = Input->type;
	//	Loop * L1;
	//	Loop * L2;

	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';


	// look at neighbouring children across each row of the loop tree
	// use energy minimization for children/branches in a multiloop

	Loop * L = RightChild;
	while (L != NUL){
		sum += L->getPartialCoaxialEnergyAll(flag);
		L = L->LeftSibling;  // go to the child of this parent

	}
	float energyToAdd = 0;
	int * sequence = Input->type;
	Loop * L1;
	Loop * L2;

	int i = 0;
	int indexMinEnergy = 0;  // index into allCoaxStacks[] array which holds the lowest energy

	int dangle0 = -1;  // free base associated with L1
	int dangle1 = -1;  // free base associated with L2
	pk_coax_features allCoaxStacks[MAXCOAXSTACK];

	//printf("Current Loop (Rightmost child): %d to %d\n", begin, end);

	// only perform this calculation for right children (since the left siblings are included in this calculation)
	if (Parent != NUL && this == Parent->RightChild)  
	{
		switch (Parent->type){
			case multi :
				i = 0;
				energyToAdd = 0;
				indexMinEnergy = 0;

				L1 = Parent;  // outer base pair
				L2 = this;  // stacks with rightmost child

				allCoaxStacks[i].index0 = i;
				allCoaxStacks[i].index1 = i+1;
				allCoaxStacks[i].loop0 = L1;
				allCoaxStacks[i].loop1 = L2;
				allCoaxStacks[i].stack_dangle = -1;
				allCoaxStacks[i].lastPairAdded = i;
				allCoaxStacks[i].lastArrayIndex = -1;

				dangle0 = -1;  // free base associated with L1
				dangle1 = -1;  // free base associated with L2

				if (L1->end == L2->end + 1)  // flush coaxial stacking
				{
					if (Input->BasePair(L2->end) - 1 <= 0)
					{
						dangle0 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and right child: dangling end 0 is in wrong place!\n");
					}
					else
						dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1 if base is paired
					if (L1->begin + 1 > Input->Size)
					{
						dangle1 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and right child: dangling end 1 is in wrong place!\n");
					}
					else
						dangle1 = (Input->BasePair(L1->begin + 1) > 0) ? -1 : L1->begin + 1;  // -1 if base is paired
					// NOTE: can't use L2->begin - 1 for dangle0, since L2 may be a pseudoloop and we want the base pair of the stem starting at i

					if (flag == CC2006a)
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), sequence);
					else
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->end)+2), sequence, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair & right child: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->begin - 1] = 1;
					if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->begin + 1] = 1;
					}
					*/
				}
				else if ( (L1->end == L2->end + 2) && ( (L1->begin + 1 <= Input->Size && Input->BasePair(L1->begin + 1) < 0) || (Input->BasePair(L2->end)- 1 > 0 && Input->BasePair(Input->BasePair(L2->end)- 1) < 0) ) )  // mismatch coaxial stacking
				{
					dangle0 = (Input->BasePair(Input->BasePair(L2->end)- 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1: this base can't act in coaxial stack since already paired
					dangle1 = (Input->BasePair(L1->begin + 1) > 0) ? -1 : L1->begin + 1;  // -1: this base can't act in coaxial stack since already paired
					// NOTE: can't use L2->begin - 1 for dangle0, since L2 may be a pseudoloop and we want the base pair of the stem starting at i

					allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->end)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair & right child: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->end, Input->BasePair(L1->end), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->begin - 1] = 1;
					if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->begin + 1] = 1;
					}
					*/
				}
				else  // coaxial stacking not possible
				{
					if (DEBUG2)
						printf("multi loop: outer pair and right child: neither flush nor mismatch\n");

					allCoaxStacks[i].coax_energy = 0;
				}

				i++;
				L1 = this;
				L2 = LeftSibling;

				while (L2 != NUL)
				{
					allCoaxStacks[i].index0 = i;
					allCoaxStacks[i].index1 = i+1;
					allCoaxStacks[i].loop0 = L1;
					allCoaxStacks[i].loop1 = L2;
					allCoaxStacks[i].stack_dangle = -1;	
					allCoaxStacks[i].lastPairAdded = i;
					allCoaxStacks[i].lastArrayIndex = -1;

					// consider coaxial stacking between children of a multiloop (not the outer base pair of the multiloop)
					if (L1->begin == L2->end + 1)  // flush coaxial stacking
					{
						if (Input->BasePair(L2->end) - 1 <= 0)
							dangle0 = -1;
						else
							dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1 if base is paired
						if (Input->BasePair(L1->begin) + 1 > Input->Size)
							dangle1 = -1;
						else
							dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1 if base is paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						if (flag == CC2006a)
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), sequence);
						else
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), sequence, no_pk_dangling_ends);

						if (DEBUG2)
							printf("multi loop: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (L2->begin - 1 >= 0 && Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally; also do array bounds check
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->begin + 1] = 1;
						}
						*/
					}
					else if ( (L1->begin == L2->end + 2) && ( (Input->BasePair(L1->begin) + 1 <= Input->Size && Input->BasePair(Input->BasePair(L1->begin) + 1) < 0) || (Input->BasePair(L2->end) - 1 > 0 && Input->BasePair(Input->BasePair(L2->end) - 1) < 0) ) )  // mismatch coaxial stacking
					{

						dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1: this base can't act in coaxial stack since already paired
						dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1: this base can't act in coaxial stack since already paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

						if (DEBUG2)
							printf("multi loop: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->end + 1] = 1;
						}
						*/
					}
					else  // coaxial stacking not possible
					{
						if (DEBUG2)
							printf("multi loop: neither flush nor mismatch\n");

						allCoaxStacks[i].coax_energy = 0;
					}

					if (allCoaxStacks[i].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
						indexMinEnergy = i;  // find min of currently calculated stacking possibilities

					i++;

					// add coaxial stacking energy of L1 and L2
					//					energyToAdd += LEcoax_stack_energy_flush_a (L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), sequence); 

					L1 = L2;		
					L2 = L2->LeftSibling;
				}
				//printf("after children are done i = %d (should be 6)", i);

				// consider coaxial stacking between leftmost child (index i) and the outer base pair of multiloop
				// L1 is now leftmost child
				L2 = Parent;  // L2 is outer base pair of multiloop

				allCoaxStacks[i].index0 = i;
				allCoaxStacks[i].index1 = 0;
				allCoaxStacks[i].loop0 = L1;
				allCoaxStacks[i].loop1 = L2;
				allCoaxStacks[i].stack_dangle = -1;
				allCoaxStacks[i].lastPairAdded = i;
				allCoaxStacks[i].lastArrayIndex = -1;

				if (L1->begin == L2->begin + 1)  // flush coaxial stacking
				{
					if (L2->end - 1 <= 0)
					{
						dangle0 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and left child: dangling end 0 is in wrong place!\n");
					}
					else
						dangle0 = (Input->BasePair(L2->end - 1) > 0) ? -1 : L2->end - 1;  // -1 if base is paired
					if (Input->BasePair(L1->begin) + 1 > Input->Size)
					{
						dangle1 = -1;
						printf("WARNING: coaxial stacking, multloop, outer pair and left child: dangling end 1 is in wrong place!\n");
					}
					else
						dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1 if base is paired
					// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip

					if (flag == CC2006a)
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), sequence);
					else
						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->begin)-2), Input->BasePair(Input->BasePair(L1->begin)+2), sequence, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair & left child: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->end - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->end - 1] = 1;
					if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->end + 1] = 1;
					}
					*/
				}
				else if ( (L1->begin == L2->begin + 2) && ( (Input->BasePair(L1->begin) + 1 <= Input->Size && Input->BasePair(Input->BasePair(L1->begin) + 1) < 0) || (L2->end - 1 > 0 && Input->BasePair(L2->end - 1) < 0) ) )  // mismatch coaxial stacking (dangling ends are free bases)
				{
					dangle0 = (Input->BasePair(L2->end - 1) > 0) ? -1 : L2->end - 1;  // -1: this base can't act in coaxial stack since already paired
					dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1: this base can't act in coaxial stack since already paired
					// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip

					allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->begin)-2), Input->BasePair(Input->BasePair(L1->begin)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

					if (DEBUG2)
						printf("multi loop: outer pair and left child: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->begin, Input->BasePair(L2->begin), L1->begin, Input->BasePair(L1->begin), COAX_MULTI, dangle0, dangle1, allCoaxStacks[i].coax_energy);
					/*
					   DONE in removeDangling()
					   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
					   {
					// make sure dangling ends are not included elsewhere (since include coaxial stacking)
					if (Input->BasePair(L2->end - 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L2->end - 1] = 1;
					if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[L1->end + 1] = 1;
					}
					*/
				}
				else  // coaxial stacking not possible
				{
					if (DEBUG2)
						printf("multi loop: outer pair and left child: neither flush nor mismatch\n");

					allCoaxStacks[i].coax_energy = 0;
				}

				if (allCoaxStacks[i].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
					indexMinEnergy = i;  // find min of currently calculated stacking possibilities

				//printf("multi: allCoaxStacks[%d].coax_energy = %d\n", i, allCoaxStacks[i].coax_energy);
				i++;

				// up to this point, allCoaxStacks[] from 0 to NumberOfChildren has been filled;
				// next, add all of the stacking energies of all pairs

				for (int j = 0; j < i; j++)
				{
					energyToAdd += allCoaxStacks[j].coax_energy;
					removeDangling(j, Parent->type, allCoaxStacks[j], allCoaxStacks[j].stack_dangle);
				}

				sum += energyToAdd;

				break;

			case external :
				i = 0;
				energyToAdd = 0;
				indexMinEnergy = -1;

				//				pk_coax_features allCoaxStacks[MAXCOAXSTACK];

				L1 = this;
				L2 = LeftSibling;

				while (L2 != NUL)
				{
					allCoaxStacks[i].index0 = i;
					allCoaxStacks[i].index1 = i+1;
					allCoaxStacks[i].loop0 = L1;
					allCoaxStacks[i].loop1 = L2;
					allCoaxStacks[i].stack_dangle = -1;	
					allCoaxStacks[i].lastPairAdded = i;
					allCoaxStacks[i].lastArrayIndex = -1;

					// consider coaxial stacking between children of a multiloop (not the outer base pair of the multiloop)
					if (L1->begin == L2->end + 1)  // flush coaxial stacking
					{
						if (Input->BasePair(L2->end) - 1 <= 0)
							dangle0 = -1;
						else
							dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1 if base is paired
						if (Input->BasePair(L1->begin) + 1 > Input->Size)
							dangle1 = -1;
						else
							dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1 if base is paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						if (flag == CC2006a)
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_a(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), sequence);
						else
							allCoaxStacks[i].coax_energy = LEcoax_stack_energy_flush_b(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), sequence, no_pk_dangling_ends);

						//						printf("base pair of Input->BasePair(L1->begin) + 1 = %d\n", Input->BasePair(Input->BasePair(L1->begin) + 1));

						if (DEBUG2)
							printf("external loop: flush stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (L2->begin - 1 >= 0 && Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally; also do array bounds check
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->begin + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->begin + 1] = 1;
						}
						*/
					}
					else if ( (L1->begin == L2->end + 2) && ( ( Input->BasePair(L1->begin) + 1 <= Input->Size && (Input->BasePair(Input->BasePair(L1->begin) + 1) < 0)) || (Input->BasePair(L2->end) - 1 > 0 && Input->BasePair(Input->BasePair(L2->end) - 1) < 0) ) )  // mismatch coaxial stacking
					{

						dangle0 = (Input->BasePair(Input->BasePair(L2->end) - 1) > 0) ? -1 : Input->BasePair(L2->end) - 1;  // -1: this base can't act in coaxial stack since already paired
						dangle1 = (Input->BasePair(Input->BasePair(L1->begin) + 1) > 0) ? -1 : Input->BasePair(L1->begin) + 1;  // -1: this base can't act in coaxial stack since already paired
						// NOTE: can't use L1->end for dangle1, since it might be a pseudoloop, in which case we want the base pair of the stem starting at ip
						// NOTE: can't use L2->begin for dangle0, since it might be a pseudoloop, in which case we want the base pair of the stem starting at i

						allCoaxStacks[i].coax_energy = LEcoax_stack_energy_mismatch(L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, sequence, structure, Input->BasePair(Input->BasePair(L2->end)-2), Input->BasePair(Input->BasePair(L1->begin)+2), allCoaxStacks[i].stack_dangle, no_pk_dangling_ends);

						if (DEBUG2)
							printf("external loop: mismatch stack(i.j=%d.%d, ip.jp=%d.%d, flag=%d, dangle_i=%d, dangle_ip=%d) = %d\n", L2->end, Input->BasePair(L2->end), L1->begin, Input->BasePair(L1->begin), COAX_OTHER, dangle0, dangle1, allCoaxStacks[i].coax_energy);
						/*
						   DONE in removeDangling()
						   if (allCoaxStacks[i].coax_energy != 0)  // coaxial stacking favourable over dangling ends
						   {
						// make sure dangling ends are not included elsewhere (since include coaxial stacking)
						if (Input->BasePair(L2->begin - 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L2->begin - 1] = 1;
						if (Input->BasePair(L1->end + 1) == 0)  // only change array if the dangling end was free originally
						Input->cannot_add_dangling[L1->end + 1] = 1;
						}
						*/
					}
					else  // coaxial stacking not possible
					{
						if (DEBUG2)
							printf("external loop: neither flush nor mismatch\n");

						allCoaxStacks[i].coax_energy = 0;
					}

					if (indexMinEnergy == -1 && allCoaxStacks[i].coax_energy <= 0)
						indexMinEnergy = i;
					else if (allCoaxStacks[i].coax_energy < allCoaxStacks[indexMinEnergy].coax_energy)
						indexMinEnergy = i;  // find min of currently calculated stacking possibilities

					if (DEBUG2)
					{
						printf("external: total coaxing energy found = %d\n", allCoaxStacks[i].coax_energy);
						printf("external: current min coaxing energy = %d\n", allCoaxStacks[indexMinEnergy].coax_energy);
					}

					i++;

					L1 = L2;		
					L2 = L2->LeftSibling;
				}

				// up to this point, allCoaxStacks[] from 0 to NumberOfChildren has been filled;
				// next, need to consider combinations of coaxial stacking pairs;
				// note that a stem can only coaxial stack once, e.g. we can't consider both 1.2 and 2.3 stacking

				if (DEBUG2)
					printf("external: starting combinations of basic pairings with NumberOfChildren = %d\n", Parent->NumberOfChildren);

				// up to this point, allCoaxStacks[] from 0 to NumberOfChildren has been filled;
				// next, add all of the stacking energies of all pairs

				for (int j = 0; j < i; j++)
				{
					energyToAdd += allCoaxStacks[j].coax_energy;
					removeDangling(j, Parent->type, allCoaxStacks[j], allCoaxStacks[j].stack_dangle);
				}

				sum += energyToAdd;

				break;

			default:
				break;
		}
	}

	finalCoaxEnergy = energyToAdd;
	return sum;
}

/**********************************************************************************
 * MUST be called after Energy(). This is because the no_dangling_array is cleared
 * at the end of this function. Returns in cal/mol. Ignores dangles if
 * no_pk_dangling_ends is 1.
 *********************************************************************************/
float Loop::EnergyDanglingViaSimfold(int model) {

	/*
	   int num_params = get_num_params_PK_DP();


	   double* counter = new double[num_params];
	   double** quadratic_matrix = new double*[num_params];
	   if (quadratic_matrix == NUL)
	   {
	   printf ("ERROR! Space could not be allocated for quadratic_matrix_known, ABORT!\n");
	   exit(1);
	   }
	   for (int i = 0; i < num_params; i++)
	   {
	   counter[i] = 0;
	   quadratic_matrix[i] = new double[num_params];
	   if (quadratic_matrix[i] == NUL)
	   {
	   printf ("ERROR! Space could not be allocated for quadratic_matrix_known[%d], ABORT!\n", i);
	   exit(1);
	   }
	   for (int j = i; j < num_params; j++)
	   quadratic_matrix[i][j] = 0;
	   }
	   */

	double free_value = 0;
	int reset_c = 0;
	int ignore_dangles = no_pk_dangling_ends;
	int ignore_AU = 0;

	float retval = 0;
	/*
	   retval = 1000*EnergyDangling(quadratic_matrix, counter, free_value, reset_c, ignore_dangles, ignore_AU);
	   */
	retval = 1000*EnergyDangling(model, NUL, NUL, free_value, reset_c, ignore_dangles, ignore_AU);

	/*
	   delete [] counter;
	   for (int i = 0; i < num_params; i++)
	   delete [] quadratic_matrix[i];
	   delete [] quadratic_matrix;
	   */
	return retval;
}

/**********************************************************************************
 * MUST be called after Energy(). This is because the no_dangling_array is cleared
 * at the end of this function. Returns in cal/mol.

Cristina: changed occurances of dangling_energy_ to use structure as input parameter
 *********************************************************************************/
float Loop::EnergyDangling() {

	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	float energy = 0;
	PARAMTYPE AUpen = 0;
	PARAMTYPE dang = 0;

	int nb_nucleotides = Input->Size + 1;

	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	for (int i=1; i < nb_nucleotides; i++)
	{
		//		printf("Consider dang = %d, Input-cannot_add_dang = %d\n", i, Input->cannot_add_dangling[i]);

		// add some AU_penalties
		// condition on if is necessary because multiloop AU penalty has already been applied
		if (i==1 && f[i].pair > i) // && Input->ClosedRegions[i] != NUL && Input->ClosedRegions[i]->type != multi)
			//		if (i==1 && f[i].pair > i && Input->ClosedRegions[i] != NUL && Input->ClosedRegions[i]->type != multi)
		{
			AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
			if (DEBUG)
				printf ("%d - AUpen1 \t- add energy %6d 10cal/mol\n", i, AUpen);
			energy += float(AUpen);  
		}            
		else if ( i > 1 && f[i].pair > i && f[i-1].pair < i-1 &&
				f[i-1].pair != 0 && !Input->cannot_add_dangling[i])
			//  )(  
		{            
			AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
			if (DEBUG)
				printf ("%d - AUpen2 \t- add energy %6d 10cal/mol\n", i, AUpen);
			energy += float(AUpen);
		}            

		// add dangling energies and AU_penalties (necessary because the dangling ends for one
		// stem in a pseudoknot is not included in multiloop dangling end calculation)
		if (f[i].pair == 0 && !Input->cannot_add_dangling[i])
		{
			if ((i == 1 || (i > 1 && f[i-1].pair == 0)) &&
					i < nb_nucleotides-1 && f[i+1].pair > i+1)
				// .( or ..(
			{
				if (no_pk_dangling_ends == 0){
					dang = MIN (0, dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]]);
				}
				AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
				/*
				   if (DEBUG)
				   {
				   printf ("%d - dangle1 \t- add energy %6d 10cal/mol\n", i, dang);
				   printf ("%d - AUpen3 \t- add energy %6d 10cal/mol\n", i, AUpen);
				   }                    
				   */
				energy += float(dang + AUpen);
			}
			else if ((i == nb_nucleotides-1 || (i < nb_nucleotides-1 && f[i+1].pair == 0)) && i > 0 && f[i-1].pair > 0 && f[i-1].pair < i-1){
				if (no_pk_dangling_ends == 0){
					dang = MIN (0, dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]]);
				}
				/*
				   if (DEBUG)                
				   printf ("%d - dangle2 \t- add energy %6d 10cal/mol\n", i, dang);
				   if (PRINT_CHAR == 0)  // double
				   energy += float(dang);
				   else if (PRINT_CHAR == 1)  // long double
				   energy += float(dang);
				   else  // int
				   energy += dang;
				   */
			} else if (i < nb_nucleotides-1 && f[i+1].pair > i+1 && f[i-1].pair < i-1 && f[i-1].pair != 0) {
				if (no_pk_dangling_ends == 0){
					dang = MIN (0, dangling_energy (sequence, structure, f[i-1].pair, i-1, i+1, f[i+1].pair));
				}
				AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
				//if (DEBUG){              
				//	printf ("%d - dangle1 \t- add energy %6d 10cal/mol\n", i, dang);
				//	printf ("%d - AUpen4 \t- add energy %6d 10cal/mol\n", i, AUpen);
				//}    
				energy += float(dang + AUpen);
			}else{
				continue;
			}
		}
	};

	// CHECK: why multiply by negative? - I think it should be without negative - changed it to be without negative

	//	printf("The total free energy without dangling ends is %f cal/mol\n", -10*GetEnergy);

	// clear the no dangling restrictions placed, in order to allow multilpe energy models to be used
	Input->clearDanglingRestriction();

	return -10*energy;  // returned result is in cal/mol (not kcal/mol)
}


// returns 1 if pos is immediately adjacent to a pseudoknot;
// returns 0 otherwise
int Loop::isAdjacentToPK(int pos, int a)
{
	// if pos not the last base, and pos+1 starts a closed region, which is a pseudoknot
	if (pos < Input->Size && Input->ClosedRegions[pos+1] != NUL && Input->ClosedRegions[pos+1]->type == pseudo)
	{
		if (DEBUG)
			printf("isAdjacentToPK = 1 for pos=%d, pos+1 starts a PK\n", pos);
		return 1;
	}

	// if pos not the first base, and pos-1 is paired, and does not form an arc across pos, and its pair
	//  does not start a closed region, then it must be a pk; in addition the base pair must fall inside [a,ap]
	if (pos > 1 && Input->BasePair(pos-1) > 0 && Input->BasePair(pos-1) < pos-1  && 
			Input->ClosedRegions[Input->BasePair(pos-1)] == NUL && Input->BasePair(pos-1) > a)
	{
		if (DEBUG)
			printf("isAdjacentToPK = 1 for pos=%d, pos-1 is part of a PK \n", pos);
		return 1;
	}

	return 0;
}


int Loop::findClosedRegionNestedIn(int pos)
	// find the first closed region in which pos is nested
	// PRE: pos >= 1
	// POST: if i = 0, it means pos is not nested in any arc or closed region
{
	int i = pos;
	for (i = pos; i >= 1; i--)
	{
		// if i is the beginning of an arc, and the arc spans across pos, and the arc starts a closed region
		// (because could be the second band of a PK based on only the first 3 conditions)

		//Doesn't for for nested in PK: //if (Input->BasePair(i) > 0 && Input->BasePair(i) > i && Input->BasePair(i) > pos && Input->ClosedRegions[i] != NUL)
		if (Input->ClosedRegions[i] != NUL and Input->ClosedRegions[i]->end >= pos)
			break;  // i represents the beginning of the region in which pos is nested
	}

	return i;
}

int Loop::isNestedInBand(int pos, int& a)
	// return 1 if base pos is nested in the band of the pseudoknot L and, if pos is paired,
	// furthermore corresponds to a closed region completely inside the band (not spanning it);
	// POST: a is the beginning of the band its nested in
{
	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	// check if pos is nested in [a,a'] or [b',b] or [c,c'] or [d',d]
	int k = 0;
	int i = begin;
	int numRegions = NumberOfBands * 2;

	for (k = 0; k < numRegions; k++){  // look at all regions between band borders, as listed above
		if (pos >= i && pos <= bandpattern->pattern[i].OtherBorder){
			// if pos is paired, return 0 if its base pair spans the band
			if (Input->BasePair(pos) > 0){
				if (Input->BasePair(pos) >= i && Input->BasePair(pos) <= bandpattern->pattern[i].OtherBorder){
					//if (DEBUG)
					//	printf("isNestedInBand = 1: found pos=%d nested in band [%d,%d]\n", pos, i, bandpattern->pattern[i].OtherBorder);
					a = i;
					return 1;
				}else{
					return 0;
				}
			}else{
				//if (DEBUG)
				//	printf("isNestedInBand = 1: found pos=%d nested in band [%d,%d]\n", pos, i, bandpattern->pattern[i].OtherBorder);
				a = i;
				return 1;
			}
		}
		i = bandpattern->pattern[i].next;
	}
	return 0;
}


int Loop::isPKDangleInMultiPseudoLoop(int i)
	// POST: return 1 if i is a dangle for the pseudoknotted child of a multiloop that spans a band
{
	T_IntList * L1 = ILoops;
	T_IntList * L2 = MLoops;
	int k_pt = 0;
	int i_pt = begin;

	pk_str_features * feat = Input->loops;

	int last_border = i_pt;
	int current_border = 0;

	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	// we want to fill in the structure between each band in turn, looking for
	// other nested pseudoknots
	// [a,a']|_|[b',b]
	// [c,c']|_|[d',d]
	for (k_pt = 0; k_pt < NumberOfBands*2; k_pt++)
	{
		// per each band
		L1 = ILoops;
		L2 = MLoops;
		last_border = i_pt;
		current_border = 0;

		// multiloops that span bands are handled separately (since has parameters not in simfold) so search for them
		while (L2 != NUL){
			// find the next multiloop that spans a band -- get energy/counts for everything above it
			// that is also inside the last multiloop  i.e. look through ILoops again for the first
			// loop that is above this multiloop but below the last one

			if (DEBUG)
				printf("Consider multiloop starting at %d in band starting at %d,%d\n", Input->looplists[L2->Num]->base1, i_pt, bandpattern->pattern[i_pt].OtherBorder);

			if (Input->looplists[L2->Num]->base1 >= i_pt 
					&& Input->looplists[L2->Num]->base1 <= bandpattern->pattern[i_pt].OtherBorder)
			{
				if (DEBUG)
					printf("This multiloop works %d\n", Input->looplists[L2->Num]->base1);

				current_border = Input->looplists[L2->Num]->base1;  // the current multiloop we are looking at

				if (current_border == last_border && current_border != i_pt)
				{
					printf("WARNING: should have moved onto a different multiloop that spans the band.\n");
				}

				L1 = ILoops;
				while (L1 != NUL && Input->looplists[L1->Num]->base1 <= current_border)  // ie. L1 outside multi
					L1 = L1->Next;
				// now L1 points to the first L1 inside the multi,
				// or is NUL (single base pair) (or is not null since there might be other ILoops members but not in this band)
				if (L1 == NUL || bandpattern->pattern[Input->looplists[L1->Num]->base1].OtherBorder >= bandpattern->pattern[Input->looplists[L2->Num]->base1].OtherBorder)
				{
					int single_base = Input->looplists[L2->Num]->base1+1;
					while (Input->BasePair(single_base) <= bandpattern->pattern[i_pt].next) // doesn't span band
						single_base++;
					// now single_base is the one that spans the band
					if (i == single_base - 1)
					{
						if (DEBUG)
							printf("isPKDangleInMultiPseudoLoop = 1 for %d with single_base = %d\n", i, single_base);
						return 1;
					}
					if (i == feat[single_base].pair + 1)
					{
						if (DEBUG)
							printf("isPKDangleInMultiPseudoLoop = 1 for %d with single_base = %d (pair)\n", i, single_base);
						return 1;
					}
				}
				else
				{
					if (i == Input->looplists[L1->Num]->base1 - 1)
					{
						if (DEBUG)
							printf("isPKDangleInMultiPseudoLoop = 1 for %d with iloops = %d\n", i, Input->looplists[L1->Num]->base1);
						return 1;
					}
					if (i == feat[Input->looplists[L1->Num]->base1].pair + 1)
					{
						if (DEBUG)
							printf("isPKDangleInMultiPseudoLoop = 1 for %d with iloops = %d\n", i, Input->looplists[L1->Num]->base1);
						return 1;
					}
				}

				last_border = current_border;  // this multiloop is done
			}
			// go to the next multiloop
			L2 = L2->Next;
		};
		// at this point, we've handled all multiloops in the band

		i_pt = bandpattern->pattern[i_pt].next;  // go to next band
	}

	return 0;
}


// FOR PARAMETER TUNING
/**********************************************************************************
 * MUST be called after Energy(). This is because the no_dangling_array is cleared
 * at the end of this function. Returns in kcal/mol.
 *********************************************************************************/
float Loop::EnergyDangling(int model, double** P_matrix, double *c, double &f_pt, int reset_c, int ignore_dangles, int ignore_AU) {

	int num_params = 0;
	if (model == DP)
		num_params = get_num_params_PK_DP();
	else if (model == CC2006b)
		num_params = get_num_params_PK_CC2006b();

	if (num_params == 0)
		printf("WARNING: Loop.cpp::EnergyDangling(int model, double P_matrix, ...) - num_params not initialized\n");

	if (reset_c == 1 && c != NUL)
	{
		f_pt = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			if (P_matrix != NUL)
				for (int j=i; j < num_params; j++) 
					P_matrix[i][j] = 0;
		}
	}

	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	PARAMTYPE energy = 0;
	PARAMTYPE AUpen = 0;
	PARAMTYPE dang = 0;

	int nb_nucleotides = Input->Size + 1;

	char paramtype[100];
	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	int index_dang = -1;  // for updating counter for dangling ends

	int band_nested = 0;  // used for some of the helper functions

	for (int i=1; i < nb_nucleotides; i++)
	{
		//		printf("Consider dang = %d, Input-cannot_add_dang = %d\n", i, Input->cannot_add_dangling[i]);

		dang = 0;
		AUpen = 0;
		// add some AU_penalties
		// condition on if is necessary because multiloop AU penalty has already been applied in simfold call
		//		if (i==1 && f[i].pair > i && Input->ClosedRegions[i] != NUL && Input->ClosedRegions[i]->type != multi)
		if (i==1 && f[i].pair > i && Input->ClosedRegions[i] != NUL && Input->ClosedRegions[i]->type != multi &&
				Input->ClosedRegions[i]->type != hairpin && Input->ClosedRegions[i]->type != stackloop &&
				Input->ClosedRegions[i]->type != interior && Input->ClosedRegions[i]->type != external)
			// Cristina: added last 4 conditions, which apply when we are using the EnergyViaSimfold and the energy with counts
		{
			/*  // not necessary: done in simfold
				if (ignore_AU == 0)
				{
				AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
				if (c != NUL)    count_AU_penalty (sequence[i], sequence[f[i].pair], c);
				}

				if (DEBUG)
				printf ("%d - AUpen1 \t- add energy %6d 10cal/mol\n", i, AUpen);
				if (PRINT_CHAR == 0)  // double
				energy += float(AUpen);
				else if (PRINT_CHAR == 1)  // long double
				energy += float(AUpen);
				else  // int
				energy += AUpen;
				*/
		}
		else if ( i > 1 && f[i].pair > i && f[i-1].pair < i-1 &&
				f[i-1].pair != 0 && !Input->cannot_add_dangling[i]) 
		{            
			/*  // not necessary: done in simfold
				int arcIndex = findClosedRegionNestedIn(i);

				if (DEBUG2)
				printf("Aupen2: arcIndex = %d\n", arcIndex);

			// don't add penalty if nested in something AND { if nested in pkfree 
			// OR nested in a multiloop with nested PKs and its not a PK
			// OR nested in the band of a pk and not a PK }
			if ( !( arcIndex > 0 && 
			( Input->ClosedRegions[arcIndex]->isPKFree() == 1                            ||
			(Input->ClosedRegions[arcIndex]->type != pseudo && Input->ClosedRegions[i]->type != pseudo) ||
			(Input->ClosedRegions[arcIndex]->type == pseudo && Input->ClosedRegions[arcIndex]->isNestedInBand(i) == 1) ) ) )
			//  )(  
			{
			if (ignore_AU == 0)
			{
			AUpen = AU_penalty (sequence[i], sequence[f[i].pair]);
			if (c != NUL)    count_AU_penalty (sequence[i], sequence[f[i].pair], c);
			}

			if (DEBUG)
			printf ("%d - AUpen2 \t- add energy %6d 10cal/mol\n", i, AUpen);
			if (PRINT_CHAR == 0)  // double
			energy += float(AUpen);
			else if (PRINT_CHAR == 1)  // long double
			energy += float(AUpen);
			else  // int
			energy += AUpen;
			}
			*/
		}	
		// add dangling energies and AU_penalties (necessary because the dangling ends for one
		// stem in a pseudoknot is not included in multiloop dangling end calculation)
		if (f[i].pair == 0 && !Input->cannot_add_dangling[i])
		{
			if ((i == 1 || (i > 1 && f[i-1].pair == 0)) &&
					i < nb_nucleotides-1 && f[i+1].pair > i+1)
				// .( or ..(
			{
				int arcIndex_i = findClosedRegionNestedIn(i);
				/* // not necessary: done in simfold
				   int arcIndex_i1 = findClosedRegionNestedIn(i+1);
				   */
				if (DEBUG2)
				{
					printf("dangle1: arcIndex_i = %d for dangle %d, Input-cannot_add_dang=%d\n", arcIndex_i, i, Input->cannot_add_dangling[i]);
					//				printf("Aupen3: arcIndex_i1 = %d\n", arcIndex_i1);
				}

				// add dangle if the array is set to 2
				// don't add dangle if nested in something AND { if nested in pkfree 
				// OR nested in a multiloop with nested PKs and the dangling end is not adjacent to a PK
				// OR nested in the band of a pk and the dangling end is not adjacent to a PK }
				if (ignore_dangles == 0 && 
						((arcIndex_i > 0 && Input->ClosedRegions[arcIndex_i]->type == pseudo && 
						  Input->ClosedRegions[arcIndex_i]->isPKDangleInMultiPseudoLoop(i) == 1) || 
						 !( arcIndex_i > 0 && ( Input->ClosedRegions[arcIndex_i]->isPKFree() == 1                            ||
								 (Input->ClosedRegions[arcIndex_i]->type != pseudo && isAdjacentToPK(i, 0) == 0) ||
								 (Input->ClosedRegions[arcIndex_i]->type == pseudo && Input->ClosedRegions[arcIndex_i]->isNestedInBand(i, band_nested) == 1 && isAdjacentToPK(i, band_nested) == 0) ))))
				{
					dang = MIN (0, dangle_bot [sequence[f[i+1].pair]] [sequence[i+1]] [sequence[i]]);
					if (dang != 0 && c != NUL)
					{
						sprintf (paramtype, "dangle_bot[%d][%d][%d]",sequence[f[i+1].pair], sequence[i+1], sequence[i]);
						index_dang = structure_type_index (paramtype);
						c[index_dang]++;
					}
				}
				/* // not necessary: done in simfold

				// don't add penalty if nested in something AND { if nested in pkfree 
				// OR nested in a multiloop with nested PKs and its not a PK
				// OR nested in the band of a pk and not a PK }
				if (ignore_AU == 0 && !( arcIndex_i1 > 0 && 
				( Input->ClosedRegions[arcIndex_i1]->isPKFree() == 1 ||
				(Input->ClosedRegions[arcIndex_i1]->type != pseudo && Input->ClosedRegions[i+1]->type != pseudo) ||
				(Input->ClosedRegions[arcIndex_i1]->type == pseudo && Input->ClosedRegions[arcIndex_i1]->isNestedInBand(i+1) == 1) )))
				{
				AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
				if (c != NUL)    count_AU_penalty (sequence[i+1], sequence[f[i+1].pair], c);
				}
				*/
				if (DEBUG)
				{
					printf ("%d - dangle1 \t- add energy %6d 10cal/mol\n", i, dang);
					printf ("%d - AUpen3 \t- add energy %6d 10cal/mol\n", i, AUpen);
				}                    
				if (PRINT_CHAR == 0)  // double
					energy += float(dang + AUpen);
				else if (PRINT_CHAR == 1)  // long double
					energy += float(dang + AUpen);
				else  // int
					energy += dang + AUpen;

			}
			else if ((i == nb_nucleotides-1 || 
						(i < nb_nucleotides-1 && f[i+1].pair == 0)) &&
					i > 0 && f[i-1].pair > 0 && f[i-1].pair < i-1)
			{
				int arcIndex = findClosedRegionNestedIn(i);

				if (DEBUG2)
					printf("dangle2: arcIndex = %d for dangle %d, Input-cannot_add_dang=%d\n", arcIndex, i, Input->cannot_add_dangling[i]);

				// don't add dangle if nested in something AND { if nested in pkfree 
				// OR nested in a multiloop with nested PKs and the dangling end is not adjacent to a PK
				// OR nested in the band of a pk and the dangling end is not adjacent to a PK }
				if ( (arcIndex > 0 && Input->ClosedRegions[arcIndex]->type == pseudo && 
							Input->ClosedRegions[arcIndex]->isPKDangleInMultiPseudoLoop(i) == 1) || 
						!( arcIndex > 0 && ( Input->ClosedRegions[arcIndex]->isPKFree() == 1                            ||
								(Input->ClosedRegions[arcIndex]->type != pseudo && isAdjacentToPK(i, 0) == 0) ||
								(Input->ClosedRegions[arcIndex]->type == pseudo && Input->ClosedRegions[arcIndex]->isNestedInBand(i, band_nested) == 1 && isAdjacentToPK(i, band_nested) == 0)
								) ) )
					// ). or )..
				{
					if (ignore_dangles == 0)
					{
						dang = MIN (0, dangle_top [sequence[i-1]] [sequence[f[i-1].pair]] [sequence[i]]);
						if (dang != 0 && c != NUL)
						{
							sprintf (paramtype, "dangle_top[%d][%d][%d]", sequence[i-1], sequence[f[i-1].pair], sequence[i]);
							index_dang = structure_type_index (paramtype);
							c[index_dang]++;
						}
					}
					if (DEBUG)                
						printf ("%d - dangle2 \t- add energy %6d 10cal/mol\n", i, dang);
					if (PRINT_CHAR == 0)  // double
						energy += float(dang);
					else if (PRINT_CHAR == 1)  // long double
						energy += float(dang);
					else  // int
						energy += dang;
				}
			}
			else if (i < nb_nucleotides-1 && f[i+1].pair > i+1 && f[i-1].pair < i-1 && f[i-1].pair != 0)
			{
				int arcIndex_i = findClosedRegionNestedIn(i);

				/* // not necessary: done in simfold
				   int arcIndex_i1 = findClosedRegionNestedIn(i+1);
				   */

				if (DEBUG2)
				{
					printf("dangle1: arcIndex_i = %d for dangle %d Input-cannot_add_dang=%d\n", arcIndex_i, i, Input->cannot_add_dangling[i]);
					//				printf("Aupen4: arcIndex_i1 = %d\n", arcIndex_i1);
				}

				// ).( 
				// don't add dangle if nested in something AND { if nested in pkfree 
				// OR nested in a multiloop with nested PKs and the dangling end is not adjacent to a PK
				// OR nested in the band of a pk and the dangling end is not adjacent to a PK }
				if (ignore_dangles == 0 && 
						((arcIndex_i > 0 && Input->ClosedRegions[arcIndex_i]->type == pseudo && 
						  Input->ClosedRegions[arcIndex_i]->isPKDangleInMultiPseudoLoop(i) == 1) || 
						 !( arcIndex_i > 0 && ( Input->ClosedRegions[arcIndex_i]->isPKFree() == 1                            ||
								 (Input->ClosedRegions[arcIndex_i]->type != pseudo && isAdjacentToPK(i, 0) == 0) ||
								 (Input->ClosedRegions[arcIndex_i]->type == pseudo && Input->ClosedRegions[arcIndex_i]->isNestedInBand(i, band_nested) == 1 && isAdjacentToPK(i, band_nested) == 0) ))))
				{
					dang = MIN (0, dangling_energy (sequence, structure, f[i-1].pair, i-1, i+1, f[i+1].pair));
					if (c != NUL)    count_LEdangling_energy (sequence, structure, -1, f[i-1].pair, i-1, i+1, f[i+1].pair, c);
				}

				/*  // not necessary: done in simfold
				// don't add penalty if nested in something AND { if nested in pkfree 
				// OR nested in a multiloop with nested PKs and its not a PK
				// OR nested in the band of a pk and not a PK }
				if (ignore_AU == 0 && !( arcIndex_i1 > 0 && 
				( Input->ClosedRegions[arcIndex_i1]->isPKFree() == 1 ||
				(Input->ClosedRegions[arcIndex_i1]->type != pseudo && Input->ClosedRegions[i+1]->type != pseudo) ||
				(Input->ClosedRegions[arcIndex_i1]->type == pseudo && Input->ClosedRegions[arcIndex_i1]->isNestedInBand(i+1) == 1) )))
				{
				AUpen = AU_penalty (sequence[i+1], sequence[f[i+1].pair]);
				if (c != NUL)    count_AU_penalty (sequence[i+1], sequence[f[i+1].pair], c);
				}
				*/
				if (DEBUG)
				{              
					printf ("%d - dangle1 \t- add energy %6d 10cal/mol\n", i, dang);
					printf ("%d - AUpen4 \t- add energy %6d 10cal/mol\n", i, AUpen);
				}    
				if (PRINT_CHAR == 0)  // double
					energy += float(dang + AUpen);
				else if (PRINT_CHAR == 1)  // long double
					energy += float(dang + AUpen);
				else  // int
					energy += dang + AUpen;
			}
			else
			{
				continue;
			}
		}
	};

	// CHECK: why multiply by negative? - I think it should be without negative - changed it to be without negative

	if (DEBUG2)
		printf("The dangling energy is %f cal/mol\n", -10*float(energy));

	// clear the no dangling restrictions placed, in order to allow multilpe energy models to be used
	Input->clearDanglingRestriction();

	int num_params_pkfree = get_num_params();
	if (DEBUG)
	{
		printf("Free Value: %f\n", f_pt);
		printf("PK Counter Values:\n");
		for (int i = num_params_pkfree; i < num_params; i++)
			if (c != NUL && c[i] != 0.0)
				printf("c[%d]=%f  ", i, c[i]);

		printf("\nNon-Zero Simfold Counter Values:\n");
		for (int i = 0; i < num_params_pkfree; i++)
			if (c != NUL && c[i] != 0.0)
				printf("c[%d]=%f  ", i, c[i]);

		printf("\nAll Non-Zero P_matrix Values:\n");
		for (int i = 0; i < num_params; i++)
		{
			for (int j = i; j < num_params; j++)
				if (P_matrix != NUL && P_matrix[i][j] != 0.0)
					printf("P[%d][%d]=%f  ", i, j, P_matrix[i][j]);
		}
	}

	return -float(energy)/100.0;  // returned result is in kcal/mol (not cal/mol)
}


///////////////////////////////////////////////////////////////////
// DIRKS & PIERCE ENERGY MODEL
///////////////////////////////////////////////////////////////////

// returns 1 is the loop has no nested pseudoknots (pk-free)
// and 0 otherwise
int Loop::isPKFree()
{
	// base case - found a pseudoknot
	if (type == pseudo){
		//if (DEBUG)
		//	printf("found a pk at (%d, %d)\n", begin, end);
		return 0;
	}

	Loop * L = RightChild;
	while (L != NUL){
		if (L->isPKFree() == 1)  // if still pk-free check siblings
			L = L->LeftSibling;
		else
			return 0;  // if find a pseudoknot, stop checking siblings
	};
	// if we didn't find a pseudoknot (in the else case above) then must be pk-free
	return 1;
}


/*********************************************************************************
getEnergy: calculating the summation of the free energy of the loops in
the secondary structure. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float Loop::getEnergyDP(){
	float sum = 0;
	float energyToAdd = 0;

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)
	Loop * L = RightChild;
	while (L != NUL){
		sum += L->getEnergyDP();
		L = L->LeftSibling;
	};

	switch (type){
		case	stackloop:
			energyToAdd = float(stackEnergy());
			//			energyTrace += "***[" + begin + ", " + end + "] Stacked Pair Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Stacked Pair Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	hairpin	:
			energyToAdd = float(hairpinEnergy());
			//			energyTrace += "***[" + begin + ", " + end + "] Hairpin Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Hairpin Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	interior:
			energyToAdd = float(interiorEnergy());
			//			energyTrace += "***[" + begin + ", " + end + "] Internal Loop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Internal Loop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	multi	:
			energyToAdd = float(multiEnergyDP());
			//			energyTrace += "***[" + begin + ", " + end + "] Multiloop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Multiloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	pseudo	:
			energyToAdd = float(pseudoEnergyDP());
			//			energyTrace += "***[" + begin + ", " + end + "] Pseudoloop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Pseudoloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		default:
			break;
	};

	finalEnergy = energyToAdd;
	return sum;
}

// FOR PARAMETER TUNING - helper function for getEnergyDP
void Loop::lookForPk(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles, float& sum){
	int num_params = get_num_params_PK_DP();
	if (reset_c == 1 && c != NUL){
		f = 0;         
		for (int i=0; i < num_params; i++){
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}
	//	printf("lookForPK: at top\n");

	if (type == pseudo){
		sum += this->getEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
		return;
	}

	Loop * L = RightChild;
	while (L != NUL){
		L->lookForPk(P_matrix, c, f, reset_c, ignore_dangles, sum);
		L = L->LeftSibling;
	}
}

// FOR PARAMETER TUNING
/*********************************************************************************
getEnergyDP: FOR PARAMETER TUNING. Same as before, but for parameter tuning. 
Uses simfold to get energy of biggest pseudoknot-free structure (does not parse 
down to basic motifs like hairpins, stacked pairs, etc). Returns energy in kcal/mol.
 **********************************************************************************/
float Loop::getEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles){

	float sum = 0;
	float energyToAdd = 0;
	float tempsum = 0;

	int num_params = get_num_params_PK_DP();
	if (reset_c == 1 && c != NUL){
		f = 0;         
		for (int i=0; i < num_params; i++){
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	Loop * L;

	if (begin == 0){
		if (RightChild == NUL){  // no loops in the structure
			return 0;
		}else{
			L = RightChild;
			while (L != NUL) {
				sum += L->getEnergyDP(P_matrix,c,f,reset_c,ignore_dangles);
				L = L->LeftSibling;
			}
			return sum;
		}
	}

	//	if (begin == 0)  // if we are at the root of the tree
	//		L = RightChild;
	//	else  // if we are not
	L = this;

	//if (DEBUG)
	//	printf("At top of getEnergyDP for loop (%d, %d)\n", L->begin, L->end);

	//	float sum = 0;
	//	float energyToAdd = 0;

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)

	// FOR PARAMETER TUNING:
	// if we find a pseudoknot-free structure (ie. children in the tree rooted at this loop 
	// have no pseudoknots), then simply call simfold's energy function without parsing the
	// tree into the lowest-level structures
	if (L->isPKFree() == 1){
		//if (DEBUG)
		//	printf("getEnergyDP: IDed as pk-free\n");
		finalEnergy = L->pkfreeEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
		sum += finalEnergy;
		return sum;
	}

	//if (DEBUG)
	//	printf("getEnergyDP: not pkfree, proceeding\n");

	// if we find a non pk-free structure but it is not a pseudoknot then there must be a nested pseudoknot
	if (L->type != pseudo){
		//if (DEBUG)
		//	printf("getEnergyDP: IDed as nested pk\n");

		// repeat for the highest level of pseudoknots -- call helper function
		tempsum = 0;
		L->lookForPk(P_matrix, c, f, reset_c, ignore_dangles, tempsum);
		sum += tempsum;

		//if (DEBUG)
		//	printf("tempsum: %f", tempsum);

		//Loop * L1 = L->RightChild;
		//while (L1 != NUL){

		//	printf("getEnergyDP: pk inside pk-free region (%d,%d)\n", L1->begin, L1->end);
		//	sum += L1->getEnergyDP(c, f, reset_c, ignore_dangles);
		//	printf("getEnergyDP: done pk inside pk-free region (%d,%d)\n", L1->begin, L1->end);
		//	L1 = L1->LeftSibling;
		//};

		//		printf("getEnergyDP: done the highest level of pks\n");
		energyToAdd = L->nestedPseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
		sum += energyToAdd;
		finalEnergy = energyToAdd;
		return sum;
	}

	if (L->type == pseudo){
		// repeat for anything nested inside
		Loop * L1 = L->RightChild;
		while (L1 != NUL){
			sum += L1->getEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
			L1 = L1->LeftSibling;
		};
		energyToAdd = L->pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
		sum += energyToAdd;
		finalEnergy = energyToAdd;
		return sum;
	}else{
		printf("WARNING: getEnergyDP(): This case should never happen.\n");
	}
	return 0;
}

/*********************************************************************************
pkfreeEnergyDP: FOR PARAMETER TUNING Calculates the free energy of a closed region that is pseudoknot 
free. Calls the SimFold function directly (without further parsing the structure 
into smaller pieces). Also gets the feature counts for parameter tuning.
Energy value is returned in kcal/mol.
 *********************************************************************************/
double Loop::pkfreeEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles){
	// P_matrix remains unchanged
	int num_params = get_num_params_PK_DP();
	if (reset_c == 1 && c != NUL){
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	int numbases = end - begin + 1;

	// create structure: a string of dot-brackets to represent the region
	//	char * structure = new char[numbases+1];
	char structure[numbases+1];
	// create sequence
	//	char * csequence = new char[numbases+1];
	char csequence[numbases+1];

	//	printf("Input:CSequence = ");

	for (int i = 0; i < numbases; i++)
	{
		//		printf("%c", Input->CSequence[begin+i]);

		csequence[i] = Input->CSequence[begin + i];
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	csequence[numbases] = '\0';
	structure[numbases] = '\0';

	// DEGUG PARAMETER TUNING
	/*
	if (DEBUG){
		printf("Parameter Tuning Input (a pk-free region):\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", csequence[i]);
		}
		printf("\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
		}
		printf("\n");
	}
	*/
	// call SimFold energy/feature counts function
	double retval = get_feature_counts_restricted (csequence, structure, c, f, reset_c, ignore_dangles, 0);
	/*
	if (DEBUG){
		printf("--> Energy from simfold get_feature_counts: %f\n", retval);
		printf("Free Value: %f\n", f);
		printf("PK Counter Values:\n");
		for (int i = get_num_params(); i < num_params; i++)
			if (c != NUL && c[i] != 0.0)
				printf("c[%d]=%f  ", i, c[i]);

		printf("\nNon-Zero Simfold Counter Values:\n");
		for (int i = 0; i < get_num_params(); i++)
			if (c != NUL && c[i] != 0.0)
				printf("c[%d]=%f  ", i, c[i]);

		printf("\nAll Non-Zero P_matrix Values:\n");
		for (int i = 0; i < num_params; i++)
		{
			for (int j = i; j < num_params; j++)
				if (P_matrix != NUL && P_matrix[i][j] != 0.0)
					printf("P[%d][%d]=%f  ", i, j, P_matrix[i][j]);
		}
	}
	*/
	return retval;
}


/*********************************************************************************
multiEnergyDP: it calculates the free energy of a multiLoop.
The code is similar to the one in simFold but it has been changed slightly
to consider the pseudoknotted tuples. Energy value is returned in 10*cal/mol.

Cristina: changed occurances of dangling_energy_ to use structure input parameter.
 *********************************************************************************/
double	Loop::multiEnergyDP(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = begin;
	PARAMTYPE AUpen;
	int h, l;
	PARAMTYPE dang;
	double misc_energy;

	dang = 0;
	misc_energy = 0;
	AUpen = 0;
	int special;
	special = 0;

	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	// add the energies/enthalpies for free bases
	for (l=i+1; l < f[i].bri[0]; l++)  // count the number of unpaired bases from i+1 to the first branch
	{
		misc_energy += pkmodelDP.c;
		//		misc_energy += misc.multi_free_base_penalty;

		g_count_c++;
	}

	// for all the branches in this multiloop...
	for (h=0; h < f[i].num_branches-1; h++)
	{
		// ...count the number of unpaired bases from the end of one branch to the beginning of the next branch
		for (l = f[f[i].bri[h]].pair + 1; l < f[i].bri[h+1]; l++)
		{
			misc_energy += pkmodelDP.c;
			//			misc_energy += misc.multi_free_base_penalty;

			g_count_c++;
		}
	}

	// count the number of unpaired bases from the end of the last branch to the end of the multiloop
	for (l = f[f[i].bri[f[i].num_branches-1]].pair + 1; l < f[i].pair; l++)
	{
		misc_energy += pkmodelDP.c;
		//		misc_energy += misc.multi_free_base_penalty;

		g_count_c++;
	}

	if (DEBUG)
		printf("[misc_energyDP]1 (unpaired bases) is %f kcal/mol\n", misc_energy);

	// add penalty for initiating a multiloop
	misc_energy += pkmodelDP.a;
	//	misc_energy += misc.multi_offset;

	g_count_a++;

	if (DEBUG)
		printf("[misc_energyDP]2 (+ initiation penalty) is %f kcal/mol\n", misc_energy);
	//		printf("[%%%%%%%%%%%%%%%%%%%%%%%%5]misc.u.tidfjfwrw = %d\n", misc.multi_offset);

	// DONECHECK: the +1 in number of branches is for the closing base pair
	// add the penalty for each base pair in the multiloop (the +1 in number of branches is for the closing base pair i,f[i].pair)
	misc_energy += pkmodelDP.b * (f[i].num_branches + 1);
	//	misc_energy += pkmodelDP.b * (f[i].pseudo_num_branches + 1);  // this is if pseudoloops count as 2 branches instead of 1
	//	misc_energy += misc.multi_helix_penalty * (f[i].pseudo_num_branches + 1);

	g_count_b += f[i].num_branches + 1;

	if (DEBUG)
		printf("[misc_energyDP]3 (+ number of branches = %d) is %f kcal/mol\n", f[i].num_branches + 1, misc_energy);

	misc_energy *= 100;  // multiply by 100 since the remaining energies in this function are in 10cal/mol


	// add AU_penalties for multi-loop
	AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);
	//	printf("MULTI: AUpen = %d\n", AUpen);	// Cristina: July 17, 2007

	for (h=0; h < f[i].num_branches; h++) {
		AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);
		//		printf("MULTI: AUpen %d = %d\n", h, AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]));	// Cristina: July 17, 2007
	}

	if (no_pk_dangling_ends == 0)
	{	
		// add dangling energies for multi-loop
		dang += dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
		if (DEBUG)
			printf("MULTI: dang left = %d\n", dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair));	// Cristina: July 17, 2007

		for (l=0; l < f[i].num_branches - 1; l++){
			dang += dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
			if (DEBUG)
				printf("MULTI: dang %d = %d: %d %d %d %d, %d %d %d %d\n", l, dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair),
						f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair,
						sequence[f[i].bri[l]], sequence[f[f[i].bri[l]].pair], sequence[f[i].bri[l+1]], sequence[f[f[i].bri[l+1]].pair]);	// Cristina: July 17, 2007
		}
		dang += dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);
		if (DEBUG)
			printf("MULTI: dang right = %d\n", dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair));	// Cristina: July 17, 2007
	}

	// add "no-dangling" restriction                                    
	for (l=0; l < f[i].num_branches; l++)
	{
		Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
		Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
	}
	Input->cannot_add_dangling [i] = 1;  // restrict start of multiloop from having an AU penalty added

	if (DEBUG)
	{
		printf ("%d - multi m DP\t- add energy %f 10cal/mol\n", i, misc_energy);
		printf ("%d - multi d DP\t- add energy %6d 10cal/mol\n", i, dang);
		printf ("%d - multi AU DP\t- add energy %6d 10cal/mol\n", i, AUpen);
	}

	if (DEBUG)            
		printf ("%d multi DP \t- add energy %f 10cal/mol\n", i, misc_energy + dang + AUpen);

	return misc_energy + float(dang + AUpen);


}

// FOR PARAMETER TUNING
void Loop::fillMultiStructure(char * structure, char * csequence, int numbases, int l_begin)
{
	//	printf("fillMultiStructure: at top of fcn\n");
	if (DEBUG)
	{
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
		}
		printf("\n");
	}

	Loop * L = RightChild;
	while (L != NUL){
		if (L->type == pseudo)
		{
			if (DEBUG)
				printf("fillMultiStructure: consider (%d,%d) as pk\n", L->begin, L->end);
			structure[L->begin - l_begin] = '<';
			structure[L->end - l_begin] = '>';
			for (int i = L->begin - l_begin + 1; i < L->end - l_begin; i++)
				structure[i] = 'x';

			//printf("fillMultiStructure: after pk\n");
			if (DEBUG)
			{
				for (int i = 0; i < numbases; i++)
				{
					printf("%c", structure[i]);
				}
				printf("\n");
			}
		}
		else
		{
			if (DEBUG)
				printf("fillMultiStructure: consider (%d,%d) as non-pk\n", L->begin, L->end);
			L->fillMultiStructure(structure, csequence, numbases, l_begin);
		}
		L = L->LeftSibling;
	};
}

/*********************************************************************************
nestedPseudoEnergyDP: FOR PARAMETER TUNING. It calculates the free energy of a 
region containing a pseudoknot. This simply calls the Simfold energy function
with the pk regions replaced by <xxx>. Energy value is returned in kcal/mol.

IF THIS FUNCTION IS CALLED, THERE IS AT LEAST ONE PK NESTED INSIDE (arbitrarily deep)
 *********************************************************************************/
double	Loop::nestedPseudoEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles)
	// P_matrix is not modified
{
	int num_params = get_num_params_PK_DP();
	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	int numbases = end - begin + 1;

	// create structure: a string of dot-brackets to represent the region
	char * structure = new char[numbases+1];
	// create sequence
	char * csequence = new char[numbases+1];

	char c_structure[MaxN];

	for (int i = 0; i < numbases; i++)
	{
		csequence[i] = Input->CSequence[begin + i];
		structure[i] = '.';
	}
	csequence[numbases] = '\0';
	structure[numbases] = '\0';
	c_structure[numbases] = '\0';

	// replace PKs with <xxx>
	fillMultiStructure(structure, csequence, numbases, begin);

	//	printf("After call to fillMultiStructure\n");
	if (DEBUG)
	{
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
			c_structure[i] = structure[i];
		}
		printf("\n");
	}

	// fill in the remaining pk-free portions
	for (int i = 0; i < numbases; i++)
	{
		if (structure[i] != '<' && structure[i] != '>' && structure[i] != 'x')
		{
			if (Input->Sequence[begin + i] <= 0)
				structure[i] = '.';
			else if (Input->Sequence[begin + i] > (begin+i))
				structure[i] = '(';
			else
				structure[i] = ')';
		}
	}

	// DEGUG PARAMETER TUNING
	if (DEBUG)
	{	
		printf("Parameter Tuning Input (a multi loop with a pk inside):\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", csequence[i]);
		}
		printf("\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
		}
		printf("\n");
	}

	// call SimFold energy/feature counts function
	float retval = get_feature_counts_restricted (csequence, structure, c, f, reset_c, ignore_dangles, 0);
	if (DEBUG)
		printf("--> Energy from simfold get_feature_counts: %f\n", retval);


	PARAMTYPE dang = 0;  // int, 10cal/mol
	// create pkfree_structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	char pkfree_structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			pkfree_structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			pkfree_structure[i] = '(';
		else
			pkfree_structure[i] = ')';
	}
	pkfree_structure[numbases] = '\0';
	// if there is a PK child, need to add dangles here
	if (ignore_dangles == 0)
	{
		pk_str_features * f = Input->loops;
		int * sequence = Input->type;
		int i_pt = begin;

		for (int j = begin; j < end; j++)
		{
			if (Input->ClosedRegions[j] != NUL && Input->ClosedRegions[j]->type == multi && Input->ClosedRegions[j]->hasPKBranches() == 1)
			{
				// add dangling energies for multi-loop
				dang += dangling_energy_left (sequence, pkfree_structure, i_pt, f[i_pt].pair, f[i_pt].bri[0], f[f[i_pt].bri[0]].pair);
				if (DEBUG)
					printf("MULTI: dang left = %d\n", dangling_energy_left (sequence, pkfree_structure, i_pt, f[i_pt].pair, f[i_pt].bri[0], f[f[i_pt].bri[0]].pair));	// Cristina: July 17, 2007
				if (c != NUL)    count_LEdangling_energy_left (sequence, pkfree_structure, -1, i_pt, f[i_pt].pair, f[i_pt].bri[0], f[f[i_pt].bri[0]].pair, c);

				for (int l=0; l < f[i_pt].num_branches - 1; l++){
					dang += dangling_energy (sequence, pkfree_structure, f[i_pt].bri[l], f[f[i_pt].bri[l]].pair, f[i_pt].bri[l+1], f[f[i_pt].bri[l+1]].pair);
					if (DEBUG)
						printf("MULTI: dang %d = %d: %d %d %d %d, %d %d %d %d\n", l, dangling_energy (sequence, pkfree_structure, f[i_pt].bri[l], f[f[i_pt].bri[l]].pair, f[i_pt].bri[l+1], f[f[i_pt].bri[l+1]].pair),
								f[i_pt].bri[l], f[f[i_pt].bri[l]].pair, f[i_pt].bri[l+1], f[f[i_pt].bri[l+1]].pair,
								sequence[f[i_pt].bri[l]], sequence[f[f[i_pt].bri[l]].pair], sequence[f[i_pt].bri[l+1]], sequence[f[f[i_pt].bri[l+1]].pair]);	// Cristina: July 17, 2007
					if (c != NUL)    count_LEdangling_energy (sequence, pkfree_structure, -1, f[i_pt].bri[l], f[f[i_pt].bri[l]].pair, f[i_pt].bri[l+1], f[f[i_pt].bri[l+1]].pair, c);
				}
				/*
				   printf("\nNon-Zero Simfold Counter Values:\n");
				   for (int i = 0; i < get_num_params(); i++)
				   if (c[i] != 0.0)
				   printf("c[%d]=%f  ", i, c[i]);
				   printf("\n");
				   */
				dang += dangling_energy_right (sequence, pkfree_structure, i_pt, f[i_pt].pair, f[i_pt].bri[f[i_pt].num_branches-1], f[f[i_pt].bri[f[i_pt].num_branches-1]].pair);
				if (DEBUG)
					printf("MULTI: dang right = %d\n", dangling_energy_right (sequence, pkfree_structure, i_pt, f[i_pt].pair, f[i_pt].bri[f[i_pt].num_branches-1], f[f[i_pt].bri[f[i_pt].num_branches-1]].pair));	// Cristina: July 17, 2007
				if (c != NUL)    count_LEdangling_energy_right (sequence, pkfree_structure, -1, i_pt, f[i_pt].pair, f[i_pt].bri[f[i_pt].num_branches-1], f[f[i_pt].bri[f[i_pt].num_branches-1]].pair, c);
				/*
				   printf("\nNon-Zero Simfold Counter Values:\n");
				   for (int i = 0; i < get_num_params(); i++)
				   if (c[i] != 0.0)
				   printf("c[%d]=%f  ", i, c[i]);
				   printf("\n");
				   */

			}
			// add "no-dangling" restriction                                    
			for (int l=0; l < f[i_pt].num_branches; l++)
			{
				Input->cannot_add_dangling [f[i_pt].bri[l] -1] = 1;
				Input->cannot_add_dangling [f[f[i_pt].bri[l]].pair + 1] = 1;
			}
			Input->cannot_add_dangling [i_pt] = 1;  // restrict start of multiloop from having an AU penalty added
		}
	}

	return retval + float(dang)/100.0;
}


/*********************************************************************************
pseudoEnergyDP: it calculates the free energy of a pseudoknotted loop.
It also calculates the free energy of interior-pseudoknotted and multi-pseudoknotted
loops regarding the pseudoknotted loop and add them to the energy
of pseudoknotted loop. Energy value is returned in 10*cal/mol.
 *********************************************************************************/

float	Loop::pseudoEnergyDP(){

	double Energy = 0;
	double initPenalty = 0;

	// Determine the number of unpaired bases in the pseudoloop
	// NOTE: NumberOfUnpaird does not give the desired result for pseudoloops! Don't use it!

	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	// the number of unpaired bases is the number of unpaired bases between:
	//  a' and c
	//  c' and b'
	//  b and d'
	// excluding any bases in closed regions nested in the three regions above
	int k = 0;
	int i = begin;
	int numRegions = NumberOfBands * 2 - 1;

	NumberOfUnpairedInPseudo = 0;  // necessary if there is more than one energy model being used

	for (k = 0; k < numRegions; k++)  // look at all regions between bands, as listed above
	{
		if (DEBUG2)
			printf("[pseudoEnergyDP]: c - a' = %d - %d\n", bandpattern->pattern[i].next, bandpattern->pattern[i].OtherBorder + 1);  
		NumberOfUnpairedInPseudo += bandpattern->pattern[i].next - (bandpattern->pattern[i].OtherBorder + 1);
		i = bandpattern->pattern[i].next;
	}
	NumberOfUnpairedInPseudo += NumberOfUnpairedInUnbandChild;  // subtract bases in unband children

	// energy in kcal/mol
	Energy = ( pkmodelDP.Pb * NumberOfBands  // number of bands
			+ pkmodelDP.Pps * NumberOfUnBandChild  // number of nested closed regions inside pseudoloop (not inside bands)
			//				+ pkmodelDP.Pup * NumberOfUnpaird );  // number of unpaired bases
			// NOTE: the number of unpaired bases in a pseudoloop is anything between inner band borders outside of nested closed regions (excluded in PseudoNestedCheck)
		+ pkmodelDP.Pup * NumberOfUnpairedInPseudo );

	if (DEBUG)
		printf("[Pseudoknot init penalty DP]:"); 

	if (Parent->type == multi || nested == inBand) {  // inBand children must be in a multiloop since the pseudoknotted base pair is a branch
		initPenalty = pkmodelDP.Psm;  // for pseudoloop inside a multiloop or multiloop that spans a band
		if (DEBUG)
			printf("Psm\n"); 
	}
	else if (nested == unBand) {  // Parent->type == pseudo
		initPenalty = pkmodelDP.Psp;  // for pseudoloop inside another pseudoloop (not in a band)
		if (DEBUG)
			printf("Psp\n"); 
	}
	else {
		initPenalty = pkmodelDP.Ps;  // for exterior pseudoloop
		if (DEBUG)
			printf("Ps\n"); 
	}

	Energy += initPenalty;  // penalty for initiating a pseudoknot

	if (DEBUG2)
		printf("[PseudoEnergyDP] without energy of loops spanning bands (10cal/mol): %f, NumBands: %d, NumUnpaired: %d, NumberOfUnpairedInPseudo = %d, NumOfBandChild = %d, NumberOfUnBandChild: %d\n", 
				100*Energy, NumberOfBands, NumberOfUnpaird, NumberOfUnpairedInPseudo, NumberOfBandChild, NumberOfUnBandChild);
	Energy *= 100;  // energy in 100*kcal/mol = 10 cal/mol

	T_IntList * L1 = ILoops;
	T_IntList * L2 = MLoops;

	// Calculating the free energy of internal loops that span a band (interior-pseudoknotted loops)

	while (L1 != NUL){
		float en = Input->looplists[L1->Num]->interiorPseudoEnergyDP();
		if (DEBUG)
			printf("[ILoopsDP]: (%d, %d) Energy: %.2f 10cal/mol\n", L1->Num, Input->BasePair(L1->Num), en); fflush(stdout);
		Energy += en;
		L1 = L1->Next;
	};

	// Calculating the free energy of multiloops that span a band (multi-pseudoknotted loops)

	while (L2 != NUL){
		float en = Input->looplists[L2->Num]->multiPseudoEnergyDP();
		if (DEBUG)
			printf("[MLoopsDP]: (%d, %d) Energy: %.2f 10cal/mol\n", L2->Num, Input->BasePair(L2->Num), en); fflush(stdout);
		Energy += en;
		if (DEBUG)
			printf("[Resulted Energy DP} : %.2f 10cal/mol\n", Energy);
		L2 = L2->Next;
	};

	if (DEBUG)
		printf("[PseudoEnergyDP] energy of [%d, %d] is %f 10cal/mol\n", begin, end, Energy);

	return Energy;

}

/*
// FOR PARAMETER TUNING -- NOT USED!!
void Loop::fillPseudoStructure(char * structure, char * csequence, int a, int ap, int bp, int b)
{
// consider only pseudoknots inside [a,ap] |_| [bp,b]
Loop * L = RightChild;
while (L != NUL){
if (L->type == pseudo)
{
// PK nested in band
if ( (L->begin > a && L->end < ap) || (L->begin > bp && L->end < b) )
{
structure[L->begin - begin] = '<';
structure[L->end - begin] = '>';
for (int i = L->begin - begin + 1; i < L->end - begin; i++)
structure[i] = 'x';
}
}
else
{
L->fillPseudoStructure(structure, csequence, a, ap, bp, b);
}
L = L->LeftSibling;
};
}
*/


// FOR PARAMETER TUNING
// used for a band [a,ap] |_| [bp,b] - assumes there can be pseudoknots inside the band 
// (ie. there are multiloops that span the band)
void Loop::fillPseudoStructure(char * csequence, char * structure, int len, int a, int ap, int bp, int b)
{
	//	if (DEBUG)
	//		printf("fillPseudoStructure: a,ap,bp,b = %d,%d,%d,%d\n",a,ap,bp,b);

	// replace PKs with <xxx>
	fillMultiStructure(structure, csequence, b-a+1, a);

	// for a band -- replace ap.bp with (xxx)
	// consider only those base pairs inside [a,ap] |_| [bp,b]
	for (int i = 0; i < len; i++)
	{
		if ((i > ap-a) && (i < bp-a))
		{
			structure[i] = 'x';
			//			printf("x placed at i = %d", i);
		}
		else
		{
			if (structure[i] != '<' && structure[i] != '>' && structure[i] != 'x')
			{
				if (Input->Sequence[a + i] <= 0)
					structure[i] = '.';
				else if (Input->Sequence[a + i] > (a+i))
					structure[i] = '(';
				else
					structure[i] = ')';
			}
		}
	}
}


// FOR PARAMETER TUNING
// used for part of a band or a band [a,ap] |_| [bp,b], but assumes there are no multiloops that span the band
void Loop::fillPseudoStructureNoMulti(char * csequence, char * structure, int len, int a, int ap, int bp, int b)
{
	//	if (DEBUG)
	//		printf("fillPseudoStructureNoMulti: a,ap,bp,b = %d,%d,%d,%d\n",a,ap,bp,b);

	// for stacks/internal loops that span a band -- replace ap.bp with (xxx)
	// consider only those base pairs inside [a,ap] |_| [bp,b]
	for (int i = 0; i < len; i++)
	{
		if ((i > ap-a) && (i < bp-a))
		{
			structure[i] = 'x';
			//			printf("x placed at i = %d", i);
		}
		else
		{
			if (Input->Sequence[a + i] <= 0)
				structure[i] = '.';
			else if (Input->Sequence[a + i] > (a+i))
				structure[i] = '(';
			else
				structure[i] = ')';
		}
	}
}

// return 1 if i_pt is the beginning of a multiloop that spans a band
int Loop::startsMultiSpanningBand(int i_pt)
{
	int startsMulti = 0;

	T_IntList * L2 = MLoops;
	while (L2 != NUL)
	{
		if (Input->looplists[L2->Num]->base1 == i_pt)  // if the start of this multi matches i_pt
		{
			startsMulti = 1;
			break;
		}
		L2 = L2->Next;
	}
	return 0;
}

// NOT  USED !!!
void Loop::setMultiPseudoLoopDangles(int a, int ap, int bp, int b)
	// PRE: a,ap |_| bp,b is a region in which the multiloop of interest resides:
	//		(a,b) is the outer pair of the multiloop; (ap,bp) is the inner band border of the band 
	//        which the multiloop spans
	// POST: set dangling restrictions on the pseudoknotted child of the multiloop
{
	int i = a;
	for (i = a; i < ap; i++)
	{
		if (Input->ClosedRegions[i] != NUL && Input->ClosedRegions[i]->end >= bp)
			break;
	}

	// now i represents the starting index of the pseudoknotted child
	Input->cannot_add_dangling[i-1] = 1;
	Input->cannot_add_dangling[Input->BasePair(i)+1] = 1;
}


/*********************************************************************************
pseudoEnergyDP: FOR PARAMETER TUNING. it calculates the free energy of a pseudoknotted loop.
It also calculates the free energy of interior-pseudoknotted and multi-pseudoknotted
loops regarding the pseudoknotted loop and add them to the energy
of pseudoknotted loop. Energy value is returned in kcal/mol.
 *********************************************************************************/

float	Loop::pseudoEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles)
	// P_matrix is modified to reflect the loops that spans bands
{
	int num_params = get_num_params_PK_DP();
	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	// NOTE: the energy value returned by the simfold function is not used here
	//       instead we use the original pseudoloop energy function

	// Calculate Energy
	double Energy = 0;
	double initPenalty = 0;

	int numbases = 0;
	// create structure: a string of dot-brackets to represent the region
	char * structure;
	// create sequence
	char * csequence;
	// create a temp counter
	double * c_temp;
	double f_temp = 0;

	//	int num_params = get_num_params_PK_DP();
	int num_params_pkfree = get_num_params();  // number of simfold parameters

	int k_pt = 0;
	int i_pt = begin;

	int last_border;
	int current_border;

	T_IntList * L1 = ILoops;
	T_IntList * L2 = MLoops;
	int a, ap, bp, b;

	float pkfree_retval = 0;

	pk_str_features * feat = Input->loops;
	int * sequence = Input->type;
	int AUpen = 0;  // in 10cal/mol


	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	// we want to fill in the structure between each band in turn, looking for
	// other nested pseudoknots
	// [a,a']|_|[b',b]
	// [c,c']|_|[d',d]
	for (k_pt = 0; k_pt < NumberOfBands*2; k_pt++)
	{
		// per each band
		L1 = ILoops;
		L2 = MLoops;
		last_border = i_pt;
		current_border = 0;

		// multiloops that span bands are handled separately (since has parameters not in simfold) so search for them
		while (L2 != NUL){
			// find the next multiloop that spans a band -- get energy/counts for everything above it
			// that is also inside the last multiloop  i.e. look through ILoops again for the first
			// loop that is above this multiloop but below the last one

			if (DEBUG)
				printf("Consider multiloop starting at %d in band starting at %d,%d\n", Input->looplists[L2->Num]->base1, i_pt, bandpattern->pattern[i_pt].OtherBorder);

			if (Input->looplists[L2->Num]->base1 >= i_pt 
					&& Input->looplists[L2->Num]->base1 <= bandpattern->pattern[i_pt].OtherBorder)
			{
				if (DEBUG)
					printf("This multiloop works %d\n", Input->looplists[L2->Num]->base1);

				current_border = Input->looplists[L2->Num]->base1;  // the current multiloop we are looking at

				if (current_border == last_border && current_border != i_pt)
				{
					printf("WARNING: should have moved onto a different multiloop that spans the band.\n");
				}

				// look through ILoops again for the first loop that is above this multiloop but below the last one
				L1 = ILoops;
				// if ILoops is null, no need to do anything since no 
				while (L1 != NUL){
					if (Input->looplists[L1->Num]->base1 >= last_border && Input->looplists[L1->Num]->base1 <= current_border){
						a = Input->looplists[L1->Num]->base1;
						ap = Input->looplists[L1->Num]->base2;
						bp = Input->BasePair(ap);
						b = Input->BasePair(a);

						numbases = b-a+1;
						structure = new char[numbases+1];
						csequence = new char[numbases+1];
						for (int i = 0; i < numbases; i++){
							csequence[i] = Input->CSequence[a + i];
						}
						csequence[numbases] = '\0';
						structure[numbases] = '\0';

						// now get energy/counts for everything above current multiloop

						fillPseudoStructureNoMulti(csequence, structure, numbases, a,ap,bp,b);

						if (c != NUL){
							c_temp = new double[num_params];
							f_temp = 0;
							pkfree_retval = get_feature_counts_restricted(csequence, structure, c_temp, f_temp, 1, ignore_dangles, 1);
						}else{
							f_temp = 0;
							pkfree_retval = get_feature_counts_restricted(csequence, structure, NUL, f_temp, 1, ignore_dangles, 1);
						}

						if (ap == a+1){  // stack spanning a band
							Energy += pkmodelDP.stP * pkfree_retval;
							if (c != NUL){
								for (int i = 0; i < num_params_pkfree; i++){
									P_matrix[i][num_params_pkfree + structure_type_index_PK( (char *) "stp") - 1] += c_temp[i];
								}
								// update the single stP entry in counter
								c[num_params_pkfree + structure_type_index_PK( (char *) "stp") - 1] += f_temp;
							}
						}else{  // internal loop spanning band
							Energy += pkmodelDP.intP * pkfree_retval;
							if (c != NUL){
								for (int i = 0; i < num_params_pkfree; i++){
									// if c_temp[i] == 0, P_matrix remains unchanged, as needed
									P_matrix[i][num_params_pkfree + structure_type_index_PK( (char *) "intp") - 1] += c_temp[i];
								}
								// update the single stP entry in counter
								c[num_params_pkfree + structure_type_index_PK( (char *) "intp") - 1] += f_temp;
							}
						}
						L1->tuning_flag = 1;  // set dirty bit flag
						// clear structure and csequence for the next spanning loop (which possibly needs these to be different length)
						delete [] structure;
						delete [] csequence;
						if (c != NUL) { delete c_temp; }
					}

					L1 = L1->Next;  // go to next possible loop within this band that is above the current multiloop
				}

				// find energy/counts for the multiloop
				float en = Input->looplists[L2->Num]->multiPseudoEnergyDP(c, f, reset_c, ignore_dangles);
				Energy += en/100;  // since Energy is in kcal/mol
				// add AU penalty for the pseudoknotted child of the multiloop and set dangling restrictions
				// this child is either the next ILoop member, or if it is a single pair then the only base
				// pair that spans the band
				L1 = ILoops;
				while (L1 != NUL && Input->looplists[L1->Num]->base1 <= current_border)  // ie. L1 outside multi
					L1 = L1->Next;
				// now L1 points to the first L1 inside the multi,
				// or is NUL (single base pair) (or is not null since there might be other ILoops members but not in this band)
				if (L1 == NUL || bandpattern->pattern[Input->looplists[L1->Num]->base1].OtherBorder >= bandpattern->pattern[Input->looplists[L2->Num]->base1].OtherBorder){
					int single_base = Input->looplists[L2->Num]->base1+1;
					while (Input->BasePair(single_base) <= bandpattern->pattern[i_pt].next) // doesn't span band
						single_base++;
					// now single_base is the one that spans the band
					AUpen += AU_penalty (sequence[single_base], sequence[feat[single_base].pair]);
					if (c != NUL)    count_AU_penalty (sequence[single_base], sequence[feat[single_base].pair], c);

				}else{
					AUpen += AU_penalty (sequence[Input->looplists[L1->Num]->base1], sequence[feat[Input->looplists[L1->Num]->base1].pair]);
					if (c != NUL)    count_AU_penalty (sequence[Input->looplists[L1->Num]->base1], sequence[feat[Input->looplists[L1->Num]->base1].pair], c);

				}

				last_border = current_border;  // this multiloop is done
			}

			// else: this loop isn't in the band we are currently looking at

			// go to the next multiloop
			L2 = L2->Next;
		};
		// at this point, we've handled the energy of multiloops in the band and some of the internal loops or stacks
		// will take care of any missed internal loops at the end

		i_pt = bandpattern->pattern[i_pt].next;  // go to next band
	}

	// take care of any missed internal loops / stacks

	// go back through ILoops and get energy of anything we missed
	L1 = ILoops;
	while (L1 != NUL){
		if (L1->tuning_flag != 1){
			a = Input->looplists[L1->Num]->base1;
			ap = Input->looplists[L1->Num]->base2;
			bp = Input->BasePair(ap);
			b = Input->BasePair(a);

			numbases = b-a+1;
			structure = new char[numbases+1];
			csequence = new char[numbases+1];
			for (int i = 0; i < numbases; i++){
				csequence[i] = Input->CSequence[a + i];
			}
			csequence[numbases] = '\0';
			structure[numbases] = '\0';

			// get energy/counts for everything above current multiloop
			fillPseudoStructureNoMulti(csequence, structure, numbases, a,ap,bp,b);

			if (c != NUL){
				c_temp = new double[num_params];
				f_temp = 0;
				pkfree_retval = get_feature_counts_restricted(csequence, structure, c_temp, f_temp, 1, ignore_dangles, 1);
			}else{
				f_temp = 0;
				pkfree_retval = get_feature_counts_restricted(csequence, structure, NUL, f_temp, 1, ignore_dangles, 1);	
			}


			if (ap == a+1){  // stack spanning a band
				Energy += pkmodelDP.stP * pkfree_retval;
				// UPDATE COUNTS HERE  -- note that to get only the counts for the stacked pair, you should probably not
				//                        reuse the same counter (which has been keeping track of ALL the stacked pairs in
				//						  the structure)
				// update the P_matrix (the stP column)
				if (c != NUL){
					for (int i = 0; i < num_params_pkfree; i++){
						// if c_temp[i] == 0, P_matrix remains unchanged, as needed
						P_matrix[i][num_params_pkfree + structure_type_index_PK( (char *) "stp") - 1] += c_temp[i];
					}
					// update the single stP entry in counter
					c[num_params_pkfree + structure_type_index_PK( (char *) "stp") - 1] += f_temp;
				}
			}else{  // internal loop spanning band
				Energy += pkmodelDP.intP * pkfree_retval;
				// UPDATE COUNTS HERE  -- note that to get only the counts for the stacked pair, you should probably not
				//                        reuse the same counter (which has been keeping track of ALL the stacked pairs in
				//						  the structure)
				// update the P_matrix (the intP column)
				if (c != NUL){
					for (int i = 0; i < num_params_pkfree; i++){
						// if c_temp[i] == 0, P_matrix remains unchanged, as needed
						P_matrix[i][num_params_pkfree + structure_type_index_PK( (char *) "intp") - 1] += c_temp[i];
					}
					// update the single stP entry in counter
					c[num_params_pkfree + structure_type_index_PK( (char *) "intp") - 1] += f_temp;
				}
			}
			L1->tuning_flag = 1;  // set dirty bit flag
			// clear structure and csequence for the next spanning loop (which possibly needs these to be different length)
			delete [] structure;
			delete [] csequence;
			if (c != NUL)  { delete c_temp; }
		}
		L1 = L1->Next;
	}
	// at this point, we've handled the energies of each individual band

	i_pt = begin;
	k_pt = 0;

	// add the appropriate AU penalties
	for (k_pt = 0; k_pt < NumberOfBands*2; k_pt++){  // look at all the bands, as listed above
		if (Input->BasePair(i_pt) > i_pt){  // if it is the [a,ap] portion of the band
			// if the outermost loop spanning the band is not a multiloop, add AU
			// the cdt is unnecessary, based on old hotknots code (energydangling function) and simfold energy of a multiloop

			AUpen += AU_penalty (sequence[i_pt], sequence[feat[i_pt].pair]);
			if (c != NUL)    count_AU_penalty (sequence[i_pt], sequence[feat[i_pt].pair], c);
		}

		i_pt = bandpattern->pattern[i_pt].next;  // next band
	}
	Energy += float(AUpen)/100;  // since AUpen is in 10cal/mol, and everything else is in kcal/mol


	// Calculate Energy/counts of pseudoknot

	// Determine the number of unpaired bases in the pseudoloop
	// NOTE: NumberOfUnpaird does not give the desired result for pseudoloops! Don't use it!

	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	// the number of unpaired bases is the number of unpaired bases between:
	//  a' and c
	//  c' and b'
	//  b and d'
	// excluding any bases in closed regions nested in the three regions above
	int k = 0;
	int i = begin;
	int numRegions = NumberOfBands * 2 - 1;

	NumberOfUnpairedInPseudo = 0;  // necessary if there is more than one energy model being used

	for (k = 0; k < numRegions; k++){  // look at all regions between bands, as listed above
		NumberOfUnpairedInPseudo += bandpattern->pattern[i].next - (bandpattern->pattern[i].OtherBorder + 1);
		i = bandpattern->pattern[i].next;
	}
	NumberOfUnpairedInPseudo += NumberOfUnpairedInUnbandChild;  // subtract bases in unband children

	// energy in kcal/mol
	Energy += ( pkmodelDP.Pb * NumberOfBands  // number of bands
			+ pkmodelDP.Pps * NumberOfUnBandChild  // number of nested closed regions inside pseudoloop (not inside bands)
			//				+ pkmodelDP.Pup * NumberOfUnpaird );  // number of unpaired bases
			// NOTE: the number of unpaired bases in a pseudoloop is anything between inner band borders outside of nested closed regions (excluded in PseudoNestedCheck)
		+ pkmodelDP.Pup * NumberOfUnpairedInPseudo );

	g_count_Pb += NumberOfBands;
	if (c != NUL)    c[num_params_pkfree + structure_type_index_PK( (char *) "pb") - 1] += float(NumberOfBands);
	g_count_Pps += NumberOfUnBandChild;
	if (c != NUL)    c[num_params_pkfree + structure_type_index_PK( (char *) "pps") - 1] += float(NumberOfUnBandChild);
	g_count_Pup += NumberOfUnpairedInPseudo;
	if (c != NUL)    c[num_params_pkfree + structure_type_index_PK( (char *) "pup") - 1] += float(NumberOfUnpairedInPseudo);


	if (Parent->type == multi || nested == inBand) {  // inBand children must be in a multiloop since the pseudoknotted base pair is a branch
		initPenalty = pkmodelDP.Psm;  // for pseudoloop inside a multiloop or multiloop that spans a band
		g_count_Psm++;
		if (c != NUL)    c[num_params_pkfree + structure_type_index_PK( (char *) "psm") - 1] += 1;

	}else if (nested == unBand) {  // Parent->type == pseudo
		initPenalty = pkmodelDP.Psp;  // for pseudoloop inside another pseudoloop (not in a band)
		g_count_Psp++;
		if (c != NUL)    c[num_params_pkfree + structure_type_index_PK( (char *) "psp") - 1] += 1;
	}else {
		initPenalty = pkmodelDP.Ps;  // for exterior pseudoloop
		g_count_Ps++;
		if (c != NUL)    c[num_params_pkfree + structure_type_index_PK( (char *) "ps") - 1] += 1;
	}

	Energy += initPenalty;  // penalty for initiating a pseudoknot

	return Energy;

}


///////////////////////////////////////////////////////////////////
// RIVAS & EDDY ENERGY MODEL
///////////////////////////////////////////////////////////////////

/**
 *  Rivas and Eddy energy model 
 */

float Loop::getEnergyRE(){
	float sum = 0;
	float energyToAdd = 0;

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)
	Loop * L = RightChild;
	while (L != NUL){
		sum += L->getEnergyRE();
		L = L->LeftSibling;
	};

	switch (type){
		case	stackloop:
			energyToAdd = stackEnergy();
			//			printf("***[%d, %d] Stacked Pair Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	hairpin	:
			energyToAdd = hairpinEnergy();
			//		printf("***[%d, %d] Hairpin Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	interior:
			energyToAdd = interiorEnergy();
			//printf("***[%d, %d] Internal Loop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	multi	:
			energyToAdd = multiEnergyRE();
			//	printf("***[%d, %d] Multiloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	pseudo	:
			energyToAdd = pseudoEnergyRE();
			//printf("***[%d, %d] Pseudoloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		default:
			break;
	};

	finalEnergy = energyToAdd;
	return sum;
}



//=================================================//
/*
 *	Rivas and Eddy energy model for multiEnergyRE() function

Cristina: added structure as a new parameter to dangling_energy functions in commonPK.cpp
*/

double	Loop::multiEnergyRE(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = begin;
	int AUpen, h, l;
	int dang;
	int misc_energy;

	dang = 0;
	misc_energy = 0;
	AUpen = 0;
	int special;
	special = 0;

	int numbases = end - begin + 1;
	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	// add the energies/enthalpies for free bases
	for (l=i+1; l < f[i].bri[0]; l++)
		misc_energy += misc.multi_free_base_penalty;

	if (DEBUG)
		printf("[misc_energyRE]1 is %d\n", misc_energy);
	for (h=0; h < f[i].num_branches-1; h++)
	{
		for (l = f[f[i].bri[h]].pair + 1; l < f[i].bri[h+1]; l++)
			misc_energy += misc.multi_free_base_penalty;
	}
	for (l = f[f[i].bri[f[i].num_branches-1]].pair + 1; l < f[i].pair; l++)
		misc_energy += misc.multi_free_base_penalty;

	if (DEBUG)
		printf("[misc_energyRE]2 is %d\n", misc_energy);
	misc_energy += misc.multi_offset;

	if (DEBUG)
		printf("[misc_energyRE]3 is %d\n", misc_energy);
	if (DEBUG)
		printf("[%%%%%%%%%%%%%%%%%%%%%%%%5]RE misc.u.tidfjfwrw = %d\n", misc.multi_offset);

	misc_energy += misc.multi_helix_penalty * (f[i].pseudo_num_branches + 1);




	// add AU_penalties for multi-loop
	AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);

	for (h=0; h < f[i].num_branches; h++)
		AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);

	// add dangling energies for multi-loop
	dang += dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);

	for (l=0; l < f[i].num_branches - 1; l++){
		dang += dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
	};
	dang += dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);

	// add "no-dangling" restriction                                    
	for (l=0; l < f[i].num_branches; l++)
	{
		Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
		Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
	}
	if (DEBUG)
	{
		printf ("%d - RE multi m\t- add energy %6d\n", i, misc_energy);            
		printf ("%d - RE multi d\t- add energy %6d\n", i, dang);
		printf ("%d - RE multi AU\t- add energy %6d\n", i, AUpen);
	}                

	if (DEBUG)            
		printf ("%d RE multi \t- add energy %6d\n", i, misc_energy + dang + AUpen);        
	return misc_energy + dang + AUpen;  


}



/*
 *	Rivas and Eddy energy model for pseudoEnergyRE() function
 */

float	Loop::pseudoEnergyRE(){


	float Energy;


	Energy = (pkmodelRE.Gw +   // starting a pseudoloop
			pkmodelRE.Gwh * (NumberOfBands - 2) +
			pkmodelRE.Q_tilda * NumberOfUnpaird +  // unpaired bases in pseduoknot
			pkmodelRE.P_tilda * 2 * NumberOfBands +  // pairs (band) in a pseudoknot
			pkmodelRE.P_i * NumberOfUnBandChild);  //nested  // CHECK_RE: should be P_tilda?

	/*
	   float P_tilda = 0.1; // pair in a pseudoknot
	   float Q_tilda = 0.2; // unpaired bases in a pseudoknot
	   float P_i = 0.1; // for E&R energy model
	   float Gw = 7; //starting a pseudoknot loop
	   float Gwh = 6; //each more band region

	   Energy = (Gw+   // Starting a pseudoloop
	   Gwh * (NumberOfBands -2 ) +
	   Q_tilda* NumberOfUnpaird + //unpaired bases part
	   P_tilda * 2 *NumberOfBands + //bands part
	   P_i * NumberOfUnBandChild);        //nested 
	   */

	if (DEBUG)
		printf("[PseudoEnergyRE] Energy (10cal/mol): %.2f, NumBands: %d,  NumUnpaired: %d,   NumberOfUnBandChild: %d\n", 100*Energy, 
				NumberOfBands, NumberOfUnpaird,  NumberOfUnBandChild);
	Energy *= 100;  // multilpy by 100 since rest of function has energy values in 10cal/mol


	//Calculating the free energy of interior-pseudoknotted loops
	T_IntList * L1 = ILoops;
	T_IntList * L2 = MLoops;

	while (L1 != NUL){
		float en = Input->looplists[L1->Num]->interiorPseudoEnergyRE();
		if (DEBUG)
			printf("[ILoops RE]: (%d, %d) Energy: %.2f \n", L1->Num, Input->BasePair(L1->Num), en); fflush(stdout);
		Energy += en;
		L1 = L1->Next;
	};

	//Calculating the free energy of multi-pseudoknotted loops


	while (L2 != NUL){
		float en = Input->looplists[L2->Num]->multiPseudoEnergyRE();
		if (DEBUG)
			printf("[MLoopsRE]: (%d, %d) Energy: %.2f\n", L2->Num, Input->BasePair(L2->Num), en); fflush(stdout);
		Energy += en;
		if (DEBUG)
			printf("[Resulted Energy RE} : %.2f \n", Energy);
		L2 = L2->Next;
	};

	if (DEBUG)            
		printf("[PseudoEnergyRE] energy of [%d, %d] is %f\n", begin, end, Energy);

	return Energy;

}


///////////////////////////////////////////////////////////////////
// CAO & CHEN 2006 (a) ENERGY MODEL
///////////////////////////////////////////////////////////////////

/*********************************************************************************
getEnergy: calculating the summation of the free energy of the loops in
the secondary structure. Energy value is returned in 10*cal/mol.
 *********************************************************************************/

float Loop::getEnergyCC2006a(){
	float sum = 0;
	float energyToAdd = 0;

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)
	Loop * L = RightChild;
	while (L != NUL){
		sum += L->getEnergyCC2006a();
		L = L->LeftSibling;
	};

	switch (type){
		case	stackloop:
			energyToAdd = stackEnergy();
			//			energyTrace += "***[" + begin + ", " + end + "] Stacked Pair Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Stacked Pair Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	hairpin	:
			energyToAdd = hairpinEnergy();
			//			energyTrace += "***[" + begin + ", " + end + "] Hairpin Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Hairpin Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	interior:
			energyToAdd = interiorEnergy();
			//			energyTrace += "***[" + begin + ", " + end + "] Internal Loop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Internal Loop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	multi	:
			energyToAdd = multiEnergyCC2006a();
			//			energyTrace += "***[" + begin + ", " + end + "] Multiloop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Multiloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	pseudo	:
			energyToAdd = pseudoEnergyCC2006a();
			//			energyTrace += "***[" + begin + ", " + end + "] Pseudoloop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Pseudoloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		default:
			break;
	};

	finalEnergy = energyToAdd;
	return sum;
}


/*********************************************************************************
getEnergy: calculating the summation of the free energy of the loops in
the secondary structure. Energy value is returned in 10*cal/mol.
 *********************************************************************************/

float Loop::getEnergyCC2006b(){
	float sum = 0;
	float energyToAdd = 0;

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)
	Loop * L = RightChild;
	while (L != NUL){
		sum += L->getEnergyCC2006b();
		L = L->LeftSibling;
	};

	switch (type){
		case	stackloop:
			energyToAdd = stackEnergy();
			//			energyTrace += "***[" + begin + ", " + end + "] Stacked Pair Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Stacked Pair Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	hairpin	:
			energyToAdd = hairpinEnergy();
			//			energyTrace += "***[" + begin + ", " + end + "] Hairpin Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Hairpin Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	interior:
			energyToAdd = interiorEnergy();
			//			energyTrace += "***[" + begin + ", " + end + "] Internal Loop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Internal Loop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	multi	:
			energyToAdd = multiEnergyCC2006a();  // NOT A MISTAKE, Versions a and b have the same multiloop function
			//			energyTrace += "***[" + begin + ", " + end + "] Multiloop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Multiloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		case	pseudo	:
			energyToAdd = pseudoEnergyCC2006b();
			//			energyTrace += "***[" + begin + ", " + end + "] Pseudoloop Energy = " + 10*energyToAdd + " cal/mol\n";
			//			printf("***[%d, %d] Pseudoloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
			sum += energyToAdd;
			break;
		default:
			break;
	};

	finalEnergy = energyToAdd;
	return sum;
}


/*********************************************************************************
getEnergy: FOR PARAMETER TUNING. calculating the summation of the free energy of the loops in
the secondary structure. Energy value is returned in kcal/mol.
 *********************************************************************************/
float Loop::getEnergyCC2006b(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles)
{
	float sum = 0;
	float energyToAdd = 0;
	float tempsum = 0;

	int num_params = get_num_params_PK_CC2006b();
	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	Loop * L;

	if (begin == 0)
	{
		if (RightChild == NUL)  // no loops in the structure
			return 0;
		else
		{
			L = RightChild;
			while (L != NUL) {
				sum += L->getEnergyCC2006b(P_matrix,c,f,reset_c,ignore_dangles);
				L = L->LeftSibling;
			}
			return sum;
		}
	}

	//	if (begin == 0)  // if we are at the root of the tree
	//		L = RightChild;
	//	else  // if we are not
	L = this;

	if (DEBUG)
		printf("At top of getEnergyCC2006b for loop (%d, %d)\n", L->begin, L->end);

	//	float sum = 0;
	//	float energyToAdd = 0;

	// determine the energy of the whole structure by summing across rows 
	// (that is, sum energies across the first row starting with the rightmost child and going left;
	// (then the energy of each child is the sum of energies of each of its children, starting from the rightmost one)

	// FOR PARAMETER TUNING:
	// if we find a pseudoknot-free structure (ie. children in the tree rooted at this loop 
	// have no pseudoknots), then simply call simfold's energy function without parsing the
	// tree into the lowest-level structures
	if (L->isPKFree() == 1)
	{
		if (DEBUG)
			printf("getEnergyCC2006b: IDed as pk-free\n");
		finalEnergy = L->pkfreeEnergyCC2006b(P_matrix, c, f, reset_c, ignore_dangles);
		sum += finalEnergy;
		return sum;
	}

	if (DEBUG)
		printf("getEnergyCC2006b: not pkfree, proceeding\n");

	// if we find a non pk-free structure but it is not a pseudoknot then there must be a nested pseudoknot
	if (L->type != pseudo)
	{
		// repeat for the highest level of pseudoknots -- call helper function
		tempsum = 0;
		L->lookForPkCC2006b(P_matrix, c, f, reset_c, ignore_dangles, tempsum);
		sum += tempsum;

		//Loop * L1 = L->RightChild;
		//while (L1 != NUL){

		//	printf("getEnergyCC2006b: pk inside pk-free region (%d,%d)\n", L1->begin, L1->end);
		//	sum += L1->getEnergyCC2006b(c, f, reset_c, ignore_dangles);
		//	printf("getEnergyCC2006b: done pk inside pk-free region (%d,%d)\n", L1->begin, L1->end);
		//	L1 = L1->LeftSibling;
		//};

		//		printf("getEnergyCC2006b: done the highest level of pks\n");
		energyToAdd = L->nestedPseudoEnergyCC2006b(P_matrix, c, f, reset_c, ignore_dangles);
		sum += energyToAdd;
		finalEnergy = energyToAdd;
		return sum;
	}


	if (L->type == pseudo)
	{
		// repeat for anything nested inside
		Loop * L1 = L->RightChild;
		while (L1 != NUL){
			sum += L1->getEnergyCC2006b(P_matrix, c, f, reset_c, ignore_dangles);
			L1 = L1->LeftSibling;
		};
		energyToAdd = L->pseudoEnergyCC2006b(P_matrix, c, f, reset_c, ignore_dangles);
		sum += energyToAdd;
		finalEnergy = energyToAdd;
		return sum;
	}
	else
	{
		printf("WARNING: getEnergyCC2006b(): This case should never happen.\n");
	}


	/*
	   switch (type){
	   case	stackloop:
	   energyToAdd = stackEnergy();
//			energyTrace += "***[" + begin + ", " + end + "] Stacked Pair Energy = " + 10*energyToAdd + " cal/mol\n";
//			printf("***[%d, %d] Stacked Pair Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
sum += energyToAdd;
break;
case	hairpin	:
energyToAdd = hairpinEnergy();
//			energyTrace += "***[" + begin + ", " + end + "] Hairpin Energy = " + 10*energyToAdd + " cal/mol\n";
//			printf("***[%d, %d] Hairpin Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
sum += energyToAdd;
break;
case	interior:
energyToAdd = interiorEnergy();
//			energyTrace += "***[" + begin + ", " + end + "] Internal Loop Energy = " + 10*energyToAdd + " cal/mol\n";
//			printf("***[%d, %d] Internal Loop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
sum += energyToAdd;
break;
case	multi	:
energyToAdd = multiEnergyCC2006a();  // NOT A MISTAKE, Versions a and b have the same multiloop function
//			energyTrace += "***[" + begin + ", " + end + "] Multiloop Energy = " + 10*energyToAdd + " cal/mol\n";
//			printf("***[%d, %d] Multiloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
sum += energyToAdd;
break;
case	pseudo	:
energyToAdd = pseudoEnergyCC2006b();
//			energyTrace += "***[" + begin + ", " + end + "] Pseudoloop Energy = " + 10*energyToAdd + " cal/mol\n";
//			printf("***[%d, %d] Pseudoloop Energy = %f cal/mol\n", begin, end, 10*energyToAdd);
sum += energyToAdd;
break;
default:
break;
};

finalEnergy = energyToAdd;
return sum;
*/
return 0;
}

// FOR PARAMETER TUNING - helper function for getEnergyCC2006b
void Loop::lookForPkCC2006b(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles, float& sum)
{
	int num_params = get_num_params_PK_CC2006b();
	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	//	printf("lookForPKCC2006b: at top\n");

	if (type == pseudo)
	{
		sum += this->getEnergyCC2006b(P_matrix, c, f, reset_c, ignore_dangles);
		return;
	}

	Loop * L = RightChild;
	while (L != NUL)
	{
		L->lookForPkCC2006b(P_matrix, c, f, reset_c, ignore_dangles, sum);
		L = L->LeftSibling;
	}
}


/*********************************************************************************
pkfreeEnergyCC2006b: FOR PARAMETER TUNING Calculates the free energy of a closed region that is pseudoknot 
free. Calls the SimFold function directly (without further parsing the structure 
into smaller pieces). Also gets the feature counts for parameter tuning.
Energy value is returned in kcal/mol.
 *********************************************************************************/
double Loop::pkfreeEnergyCC2006b(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles)
	// P_matrix remains unchanged
{
	int num_params = get_num_params_PK_CC2006b();
	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	int numbases = end - begin + 1;

	// create structure: a string of dot-brackets to represent the region
	//	char * structure = new char[numbases+1];
	char structure[numbases+1];
	// create sequence
	//	char * csequence = new char[numbases+1];
	char csequence[numbases+1];

	//	printf("Input:CSequence = ");

	for (int i = 0; i < numbases; i++)
	{
		//		printf("%c", Input->CSequence[begin+i]);

		csequence[i] = Input->CSequence[begin + i];
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	csequence[numbases] = '\0';
	structure[numbases] = '\0';

	// DEGUG PARAMETER TUNING
	if (DEBUG)
	{
		printf("Parameter Tuning Input (a pk-free region):\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", csequence[i]);
		}
		printf("\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
		}
		printf("\n");
	}

	// call SimFold energy/feature counts function
	double retval = get_feature_counts_restricted (csequence, structure, c, f, reset_c, ignore_dangles, 0);
	if (DEBUG)
		printf("--> Energy from simfold get_feature_counts: %f\n", retval);

	if (DEBUG)
	{
		printf("Free Value: %f\n", f);
		printf("PK Counter Values:\n");
		for (int i = get_num_params(); i < num_params; i++)
			if (c != NUL && c[i] != 0.0)
				printf("c[%d]=%f  ", i, c[i]);

		printf("\nNon-Zero Simfold Counter Values:\n");
		for (int i = 0; i < get_num_params(); i++)
			if (c != NUL && c[i] != 0.0)
				printf("c[%d]=%f  ", i, c[i]);

		printf("\nAll Non-Zero P_matrix Values:\n");
		for (int i = 0; i < num_params; i++)
		{
			for (int j = i; j < num_params; j++)
				if (P_matrix != NUL && P_matrix[i][j] != 0.0)
					printf("P[%d][%d]=%f  ", i, j, P_matrix[i][j]);
		}
	}

	return retval;
}

int Loop::hasPKBranches()
	// PRE: Loop is a multiloop
	// POST: return 1 if it has PK branches as children
{
	Loop * L = RightChild;
	while (L != NUL)
	{
		if (L->type == pseudo)
		{
			return 1;
		}
		L = L->LeftSibling;
	}
	return 0;
}

/*********************************************************************************
nestedPseudoEnergyCC2006b: FOR PARAMETER TUNING. It calculates the free energy of a 
region containing a pseudoknot. This simply calls the Simfold energy function
with the pk regions replaced by <xxx>. Energy value is returned in kcal/mol.

IF THIS FUNCTION IS CALLED, THERE IS AT LEAST ONE PK NESTED INSIDE (arbitrarily deep)
 *********************************************************************************/
double	Loop::nestedPseudoEnergyCC2006b(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles)
	// P_matrix is not modified
{
	int num_params = get_num_params_PK_CC2006b();
	if (reset_c == 1 && c != NUL)
	{
		f = 0;         
		for (int i=0; i < num_params; i++)
		{
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	int numbases = end - begin + 1;

	// create structure: a string of dot-brackets to represent the region
	char * structure = new char[numbases+1];
	// create sequence
	char * csequence = new char[numbases+1];

	char c_structure[MaxN];

	for (int i = 0; i < numbases; i++)
	{
		csequence[i] = Input->CSequence[begin + i];
		structure[i] = '.';
	}
	csequence[numbases] = '\0';
	structure[numbases] = '\0';
	c_structure[numbases] = '\0';

	// replace PKs with <xxx>
	fillMultiStructure(structure, csequence, numbases, begin);

	//	printf("After call to fillMultiStructure\n");
	if (DEBUG)
	{
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
			c_structure[i] = structure[i];
		}
		printf("\n");
	}

	// fill in the remaining pk-free portions
	for (int i = 0; i < numbases; i++)
	{
		if (structure[i] != '<' && structure[i] != '>' && structure[i] != 'x')
		{
			if (Input->Sequence[begin + i] <= 0)
				structure[i] = '.';
			else if (Input->Sequence[begin + i] > (begin+i))
				structure[i] = '(';
			else
				structure[i] = ')';
		}
	}

	// DEGUG PARAMETER TUNING
	if (DEBUG)
	{	
		printf("Parameter Tuning Input (a multi loop with a pk inside):\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", csequence[i]);
		}
		printf("\n");
		for (int i = 0; i < numbases; i++)
		{
			printf("%c", structure[i]);
		}
		printf("\n");
	}

	// call SimFold energy/feature counts function
	float retval = get_feature_counts_restricted (csequence, structure, c, f, reset_c, ignore_dangles, 0);
	if (DEBUG)
		printf("--> Energy from simfold get_feature_counts: %f\n", retval);


	PARAMTYPE dang = 0;  // int, 10cal/mol
	// create pkfree_structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	char pkfree_structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			pkfree_structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			pkfree_structure[i] = '(';
		else
			pkfree_structure[i] = ')';
	}
	pkfree_structure[numbases] = '\0';
	// if there is a PK child, need to add dangles here
	if (ignore_dangles == 0){
		pk_str_features * f = Input->loops;
		int * sequence = Input->type;
		int i_pt = begin;

		for (int j = begin; j < end; j++){
			if (Input->ClosedRegions[j] != NUL && Input->ClosedRegions[j]->type == multi && Input->ClosedRegions[j]->hasPKBranches() == 1)
			{
				// add dangling energies for multi-loop
				dang += dangling_energy_left (sequence, pkfree_structure, i_pt, f[i_pt].pair, f[i_pt].bri[0], f[f[i_pt].bri[0]].pair);
				if (c != NUL)    count_LEdangling_energy_left (sequence, pkfree_structure, -1, i_pt, f[i_pt].pair, f[i_pt].bri[0], f[f[i_pt].bri[0]].pair, c);

				for (int l=0; l < f[i_pt].num_branches - 1; l++){
					dang += dangling_energy (sequence, pkfree_structure, f[i_pt].bri[l], f[f[i_pt].bri[l]].pair, f[i_pt].bri[l+1], f[f[i_pt].bri[l+1]].pair);
					if (c != NUL)    count_LEdangling_energy (sequence, pkfree_structure, -1, f[i_pt].bri[l], f[f[i_pt].bri[l]].pair, f[i_pt].bri[l+1], f[f[i_pt].bri[l+1]].pair, c);
				}
				dang += dangling_energy_right(sequence, pkfree_structure, i_pt, f[i_pt].pair, f[i_pt].bri[f[i_pt].num_branches-1], f[f[i_pt].bri[f[i_pt].num_branches-1]].pair);
				if (c != NUL)    count_LEdangling_energy_right (sequence, pkfree_structure, -1, i_pt, f[i_pt].pair, f[i_pt].bri[f[i_pt].num_branches-1], f[f[i_pt].bri[f[i_pt].num_branches-1]].pair, c);

			}
			// add "no-dangling" restriction                                    
			for (int l=0; l < f[i_pt].num_branches; l++){
				Input->cannot_add_dangling [f[i_pt].bri[l] -1] = 1;
				Input->cannot_add_dangling [f[f[i_pt].bri[l]].pair + 1] = 1;
			}
			Input->cannot_add_dangling [i_pt] = 1;  // restrict start of multiloop from having an AU penalty added
		}
	}

	return retval + float(dang)/100.0;
}

/*********************************************************************************
multiEnergyCC2006: it calculates the free energy of a multiLoop.
The code is similar to the one in simFold but it has been changed slightly
to consider the pseudoknotted tuples. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
double	Loop::multiEnergyCC2006a(){
	if (no_coax_in_pkfree == 1)
		return multiEnergyDP();

	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = begin;
	int AUpen, h, l;
	int dang;
	double misc_energy;

	dang = 0;
	misc_energy = 0;
	AUpen = 0;
	int special;
	special = 0;

	// add the energies/enthalpies for free bases
	for (l=i+1; l < f[i].bri[0]; l++)  // count the number of unpaired bases from i+1 to the first branch
		misc_energy += pkmodelDP.c;
	//		misc_energy += misc.multi_free_base_penalty;

	// for all the branches in this multiloop...
	for (h=0; h < f[i].num_branches-1; h++){
		// ...count the number of unpaired bases from the end of one branch to the beginning of the next branch
		for (l = f[f[i].bri[h]].pair + 1; l < f[i].bri[h+1]; l++)
			misc_energy += pkmodelDP.c;
		//			misc_energy += misc.multi_free_base_penalty;
	}

	// count the number of unpaired bases from the end of the last branch to the end of the multiloop
	for (l = f[f[i].bri[f[i].num_branches-1]].pair + 1; l < f[i].pair; l++)
		misc_energy += pkmodelDP.c;
	//		misc_energy += misc.multi_free_base_penalty;

	// add penalty for initiating a multiloop
	misc_energy += pkmodelDP.a;
	//	misc_energy += misc.multi_offset;

	// DONECHECK: the +1 in number of branches is for the closing base pair
	// add the penalty for each base pair in the multiloop (the +1 in number of branches is for the closing base pair i,f[i].pair)
	misc_energy += pkmodelDP.b * (f[i].num_branches + 1);
	//	misc_energy += pkmodelDP.b * (f[i].pseudo_num_branches + 1);  // this is if pseudoloops count as 2 branches instead of 1
	//	misc_energy += misc.multi_helix_penalty * (f[i].pseudo_num_branches + 1);

	misc_energy *= 100;  // multiply by 100 since the remaining energies in this function are in 10cal/mol


	// add AU_penalties for multi-loop
	AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);
	//	printf("MULTI: %d.%d, AUpen = %d (%d, %d)\n", i,f[i].pair, AUpen, sequence[i], sequence[f[i].pair]);	// Cristina: July 17, 2007

	for (h=0; h < f[i].num_branches; h++) {
		AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);
		//		printf("MULTI: %d AUpen = %d\n", f[i].bri[h], AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]));	// Cristina: July 17, 2007
	}

	// NOTE: Cristina: This uses a variation of dangling_energy_ simfold functions that takes two extra
	//       parameters indicating whether or not the dangling ends should be included in the calculation.
	//       This restriction is set in partial_coaxial_energy, and prevents ends involved in coaxial
	//       stacking from being included again as dangles.

	// add dangling energies for multi-loop
	dang += dangling_energy_left_res (sequence, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair, Input->cannot_add_dangling[i+1], Input->cannot_add_dangling[f[i].bri[0]-1]); //d_12=i1+1; d_34=i3-1
	//	printf("MULTI: dang left = %d\n", dangling_energy_left_res (sequence, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair, Input->cannot_add_dangling[i+1], Input->cannot_add_dangling[f[i].bri[0]-1]));	// Cristina: July 17, 2007

	for (l=0; l < f[i].num_branches - 1; l++){
		dang += dangling_energy_res (sequence, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair, Input->cannot_add_dangling[f[f[i].bri[l]].pair+1], Input->cannot_add_dangling[f[i].bri[l+1]-1]); //d_12=i2+1; d_34=i3-1

		//		printf("MULTI: dang %d = %d: %d %d %d %d, %d %d %d %d\n", l, dangling_energy_res (sequence, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair, Input->cannot_add_dangling[f[f[i].bri[l]].pair+1], Input->cannot_add_dangling[f[i].bri[l+1]-1]),
		//		 	f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair,
		//			sequence[f[i].bri[l]], sequence[f[f[i].bri[l]].pair], sequence[f[i].bri[l+1]], sequence[f[f[i].bri[l+1]].pair]);	// Cristina: July 17, 2007
	}
	dang += dangling_energy_right_res (sequence, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair, Input->cannot_add_dangling[f[i].pair-1], Input->cannot_add_dangling[f[f[i].bri[f[i].num_branches-1]].pair+1]); //d_12=i2-1; d_34=i4+1
	//	printf("MULTI: dang right = %d\n", dangling_energy_right_res (sequence, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair, Input->cannot_add_dangling[f[i].pair-1], Input->cannot_add_dangling[f[f[i].bri[f[i].num_branches-1]].pair+1]));	// Cristina: July 17, 2007

	// add "no-dangling" restriction                                    
	for (l=0; l < f[i].num_branches; l++){
		Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
		Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
	}
	Input->cannot_add_dangling [i] = 1;  // restrict start of multiloop from having an AU penalty added

	return misc_energy + dang + AUpen;
}


/*********************************************************************************
pseudoEnergyCC2006b: it calculates the free energy of a pseudoknotted loop.
It also calculates the free energy of interior-pseudoknotted and multi-pseudoknotted
loops regarding the pseudoknotted loop and add them to the energy
of pseudoknotted loop. Energy value is returned in 10*cal/mol.
 *********************************************************************************/

float	Loop::pseudoEnergyCC2006b(){

	double Energy = 0;
	int * sequence = Input->type;

	int numbases = end - begin + 1;
	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	char structure_all[numbases+1];
	for (int i = 0; i < numbases; i++){
		if (Input->Sequence[begin + i] <= 0)
			structure_all[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure_all[i] = '(';
		else
			structure_all[i] = ')';
	}
	structure_all[numbases] = '\0';

	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	int k = 0;
	int i = begin;

	int internalLoop = 0;  // does this type of pseudoloop contain an internal loop? (handled by CC2006: 1 = yes; 0 = no)
	int multi = 0;  // does this type of pseudoloop contain a multi loop? (not handled by CC2006: 1 = yes; 0 = no)
	int validBasepair = Input->BasePair(i);  // first valid base pair in first band is "b" above i.e. BasePair(a)
	int valid = 0;  // is this free of kissing hairpins or chains of bands? (not handled by CC2006: 1 = yes; 0 = no)

	if (NumberOfBands == 2)
		valid = 1;

	if (valid == 0)
		return pseudoEnergyDP();

	// check if it's an H-type pseudoloop
	// for [a,a']
	for (k = i; k < bandpattern->pattern[i].OtherBorder; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			}else if (Input->BasePair(k) < bandpattern->pattern[i].OtherBorder){ // if base pair is less than a' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}else{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}
	// for [c,c']
	validBasepair = Input->BasePair(bandpattern->pattern[i].next);  // first valid base pair in second band is "d" above i.e. BasePair(c)
	for (k = bandpattern->pattern[i].next; k < bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			}else if (Input->BasePair(k) < bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder){ // if base pair is less than c' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}else{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}
	// for [b',b]  // added: August 11 - to check if there are nested children on the other side of the band
	int border_bp = bandpattern->pattern[bandpattern->pattern[i].next].next;
	int border_b = bandpattern->pattern[border_bp].OtherBorder;
	int border_dp = bandpattern->pattern[border_bp].next;
	int border_d = bandpattern->pattern[border_dp].OtherBorder;

	validBasepair = Input->BasePair(border_bp);  // first valid base pair in second band is "ap" above i.e. BasePair(bp)
	for (k = border_bp; k < border_b; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			}else if (Input->BasePair(k) > border_bp){ // if base pair is less than c' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}else{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}
	// for [d',d]  // added: August 11 - to check if there are nested children on the other side of the band
	validBasepair = Input->BasePair(border_dp);  // first valid base pair in second band is "c" above i.e. BasePair(d)
	for (k = border_dp; k < border_d; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			} else if (Input->BasePair(k) > border_dp) { // if base pair is less than c' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}else{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}

	int stem1 = 0;  // size of S1 (first stem from 5' end)
	int stem2 = 0;  // size of S2
	int loop1 = 0;  // size of L1 (first loop from 5' end)
	int loop2 = 0;  // size of L2

	if (valid == 1 && multi == 0 && internalLoop == 0){  // H-type pseudoknot: Use Cao&Chen model directly
		// i.e. a'-a+1
		stem1 = bandpattern->pattern[i].OtherBorder - i + 1;
		// i.e. c'-c+1
		stem2 = bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder - bandpattern->pattern[i].next + 1;
		// i.e. c-a'+1
		loop1 = bandpattern->pattern[i].next - bandpattern->pattern[i].OtherBorder - 1;
		// i.e. d'-b+1 
		loop2 = Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - Input->BasePair(i) - 1;

	}else if (valid == 1 && multi == 0 && internalLoop == 1){  // internal loop inside stem: User Cao&Chen model with stem length = number of basepairs
		stem1 = bandpattern->pattern[i].OtherBorder - i + 1 - UnpairedBases(Input, i, bandpattern->pattern[i].OtherBorder);  // i.e. a'-a+1
		stem2 = bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder - bandpattern->pattern[i].next + 1 - UnpairedBases(Input, bandpattern->pattern[i].next, bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder);  // i.e. c'-c+1
		loop1 = bandpattern->pattern[i].next - bandpattern->pattern[i].OtherBorder - 1; // i.e. c-a'+1
		loop2 = Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - Input->BasePair(i) - 1; // i.e. d'-b+1 

	}else{ // Multiloop or chain exists: Use Dirks&Pierce instead
		return pseudoEnergyDP();
	}

	// added: August 11 - conditions on loop1 and loop2
	if (stem1 > 12 || stem2 > 12 || stem1 <= 1 || stem2 <= 1 || loop1 == 0 || loop2 == 0){  // not handled by Cao & Chen tables
		return pseudoEnergyDP();
	}

	float TdeltaS_L1 = 0;
	float TdeltaS_L2 = 0;
	float ln_omega_coil_L1 = 0;
	float ln_omega_L1 = 0;
	float ln_omega_coil_L2 = 0;
	float ln_omega_L2 = 0;


	// Not in the CC tables, but possible pseudoknot conformations:
	if (loop1 <= 12 && cc2006_s2_l1[stem2-2][loop1-1] >= INF){
		return pseudoEnergyDP();  // not possible using CC, but may occur in real life
	}

	if (loop2 <= 12 && cc2006_s1_l2[stem1-2][loop2-1] >= INF){
		return pseudoEnergyDP();  // not possible using CC, but may occur in real life
	}


	if (loop1 > 12){
		ln_omega_coil_L1 = 2.14 * loop1 + 0.10;
		// a * ln(loop1 - l_min + 1) + b * (loop1 - l_min + 1) + c;
		ln_omega_L1 = cc2006_s2_formula[1][stem2-2] * log(loop1 - cc2006_s2_formula[0][stem2-2] + 1) + cc2006_s2_formula[2][stem2-2] * (loop1 - cc2006_s2_formula[0][stem2-2] + 1) + cc2006_s2_formula[3][stem2-2];
		// all other energies are in 10cal/mol
		TdeltaS_L1 = -KB * (ln_omega_coil_L1 - ln_omega_L1) * pkmodelCC2006.temp * 100; 
	}else{
		TdeltaS_L1 = -cc2006_s2_l1[stem2-2][loop1-1];

		// Cristina: added this since cases can occur where a pseudoknot of dimensions not
		// included in the CC table do occur in real life
		if (cc2006_s2_l1[stem2-2][loop1-1] >= INF)
			return pseudoEnergyDP();  // not possible using CC, but may occur in real life
	}

	if (loop2 > 12){
		ln_omega_coil_L2 = 2.14 * loop2 + 0.10;
		ln_omega_L2 = cc2006_s1_formula[1][stem1-2] * log(loop2 - cc2006_s1_formula[0][stem1-2] + 1) + cc2006_s1_formula[2][stem1-2] * (loop2 - cc2006_s1_formula[0][stem1-2] + 1) + cc2006_s1_formula[3][stem1-2];  // a * ln(loop1 - l_min + 1) + b * (loop1 - l_min + 1) + c;
		TdeltaS_L2 = -KB * (ln_omega_coil_L2 - ln_omega_L2) * pkmodelCC2006.temp * 100;  // all other energies are in 10cal/mol
	}else{
		TdeltaS_L2 = -cc2006_s1_l2[stem1-2][loop2-1];

		// Cristina: added this since cases can occur where a pseudoknot of dimensions not 
		// included in the CC table do occur in real life
		if (cc2006_s2_l1[stem2-2][loop1-1] >= INF)
			return pseudoEnergyDP();  // not possible using CC, but may occur in real life
	}

	// energy in 10cal/mol
	Energy = (pkmodelCC2006.deltaG_assemble	// pseudoknot initiation penalty
			- TdeltaS_L1	// energy of stem 2, loop 1
			- TdeltaS_L2);  // energy of stem 1, loop 2

	// Add energy values for stacked pairs in the stem below (similar to D&P process):

	// CHECK: change the rest of this function to get the energy of the stacked pairs in the H-type pseudoloop

	T_IntList * L1 = ILoops;

	// Calculating the free energy of internal loops that span a band (interior-pseudoknotted loops)

	while (L1 != NUL){
		float en = Input->looplists[L1->Num]->interiorPseudoEnergyCC2006a();  // same function for Cao&Chen a and b
		Energy += en;
		L1 = L1->Next;
	};

	// Add coaxial stacking term...
	int CoaxialStacking = 0;
	int stack_dangle = 0;
	int dangle0 = -1;
	// use dangle only if allowed (i.e. not involved in coaxial stacking elsewhere)
	if (Input->cannot_add_dangling[bandpattern->pattern[i].OtherBorder + 1] != 1)
		// -1 if already paired
		dangle0 = (Input->BasePair(bandpattern->pattern[i].OtherBorder + 1) > 0) ? -1 : bandpattern->pattern[i].OtherBorder + 1;
	int dangle1 = -1;
	// use dangle only if allowed (i.e. not involved in coaxial stacking elsewhere)
	if (Input->cannot_add_dangling[Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1] != 1) 
		// -1 if already paired
		dangle1 = (Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1) > 0) ? -1 : Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1; 
	//... as long as c'+1=b'
	if ( (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder + 1) == Input->BasePair(bandpattern->pattern[i].OtherBorder) ){
		// coaxial stacking energy of (c', bp(c'), b', bp(b')=a')
		CoaxialStacking = LEcoax_stack_energy_flush_b (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, dangle0, dangle1, 
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder)-2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				sequence, no_pk_dangling_ends);

		if (CoaxialStacking != 0){  // coaxial stacking favourable over dangling ends
			// make sure dangling ends are not included elsewhere (since include coaxial stacking)
			if (dangle0 != -1)  // only change array if the dangling end was free originally
				Input->cannot_add_dangling[dangle0] = 1;
			if (dangle1 != -1)  // only change array if the dangling end was free originally
				Input->cannot_add_dangling[dangle1] = 1;
		}
	//... or as long as c'+2 = b' (mismatch coaxial stacking)
	} else if ( (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder + 2) == Input->BasePair(bandpattern->pattern[i].OtherBorder) ){
		// coaxial stacking energy of (c', bp(c'), b', bp(b')=a', dangle_a'=a'+1, dangle_bp(c')=bp(c')-1)
		CoaxialStacking = LEcoax_stack_energy_mismatch (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, 
				dangle0, dangle1, sequence, structure_all,
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				stack_dangle, no_pk_dangling_ends);

		if (CoaxialStacking != 0)  // coaxial stacking favourable over dangling ends
		{
			if (RES_STACK_DANGLE == 1)
			{
				// make sure the dangling end not involved in stacking (the opposite from stack_dangle) is not included elsewhere
				if (stack_dangle == 0 && dangle0 != -1)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[dangle0] = 1;
				if (stack_dangle == 1 && dangle1 != -1)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[dangle1] = 1;
				//if (stack_dangle != 0 || stack_dangle != 1)
				printf("ERROR: stack_dangle = %d not set properly for mismatch pseudoloop; must be 0 or 1\n", stack_dangle); 
			}else{
				Input->cannot_add_dangling[dangle0] = 1;
				Input->cannot_add_dangling[dangle1] = 1;
			}
		}
	}

	Energy += CoaxialStacking;

	return Energy;

}


/*********************************************************************************
pseudoEnergyCC2006b: FOR PARAMETER TUNING. it calculates the free energy of a pseudoknotted loop.
It also calculates the free energy of interior-pseudoknotted and multi-pseudoknotted
loops regarding the pseudoknotted loop and add them to the energy
of pseudoknotted loop. Energy value is returned in kcal/mol.  Assumes temp is 37.0 degree Celsius.
 *********************************************************************************/
float	Loop::pseudoEnergyCC2006b(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles)
{

	// TODO: globalize this?
	float temp = 37.0 + 273.15;  // TODO

	int num_params_pkfree = get_num_params();
	int num_params = get_num_params_PK_CC2006b();
	int num_params_DP_and_pkfree = get_num_params_PK_DP();
	if (reset_c == 1 && c != NUL){
		f = 0;         
		for (int i=0; i < num_params; i++){
			c[i] = 0;
			for (int j=i; j < num_params; j++) 
				P_matrix[i][j] = 0;
		}
	}

	double Energy = 0;
	int * sequence = Input->type;

	int numbases = end - begin + 1;
	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	char structure_all[numbases+1];
	for (int i = 0; i < numbases; i++){
		if (Input->Sequence[begin + i] <= 0)
			structure_all[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure_all[i] = '(';
		else
			structure_all[i] = ')';
	}
	structure_all[numbases] = '\0';


	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	int k = 0;
	int i = begin;

	int internalLoop = 0;  // does this type of pseudoloop contain an internal loop? (handled by CC2006: 1 = yes; 0 = no)
	int multi = 0;  // does this type of pseudoloop contain a multi loop? (not handled by CC2006: 1 = yes; 0 = no)
	int validBasepair = Input->BasePair(i);  // first valid base pair in first band is "b" above i.e. BasePair(a)
	int valid = 0;  // is this free of kissing hairpins or chains of bands? (not handled by CC2006: 1 = yes; 0 = no)

	if (NumberOfBands == 2)
		valid = 1;

	if (valid == 0)
		return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);

	// check if it's an H-type pseudoloop
	// for [a,a']
	for (k = i; k < bandpattern->pattern[i].OtherBorder; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			} else if (Input->BasePair(k) < bandpattern->pattern[i].OtherBorder){ // if base pair is less than a' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			} else {
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}
	// for [c,c']
	validBasepair = Input->BasePair(bandpattern->pattern[i].next);  // first valid base pair in second band is "d" above i.e. BasePair(c)
	for (k = bandpattern->pattern[i].next; k < bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			} else if (Input->BasePair(k) < bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder){ // if base pair is less than c' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			} else {
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}
	// for [b',b]  // added: August 11 - to check if there are nested children on the other side of the band
	int border_bp = bandpattern->pattern[bandpattern->pattern[i].next].next;
	int border_b = bandpattern->pattern[border_bp].OtherBorder;
	int border_dp = bandpattern->pattern[border_bp].next;
	int border_d = bandpattern->pattern[border_dp].OtherBorder;

	validBasepair = Input->BasePair(border_bp);  // first valid base pair in second band is "ap" above i.e. BasePair(bp)
	for (k = border_bp; k < border_b; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			} else if (Input->BasePair(k) > border_bp) { // if base pair is less than c' above
				multi = 1;
				break;  // break since not handled by Cao & Chen
			} else {
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}
	// for [d',d]  // added: August 11 - to check if there are nested children on the other side of the band
	validBasepair = Input->BasePair(border_dp);  // first valid base pair in second band is "c" above i.e. BasePair(d)
	for (k = border_dp; k < border_d; k++){
		if (Input->BasePair(k) != validBasepair){
			if (Input->BasePair(k) <= 0) {
				internalLoop = 1;
			}
			else if (Input->BasePair(k) > border_dp) // if base pair is less than c' above
			{
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}
			else
			{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}

	int stem1 = 0;  // size of S1 (first stem from 5' end)
	int stem2 = 0;  // size of S2
	int loop1 = 0;  // size of L1 (first loop from 5' end)
	int loop2 = 0;  // size of L2

	if (valid == 1 && multi == 0 && internalLoop == 0)  // H-type pseudoknot: Use Cao&Chen model directly
	{
		stem1 = bandpattern->pattern[i].OtherBorder - i + 1;  // i.e. a'-a+1
		stem2 = bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder - bandpattern->pattern[i].next + 1;  // i.e. c'-c+1
		loop1 = bandpattern->pattern[i].next - bandpattern->pattern[i].OtherBorder - 1; // i.e. c-a'+1
		loop2 = Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - Input->BasePair(i) - 1; // i.e. d'-b+1 

	}
	else if (valid == 1 && multi == 0 && internalLoop == 1)  // internal loop inside stem: User Cao&Chen model with stem length = number of basepairs
	{
		stem1 = bandpattern->pattern[i].OtherBorder - i + 1 - UnpairedBases(Input, i, bandpattern->pattern[i].OtherBorder);  // i.e. a'-a+1
		stem2 = bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder - bandpattern->pattern[i].next + 1 - UnpairedBases(Input, bandpattern->pattern[i].next, bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder);  // i.e. c'-c+1
		loop1 = bandpattern->pattern[i].next - bandpattern->pattern[i].OtherBorder - 1; // i.e. c-a'+1
		loop2 = Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - Input->BasePair(i) - 1; // i.e. d'-b+1 

	}
	else // Multiloop or chain exists: Use Dirks&Pierce instead
	{
		return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
	}

	// added: August 11 - conditions on loop1 and loop2
	if (stem1 > 12 || stem2 > 12 || stem1 <= 1 || stem2 <= 1 || loop1 == 0 || loop2 == 0)  // not handled by Cao & Chen tables
	{
		return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);
	}

	float TdeltaS_L1 = 0;
	float TdeltaS_L2 = 0;
	float ln_omega_coil_L1 = 0;
	float ln_omega_L1 = 0;
	float ln_omega_coil_L2 = 0;
	float ln_omega_L2 = 0;


	char paramtype[100];

	// Not in the CC tables, but possible pseudoknot conformations:
	if (loop1 <= 12 && cc2006_s2_l1[stem2-2][loop1-1] >= INF)
	{
		return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);  // not possible using CC, but may occur in real life
	}

	if (loop2 <= 12 && cc2006_s1_l2[stem1-2][loop2-1] >= INF)
	{
		return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);  // not possible using CC, but may occur in real life
	}


	if (loop1 > 12)
	{
		ln_omega_coil_L1 = 2.14 * loop1 + 0.10;
		ln_omega_L1 = cc2006_s2_formula[1][stem2-2] * log(loop1 - cc2006_s2_formula[0][stem2-2] + 1) + cc2006_s2_formula[2][stem2-2] * (loop1 - cc2006_s2_formula[0][stem2-2] + 1) + cc2006_s2_formula[3][stem2-2];  // a * ln(loop1 - l_min + 1) + b * (loop1 - l_min + 1) + c;
		TdeltaS_L1 = -KB * (ln_omega_coil_L1 - ln_omega_L1) * pkmodelCC2006.temp * 100;  // all other energies are in 10cal/mol
		//		printf("ln_omega_coil_L1 = %f, ln_omega_L1 = %f, cc2006_s2_formula[1][stem2-2] = %f, log = %f", ln_omega_coil_L1, ln_omega_L1, cc2006_s2_formula[1][stem2-2], log(loop1 - cc2006_s2_formula[0][stem2-2] + 1));

		// TODO
		// NOTE: the "-=" comes from below where we do Energy += -TdeltaS
		// update counter for cc2006_s2_formula[1][stem2-2], cc2006_s2_formula[2][stem2-2], and cc2006_s2_formula[3][stem2-2]
		sprintf(paramtype, "cc2006_s2_formula[%d][%d]", 1, stem2-2);
		if (c!= NUL)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= log(loop1 - cc2006_s2_formula[0][stem2-2] + 1)*KB*pkmodelCC2006.temp;
		sprintf(paramtype, "cc2006_s2_formula[%d][%d]", 2, stem2-2);
		if (c!= NUL)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= (loop1 - cc2006_s2_formula[0][stem2-2] + 1)*KB*pkmodelCC2006.temp;
		sprintf(paramtype, "cc2006_s2_formula[%d][%d]", 3, stem2-2);
		if (c!= NUL)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= KB*pkmodelCC2006.temp;
		f -= -KB * ln_omega_coil_L1 * pkmodelCC2006.temp;
	}
	else
	{
		TdeltaS_L1 = -cc2006_s2_l1[stem2-2][loop1-1];

		if (cc2006_s2_l1[stem2-2][loop1-1] >= INF)
		{
			return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);  // not possible using CC, but may occur in real life
		}
		// TODO
		// NOTE: the "-=" comes from below where we do Energy += -TdeltaS
		// update counter for cc2006_s2_l1[stem2-2][loop1-1];
		sprintf(paramtype, "cc2006_s2_l1[%d][%d]", stem2-2, loop1-1);
		if (c!= NUL && cc2006_s2_l1[stem2-2][loop1-1] < INF)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= -1*KB*temp;
	}

	if (loop2 > 12)
	{
		ln_omega_coil_L2 = 2.14 * loop2 + 0.10;
		ln_omega_L2 = cc2006_s1_formula[1][stem1-2] * log(loop2 - cc2006_s1_formula[0][stem1-2] + 1) + cc2006_s1_formula[2][stem1-2] * (loop2 - cc2006_s1_formula[0][stem1-2] + 1) + cc2006_s1_formula[3][stem1-2];  // a * ln(loop1 - l_min + 1) + b * (loop1 - l_min + 1) + c;
		TdeltaS_L2 = -KB * (ln_omega_coil_L2 - ln_omega_L2) * pkmodelCC2006.temp * 100;  // all other energies are in 10cal/mol

		// TODO
		// NOTE: the "-=" comes from below where we do Energy += -TdeltaS
		// update counter for cc2006_s1_formula[1][stem1-2], cc2006_s1_formula[2][stem1-2], and cc2006_s1_formula[3][stem1-2]
		sprintf(paramtype, "cc2006_s1_formula[%d][%d]", 1, stem1-2);
		if (c!= NUL)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= log(loop2 - cc2006_s1_formula[0][stem1-2] + 1)*KB*pkmodelCC2006.temp;
		sprintf(paramtype, "cc2006_s1_formula[%d][%d]", 2, stem1-2);
		if (c!= NUL)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= (loop2 - cc2006_s1_formula[0][stem1-2] + 1)*KB*pkmodelCC2006.temp;
		sprintf(paramtype, "cc2006_s1_formula[%d][%d]", 3, stem1-2);
		if (c!= NUL)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= KB*pkmodelCC2006.temp;
		f -= -KB * ln_omega_coil_L2 * pkmodelCC2006.temp;
	}
	else
	{
		TdeltaS_L2 = -cc2006_s1_l2[stem1-2][loop2-1];

		if (cc2006_s1_l2[stem1-2][loop2-1] >= INF)
		{
			return pseudoEnergyDP(P_matrix, c, f, reset_c, ignore_dangles);  // not possible using CC, but may occur in real life
		}

		// TODO
		// NOTE: the "-=" comes from below where we do Energy += -TdeltaS
		// update counter for cc2006_s1_l2[stem1-2][loop2-1];
		sprintf(paramtype, "cc2006_s1_l2[%d][%d]", stem1-2, loop2-1);
		if (c!= NUL && cc2006_s1_l2[stem1-2][loop2-1] < INF)   c[num_params_DP_and_pkfree + structure_type_index_PK_CC(paramtype) - 1] -= -1*KB*temp;
	}

	// energy in 10cal/mol
	Energy += (pkmodelCC2006.deltaG_assemble	// pseudoknot initiation penalty
			- TdeltaS_L1	// energy of stem 2, loop 1
			- TdeltaS_L2);  // energy of stem 1, loop 2

	// update counts:
	f += pkmodelCC2006.deltaG_assemble/100.0;

	// Energy so far is in 10cal/mol. Change to kcal/mol, which matches the rest of this function
	Energy = Energy / 100.0;

	// Add energy values for stacked pairs in the stem below (similar to D&P process):

	// Calculating the free energy of internal loops that span a band (interior-pseudoknotted loops)
	// THIS IS DONE BELOW
	/*
	   while (L1 != NUL){
	   float en = Input->looplists[L1->Num]->interiorPseudoEnergyCC2006a();  // same function for Cao&Chen a and b
	   Energy += en;
	   L1 = L1->Next;
	   };
	 */

	T_IntList * L1 = ILoops;
	int a, ap, bp, b;
	int k_pt = 0;
	int i_pt = begin;
	numbases = 0;
	// create structure: a string of dot-brackets to represent the region
	char * structure;
	// create sequence
	char * csequence;

	float pkfree_retval = 0;

	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	// we want to send each band structure to simfold

	// NOTE: unlike in the similar loop for pseudoEnergyDP(), we want only NumberOfBands, since there
	//       are only two anyway. We need NumberOfBands*2 if there were more than 2 bands.
	for (k_pt = 0; k_pt < NumberOfBands; k_pt++){
		if (k_pt >= 2)
			printf("WARNING: Loop.cpp::pseudoEnergyCC2006b - k_pt is %d >= 2 which shouldn't happen for CC model!\n", k_pt);

		// only find energy if the band has stacked pairs / internal loops - otherwise, no point
		if (L1 != NUL){
			a = i_pt;
			ap = bandpattern->pattern[i_pt].OtherBorder;
			bp = Input->BasePair(ap);
			b = Input->BasePair(a);

			numbases = b-a+1;
			structure = new char[numbases+1];
			csequence = new char[numbases+1];
			for (int i = 0; i < numbases; i++){
				csequence[i] = Input->CSequence[a + i];
			}
			csequence[numbases] = '\0';
			structure[numbases] = '\0';

			// now get energy/counts for this band

			// this is ok since we expect there to be no pseudoknots in this band, since otherwise would have been passed
			// onto DP model instead (pseudoknots in band --> multiloop in band --> not handled by CC)
			fillPseudoStructureNoMulti(csequence, structure, numbases, a,ap,bp,b);

			if (c != NUL){
				pkfree_retval = get_feature_counts_restricted(csequence, structure, c, f, 0, ignore_dangles, 1);
			}else{
				pkfree_retval = get_feature_counts_restricted(csequence, structure, NUL, f, 0, ignore_dangles, 1);
			}


			Energy += pkfree_retval;

			// clear structure and csequence for the next spanning loop (which possibly needs these to be different length)
			delete [] structure;
			delete [] csequence;
		}
		i_pt = bandpattern->pattern[i_pt].next;
	}


	pk_str_features * feat = Input->loops;
	int AUpen = 0;
	i_pt = begin;
	k_pt = 0;


	// add the appropriate AU penalties
	for (k_pt = 0; k_pt < NumberOfBands*2; k_pt++)  // look at all the bands, as listed above
	{
		if (Input->BasePair(i_pt) > i_pt)  // if it is the [a,ap] portion of the band
		{
			// if the outermost loop spanning the band is not a multiloop, add AU
			// the cdt is unnecessary, based on old hotknots code (energydangling function) and simfold energy of a multiloop

			AUpen += AU_penalty (sequence[i_pt], sequence[feat[i_pt].pair]);
			if (c != NUL)    count_AU_penalty (sequence[i_pt], sequence[feat[i_pt].pair], c);
		}

		i_pt = bandpattern->pattern[i_pt].next;  // next band
	}
	Energy += float(AUpen)/100;  // since AUpen is in 10cal/mol, and everything else is in kcal/mol


	// TODO - get rid of dangles added by simfold
	// add counts for dangles and coax stacking terms

	// Add coaxial stacking term...
	int CoaxialStacking = 0;
	int stack_dangle = 0;
	int dangle0 = -1;
	if (Input->cannot_add_dangling[bandpattern->pattern[i].OtherBorder + 1] != 1)  // use dangle only if allowed (i.e. not involved in coaxial stacking elsewhere)
		dangle0 = (Input->BasePair(bandpattern->pattern[i].OtherBorder + 1) > 0) ? -1 : bandpattern->pattern[i].OtherBorder + 1;  // -1 if already paired
	int dangle1 = -1;
	if (Input->cannot_add_dangling[Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1] != 1)  // use dangle only if allowed (i.e. not involved in coaxial stacking elsewhere)
		dangle1 = (Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1) > 0) ? -1 : Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1;  // -1 if already paired

	//... as long as c'+1=b'
	if ( (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder + 1) == Input->BasePair(bandpattern->pattern[i].OtherBorder) ){
		// coaxial stacking energy of (c', bp(c'), b', bp(b')=a')
		CoaxialStacking = LEcoax_stack_energy_flush_b (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, dangle0, dangle1, 
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder)-2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				sequence, ignore_dangles);

		if (CoaxialStacking != 0)  // coaxial stacking favourable over dangling ends
		{
			// make sure dangling ends are not included elsewhere (since include coaxial stacking)
			if (dangle0 != -1)  // only change array if the dangling end was free originally
				Input->cannot_add_dangling[dangle0] = 1;
			if (dangle1 != -1)  // only change array if the dangling end was free originally
				Input->cannot_add_dangling[dangle1] = 1;
		}

		// TODO update counts for dangling ends
		// NOTE: dangling ends included here will not conflict with those from simfold function
		// If no coaxial stacking is used, dangling end energies will be added in EnergyDangling()
		// get the counts for the coaxial stacking parameters
		if (c!= NUL)   count_LEcoax_stack_energy_flush_b(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, dangle0, dangle1, 
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder)-2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				sequence, c, ignore_dangles);

	//... or as long as c'+2 = b' (mismatch coaxial stacking)
	} else if ( (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder + 2) == Input->BasePair(bandpattern->pattern[i].OtherBorder) ) {
		// coaxial stacking energy of (c', bp(c'), b', bp(b')=a', dangle_a'=a'+1, dangle_bp(c')=bp(c')-1)
		CoaxialStacking = LEcoax_stack_energy_mismatch (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, 
				dangle0, dangle1, sequence, structure_all, 
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				stack_dangle, ignore_dangles);

		if (CoaxialStacking != 0){  // coaxial stacking favourable over dangling ends
			if (RES_STACK_DANGLE == 1){
				// make sure the dangling end not involved in stacking (the opposite from stack_dangle) is not included elsewhere
				if (stack_dangle == 0 && dangle0 != -1)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[dangle0] = 1;
				if (stack_dangle == 1 && dangle1 != -1)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[dangle1] = 1;
				//if (stack_dangle != 0 || stack_dangle != 1)
				printf("ERROR: stack_dangle = %d not set properly for mismatch pseudoloop; must be 0 or 1\n", stack_dangle); 
			}else{
				Input->cannot_add_dangling[dangle0] = 1;
				Input->cannot_add_dangling[dangle1] = 1;
			}
		}

		// TODO update counts for dangling ends
		// NOTE: dangling ends included here will not conflict with those from simfold function
		// If no coaxial stacking is used, dangling end energies will be added in EnergyDangling()
		// get the counts for the coaxial stacking parameters
		if (c!= NUL)   count_LEcoax_stack_energy_mismatch (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, 
				dangle0, dangle1, sequence, structure_all, 
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				stack_dangle, c, ignore_dangles);

	}


	Energy += float(CoaxialStacking)/100.0;

	return Energy;

}


/*********************************************************************************
pseudoEnergyCC2006a: it calculates the free energy of a pseudoknotted loop.
It also calculates the free energy of interior-pseudoknotted and multi-pseudoknotted
loops regarding the pseudoknotted loop and add them to the energy
of pseudoknotted loop. Energy value is returned in 10*cal/mol.
 *********************************************************************************/

float	Loop::pseudoEnergyCC2006a(){


	double Energy = 0;
	int * sequence = Input->type;

	int numbases = end - begin + 1;
	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	char structure_all[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure_all[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure_all[i] = '(';
		else
			structure_all[i] = ')';
	}
	structure_all[numbases] = '\0';

	// Given pseudoloop defined by 2 bands: [a,a']|_|[b',b] and [c,c']|_|[d',d]
	int k = 0;
	int i = begin;

	int internalLoop = 0;  // does this type of pseudoloop contain an internal loop? (handled by CC2006: 1 = yes; 0 = no)
	int multi = 0;  // does this type of pseudoloop contain a multi loop? (not handled by CC2006: 1 = yes; 0 = no)
	int validBasepair = Input->BasePair(i);  // first valid base pair in first band is "b" above i.e. BasePair(a)
	int valid = 0;  // is this free of kissing hairpins or chains of bands? (not handled by CC2006: 1 = yes; 0 = no)

	if (NumberOfBands == 2)
		valid = 1;

	if (valid == 0)
		return pseudoEnergyDP();

	// check if there are only stems in the pseudoloop: it's an H-type pseudoloop
	for (k = i; k < bandpattern->pattern[i].OtherBorder; k++)
	{
		if (Input->BasePair(k) != validBasepair)
		{
			if (Input->BasePair(k) < bandpattern->pattern[i].OtherBorder) // if base pair is less than a' above
			{
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}
			else
			{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the first band are stacked together
	}

	validBasepair = Input->BasePair(bandpattern->pattern[i].next);  // first valid base pair in second band is "d" above i.e. BasePair(c)
	for (k = bandpattern->pattern[i].next; k < bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder; k++)
	{
		if (Input->BasePair(k) != validBasepair)
		{
			if (Input->BasePair(k) < bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) // if base pair is less than c' above
			{
				multi = 1;
				break;  // break since not handled by Cao & Chen
			}
			else
			{
				internalLoop = 1;  // don't break in case there is a multiloop left
			}
		}
		validBasepair--;  // make sure all basepairs in the second band are stacked together
	}

	int stem1 = 0;  // size of S1 (first stem from 5' end)
	int stem2 = 0;  // size of S2
	int loop1 = 0;  // size of L1 (first loop from 5' end)
	int loop2 = 0;  // size of L2

	if (valid == 1 && multi == 0 && internalLoop == 0)  // H-type pseudoknot: Use Cao&Chen model directly
	{
		stem1 = bandpattern->pattern[i].OtherBorder - i + 1;  // i.e. a'-a+1
		stem2 = bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder - bandpattern->pattern[i].next + 1;  // i.e. c'-c+1
		loop1 = bandpattern->pattern[i].next - bandpattern->pattern[i].OtherBorder - 1; // i.e. c-a'+1
		loop2 = Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - Input->BasePair(i) - 1; // i.e. d'-b+1 

	}	
	else // Multiloop or internal loop or chain exists: Use Dirks&Pierce instead
	{
		return pseudoEnergyDP();
	}

	if (stem1 > 12 || stem2 > 12 || stem1 <= 1 || stem2 <= 1)  // not handled by Cao & Chen tables
	{
		return pseudoEnergyDP();
	}

	float TdeltaS_L1 = 0;
	float TdeltaS_L2 = 0;
	float ln_omega_coil_L1 = 0;
	float ln_omega_L1 = 0;
	float ln_omega_coil_L2 = 0;
	float ln_omega_L2 = 0;

	if (loop1 > 12)
	{
		ln_omega_coil_L1 = 2.14 * loop1 + 0.10;
		ln_omega_L1 = cc2006_s2_formula[1][stem2-2] * log(loop1 - cc2006_s2_formula[0][stem2-2] + 1) + cc2006_s2_formula[2][stem2-2] * (loop1 - cc2006_s2_formula[0][stem2-2] + 1) + cc2006_s2_formula[3][stem2-2];  // a * ln(loop1 - l_min + 1) + b * (loop1 - l_min + 1) + c;
		TdeltaS_L1 = -KB * (ln_omega_coil_L1 - ln_omega_L1) * pkmodelCC2006.temp * 100;  // all other energies are in 10cal/mol
		//		printf("ln_omega_coil_L1 = %f, ln_omega_L1 = %f, cc2006_s2_formula[1][stem2-2] = %f, log = %f", ln_omega_coil_L1, ln_omega_L1, cc2006_s2_formula[1][stem2-2], log(loop1 - cc2006_s2_formula[0][stem2-2] + 1));
	}
	else
	{
		TdeltaS_L1 = -cc2006_s2_l1[stem2-2][loop1-1];
	}

	if (loop2 > 12)
	{
		ln_omega_coil_L2 = 2.14 * loop2 + 0.10;
		ln_omega_L2 = cc2006_s1_formula[1][stem1-2] * log(loop2 - cc2006_s1_formula[0][stem1-2] + 1) + cc2006_s1_formula[2][stem1-2] * (loop2 - cc2006_s1_formula[0][stem1-2] + 1) + cc2006_s1_formula[3][stem1-2];  // a * ln(loop1 - l_min + 1) + b * (loop1 - l_min + 1) + c;
		TdeltaS_L2 = -KB * (ln_omega_coil_L2 - ln_omega_L2) * pkmodelCC2006.temp * 100;  // all other energies are in 10cal/mol
	}
	else
	{
		TdeltaS_L2 = -cc2006_s1_l2[stem1-2][loop2-1];
	}

	// energy in 10cal/mol
	Energy = (pkmodelCC2006.deltaG_assemble	// pseudoknot initiation penalty
			- TdeltaS_L1	// energy of stem 2, loop 1
			- TdeltaS_L2);  // energy of stem 1, loop 2

	// Add energy values for stacked pairs in the stem below (similar to D&P process):

	// CHECK: change the rest of this function to get the energy of the stacked pairs in the H-type pseudoloop

	T_IntList * L1 = ILoops;


	// Calculating the free energy of internal loops that span a band (interior-pseudoknotted loops)


	while (L1 != NUL){
		float en = Input->looplists[L1->Num]->interiorPseudoEnergyCC2006a();  // same function for Cao&Chen a and b
		Energy += en;
		L1 = L1->Next;
	};


	// Add coaxial stacking term whether of not pseudoknot contains a multi or internal loop...
	int CoaxialStacking = 0;
	int stack_dangle = 0;
	int dangle0 = -1;
	if (Input->cannot_add_dangling[bandpattern->pattern[i].OtherBorder + 1] != 1)  // use dangle only if allowed (i.e. not involved in coaxial stacking elsewhere)
		dangle0 = (Input->BasePair(bandpattern->pattern[i].OtherBorder + 1) > 0) ? -1 : bandpattern->pattern[i].OtherBorder + 1;  // -1 if already paired
	int dangle1 = -1;
	if (Input->cannot_add_dangling[Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1] != 1)  // use dangle only if allowed (i.e. not involved in coaxial stacking elsewhere)
		dangle1 = (Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1) > 0) ? -1 : Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 1;  // -1 if already paired

	//... as long as c'+1=b'
	if ( (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder + 1) == Input->BasePair(bandpattern->pattern[i].OtherBorder) )
	{
		// coaxial stacking energy of (c', bp(c'), b', bp(b')=a')
		CoaxialStacking = LEcoax_stack_energy_flush_a (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, sequence);
	}
	//... or as long as c'+2 = b' (mismatch coaxial stacking)
	else if ( (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder + 2) == Input->BasePair(bandpattern->pattern[i].OtherBorder) )
	{
		// coaxial stacking energy of (c', bp(c'), b', bp(b')=a', dangle_a'=a'+1, dangle_bp(c')=bp(c')-1)
		CoaxialStacking = LEcoax_stack_energy_mismatch (bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder, 
				Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder),
				Input->BasePair(bandpattern->pattern[i].OtherBorder),
				bandpattern->pattern[i].OtherBorder, COAX_PSEUDO, 
				dangle0, dangle1, sequence, structure_all,
				Input->BasePair(Input->BasePair(bandpattern->pattern[bandpattern->pattern[i].next].OtherBorder) - 2),
				Input->BasePair(bandpattern->pattern[i].OtherBorder+2),
				stack_dangle, no_pk_dangling_ends);

		if (CoaxialStacking != 0)  // coaxial stacking favourable over dangling ends
		{
			if (RES_STACK_DANGLE == 1)
			{
				// make sure the dangling end not involved in stacking (the opposite from stack_dangle) is not included elsewhere
				if (stack_dangle == 0 && dangle0 != -1)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[dangle0] = 1;
				if (stack_dangle == 1 && dangle1 != -1)  // only change array if the dangling end was free originally
					Input->cannot_add_dangling[dangle1] = 1;
				//if (stack_dangle != 0 || stack_dangle != 1)
				printf("ERROR: stack_dangle = %d not set properly for mismatch pseudoloop; must be 0 or 1\n", stack_dangle); 
			}
			else
			{
				Input->cannot_add_dangling[dangle0] = 1;
				Input->cannot_add_dangling[dangle1] = 1;
			}
		}
	}

	Energy += CoaxialStacking;

	return Energy;

}
