#ifndef INPUT_H
#define INPUT_H		

#include "Defines.h"
#include <ctype.h>

class Loop;
class LoopList;


class ReadInput {
private:

public:
	
//FUNCTIONS

    void AddPair(int m, int n, char c);     // add a base pair: pair[m] = n and type[m] = c
	void RemoveUnpaired();                  //Remove unpaired bases from the secondary structure

	ReadInput(char* fbpseq);  // ADDED: added this - not present in computeEnergy
	ReadInput(char* fseq, char* fbpseq);
	ReadInput(int size, char* baseSequence, short* pairRefSequence); 		
    ~ReadInput();

	int BasePair(int a);			//The base pair of an element. Is -1 if the element is not paired

//ATTRIBUTES

	LoopList *looplists[MaxN];      //For pseudoknotted loops, contains the list of interior-pseduoknotted
	                                //and multi-pseudoknotted nested in them!
	                                
    pk_str_features loops[MaxN];       //Loops in input structure where loops[i] is the loop starting at i;
    								//Contains some essential information about closed region-loops
                                    //such as closing pair, type of the loop ...
	
	int cannot_add_dangling[MaxN];  // 1 if can't add dangling (no-dangling restriction), 0 if can
//	int must_add_dangling[MaxN];    // 1 if must add dangling
	int Size;                       //Size of the RNA strand


	char CSequence[MaxN];           //RNA primary structure
	int  Sequence[MaxN];            //RNA secondary structure
	Loop * ClosedRegions[MaxN];     //ClosedRegion[i] is the pointer to the ClosedRegion/Loop starting at i
	
	int  type[MaxN];			    //whether it is  A, C, G, or T(U)
	
	int  Next[MaxN];		//Next paired Base
	int Prev[MaxN];		    //Previous paired Base

	void clearDanglingRestriction();  // clears the cannot_add_dangling array
};




#endif
