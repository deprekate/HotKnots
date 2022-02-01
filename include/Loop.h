//#include "Sequence.h"
#ifndef LOOP_H
#define LOOP_H

#include <string.h>

#include "Input.h"
#include "Bands.h"
#include "LoopList.h"

//#include "Defines.h"  // July 16 - put in .cpp file instead
//#include "common.h"  // July 16 - put in .cpp file instead
//#include "commonPK.h"  // July 16 - put in .cpp file instead

// Struct used in Coaxial Stacking
typedef struct pk_coax_features
{
	int index0;  // index of first loop in coaxial stacks
	int index1;  // index of last loop in coaxial stacks; if there is more than one pair, this might not refer to loop that stacks with index0
	Loop * loop0;  // pointer to first loop in the coaxial stacking; only set if this stacking involves one pair of loops
	Loop * loop1;  // pointer to first loop in the coaxial stacking; only set if this stacking involves one pair of loops
	int lastPairAdded;  // indexes the array pk_coax_features and points to the index of the last pair added in the coaxial stacking (between 0 and NumberOfChildren)
	int lastArrayIndex; // indexes the array pk_coax_features and points to the index of the combination of pairs without lastPairAdded
	int coax_energy;  // energy of the coaxial stacking of all combinations of pairs in this struct
	int stack_dangle;  // 0 if the dangle involved in stacking is associated with i.j and 1 if with ip.jp
					   // set only in basic pairs, -1 or undefined otherwise
} pk_coax_features;


class Loop {

public:

//FUNCTIONS


    Loop(int b, int e, ReadInput * R, Bands * B, Stack * S);				
	~Loop();

	//void clearLoop();  // clears the loop structure
	
    void addLoop(int begin, int end);		//Adds a loop [begin,end] to the tree
    void Print(int i);						//Prints the loop tree in a hierarchical structure

	void loopSetType();						
	void setTypes();						//Set the type of each Loop in the tree
											//(interior, hairin, pseudo, ...)
	void LoopFindBands();					//Finds the bands(band regions) of a pseudoknotted region
	
	void FindInnerLoops(int border1, int border2);		//Figuring out the span-Band loops regarding a 														
														//pseudoknotted closed region																		
	void PseudoNestedCheck();
	void WhereLocated();	//Find the locatin status of the children-loops of a pseudoknotted closed region

	//Free energy calculatin functions
	float Energy(int model);			                //Calculating the free energy of the secondary structure
	float EnergyDangling();
    float getEnergyDP();						//Computes the free energy of the loops
    float getEnergyRE();						//Computes the free energy of the loops
    float getEnergyCC2006a();					//Computes the free energy of the loops
    float getEnergyCC2006b();					//Computes the free energy of the loops
	float	pseudoEnergyDP();
	float	pseudoEnergyRE();
	float	pseudoEnergyCC2006a();
	float	pseudoEnergyCC2006b();
	PARAMTYPE	hairpinEnergy();
	PARAMTYPE	stackEnergy();
	PARAMTYPE	interiorEnergy();
	double multiEnergyDP();
	double multiEnergyRE();
	double multiEnergyCC2006a();
	
	// FOR PARAMETER TUNING
	int isPKFree();  // determines if a closed region is pseudoknot free (down to its lowest level structures)
					 // returns 1 if the loop has no nested pseudoknots (pk-free); 0 otherwise
	void fillMultiStructure(char * structure, char * csequence, int numbases, int l_begin);
	void fillPseudoStructure(char * csequence, char * structure, int len, int a, int ap, int bp, int b);

	float EnergyViaSimfold(int model);  // wrapper for Energy funcion below
	float Energy(int model, double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles); //Calculating the free energy of the secondary structure
	float EnergyDanglingViaSimfold(int model);  // wrapper for EnergyDangling function below
	float EnergyDangling(int model, double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles, int ignore_AU);
	void lookForPk(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles, float& sum);

	float getEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);						//Computes the free energy of the loops
	float pseudoEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);
	double nestedPseudoEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);
	double pkfreeEnergyDP(double** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);  // finds energy of a pseudoknot free region

	void fillPseudoStructureNoMulti(char * structure, char * csequence, int len, int a, int ap, int bp, int b);

	// for CC model
	float getEnergyCC2006b(double ** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);
	void lookForPkCC2006b(double ** P_matrix, double *c, double &f, int reset_c, int ignore_dangles, float& sum);
	double pkfreeEnergyCC2006b(double ** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);
	double nestedPseudoEnergyCC2006b(double ** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);
	float pseudoEnergyCC2006b(double ** P_matrix, double *c, double &f, int reset_c, int ignore_dangles);

	// for EnergyDangling
	int isAdjacentToPK(int pos, int a);
	int findClosedRegionNestedIn(int pos);
	int isNestedInBand(int pos, int& a);  // return 1 if base pos is nested in the band of this loop

	int startsMultiSpanningBand(int i_pt);  // return 1 if i_pt is the beginning of a multiloop that spans a band
	void setMultiPseudoLoopDangles(int a, int ap, int bp, int b);
	int isPKDangleInMultiPseudoLoop(int i);
	int hasPKBranches();


	// END PARAMETER TUNING FUNCTIONS

	void FindSpanBandLoops();  // determine the types of loops that span each band
	
	void setDotParanthStruct(char * structure);

//ATTRIBUTES	
	T_IntList *ILoops, *MLoops; //interior-pseudoknotted loops list and multi-pseudoknotted loops list
	Stack * St;

    PseudoNestedType nested;   //location status of a loop
	int NumberOfUnBandChild;
	int NumberOfBandChild;
	int CurrentBandRegion;
    int NumberOfBands;
    int NumberOfChildren;				
	int NumberOfUnpaird;		//Number of unpaired bases within the loop (NOT applicable to pseudoloops)
	int Pseudo_NumBranches;		//Number of Branches by considering pseudoknot children as 2.
	int NumberOfUnpairedInPseudo;  // Number of unpaired bases in a pseudoloop (is only set to true value after pseudoEnergy() call)
	int NumberOfUnpairedInUnbandChild;

	Loop* RightChild, *LeftSibling, *Parent; 
	int begin, end;				//The closed region borders
	LoopType type;
	
	ReadInput * Input;
	Bands  * bandpattern;	  //List of the band regions, sorted by the left border

//	SpanBand * spanBandLoops;  // List of loops that span the bands of closed region [begin,end]

	float getPartialCoaxialEnergy(int flag);
//	float getPartialCoaxialEnergyHelper(int flag);
	float getPartialCoaxialEnergyAll(int flag);
	void removeDangling(int i, LoopType loop_type, pk_coax_features oneCoaxStack, int stack_dangling);
	void countNumberOfChildren();

	void printEnergyTrace();
	float finalEnergy;
	float finalCoaxEnergy;

	char* DotParanthStructure;  // represents the loop tree (secondary structure) in dot-paranthesis format
};


#endif
