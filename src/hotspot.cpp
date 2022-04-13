/*****************************************************************
HotKnot: A heuristic algorithm for RNA secondary 
structure prediction including pseudoknots
File: hotspot.c
Description: 
This file is used to generate a tree of secondary structures, 
using HotKnot heuristic algorithm.

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



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "simfold.h"
#include "hotspot.h"
#include "utils.h"
#include "goodStem.h"
#include "sc.h"
#include "constantsPK.h"

short  *pair_table;
extern short  *S, *S1;
extern struct Hotspots* hotspots;   //global variable
extern int noPS;  /* do not output PS files for RNA structure*/
extern int TRACE;
extern float THR;
extern float T_RATIO;
extern struct Node* listOfNodes[50]; 
//an array that stores node pointers according to their scores. Highest score first.
extern int numRnaStruct;  //total number of different Rna structures
extern int count;
extern int MaxSubOpt; 

void usage(void);

long parentScore=0;

void InsertRna(struct Node* currentNode);
void InitHotspots(int MaxHotspots, int length);
void ClearHotspots(int MaxHotspots);
void ChildNodeConstraint(int length, struct Node *currentNode, struct Node *rootNode);
int FindHotspots(char* sequence, char* structure, int MaxHotspots);
int IsAlreadyExist(char *constraint, struct Node *rootNode, struct Node *currentNode);
int secondaryStruct(char *sequence, int length, struct Node *currentNode, 
struct Node *rootNode, int MaxHotspots, int energyModel);

void InitHotspots(int MaxHotspots, int length)      /*initialization, used before FindHotspots*/
{
	int i;
	hotspots = (struct Hotspots *)malloc(sizeof(struct Hotspots));
	hotspots->substruct = (struct Substruct *)malloc(sizeof(struct Substruct)*MaxHotspots);

	for (i = 0; i < MaxHotspots; i++) 
		hotspots->substruct[i].pairTable = (short *) malloc(sizeof(short)*2*length);
	hotspots->numHotspots = 0;
}

void ClearHotspots(int MaxHotspots) /*housekeeping*/
{
	int i;

	//if (DEBUGH)
	//	printf("In total, %d hotspots found. \n", hotspots->numHotspots);

	for (i = 0; i<MaxHotspots; i++) {
		free(hotspots->substruct[i].pairTable);
	}
	free(hotspots->substruct);
	free(hotspots);
}

int FindHotspots(char* sequence, char* structure, int MaxHotspots)
	/*char* structure points to the output of fold.c, 
	  note: it's the current output of fold.c, not the cumulative global structure
	  function returns the number of hotspots if it finds substructures with energy lower than THR, 
	  THR is a global variable, which is threshold
	  otherwise, it returns 0*/

{
	int i, j, p, q;
	int ee = 0;
	int k,mm;
	int length, startbase;
	int si, sj;
	int nh=0;
	int numHotspots = hotspots->numHotspots;
	int min_en;
	int replace;
	int num;
	int flag = 0;

	int len = strlen(structure); 

	//if (TRACE == 1) printf("%s\n%s\n",sequence, structure);

	pair_table = make_pair_table(structure);
	i = 1;
	while (i <= len) {
		ee = 0;
		while((i <= len) && (pair_table[i] == 0)) i++;

		while (i <= len && i > pair_table[i] && i <= (int)strlen(structure)) i++;   
		/*only consider pair to the right*/
		if (i >= (int)strlen(structure)) break;

		j = pair_table[i];
		si = i;
		sj = j;	    
		p = i;
		q = j;
		while(p<q) {
			p++;
			q--;
			if ( pair_table[p]!=(short)q)  { 
				/*p is not paired with q, check whether q is paired with p+1, if not, check whether
				  p is paired with q - 1, if not, check whether p+1 is paired with q-1*/
				if (pair_table[p+1] == (short)q ) p++;
				else if(pair_table[p+1] == (short)(q - 1)&& pair_table[p+2] == (short)(q-2) ) {
					p++;
					q--;
				}
				else if(pair_table[p] == (short)(q-1))  q--;
				else {
					p--;
					q++;         /*back to the closing pair*/
					break;  /*there is internal loop with size larger than 2 or there is multiloop*/   	
				}
			}

			ee += -F2(len, sequence, si,p,q,sj); 

			si = p;
			sj = q;
		}

		/*we find a structure  with closing pairs i.j and p.q*/
		//if (TRACE == 1) printf("substructure %d-%d, %d-%d with energy %d  THR  %f\n", i,p,q, j, ee, THR);

		if(ee < THR) {
			i = p + 1;
			continue;
		}

		/*hotspots->substruct[nh].energy always store the energy of the substructure, 
		  initialized with INF*/
		/*hotspots->substruct[nh].length always store the length of the table, 
		  note the size of the table always 2*length*/

		if(p-i-1 < j-q-1) {
			length =   j-q + 1;
			startbase = q;
		}
		else {
			length = p - i + 1;
			startbase = i;
		}
		//numHotspots: current number of hotspots (not including the new one)
		if (numHotspots < MaxHotspots) { 
			flag = 0;
			for (num = 0; num < numHotspots; num++) {
				if (ee == hotspots->substruct[num].energy && length == hotspots->substruct[num].length&& 
						startbase == hotspots->substruct[num].pairTable[0]&&
						pair_table[startbase] == hotspots->substruct[num].pairTable[length]) {
					flag = 1; 
					//if (TRACE == 1) printf("This hotspot already exists\n");
					break; 
				}
			}
			if (flag != 1) { 
				hotspots->substruct[numHotspots].energy = (int) ee;
				hotspots->substruct[numHotspots].length = length;
				for ( mm = 0; mm < length; mm++){
					hotspots->substruct[numHotspots].pairTable[mm] = startbase + mm;
					hotspots->substruct[numHotspots].pairTable[mm+length] = pair_table[startbase + mm];
				}
				numHotspots++;
			}
		}
		else {
			//search for the highest score(negative)
			min_en = -1000000000;
			for (k = 0; k < numHotspots-1; k++) {
				if (hotspots->substruct[k].energy > min_en) replace = k;
			}
			if (ee < min_en) {  //replace
				nh = replace;
				//if (TRACE == 1) printf("hotspot list full, replace a hotspot with the newly found better one\n");
				for ( mm = 0; mm < length; mm++){
					hotspots->substruct[nh].pairTable[mm] = startbase + mm;
					hotspots->substruct[nh].pairTable[mm+length] = pair_table[startbase + mm];
				}
			}
		}
		i = p+1;
	}  /*end of the first while loop => finsh search one nested structure*/
	hotspots->numHotspots = numHotspots;
	free(pair_table);
	return(numHotspots);

}


void ChildNodeConstraint(int length, struct Node *currentNode, struct Node *rootNode){
	//create childNodes for currentNode according to hotspots
	short i,j,k,hotspotLength,numChild=0,temp1, temp2;
	char tempConstraint[length+1];
	extern struct Hotspots *hotspots;
	short exist[hotspots->numHotspots];
	int overlap;
	currentNode->numChild = 0;

	tempConstraint[length] = 0;

	for (i=0; i<hotspots->numHotspots; i++) {
		exist[i] = 1;
		overlap = 0;
		strncpy(tempConstraint, currentNode->constraint, length);
		hotspotLength = hotspots->substruct[i].length;

		for(j=0; j<hotspotLength; j++) {
			//check whether the hotspot is overlapping with the current constraint
			if (((hotspots->substruct[i].pairTable[j] < 0) || 
						(tempConstraint[hotspots->substruct[i].pairTable[j]-1] == 'x')) ||
					(hotspots->substruct[i].pairTable[j+hotspotLength] > 0 && 
					 tempConstraint[hotspots->substruct[i].pairTable[j+hotspotLength]-1] == 'x')) {
				overlap = 1;
				break;
			}
			if (hotspots->substruct[i].pairTable[j] > 0) 
				tempConstraint[hotspots->substruct[i].pairTable[j]-1] = 'x';
			if (hotspots->substruct[i].pairTable[j+hotspotLength] > 0) 
				tempConstraint[hotspots->substruct[i].pairTable[j+hotspotLength]-1] = 'x';
		}
		if (overlap == 0 && IsAlreadyExist(tempConstraint, rootNode, currentNode)==0) {
			numChild++;
			exist[i] = 0;
		}else {
			exist[i] = 1;
		}

	}

	currentNode->numChild = numChild;
	currentNode->children = (struct Node**)malloc(sizeof(struct Node*)*numChild);
	k = 0;
	for (i=0; i<hotspots->numHotspots; i++) {
		if(exist[i] == 0) {   //new constraints, create child node
			currentNode->children[k] = (struct Node*)malloc(sizeof(struct Node));
			currentNode->children[k]->constraint = (char*)malloc(sizeof(char)*(length+1));
			//memset(currentNode->children[k]->constraint, 0, sizeof(currentNode->children[k]->constraint));
			currentNode->children[k]->fixedPairs = (short*)malloc(sizeof(short)*length);
			currentNode->children[k]->length = length;
			currentNode->children[k]->numChild = 0;
			currentNode->children[k]->secStructure = (short *)malloc(sizeof(short)*length);
			//memset(currentNode->children[k]->secStructure, 0, sizeof(currentNode->children[k]->secStructure));
			currentNode->children[k]->children = NULL;
			strncpy(currentNode->children[k]->constraint, currentNode->constraint,length);
			currentNode->children[k]->constraint[length] = 0;
			hotspotLength = hotspots->substruct[i].length;
			for(j=0; j<hotspotLength; j++) {          
				if (hotspots->substruct[i].pairTable[j] > 0) currentNode->children[k]->constraint[hotspots->substruct[i].pairTable[j]-1] = 'x';
				if (hotspots->substruct[i].pairTable[j+hotspotLength] > 0) currentNode->children[k]->constraint[hotspots->substruct[i].pairTable[j+hotspotLength]-1]='x';
			}
			for (j = 0; j<length; j++) {
				currentNode->children[k]->fixedPairs[j] = currentNode->fixedPairs[j];
			}
			for(j = 0; j<hotspotLength; j++) {
				temp1 = hotspots->substruct[i].pairTable[j];
				temp2 = hotspots->substruct[i].pairTable[j+hotspotLength];
				if (temp1 > 0) currentNode->children[k]->fixedPairs[temp1-1] = temp2;
				if (temp2 > 0) currentNode->children[k]->fixedPairs[temp2-1] = temp1;
			}
			k++;
		}
	}
	if (currentNode->numChild != k) exit(1);
}

int secondaryStruct(char *sequence, int length, struct Node *currentNode, struct Node *rootNode, int MaxHotspots, int energyModel)
{
	int i,j;
	int start_piece, end_piece, length_piece;
	short *temp_Pairs;
	float min_en, energy;
	char *temp_string, *temp_constraint, *temp_piece, *temp_constraint_piece;

	count++;  

	if (count >= 10000) {
		return 0;
	}			
	temp_string = (char*)malloc(sizeof(char)*(length+1));
	temp_constraint = (char*)malloc(sizeof(char)*(length+1));

	strncpy(temp_constraint,currentNode->constraint,length);
	temp_constraint[length] = 0;
	//-----fold current sequence=>2nd struct. in ...(...)... format at => *temp_constraint

	i = 0;
	while( i < length) {
		while((i < length) && currentNode->fixedPairs[i] != 0) i++;
		start_piece = i;   //starting base of a piece
		while((i<length) && (currentNode->fixedPairs[i] == 0)) i++;  
		end_piece = i-1;     //find the end base of a piece
		length_piece = end_piece - start_piece + 1;
		if (length_piece > 3) {
			temp_piece = (char*)malloc(sizeof(char)*(length_piece+1));
			temp_constraint_piece = (char*)malloc(sizeof(char)*(length_piece+1));
			for (j = 0; j < length_piece; j++) {
				temp_piece[j] = sequence[start_piece + j];
				temp_constraint_piece[j] = temp_constraint[start_piece + j];
			}
			temp_piece[length_piece] = 0;
			temp_constraint_piece[length_piece] = 0; 

			min_en = simfold(temp_piece, temp_constraint_piece);

			//paste structure of one piece back to the whole structure
			for (j = 0; j < length_piece; j++) {
				temp_constraint[start_piece + j] = temp_constraint_piece[j];
			}
			free(temp_piece);
			free(temp_constraint_piece);
		}
	}

	//-----transform 2nd struct. to number format
	temp_Pairs = make_pair_table(temp_constraint);

	//-----combine new structure(temp_Pairs) with the fix ones from parent.  
	for(i=0;i<length;i++){
		if(currentNode->fixedPairs[i]==0)
			currentNode->secStructure[i] = temp_Pairs[i+1];
		else
			currentNode->secStructure[i] = currentNode->fixedPairs[i];
	}
	free(temp_Pairs);

	//-----call score() to calculate the score of current 2nd struct and choose next step

	currentNode->score = score(length,sequence,currentNode->secStructure,TRACE,energyModel);


	int fl = (currentNode->score>rootNode->score*T_RATIO  || (rootNode->score - currentNode->score) < 4000) && IsAlreadyExist(NULL, rootNode, currentNode)==0 ; 
	if(fl) {
		InsertRna(currentNode);
	} else {
	}//---need future coding
	if(currentNode->score < rootNode->score*T_RATIO && (rootNode->score - currentNode->score) > 5000) {
		// hotspots are no longer good. Stop searching along this branch.
		free(temp_string); free(temp_constraint);
		return currentNode->score;
	} else {
		for(i=0;i<length;i++) {
			temp_string[i]=currentNode->constraint[i];
		}
		temp_string[length] = 0;

		min_en = simfold(sequence,temp_string);

		FindHotspots(sequence,temp_string,MaxHotspots);

		//clean up before recursive calls
		free(temp_string);
		free(temp_constraint);

		ChildNodeConstraint(length, currentNode, rootNode);

		//hotspots->numHotspots = 0;

		min_en = currentNode->score;
		for (i = 0; i<currentNode->numChild; i++) {
			energy = secondaryStruct(sequence, length, currentNode->children[i], rootNode, MaxHotspots, energyModel);
			if (min_en < energy) min_en = energy;
		}
	}
	return(min_en);
}

int IsAlreadyExist(char *constraint, struct Node *rootNode, struct Node *currentNode){
	int i, notsame;
	
	notsame = 0;

	if (constraint != NULL) {
		if(rootNode == currentNode) return(0);
		if (strcmp(constraint, rootNode->constraint)==0)  {
			return(1);
		} else { 
			for(i = 0; i < rootNode->numChild; i++) {
				if (IsAlreadyExist(constraint, rootNode->children[i],currentNode) == 1) return(1);
			}
		}
		return(0);
	} else {
		if(rootNode == currentNode) return(0);
		notsame = 0;
		for(i = 0; i < (int) strlen(currentNode->constraint); i++) {
			if(currentNode->secStructure[i] != rootNode->secStructure[i]) {
				notsame++;
			}
		}
		if(notsame <= 2) {   //allow one base pair difference
			return(1);
		} else {
			for(i = 0; i < rootNode->numChild; i++) {
				if (IsAlreadyExist(NULL, rootNode->children[i],currentNode) == 1) return(1);
			}
		}
		return(0);
	} 
}



void InsertRna(struct Node* currentNode){
	//insert the node into the array according to its score.
	// by divide-and-conquer
	int i,j, seg;

	if (numRnaStruct == MaxSubOpt) {
		if (listOfNodes[numRnaStruct-1]->score < currentNode->score) numRnaStruct--; 
		//get rid of the smallest score one
		else return;
	}
	if (numRnaStruct == 0) {
		listOfNodes[0] = currentNode;
		numRnaStruct++;
		return;
	}
	if (currentNode->score >= listOfNodes[0]->score) {
		for (i = numRnaStruct-1; i >= 0; i--) {
			listOfNodes[i+1] = listOfNodes[i];
		}
		listOfNodes[0] = currentNode;
		numRnaStruct++;
		return;
	}
	if (currentNode->score <= listOfNodes[numRnaStruct-1]->score) {
		listOfNodes[numRnaStruct] = currentNode;
		numRnaStruct++;
		return;
	}
	seg = numRnaStruct;
	i = numRnaStruct/2;
	while (i >= 0 && i <= numRnaStruct-1) {
		if (i == numRnaStruct-1) {
			listOfNodes[numRnaStruct] = listOfNodes[numRnaStruct-1];
			listOfNodes[numRnaStruct-1] = currentNode;
			numRnaStruct++;
			return;
		}
		if (listOfNodes[i]->score >= currentNode->score && listOfNodes[i+1]->score <= currentNode->score) {
			for (j = numRnaStruct-1; j > i; j--) {
				listOfNodes[j+1] = listOfNodes[j];
			}
			listOfNodes[i+1] = currentNode;
			numRnaStruct++;
			return;
		}
		if (currentNode->score >= listOfNodes[i]->score ){ 
			i -= (seg/4 > 1)?(seg/4):1;
		} else {
			i  += (seg/4 > 1)?(seg/4):1;
		}
		seg = seg/2?seg/2:1;
	}
	return;
}





