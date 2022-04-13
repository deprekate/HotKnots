/*****************************************************************
HotKnot: A heuristic algorithm for RNA secondary 
structure prediction including pseudoknots
File: goodStem.c
Description: 
This file is used to generate a list of hotspots, 
using local sequence alignment between the original 
sequence and the reverse complement of the sequence.

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
#include "hotspot.h"
#include "constantsPK.h"

extern int MaxSubOpt; 
extern float THR;    //threshold for a stem being good.
extern int TRACE; 
extern struct Hotspots* hotspots;   //global variable
extern int sPair(int len, char* s0, int i, int j);
extern long score(short length, char *s0, short *p, int TRACE);
extern int F2(int len, char* s,
		int i,int k,int l,int j);
int bPenalty = -1000;  // CHECK: this is used later -- what does it mean?
int intPenalty = -9999;  // CHECK: this is used later -- what does it mean?
// Mirela: according to the Hotknots 05 paper, these penalties are added 
//  in order to keep the total list of hotspots from growing too large.
// These penalties are in cal/mol
// Mirela: tried changing thesepenalties to 0, but I got a higher miss rate


// Mirela: I'm trying to adjust these parameters so that I get a lower lower bound of the miss rate.

void revComp(int length, char* sequence, char* rcSequence); 
//This function is used to generate reverse complement
void GenerateStemList(int length, char* sequence, char* constraint); 
//This function generates a list of hotspots and stores it in hotspots.
void printMatrix(int** m, int length);
void trace(int length, char* seq, int si, int sj,  int** V, int** G, int** F, int**E, int**R, int** flag, short record);
void traceback(int length, char* seq, int** V, int** G, int** F, int**E, int**R, int** flag);
extern void InitHotspots(int MaxHotspots, int length);

/*
   int main()
   {
   char Seq[] = "CGAGGGGCGGUUGGCCUCGUAAAAAGCCGC";
   char constraints[200];
   for (int i=0; i < strlen(Seq); i++)   constraints[i]='.';
   constraints[strlen(Seq)] = '\0';
   GenerateStemList(strlen(Seq), Seq, constraints);
// InitHotspots(10, 35); 
return 0;
}
*/

void GenerateStemList(int length, char* sequence, char* constraint)
{

	char* rvSeq;
	int **V, **G, **F, **E, **R;
	int **flag; //used for traceback
	int i,j,jo,delta;
	int infty = -9999999;

	length = length + 1;

	rvSeq = (char*)malloc(length*sizeof(char));
	V = (int**)malloc(length*sizeof(int*));  //optimal alignment
	G = (int**)malloc(length*sizeof(int*));  //(s[i],T[j]) matches
	F = (int**)malloc(length*sizeof(int*));  //s[i] is unmatched, s[i-1] matched with T[j]  
	E = (int**)malloc(length*sizeof(int*));  //T[j] is unmatched, s[i] matched with T[j-1]
	R = (int**)malloc(length*sizeof(int*));  //s[i-1] matched with T[j-1], s[i],T[j] both unmatched.
	flag = (int**)malloc(length*sizeof(int*)); //used to mark whether the trace already visited or not

	for (i = 0; i < length; i++) {
		V[i] = (int *) malloc(length*sizeof(int));
		G[i] = (int *) malloc(length*sizeof(int));
		F[i] = (int *) malloc(length*sizeof(int));
		E[i] = (int *) malloc(length*sizeof(int));
		R[i] = (int *) malloc(length*sizeof(int));
		flag[i] = (int *) malloc(length*sizeof(int));
	}


	//generate its reverse complement
	revComp(length-1, sequence, rvSeq);

	//start dynamic programming, fill in matrices V, F, E and R.

	//initialize row 0 and column 0
	for (i = 0; i < length; i++) {
		V[0][i] = 0; V[i][0] = 0;
		G[0][i] = 0; G[i][0] = 0;
		E[0][i] = 0; E[i][0] = 0;
		F[0][i] = 0; F[i][0] = 0;
		R[0][i] = 0; R[i][0] = 0;
		for (j = 0; j < length; j++) {
			flag[i][j] = 0;
		}
	}

	for (i = 1; i < length; i++) {
		for (j = 1; j < length - i - 3; j++) {
			//initialize flag matrix to false
			flag[i][j] = 0;

			//recurrence
			// Mirela: try to skip adding this sPair
			if (V[i-1][j-1] == 0) {
				jo = length-j;
				int c = ((sequence[i-1] == 'A' && sequence[jo-1] == 'U') ||
						(sequence[i-1] == 'C' && sequence[jo-1] == 'G') ||
						(sequence[i-1] == 'G' && sequence[jo-1] == 'C') ||
						(sequence[i-1] == 'G' && sequence[jo-1] == 'U') ||
						(sequence[i-1] == 'U' && sequence[jo-1] == 'A') ||
						(sequence[i-1] == 'U' && sequence[jo-1] == 'G'));	
				//printf ("i=%d, jo=%d, c=%d\n", i, jo, c);
				if (constraint[i-1] == 'x' || constraint[jo-1] == 'x') {   
					// the -1 is because it's shifted: V starts from 1, but sequence/constraint start from 0
					delta = infty;
				}	
				else delta = c*100;
				//else delta = -sPair(length-1, sequence, i, jo);  // Mirela:bug?
				//else delta = -sPair(length-1, sequence, i-1, jo+1);    // Mirela: changed these indeces
				//when calculating energy,the index of the first letter is assumed to be 1.
				//the index of the first letter in sequence and constraint is 0.
			}
			else if (V[i-1][j-1] == G[i-1][j-1]) { //s[i-1] matched with s[j-1]
				jo = length-j; 
				if (constraint[i-1] == 'x' || constraint[jo-1] == 'x')
					delta = infty;
				else 
					delta = -F2(length-1, sequence, i-1, i, jo, jo+1); // stack base pairs (i-1.jo+1) with (i,jo)

			}
			else if (V[i-1][j-1] == F[i-1][j-1]) { //s[i-1] unmatched, s[i-2] matched with s[j-1]
				jo = length-j;
				if (constraint[i-1] == 'x' || constraint[jo-1] == 'x')
					delta = infty;
				else
					delta = -F2(length-1, sequence, i-2,i,jo,jo+1);  //bulge of length 1
			}
			else if (V[i-1][j-1] == E[i-1][j-1]) { //T[j-1] unmatched, s[i-1] matched with T[j-2]
				jo = length-j;
				if (constraint[i-1] == 'x' || constraint[jo-1] == 'x')
					delta = infty;
				else
					delta = -F2(length-1, sequence, i-1, i, jo, jo+2); //bulge of length 1
			}
			else if (V[i-1][j-1] == R[i-1][j-1]) { //s[i-1],T[j-1] both unmatched, s[i-2]matched with T[j-2]
				jo = length-j;
				if (constraint[i-1] == 'x' || constraint[jo-1] == 'x')
					delta = infty;
				else 
					delta = -F2(length-1,sequence, i-2,i,jo,jo+2); //interloop of length 1
			}


			G[i][j] = V[i-1][j-1] + delta;  //s[i] matched with T[j]
			F[i][j] = G[i-1][j] + bPenalty; //s[i] unmatched, thus s[i-1] has to be matched with T[j]
			E[i][j] = G[i][j-1] + bPenalty; 
			R[i][j] = G[i-1][j-1] + intPenalty; //both s[i],T[j] unmatched.


			V[i][j] = 0;
			if (G[i][j] >= V[i][j]) V[i][j] = G[i][j];
			if (F[i][j] >= V[i][j]) V[i][j] = F[i][j];
			if (E[i][j] >= V[i][j]) V[i][j] = E[i][j];
			if (R[i][j] >= V[i][j]) V[i][j] = R[i][j];
		}
	}

	traceback(length, sequence, V, G, F, E, R, flag);


	for (i = 0; i < length; i++) {
		free(V[i]);
		free(G[i]);
		free(F[i]);
		free(E[i]);
		free(R[i]);
		free(flag[i]);
	}
	free(rvSeq);
	free(V);
	free(G);
	free(F);
	free(E);
	free(R);
	free(flag);
}

void traceback(int length, char* seq, int** V, int** G, int** F, int**E, int**R, int** flag)
{
	int i,j,num=0;

	int currentBest = -1, previousBest = 999999, cBi, cBj;

	while(1) {
		for (i = 1; i < length; i++){
			for (j = 1; j < length - i - 3; j++) {
				if ((V[i][j] > currentBest) && (V[i][j] <= previousBest) && (V[i][j] == G[i][j]) && (flag[i][j]==0)) {

					cBi = i;
					cBj = j;
					currentBest = V[i][j]; 
				}
			}
		}
		if (currentBest <= THR || num == 20) {
			// Done
			break;
		}
		else {
			//trace back
			//printMatrix(flag, length);
			trace(length, seq, cBi, cBj, V, G, F, E, R, flag,1); //start tracing from si, sj and recording
			previousBest = currentBest; currentBest = 0;
			num++;
		}
	}
}

void trace(int length, char* seq, int si, int sj,  int** V, int** G, int** F, int**E, int**R, int** flag, short record)
	//trace strating at si, sj (si and sj should be matched).
	//if record = 0, only marking the trace, don't record
	//if record = 1, marking the trace, and also record the stem into hotspot list.
{
	int i, j, k,l;
	int temp[length];
	int ti, tj, startbase,mm,hl;

	for (i = 0; i < length; i++) {temp[i] = 0;}

	i = si;
	j = sj; 

	while (V[i][j] != 0) {
		if (V[i][j] == G[i][j]) //match
		{ 
			if (record == 1) {temp[i] = length - j; temp[length-j] = i;}  
			flag[i][j] = 1;  //marking visited.
			ti = i-1;tj = j-1;
		}

		if (V[i][j] == F[i][j]) {
			if (V[i][j]==G[i][j])  //there are several ways to trace back
			{
				trace(length, seq, i-1, j, V, G, F, E, R, flag,0);
			}
			else {
				if (record == 1) temp[i] = 0;
				flag[i][j] = 1;
				ti = i-1; tj = j;
			}
		}
		if (V[i][j] == E[i][j]) {
			if (V[i][j] == G[i][j] || V[i][j] == F[i][j]) //there are several ways to trace back
			{
				trace(length, seq, i, j-1, V, G, F, E, R, flag,0);
			}
			else {
				if (record == 1) temp[length-j] = 0;
				flag[i][j] = 1;
				ti = i; tj = j-1;
			}
		}
		if (V[i][j] == R[i][j]) {
			if (V[i][j] == G[i][j] || V[i][j] == F[i][j] || V[i][j] == E[i][j]) //there are several ways to trace back
			{
				trace(length, seq, i-1, j-1, V, G, F, E, R, flag,0);
			}
			else {
				if (record == 1) {temp[length-j] = 0; temp[i] = 0;}
				flag[i][j] = 1;
				ti = i-1; tj = j-1;
			}
		}
		if (V[ti][tj] == 0) break;
		i = ti; j = tj;

	}
	ti = i; tj = j;
	if (V[i][j] != G[i][j]) printf("The first one must be a match!\n");
	//till here, we know that the hotspot is a stem between base pairs [ti, length-tj][si,length-sj] 
	//si > ti, sj > tj
	i = ti; k = si; l = length - sj; j = length-tj;
	if (record == 1) {

		hl = (k-i >= j-l)?(k-i+1):(j-l+1);
		startbase = (k-i >= j-l)?i:l;
		if (hl > 2 && V[si][sj] >= THR)   //at least two base pairs
		{ 
			hotspots->substruct[hotspots->numHotspots].energy = V[si][sj];
			hotspots->substruct[hotspots->numHotspots].length = hl;
			for ( mm = 0; mm < hl; mm++){
				hotspots->substruct[hotspots->numHotspots].pairTable[mm] = startbase + mm;
				hotspots->substruct[hotspots->numHotspots].pairTable[mm+hl] = temp[startbase + mm];

			}

			hotspots->numHotspots++;
		} 
	}
}



void printMatrix(int** m, int length)
{
	int i, j;

	for (i = 0; i < length; i++) {
		for (j = 0; j < length; j++) {
			printf(" %3.0f ", (float)m[i][j]);
		}
		printf("\n");
	}
}

void revComp(int length, char* seq, char* rcSeq)
{
	int i;
	char c;

	for (i = 0; i <= length; i++){rcSeq[i] = '\0';}

	for (i = 0; i < length; i++) {
		switch (seq[i]){
			case 'A':  
				c = 'U';
				break;
			case 'U':
				c = 'A';
				break;
			case 'T':
				c = 'A';
				break;
			case 'G':
				c = 'C';
				break;
			case 'C':
				c = 'G';
				break;
			case 'a':  
				c = 'U';
				break;
			case 'u':
				c = 'A';
				break;
			case 't':
				c = 'A';
				break;
			case 'g':
				c = 'C';
				break;
			case 'c':
				c = 'G';
				break;
			default:
				printf("Wrong RNA letters: seq[%d] = %c\n", i, seq[i]);
				exit(1);
		}
		rcSeq[length-1-i] = c;
	}
}



