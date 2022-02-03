/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots
         File: defines.h
         
         Includes all globals and structs used. See structs.h, global.h, constants.h
         for defines used by Simfold functionality.
 
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

#ifndef DEFINES_H
#define DEFINES_H
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "constants.h"
#include "constantsPK.h"
#include "structs.h"
#include "structsPK.h"
#include "externs.h"
#include "externsPK.h"
#include "common.h"
#include "commonPK.h"

// Structs
struct T_IntList{

	int Num;	T_IntList * Next;
	int tuning_flag;  //dirty bit flag - used in Loop.cpp parameter tuning code 
					  //(set to 1 if this loop's energy was already counted) 
};
/******************************************
//elements of the stack -- old usage, no longer necessary
*******************************************/
/*
struct T_stackelem { 
    int MarkB , MarkP , MarkN , next, prev, nextM, prevM;
    //$next and $prev keep track of the next and the previous elements on the stack
	//nextM and pervM are the same as next and prev with skipping those nodes who have markN&markP&!markB
};
*/  // July 16 - remove
/******************************************
//elements of the stack -- Hosna, March 8th, 2007
*******************************************/

typedef struct region {
	int begin;
	int end;
	region()
	{
		begin = -1;
		end = -1;
	}
} region;

/*
typedef struct {
	int begin;
	int end;
} region;
region reg = {
	.begin = -1,
	.end = -1
};
*/
/******************************************
//possible types of the loops
*******************************************/
enum LoopType{      stackloop, hairpin,	interior,	multi,	external,	pseudo
};
/******************************************
//possible location status for the loops
*******************************************/
enum PseudoNestedType{	nothing, inBand, unBand, inMulti
};
/******************************************
*******************************************/
struct T_AllLoops{	int base1, base2;
};

/******************************************
*******************************************/
struct B_pattern {

	bool isLeftBorder; // true for the indexs if it is recognized to be the left border of a band region
	int next, prev; // next and prev elements in the pattern	
	int OtherBorder; // right border of the band region
};


//enum SpanBandType{
//	nothing, internal, multi
//};

//struct SpanBand {
//	int next, prev;
//	SpanBandType spanType;  // type of span band loop it is
//};



#endif
