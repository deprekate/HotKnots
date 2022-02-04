/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots
         File: input.cpp
         Description: 
             Main part of parsing algorithm. It Contains functions for identifying 
     	  	 closed regions in the secondary structure.	
 
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

// This file was changed by Hosna to match Baharak's parsing paper


#include "Stack.h"

/*********************************************************************************
*********************************************************************************/
Stack::Stack(ReadInput * R){
//Initialize
//    Number = 0;
//	NumMarked = 0;
//	NumberM = 0;
////    elements = new T_stackelem[MaxN];
//	elem = new stack_elem[MaxN];
	PrevInStack = new int[MaxN];
	memset(PrevInStack, 0, MaxN*sizeof(int));
	Input = R;
	top = 0;
};

/*********************************************************************************
*********************************************************************************/
Stack::~Stack(){
//Initialize
//    delete []elements;
//	delete [] elem;
//	delete [] PrevInStack;
};

region Stack::pop(){
	if (top <= 0)
    {
        fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
        exit (1);
    }
    region result = elem[top];
    top = top -1 ;
    return result;	
	
}

void Stack::push(region r){
	top = top +1;
    elem[top] = r;	
}

region Stack::Top(){
	if (top <= 0){
		region result;
		result.begin = -1;
		result.end = -1;
		return result;	
	}
	return elem[top];
}

void Stack::printPrevStack(int i) {
	printf("PrevInStack[%d]: %d\n", i, PrevInStack[i]);
}

int Stack::Add(int a, int&b, int&e){
	
	PrevInStack[a] = Top().begin;
	if (a < Input->BasePair(a)){  // potentially closed region [a,bp(a)] added to stack
		region toPush;
		toPush.begin = a;
		toPush.end = Input->BasePair(a);
		push(toPush);
		if (DEBUG2){
			printf("Pushed [%d,%d] into the stack \n", toPush.begin, toPush.end);	
		}
		return 0;
	}
	if (Input->BasePair(a) <= 0 ){  // a is unpaired; do nothing
		return 0;
	}
	if (Input->BasePair(a) < a){  // 
		int E = a;
		while(Top().begin > Input->BasePair(a))	{
			region popped = pop();
			if (DEBUG2){
				printf("Popped [%d,%d] from the stack \n", popped.begin, popped.end);	
				printf("E = MAX(%d,%d) \n", E, popped.end);
			}
			E = MAX(E,popped.end);
			
		}
//		if (DEBUG){
//			printf("Top().end \n", popped.begin, popped.end);	
//		}		
		elem[top].end = MAX(E,Top().end);
	}
	if (a == Top().end){  // closed region found; return 1 indicating that can add [b,e] to tree of closed regions
		b = Top().begin;
		e = Top().end;
		pop();
		return 1;
	}
	return 0;
}










/*********************************************************************************
Remove: removes an element from the stack and updates Number
*********************************************************************************/
//void Stack::Remove(int a){    
////    elements[elements[a].prev].next = elements[a].next;
////    elements[elements[a].next].prev = elements[a].prev;
//	elem[elem[a].prev].next = elem[a].next;
//	elem[elem[a].next].prev = elem[a].prev;
//
////	RemoveM(a);
//
//    if (Number == a)
////        Number = elements[a].prev;
//		Number = elem[a].prev;
//};

/*********************************************************************************
RemoveM: removes an element from the stack in case it has P&N mark (and not B mark)
and updates NumberM
*********************************************************************************/
//void Stack::RemoveM(int a){
//    elements[elements[a].prevM].nextM = elements[a].nextM;
//    elements[elements[a].nextM].prevM = elements[a].prevM;
//
//    if (NumberM == a)
//        NumberM = elements[a].prevM;
//};

/*********************************************************************************
Push: Pushes an element on the stack
*********************************************************************************/
//void Stack::Push(int a){
//Push element 'a' in the stack
//Number start form 1 to $Number
//    elements[a].prev = Number;
//    elements[a].next = 0;
//	elements[a].prevM = NumberM;
//	elements[a].nextM = 0;
//    elements[a].MarkB = 0;
//    elements[a].MarkP = 0;
//    elements[a].MarkN = 0;
//
//    elements[Number].next = a;
//    Number = a;
//
//    elements[NumberM].nextM = a;
//    NumberM = a;
//};

/*********************************************************************************
Top: Returns The element on top of the stack
*********************************************************************************/
//int Stack::Top(){
//  return Number;
//
//};

/*********************************************************************************
putMarkB: Marks an element as 'B'
*********************************************************************************/
//void Stack::putMarkB(int a){
//    elements[a].MarkB = 1;
//	NumMarked++;
//
//};

/*********************************************************************************
putMarkP: Marks an element as 'P'
*********************************************************************************/
//void Stack::putMarkP(int a){
//    elements[a].MarkP = 1;
//	NumMarked++;
//};

/*********************************************************************************
putMarkN: Marks an element as 'N'
*********************************************************************************/
//void Stack::putMarkN(int a){
//    elements[a].MarkN = 1;
//	NumMarked++;
//};

/*********************************************************************************
markB: Has the element mark B?
*********************************************************************************/
//int Stack::markB(int a){
//    return elements[a].MarkB;
//};

/*********************************************************************************
markP: Has the element mark P?
*********************************************************************************/
//int Stack::markP(int a){
//    return elements[a].MarkP;
//};

/*********************************************************************************
markN: Has the element mark N?
*********************************************************************************/
//int Stack::markN(int a){
//    return elements[a].MarkN;
//};

/*********************************************************************************
*********************************************************************************/
//int Stack::Pred(int a){    
////    return elements[a].prev;
//	return elem[a].prev;
//}

/*********************************************************************************
*********************************************************************************/
//void Stack::Add(int a){
//    int b, e;
//    if (Add(a, b, e))
//        printf("Found Close Region [%d, %d]\n",b,e); fflush(stdout);
//}

/*********************************************************************************
Add: Adds an element to the stack and if a closed region is found the function returs
1 and [b,e] (borders of the closed region)
*********************************************************************************/
//int Stack::Add(int a, int&b, int &e){
//
//	PrevInStack[a] = Top();
//    
//    //the element is the left base of a base pair
//    if (a < Input->BasePair(a)){ 
//        Push(a);
//		return 0;
//    };
//
//    //the element is unpaired
//    if (Input->BasePair(a) <= 0)  
//        return 0;
//                
//    //the element is the right base of a base pair
//	int pa = Input->BasePair(a);
//    putMarkB(pa);
//    if (Top() == pa)	       //bp(a) is on top of the stack
//		if (!markP(pa)){       // [pa, a] is a closed region
//        Remove(pa);
//        //closed region [h,a] is identified
//        //which is an unpseudoknotted closed region
//        b = pa; e = a;
//        return 1;
//    	}
//    else{
//		if ( markP(Pred(pa)) && markB(Pred(pa)) && !markN(Pred(pa)) ) {
//			Remove(pa);
//			int h = Top();
//			Remove(h);
//			//closed region [h,a] is identified
//			//which is a pseudoknotted closed region
//			b = h; e = a;
//			return 1;
//
//		};
//	}
//
//    //pseudoknotted base pair is found
//    //Marking elements which are now identified as pseudoknotted ones
//    putMarkP(pa);
//    int top = NumberM;
//	while  ( top > pa) {
//        putMarkP(top);
//        putMarkN(top);
//        if (markP(top) && markB(top) && markN(top))
//            Remove(top);
//
//		else if (markN(top))
//			RemoveM(top);
//
//		top = NumberM;
//    };
//
//	if (markP(pa) && markB(pa) && markN(pa))
//        Remove(pa);
//
//	
//
//    return 0;
//}



