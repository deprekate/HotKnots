#ifndef STACK_H
#define STACK_H


#include "Defines.h"
#include "Input.h"

class Stack{ 

public:  
//FUNCTIONS
	Stack(ReadInput * R);
	~Stack();
	
	int Add(int a, int&b, int &e);		//Adds an element to the stack 
										//and if a closed region is found the function returs 1 and [b,e] )
	region pop(); 						// pops the top element of the stack and returns it
	region Top();						// returns the top elements on the stack without removing it
	void push(region r);
	
	void printPrevStack(int i);
	
//	void Add(int a);					//Adds an element to the stack
//	int Add(int a, int&b, int &e);		//Adds an element to the stack 
//										//and if a closed region is found the function returs 1 and [b,e] 
//										
//	void Remove(int a);					//Removes an element from the stack
//	void RemoveM(int a);				//Updates prevM and nextM
//	void Push(int a);					//Pushes an element on the stack
//	int Top();							// Returns The element on top of the stack
//	void putMarkB(int a);				// Marks an element as 'B'
//	void putMarkP(int a);				// Marks an element as 'P'
//	void putMarkN(int a);				// Marks an element as 'N'
//	int markB(int a);					// Has the element mark B?
//	int markP(int a); 					// Has the element mark P?
//	int markN(int a); 					// Has the element mark N?
//	int Pred(int a);					// The previous element before $a$ on the stack

//ATTRIBUTES
	ReadInput * Input;
	int * PrevInStack;	//The top element on the stack when you put this element
//	int NumMarked; 

private: 
//	T_stackelem* elements; 
	int top;
	region elem[MaxN];
//	int Number;  //The Stack Pointer 
//	int NumberM; //The stack Pointer for those elements that have P&N&!B


};
#endif

