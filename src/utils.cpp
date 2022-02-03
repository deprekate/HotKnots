#include <stdio.h>
#include <string.h>
#include <stdlib.h>

short *make_pair_table(char *structure)
{
	/* returns array representation of structure.
	   table[i] is 0 if unpaired or j if (i.j) pair.  */
	int i,j,hx;
	int length;
	short *stack;
	short *table;

	length = strlen(structure);
	stack = (short *) malloc(sizeof(short)*(length+1));
	table = (short *) malloc(sizeof(short)*(length+2));
	table[0] = length;

	for (hx=0, i=1; i<=length; i++) {
		switch (structure[i-1]) {
			case '(': 
				stack[hx++]=i;
				break;
			case ')':
				j = stack[--hx];
				if (hx<0) {
					fprintf(stderr, "unbalanced brackets in %s\n", structure);
					free(stack); free(table); return NULL;
				}
				table[i]=j;
				table[j]=i;
				break;
			default:   /* unpaired base, usually '.' */
				table[i]= 0;
				break;
		}
	}
	free(stack);
	if (hx!=0) {
		fprintf(stderr, "unbalanced brackets %s\n", structure);
		free(table);
		return NULL;
	}
	return(table);
}

void nrerror(char *message)       /* output message upon error */
{
	fprintf(stderr, "\n%s\n", message);
	exit(0);
}

void printRnaStruct(short* structure, int length) 
	//print out secondary struct in the form: 1-4;13-16, 6-10; 30-34...
{
	int i;
	int st, stp;

	i = 0;
	while(i < length) {
		while(i < length && (structure[i] == 0 || structure[i] < i)) i++;
		st = i;    
		while(i < length && (structure[i]-1 == structure[i+1])) i++; 
		stp = i;
		if (i >= length) break;
		printf("%d-%d; %d-%d,  ",st+1, stp+1, structure[stp],structure[st]);
		i++;
	}
}



//////////////////////////////////////////////////////////////////
////////////// Added by Mirela on August 2, 2008 /////////////////
//////////////////////////////////////////////////////////////////

typedef struct
{
	int top;
	int elem[5000];
} stack_b2d;


void init (stack_b2d *st)
	// PRE:  None
	// POST: Initialize the stack st
{
	st->top = 0;
}

void push (stack_b2d *st, int el)
	// PRE:  st is a valid stack
	// POST: Push an element to the stack
{
	st->elem[st->top++] = el;
}

int pop (stack_b2d *st)
	// PRE:  st is a valid stack, that is not empty
	// POST: pop an element from the stack and return it; decrease top
{
	if (st->top <= 0)
	{
		fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
		exit (1);
	}
	return st->elem[--st->top];
}

int peak (stack_b2d *st)
	// PRE:  st is a valid stack, that is not empty
	// POST: pop an element from the stack and return it; don't decrease top
{
	if (st->top <= 0)
	{
		fprintf (stderr, "The given structure is not valid: more right parentheses than left parentheses\n");
		exit (1);
	}
	return st->elem[st->top-1];
}


int size (stack_b2d *st)
{
	return st->top;
}


void bpseq2dp (int length, short int *pair, char *structure)
	// Function that transforms from "bpseq" format to dot-parentheses format
	// Input parameters:
	//      length = the length of the sequence
	//      pair = an int array with the pairings. unpaired is 0 (not -1!), so all indeces are shifted by 1
	// Output: structure
	// It support 30 types of symbols (see symbol_left and symbol_right below.
	// It uses the first available symbol, and when needed it goes to the next symbol.
	// Note: there's not relationship between the symbols used and the
	//      minimum number of base pairs to break in order to remove the pseudoknots!!
	// Written by Mirela on August 2, 2008.
{        
#define MAXST 30
	char symbol_left[] =  "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	char symbol_right[] = ")]}>abcdefghijklmnopqrstuvwxyz";
	int i, j;
	int *stackindex_pop = new int[length];    // which stack I should pop from
	int stackindex_push;

	// make a local copy of pair, otherwise it modifies the good one
	short int ppair[length];

	stack_b2d stack[MAXST];     // assume max 10 stacks
	for (i=0; i < MAXST; i++)
	{
		init (&stack[i]);
	}
	for (i=0; i < length; i++)
	{
		//printf ("%d\t%d\n", i, pair[i]);
		//pair[i]--;
		ppair[i] = pair[i]-1;
		structure[i] = '.';
		stackindex_pop[i] = -1;
	}
	structure[length] = '\0';

	for (i=0; i < length; i++)
	{
		if (ppair[i] < i && ppair[i] >= 0)    // pop
		{
			//printf ("Stack to pop from for %d is %d\n", i, stackindex_pop[i]);
			pop (&stack[stackindex_pop[i]]);
		}
		if (ppair[i] > i)    // push
		{
			// figure out which stack to push on
			for (j=0; j < MAXST; j++)
			{
				//printf ("i=%d, ppair=%d, j=%d, size=%d\n", i, ppair[i], j, size(&stack[j]));
				// if it's 0, then that's the one
				if (size (&stack[j]) == 0)
				{
					stackindex_push = j;
					break;
				}
				else
				{
					//printf ("i=%d, j=%d, peak=%d\n", i, j, peak (&stack[j]));
					if (peak (&stack[j]) > ppair[i]) // that's the one
					{
						stackindex_push = j;
						break;
					}
				}
			}
			//printf ("Stack for i=%d, ppair=%d is %d\n", i, ppair[i], stackindex_push);
			push (&stack[stackindex_push], ppair[i]);
			stackindex_pop[ppair[i]] = stackindex_push;  // to know where to pop from
			structure[i] = symbol_left[stackindex_push];
			structure[ppair[i]] = symbol_right[stackindex_push];
		}
	}

	//printf ("%s\n", structure);
	delete [] stackindex_pop;
}
