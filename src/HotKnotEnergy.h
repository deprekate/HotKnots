#ifndef HOTKNOTE_H_
#define HOTKNOTE_H_

// Prototype: double hotknots_DP (char *sequence, char *structure)
//      - sequence is the given RNA sequence to be folded (input parameter)
//      - structure is the predicted secondary structure (output parameter)
//      - the function may return the energy
// Uses the Dirks&Pierce model
// Defaults: no PS is output, GU is included (if want to change this, set it manually in ? TODO )
double hotknots_DP(char *sequence, char *structure);

double hotknots_DP_name(char *sequence, char *structure, char* fileName);

double hotknots_DP_suboptimals(char *seq, char **structures, double* energies, int& num_suboptimals);

double hotknots_CC2006b(char *sequence, char *structure);

double hotknots_CC2006b_name(char *sequence, char *structure, char* fileName);

double hotknots_CC2006b_suboptimals(char *seq, char **structures, double* energies, int& num_suboptimals);

extern struct Hotspot* hotspots;
extern int noPS;
extern int noGU;
extern int TRACE; 
extern int *sp;
extern float THR;//threshold for sub-sequence
extern float T_RATIO;
/*Consider only subsequent 2ndary structures that have energy lower than
  T_RATIO*energy of the best non-pseudoknotted structure */
extern int MaxSubOpt;
extern struct Node* listOfNodes[50];   // for Rivas&Eddy
extern struct Node* listOfNodes2[50];  // for Dirks&Pierce
extern struct Node* listOfNodes3[50];  // for Cao&Chen (a)
extern struct Node* listOfNodes4[50];  // for Cao&Chen (b)
extern struct Node* listOfNodes5[50];  // for Cao&Chen (c)
// *** Add new energy model code here (add another listOfNodes)

extern int count;//number of nodes
extern int numRnaStruct;  //total number of different Rna structures

#endif /*HOTKNOTE_H_*/
