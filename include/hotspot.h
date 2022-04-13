#ifndef HOTSPOT_H
#define HOTSPOT_H

struct Substruct {
  int energy;
  int length;
  short* pairTable;
};

struct Hotspots {
  struct Substruct* substruct;
  int numHotspots;
};

struct Node {
    short *secStructure;  //the secondary strucuter of current node.
                            // format is 10,11,12,0,0,0,...1,2,3,...
    short length;             //length of current 2nd struct,i.e.length of sequence
    char *constraint;       //constrains for the current node.
    short *fixedPairs;       //contains fixed pairs of current node.
    short numChild;
    struct Node **children;
    float score;
};

struct Fold {
	char *structure;
	float score;
};

extern void InitHotspots(int MaxHotspots, int length);
extern int  FindHotspots(char* sequence, char* structure,int MaxHotspots);
extern void ClearHotspots(int MaxHotspots);
extern int secondaryStruct(char *sequence, int length, struct Node *currentNode, struct Node *rootNode, int MaxHotspots, int energyModel);
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */

extern void printRnaStruct(short* structure, int length); 

 
#endif
