#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <iostream>

#include "goodStem.h"
#include "simfold.h"
#include "hotspot.h"
#include "utils.h"
#include "initPK.h"  // for init_dataPK()
#include "params.h"  // for init_dataPK()
#include "init.h"  // for init_dataPK()
#include "constantsPK.h"
//#include "externs.h"
#include "HotKnotEnergy.h"
#include "paramsPK.h"

//using namespace std;


extern int count;

struct Hotspot* hotspots;
int noGU = 0;
int TRACE=0;    // Mirela: to change back to 0
//int *sp;
float THR=400.0;  // Mirela: was 400.0. //threshold for sub-sequence
float T_RATIO = 0.8; //Mirela: was 0.8
int MaxSubOpt = 20;
int count;//number of nodes
int numRnaStruct;  //total number of different Rna structures
struct Node* listOfNodes[50];   // for Rivas&Eddy

/*--------------------------------------------------------------------------*/

int initiate(char *currentModel, char *paramsFile, char *config_file, char *config_filePK){
	//    char config_file[200] = "params/pairfold.conf";  // CHANGED to lines below
	//char config_file[200] = "./params/multirnafold.conf";
	//char config_filePK[200] = "./params/pkenergy.conf";

	int dna_or_rna = RNA;
	double temperature = 37;

	char *temp;
	temp = (char *) malloc((3)*sizeof(char));
	temp[0] = '.';
	temp[1] = '/';
	temp[2] = '\0';

	//   init_data (config_file, dna_or_rna, temperature);  // CHANGED to lines below
	init_data(config_file, dna_or_rna, temperature);
	init_dataPK(temp, config_filePK, dna_or_rna, temperature);
	// added for easy command line usage
	if (strcmp(currentModel,"DP") == 0){
		fill_data_structures_with_new_parameters_PK_DP(paramsFile);
	} else if (strcmp(currentModel,"CC") == 0){
		fill_data_structures_with_new_parameters_PK_CC2006b(paramsFile); 
	}else {
		fill_data_structures_with_new_parameters(  (char *) "./params/turner_parameters_fm363_constrdangles.txt");
	}
	return 1;
}

struct Node* best( char *sequence, char *currentModel){
	char *string=NULL;
	char *structure=NULL;

	int   i, length, l;
	double min_en;
	int endFlag = 0;
	// *** Add new energy model code here

	struct Node *rootNode;
	// *** Add new energy model code here

	int MaxHotspots = 200;
	length = (int) strlen(sequence);
	string = (char *) malloc((length+1)*sizeof(char));
	structure = (char *) malloc((length+1)*sizeof(char));

	strcpy(string, sequence);

	for (l = 0; l < length; l++) 
		structure[l] = '.';
	string[length] = '\0';
	structure[length] = '\0';

	for (l = 0; l < length; l++) {
		if (structure[l] != '.' && structure[l] != 'x') {
			structure[l] = '.';
		}
	}
	/////////////////////////// RIVAS and EDDY ENERGY MODEL /////////////////////////
	//-----initialization of rootNode
	rootNode=(struct Node *)malloc(sizeof(struct Node));
	rootNode->secStructure=(short *)malloc((length)*sizeof(short));
	rootNode->fixedPairs=(short *)malloc((length)*sizeof(short));
	rootNode->constraint=(char *)malloc((length+1)*sizeof(char));
	rootNode->constraint[length] = 0;
	rootNode->numChild = 0;
	rootNode->length = length;
	rootNode->score = 0;
	for(i=0;i<length;i++){
		rootNode->secStructure[i]=0;
		rootNode->fixedPairs[i]=0;
		//rootNode->constraint[i]='.';
	}
	strcpy(rootNode->constraint, structure);
	count = 0;
	numRnaStruct = 0;
	//===end of initialization    


	InitHotspots(MaxHotspots,length);
	GenerateStemList(length, string, structure);
	if (strcmp(currentModel,"RE")==0) {
		endFlag = secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, RE);  // ADDED: for RE energy model
	} else if (strcmp(currentModel,"DP")==0) {
		endFlag = secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, DP);  // ADDED: for DP energy model
	} else if (strcmp(currentModel,"CC")==0) {
		endFlag = secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, CC2006b);  // ADDED: for CC2006withDP energy model
	}
	ClearHotspots(MaxHotspots);
	min_en=endFlag;
	free(structure);
	return listOfNodes[0];
}

