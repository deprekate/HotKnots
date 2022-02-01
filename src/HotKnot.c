/*****************************************************************
HotKnot: A heuristic algorithm for RNA secondary 
structure prediction including pseudoknots 
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

//#define INF 10000;

// *** Add new energy model code here (MUST match the definitions in Defines.h)
//const int RE = 0;  // id for Rivas&Eddy energy computation in Loop.cpp
//const int DP = 1;  // id for Dirks&Pierce energy computation in Loop.cpp
//const int CC2006 = 2;  // id for Cao&Chen 2006 energy computation in Loop.cpp
//const int CC2006withDP = 3;  // id for Cao&Chen 2006 (b) energy computation in Loop.cpp

/*
   extern struct Hotspot* hotspots;
   extern int noPS=0;
   extern int noGU = 0;
   extern int TRACE=0;
   extern int *sp;
   extern float THR=400.0;//threshold for sub-sequence
   extern float T_RATIO = 0.8; 
   */
/*Consider only subsequent 2ndary structures that have energy lower than 
  T_RATIO*energy of the best non-pseudoknotted structure */
/*
   extern int MaxSubOpt = 20;
   extern struct Node* listOfNodes[50];   // for Rivas&Eddy
   extern struct Node* listOfNodes2[50];  // for Dirks&Pierce
   extern struct Node* listOfNodes3[50];  // for Cao&Chen (a)
   extern struct Node* listOfNodes4[50];  // for Cao&Chen (b)
   extern struct Node* listOfNodes5[50];  // for Cao&Chen (c)
// *** Add new energy model code here (add another listOfNodes)

extern int count;//number of nodes
extern int numRnaStruct;  //total number of different Rna structures
*/

//extern void PlotRna(char* seqName, char *sequence, short *structure, char *filename, float score);


#define PRIVATE static

//static char  scale1[] = "....,....1....,....2....,....3....,....4";
//static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);

/*--------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

	char *string; 
	char line[5000];
	char *structure=NULL, *cstruc=NULL;
	char outputStructure[5000]; 
	char  fname[13], ffname[20], gfname[20];

	char  *ns_bases=NULL, *c;
	int   i, length, l, sym, r;
	double energy, min_en;
	double kT, sfact=1.07;
	int   pf=0, istty;
	int noconv=0;
	int InputFile = 1;
	int StrCon = 0;    
	int endFlag;
	int useParamsFlag = 0;  // added for easy command line usage
	int sequenceDirect = 0;  // Mirela: added for easy command line usage
	char outPSFile[100], bpseqFile[100], seqName[1000], inFile[100], ctFile[100];
	char ctFile2[100], outPSFile2[100], bpseqFile2[100];
	char ctFile3[100], outPSFile3[100], bpseqFile3[100];
	char ctFile4[100], outPSFile4[100], bpseqFile4[100];
	char ctFile5[100], outPSFile5[100], bpseqFile5[100];
	char paramsFile[1000];  // added for easy command line usage
	char currentModel[5];   // added for easy command line usage
	// *** Add new energy model code here

	struct Node *rootNode;
	struct Node *rootNode2;
	struct Node *rootNode3;
	struct Node *rootNode4;
	struct Node *rootNode5;
	// *** Add new energy model code here

	int len = 0;
	int index = -1;  
	char separator;
	char fileName[200];  // filename of input file sequence, without the path

	// *** Add new energy model code here (also change output path as desired)
	char outpath[30] = "output/";	//"TestSeq/";  // starting part of path for output files assuming directory of HotKnot executable is HotKnots/bin/
	char outpath1[20] = "";	//"RivasEddy/";  // directory to put files generated using Rivas&Eddy (preceeded by outpath)
	char outpath2[20] = "";	//"DirksPierce/";  // see above
	char outpath3[20] = "";	//"CaoChen_a/";  // see above
	char outpath4[20] = "";	//"CaoChen_b/";  // see above
	char outpath5[20] = "";	//"CaoChen_c/";  // see above


	char tempFile[10]="hello";
	FILE* input_file;
	int MaxHotspots = 200;
	int first = 0;

	string=NULL;

	strcpy(currentModel,"DP");  // added for each command line usage - default model is DP 
	for (i=1; i<argc; i++) {
		if (argv[i][0]=='-') 
			switch ( argv[i][1] )
			{
				case 'n':
					if ( strcmp(argv[i], "-noGU")==0) noGU = 1;
					break;
				case 'I':
					if (i==argc-1) usage();
					strcpy(seqName, argv[++i]);
					break;
				case 'c':
					StrCon = 1;
					break;
				case 't':
					TRACE = 1;
					break;	
				case 'b':
					first = 1;
					break;
				case 'm':  // added for easy command line usage
					if (i==argc-1) usage();
					strcpy(currentModel, argv[++i]);
					break;
				case 'p':  // added for easy command line usage 
					if (i==argc-1) usage();
					useParamsFlag = 1;
					strcpy(paramsFile, argv[++i]);
					break;
				case 's':  // Mirela: added  for easy command line usage 
					if (i==argc-1) usage();
					sequenceDirect = 1;
					i++;
					string = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
					strcpy(string, argv[i]);
					break;
				default: usage();
			} 
	}

	//    char config_file[200] = "params/pairfold.conf";  // CHANGED to lines below
	char config_file[200] = "./params/multirnafold.conf";
	char config_filePK[200] = "./params/pkenergy.conf";

	int dna_or_rna = RNA;
	double temperature = 37;

	//   init_data (config_file, dna_or_rna, temperature);  // CHANGED to lines below
	init_data(argv[0], config_file, dna_or_rna, temperature);
	fill_data_structures_with_new_parameters ("./params/turner_parameters_fm363_constrdangles.txt");
	init_dataPK(argv[0], config_filePK, dna_or_rna, temperature);

	// added for easy command line usage
	if (useParamsFlag == 1) {
		if (strcmp(currentModel,"DP") == 0)
			fill_data_structures_with_new_parameters_PK_DP (paramsFile);
		else if (strcmp(currentModel,"CC") == 0)
			fill_data_structures_with_new_parameters_PK_CC2006b (paramsFile); 
		else {
		}
	}

	if (sequenceDirect == 1) {
		strcpy(fileName,"commandline");
	}

	if (!sequenceDirect)    // Mirela: added  for easy command line usage 
	{
		// determine the name of the input file (without preceeding path)
		len = strlen(seqName);
		separator = '/';

		strcpy (fileName, "");

		for (i = len; i >=0; i--){
			// make it work on Linux 
			if (seqName[i] == '/'){
				separator = '/';
				index = i;
				break;
			}
		}
		for (i=index+1; i < len; i++){
			fileName[i-(index+1)] = seqName[i];
		}
		fileName[i-(index+1)] = '\0';


		// open input file and extract necessary contents       
		strcpy(inFile, seqName);
		strcat(inFile,".seq");
		if ((input_file = fopen(inFile, "r")) == NULL) {
			fprintf(stderr, "Cannot open %s\n", inFile);
			fprintf(stderr, "please only provide the sequence name\n");
			usage();
			return 0;
		}
		fscanf(input_file,"%s", line);
		string = (char *) malloc((strlen(line)+1)*sizeof(char));
		strcpy(string, line);
	}   // end if sequenceDirect

	length = (int) strlen(string);
	structure = (char *) malloc((length+1)*sizeof(char));

	for (l = 0; l < length; l++) 
		structure[l] = '.';
	structure[length] = 0;

	if (StrCon == 1) {
		fscanf(input_file, "%s", line);

		int sclen; 
		sclen = strlen(line);
		if (strlen(line) > length) {
			fprintf(stderr, "--------------------WARNING-----------------------\n");
			fprintf(stderr, "The length of the structure constraint is %d\n", sclen);
			fprintf(stderr, "The length of the sequence is %d\n", length);
			fprintf(stderr, "Extra constraints are ignored!\n");
			line[length] = 0;
			strcpy(structure, line);
		}
		else if (strlen(line) < length) {
			fprintf(stderr, "--------------------WARNING----------------------\n");
			fprintf(stderr, "The length of the structure constraint is %d\n", sclen);
			fprintf(stderr, "The length of the sequence is %d\n", length);
			fprintf(stderr, "The rest of the sequence is assumed to be unconstrained!\n");
			strcpy(structure, line);
			structure[sclen] = '.';
		}
		else
			strcpy(structure, line);
		printf("\nThe following bases are forced to be single stranded (the first base has index 1): \n");
		for (l = 0; l < length; l++){ 
			if (structure[l] == 'x') printf("%d   %c\n", l+1, string[l]); 
		}
	}

	for (l = 0; l < length; l++) {
		string[l] = toupper(string[l]);
		if (string[l] == 'T') string[l] = 'U';
		if (structure[l] != '.' && structure[l] != 'x') {
			fprintf(stderr, "-------------------WARNING----------------------\n");
			fprintf(stderr, "There are letters other than . and x in the structure constraint, treated as .\n");
			structure[l] = '.';
		}
	}



	if (strcmp(currentModel,"RE")==0) {
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
		//===end of initialization    

		printf("LENGTH OF THE RNA IS %d.\n",length);
		printf("%s\n", string);    

		numRnaStruct = 0;
		InitHotspots(MaxHotspots,length);
		GenerateStemList(length, string, structure);

		endFlag=secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, RE);  // ADDED: for RE energy model
		ClearHotspots(MaxHotspots);
		min_en=endFlag;
		//printf("In total, %d nodes created. \n", count);
		printf("Total number of RNA structures: %d \n", numRnaStruct);
		printf ("Seq: %s\n", string);

		for (i=0; i < numRnaStruct; i++) {
			//printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			bpseq2dp (length, listOfNodes[i]->secStructure, outputStructure);
			printf ("S%d:  %s\t%.2lf\n", i, outputStructure, -listOfNodes[i]->score /1000.0);
			//printRnaStruct(listOfNodes[i]->secStructure,length);
			//printf("\n");
		}

	} else if (strcmp(currentModel,"DP")==0) {
		/////////////////////////// DIRKS and PIERCE ENERGY MODEL /////////////////////////
		//-----initialization of rootNode
		rootNode2=(struct Node *)malloc(sizeof(struct Node));
		rootNode2->secStructure=(short *)malloc((length)*sizeof(short));
		rootNode2->fixedPairs=(short *)malloc((length)*sizeof(short));
		rootNode2->constraint=(char *)malloc((length+1)*sizeof(char));
		rootNode2->constraint[length] = 0;
		rootNode2->numChild = 0;
		rootNode2->length = length;
		rootNode2->score = 0;
		for(i=0;i<length;i++){
			rootNode2->secStructure[i]=0;
			rootNode2->fixedPairs[i]=0;
			//rootNode2->constraint[i]='.';
		}
		strcpy(rootNode2->constraint, structure);
		//===end of initialization    

		numRnaStruct = 0;
		InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
		GenerateStemList(length, string, structure);  // depends only on hairpin, stacked, internal loop parameters

		endFlag=secondaryStruct(string,length,rootNode2,rootNode2, MaxHotspots, DP);  // ADDED: for DP energy model
		ClearHotspots(MaxHotspots);
		min_en=endFlag;
		//printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
		printf("Total number of RNA structures: %d \n", numRnaStruct);

		printf ("Seq: %s\n", string);
		for (i=0; i < numRnaStruct; i++) {
			//printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			bpseq2dp (length, listOfNodes[i]->secStructure, outputStructure);
			printf ("S%d:  %s\t%.2lf\n", i, outputStructure, -listOfNodes[i]->score /1000.0);
			//printRnaStruct(listOfNodes[i]->secStructure,length);
			//printf("\n");
		}
	} else if (strcmp(currentModel,"CC")==0) {
		/////////////////////////// CAO and CHEN (b) ENERGY MODEL /////////////////////////
		//-----initialization of rootNode
		rootNode4=(struct Node *)malloc(sizeof(struct Node));
		rootNode4->secStructure=(short *)malloc((length)*sizeof(short));
		rootNode4->fixedPairs=(short *)malloc((length)*sizeof(short));
		rootNode4->constraint=(char *)malloc((length+1)*sizeof(char));
		rootNode4->constraint[length] = 0;
		rootNode4->numChild = 0;
		rootNode4->length = length;
		rootNode4->score = 0;
		for(i=0;i<length;i++) {
			rootNode4->secStructure[i]=0;
			rootNode4->fixedPairs[i]=0;
			//rootNode4->constraint[i]='.';
		}
		strcpy(rootNode4->constraint, structure);
		//===end of initialization    

		int ii;
		printf("string: %s\n", string);
		for (ii=0 ; ii < 2; ii++){
		string[ii+2] = 'C';
		numRnaStruct = 0;
		InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
		GenerateStemList(length, string, structure);  // depends only on hairpin, stacked, internal loop parameters

		endFlag=secondaryStruct(string,length,rootNode4,rootNode4, MaxHotspots, CC2006b);  // ADDED: for CC2006withDP energy model
		ClearHotspots(MaxHotspots);
		min_en=endFlag;
		//printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
		printf("Total number of RNA structures: %d \n", numRnaStruct);
		printf ("Seq: %s\n", string);
		for (i=0; i < numRnaStruct; i++) {
			//printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			bpseq2dp (length, listOfNodes[i]->secStructure, outputStructure);
			printf ("S%d:  %s\t%.2lf\n", i, outputStructure, -listOfNodes[i]->score /1000.0);
			//printRnaStruct(listOfNodes[i]->secStructure,length);
			//printf("\n");
		}
		}
	} else {
		usage();
		exit(0);
	}

	(void) fflush(stdout);
	free(string);
	free(structure);
	return 0;   //the following code was from the original RNAfold.c, not used here.

}

PRIVATE void usage(void)
{
	nrerror((char *) "Usage:\n"
			"HotKnots { -s sequence | -I filename } <options>\n"
			"   where the input sequence can be given as a string, or in a file called filename.seq\n"
			"   Output:\n"
			"     The optimal and a set of suboptimal structures are displayed on the screen.\n"
			"     In addition, a set of files are writen in the \"output\" directory.\n"
			"     Previous files with the same name are overwriten.\n\n"
			"   Options are:\n"
			"   -m energyModel (RE=Rivas&Eddy model, DP=Dirks&Pierce model (default), CC=Cao&Chen model)"
			"   -p parameterFilename (filename of parameters corresponding to the energy model)\n"
			"\tIf no file is given, the Dirks&Pierce parameters 2003 are used, like in HotKnots 1.0\n"
			"   -noGU (do not allow GU pair)\n"
			"   -t (trace)\n"
			"   -b (output only the ps, bpseq file for the best structure)\n"
			"\nExamples:\n"
			"   ./HotKnots -s GGCGCGGCACCGUCCGCGGAACAAACGG -m DP -p params/parameters_DP09.txt\n"
			"   ./HotKnots -I myseq -m CC -p params/parameters_CC09.txt\n" 

		   );
}


