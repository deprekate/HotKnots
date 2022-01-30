/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots
         File: computeEnergy.cpp
         Description: 
			Given a seq file(contains RNA primary structure) and its bpseq file 
			(base pair file, contains RNA secondary structure, see examples in TestSeq/RealStruct),
			compute its free energy value and plot the arc diagram for its structure.
			
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

#include <iostream>
#include <iomanip>

#include "Stack.h"
#include "Loop.h"
#include "Bands.h"
#include "Input.h"
#include "initPK.h"  // for parameter_init()
#include "simfold.h"
#include "paramsPK.h"
#include "commonPK.h"
//#include "init.h"        // for parameter_init()
//#include "params.h"		 // for parameter_init()
//#include "Defines.h"  // July 16 - removed

using namespace std;


extern void PlotRna(char* seqName, char *sequence, short *structure, char *filename,float score);
void usage(void);


/******************************************************************
parameter_init: initializing using simfold program (written by Mirela Andronescu)
*******************************************************************/
void parameter_init(char* argv0)
{
	// argv0 is the path of the executable
    // *** configuration file

//    char config_file[200] = "/cs/beta/People/Cpop/2007/PKEnergy/Version1.0/params/pairfold.conf";
    char config_file[200] = "../bin/params/multirnafold.conf";
    char config_filePK[200] = "../bin/params/pkenergy.conf";

    // what to fold: RNA or DNA
    int dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37;

    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again (see below)
    init_data(argv0, config_file, dna_or_rna, temperature);
    fill_data_structures_with_new_parameters ("../bin/params/turner_parameters_fm363_constrdangles.txt");

    init_dataPK(argv0, config_filePK, dna_or_rna, temperature);

}


/******************************************************************
main: taking input files, call appropriate functions for: storing
essential information in the structures which will be used by the 
program, identifying closed regions, adding the closed regions to 
the tree, computing the free energy of the secondary structure
and drawing the plots.
*******************************************************************/
int main(int argc, char ** argv){

	parameter_init(argv[0]);

	if (DEBUG2) {
		printf("Read these PK energy parameters for DP (kcal/mol):\n");
		printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
		       pkmodelDP.Ps, pkmodelDP.Psm, pkmodelDP.Psp, pkmodelDP.Pb, pkmodelDP.Pup,
		       pkmodelDP.Pps, pkmodelDP.stP, pkmodelDP.intP, pkmodelDP.a, pkmodelDP.a_p,
		       pkmodelDP.b, pkmodelDP.b_p, pkmodelDP.c, pkmodelDP.c_p);

		printf("Read these PK energy parameters for RE (kcal/mol):\n");
		printf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
		       pkmodelRE.g_interiorPseudo, pkmodelRE.P_tilda, pkmodelRE.P_i, pkmodelRE.Q_tilda, pkmodelRE.M_tilda,
		       pkmodelRE.Gw, pkmodelRE.Gwh, pkmodelRE.Gwi,
		       pkmodelRE.q_unpairedMultiPseudo, pkmodelRE.p_pairedMultiPseudo);

	}


	printf("-------------------------------\nReading the Input\n");

	int mode = 1;
	int printTrace = 0;  // by default, don't print energy trace
	int	printDiagram = 0;  // by default, don't make the PS file with the arc diagram

    char fileSeq[100], fileStruct[100];
	char outPSFile[100];  // RE = "ArcDiagram_RE.ps";
    char outPSFileDP[100];  // DP = "ArcDiagram_DP.ps";
    char outPSFileCC2006[100];  // CC2006 = "ArcDiagram_CCa.ps";
    char outPSFileCC2006withDP[100];  // CC2006withDP = "ArcDiagram_CCb.ps";
    // *** Add new energy model code here

        char paramsFile[1000];  // added for easy command line usage
        char currentModel[5];   // added for easy command line usage

	char prefix[100];

	ReadInput * R;

	if (argc < 2){
		printf("argc = %d < 2 \n", argc);
		usage();
        return 0;
    }

    int useParamsFlag = 0;  // added for easy command line usage
    strcpy(currentModel,"DP");  // added for each command line usage - default model is DP
	int sequenceDirect = 0;
	char* string = NULL;
	char* structure = NULL;
	char* filename = NULL;

    // identifying if the secondary structure input file is in bpseq format or stem format
    // identifying if user wishes to view energy values of components (print energy trace)

        for (int i=1; i<argc; i++) {
            if (argv[i][0]=='-')
                        switch ( argv[i][1] )
                        {
                                case 'b':
					mode = 0;
	                                  i++;
        	                          filename = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
                	                  strcpy(filename, argv[i]);
					break;
				case 'c':
					mode = 1;
	                                  i++;
        	                          filename = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
                	                  strcpy(filename, argv[i]);
					break;
				case 'd':
					mode = 2;
					break;
			    case 's':
                                  if (i==argc-1) usage();
                                  sequenceDirect = 1;
                                  i++;
                                  string = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
                                  strcpy(string, argv[i]);
				  i++;
				  structure = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
                                  strcpy(structure, argv[i]);
                                  break;
			    case 't':
		    		printTrace = 1;
			    	break;
				case 'P':
					printDiagram = 1;
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
			    default: usage();
                        }
        }


/*** Cristina: old parsing code
	if (argv[1][0] == '-')
	{
		switch (argv[1][1] )
    	{
		    case 'b':
				mode = 0;
				break;
		    case 'c':
				mode = 1;
		     	break;
			case 'd':
				mode = 2;
				break;
	}
	else {
		printf("argv[1][0] = %c  \n", argv[1][0]);
		usage();
	}

	if (argc > 3)
	{
		if (argv[2][0] == '-')
		{
			switch (argv[2][1] )
	    	{
			    case 'b':
					mode = 0;
					break;
			    case 'c':
					mode = 1;
			     	break;
				case 's':
					sequenceDirect = 1;
					break;
			    case 't':
			    	printTrace = 1;
			    	break;
				case 'P':
					printDiagram = 1;
					break;
	                    case 'm':  // added for easy command line usage
                        	if (argc<=3) usage();
                	        strcpy(currentModel, argv[3]);
        	                break;
	                    case 'p':  // added for easy command line usage
                        	if (argc<=3) usage();
	                        useParamsFlag = 1;
        	                strcpy(paramsFile, argv[3]);
                	        break;
			    default: usage();
		    }
		}
		else {
			if (mode != 2){
				printf("argv[2][0] = %c  \n", argv[2][0]);
				usage();
			}
		}
	}
*****/

	if (mode == 0 || mode == 1){
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
			strcpy(prefix,"commandline");

		        int size = strlen (structure);
		        short pairseq[size+1];
			detect_original_PKed_pairs_many (structure, pairseq);
			R = new ReadInput(size,string,pairseq);
		}

		if (!sequenceDirect) {
//		strcpy(prefix, argv[argc-1]);  // get the filename
		strcpy(prefix, filename);  // get the filename
		
		strcpy(fileStruct, prefix);
		strcpy(fileSeq, prefix);

		strcat(fileSeq,".seq");

		if (mode == 0) {
		  strcat(fileStruct, ".bpseq");
		}
		else {
		  strcat(fileStruct, ".stem");
		}
		R = new ReadInput(fileSeq, fileStruct);

		}  // end if sequenceDirect

		strcpy(outPSFile, prefix);
		strcpy(outPSFileDP, prefix);
		strcpy(outPSFileCC2006, prefix);
		strcpy(outPSFileCC2006withDP, prefix);


		strcat(outPSFile, "_ArcDiagram_RE.ps");
		strcat(outPSFileDP, "_ArcDiagram_DP.ps");
		strcat(outPSFileCC2006, "_ArcDiagram_CCa.ps");
		strcat(outPSFileCC2006withDP, "_ArcDiagram_CCb.ps");
		// *** Add new energy model code here

		if (DEBUG)
			printf("DEBUG = After readinput\n");

		Stack * s = new Stack(R);
		Bands * B = new Bands(R, s);

		printf("Seq: %s \n", R->CSequence);
		printf("Size: %d \n", R->Size);
	    for (int i = 1; i <= R->Size; i++) {
			printf("%d ", R->Sequence[i]);
		}
		printf("\n-------------------------------\n Making the Loop Tree\n");
		Loop * L = new Loop(0, MaxN+1, R, B, s);

		int a, b; //will store the borders of a closed regoin
		for (int i = 1; i <= R->Size; i++){
			if (R->BasePair(i)>= 0){
			  if (s->Add(i, a, b)){
			  //If a closed region is identifed add it to the tree by calling addLoop
			    L->addLoop(a,b);
			  };
			};
		};

		L->countNumberOfChildren();  // set number of children for the top loop

		L->Print(-1);
		printf("-------------------------------\n");

		if (DEBUG2)
		{
			for (int i = 1; i <= R->Size; i++){
				if (R->BasePair(i)>= 0){
					s->printPrevStack(i);
				}
			}
			printf("\n");

			printf("L->NumberOfUnpairedInPseudo = %d\n", L->RightChild->NumberOfUnpairedInPseudo);
		}

		short* secstructure = new short[R->Size+1];
		char *sequence = new char[R->Size+2];
	        sequence[R->Size+1] = '\0';
		for (int i = 0; i < R->Size+1; i++) {
			secstructure[i] = (short)(R->Sequence[i]);
			if (secstructure[i] == -1) secstructure[i] = 0;
		  	sequence[i] = R->CSequence[i];
	     	printf("%d %c %d \n", i, sequence[i], secstructure[i]);
		}

		if (printDiagram == 1)
		{
			// ADDED: called L->Energy() with RE and DP as parameters
			if (strcmp(currentModel,"RE")==0)
//			if (RE_FLAG)
				PlotRna(prefix, &sequence[1], &secstructure[1], outPSFile, L->EnergyViaSimfold(RE));
			if (strcmp(currentModel,"DP")==0)
				PlotRna(prefix, &sequence[1], &secstructure[1], outPSFileDP, L->EnergyViaSimfold(DP) + L->EnergyDanglingViaSimfold(DP));
//			if (CC2006a_FLAG)
//				PlotRna(prefix, &sequence[1], &secstructure[1], outPSFileCC2006, L->EnergyViaSimfold(CC2006a));
			if (strcmp(currentModel,"CC")==0)
				PlotRna(prefix, &sequence[1], &secstructure[1], outPSFileCC2006withDP, L->EnergyViaSimfold(CC2006b) + L->EnergyDanglingViaSimfold(CC2006b));

			printf("Arc Diagrams printed");
		}


		float totalEnergy = 0;

		cout << setw(15) << left << "Energy Model" << setw(25) << left << "Free Energy (kcal/mol)" << setw(40) << left << "Free Energy without Dangling (kcal/mol)" << endl;
		printf("--------------------------------------------------------------\n");
//		if (RE_FLAG)
		if (strcmp(currentModel,"RE")==0)
		{
			totalEnergy = -L->EnergyViaSimfold(RE);

			cout << setw(15) << left << "Rivas&Eddy" << setw(25) << left << (totalEnergy - L->EnergyDanglingViaSimfold(RE))/1000 << setw(40) << left << totalEnergy/1000 << endl;

			if (printTrace)
				L->printEnergyTrace();
			cout << endl;
		}
		if (strcmp(currentModel,"DP")==0)
		{
			totalEnergy = -L->EnergyViaSimfold(DP);
			cout << setw(15) << left << "Dirks&Pierce" << setw(25) << left << (totalEnergy - L->EnergyDanglingViaSimfold(DP))/1000 << setw(40) << left << totalEnergy/1000 << endl;

			if (printTrace)
				L->printEnergyTrace();
		}
/*
		if (CC2006a_FLAG)
		{
			totalEnergy = -L->EnergyViaSimfold(CC2006a);
			cout << setw(15) << left << "Cao&Chen(a)" << setw(25) << left << (totalEnergy - L->EnergyDanglingViaSimfold(CC2006a))/1000 << setw(40) << left << totalEnergy/1000 << endl;
			if (printTrace)
				L->printEnergyTrace();
		}
*/
		if (strcmp(currentModel,"CC")==0)
		{
			totalEnergy = -L->EnergyViaSimfold(CC2006b);
			cout << setw(15) << left << "Cao&Chen(b)" << setw(25) << left << (totalEnergy - L->EnergyDanglingViaSimfold(CC2006b))/1000 << setw(40) << left << totalEnergy/1000 << endl;
			if (printTrace)
				L->printEnergyTrace();
		}

		// *** Add new energy model code here
/*
		printf("\n");
		printf("PARAMETER TUNING\n");
		printf("Ps Psm Psp Pb Pup Pps a b c stP intP a_p b_p c_p\n");
		cout << g_count_Ps << " " << g_count_Psm << " " << g_count_Psp << " " << g_count_Pb << " " << g_count_Pup << " " <<
				g_count_Pps << " " << g_count_a << " " << g_count_b << " " << g_count_c << " " << g_count_stP << " " <<
				g_count_intP << " " << g_count_a_p << " " << g_count_b_p << " " << g_count_c_p << endl;
*/

	} else if (mode == 2){
		int num_params = get_num_params_PK_DP();
		int num_params_pkfree = get_num_params();
		char* sequence = argv[2];
		char* structure = argv[3];

// density-30 pk with some nested structures
//		char* sequence  = "ccaaagccccccccccccccccccccccccccccccaacgaaggggggggggggggggggggggggggggggg";
//		char* structure = "((...)[{<ABCDEFGHIJKLMNOPQRSTUVWXYZ(..[)..])]}>abcdefghijklmnopqrstuvwxyz";

// density-30 pk
//		char* sequence  = "ccccccccccccccccccccccccccccccgggggggggggggggggggggggggggggg";
//		char* structure = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ)]}>abcdefghijklmnopqrstuvwxyz";

// density-5 pk
//		char* sequence  = "cccccggggg";
//		char* structure = "([{<A)]}>a";

// density-3 pk 
//		char* sequence  = "cccaggg";
//		char* structure = "([{.)]}";

// three nested pks (case 7 in notes - pktest9)
//		char* sequence  = "ccaccggcggccgagg";
//		char* structure = "((.[[))(]][[).]]";

// three nested pks (case 5c in notes - pktest7)
//		char* sequence  = "cccaacgaaggcaacgaagcaacgaagg";
//		char* structure = "(((..[)..])(..[)..](..[)..])";

// nested pk (case 5 in notes - pktest6)
//		char* sequence  = "cccaaaggacaacgaagg";
//		char* structure = "(((...)).(..[)..])";

// two multiloops inside pk, one with internal loop around it (case 3 in notes - pktest5)
//		char* sequence  = "ccaaagcaccaaagcaacggggaaag";
//		char* structure = "((...)(.((...)(..[))))...]";

// two multiloops inside pk, one with stacked pair and internal loop around it (case 3b in notes - pktest5b)
//		char* sequence  = "ccaaagcacccaaagcaacgggggaaag";
//		char* structure = "((...)(.(((...)(..[)))))...]";

// two multiloops inside pk (case 2 in notes - pktest4)
//		char* sequence  = "ccaaagccaaagcaacgggaaag";
//		char* structure = "((...)((...)(..[)))...]";

// one multiloop inside pk (case 1 in notes - pktest3)
//		char* sequence  = "ccaaagcaacggaaag";
//		char* structure = "((...)(..[))...]";

// single pk (call to simfold for stems) (case 4 in notes - pktest2)
//		char* sequence = "ccaaaccggaaagg";
//		char* structure = "((...[[))...]]";

// single pk (no call to simfold)
//		char* sequence = "caaacgaaag";
//		char* structure = "(...[)...]";

// pkfree (same as simfold result)
//		char* sequence = "cccaaaggg";
//		char* structure = "(((...)))";

		double* counter = new double[num_params];
		double** quadratic_matrix = new double*[num_params];
		if (quadratic_matrix == NULL)
		{
			printf ("ERROR! Space could not be allocated for quadratic_matrix_known, ABORT!\n");
			exit(1);
		}
		for (int i = 0; i < num_params; i++)
		{
			counter[i] = 0;
			quadratic_matrix[i] = new double[num_params];
			if (quadratic_matrix[i] == NULL)
			{
				printf ("ERROR! Space could not be allocated for quadratic_matrix_known[%d], ABORT!\n", i);
				exit(1);
			}
			for (int j = i; j < num_params; j++)
				quadratic_matrix[i][j] = 0;
		}

		double free_value = 0;
		printf(" seq = %s, str = %s \n",sequence, structure);

//		get_feature_counts (sequence, structure, counter);
		double retval1 = get_feature_counts_quadratic_PK_DP (sequence, structure, quadratic_matrix, counter, free_value);

		double retval2 = free_energy_PK_DP(sequence, structure);

		printf("Energy1 = %f Energy2 = %f \n", retval1, retval2);
		printf("Free Value: %f\n", free_value);
		printf("PK Counter Values:\n");
		for (int i = num_params_pkfree; i < num_params; i++)
			printf("c[%d]=%f  ", i, counter[i]);

		printf("\nNon-Zero Simfold Counter Values:\n");
		for (int i = 0; i < num_params_pkfree; i++)
			if (counter[i] != 0.0)
				printf("c[%d]=%f  ", i, counter[i]);

		printf("\nAll Non-Zero P_matrix Values:\n");
		for (int i = 0; i < num_params; i++)
		{
			for (int j = i; j < num_params; j++)
				if (quadratic_matrix[i][j] != 0.0)
					printf("P[%d][%d]=%f  ", i, j, quadratic_matrix[i][j]);
		}
		printf("\n");

	}

}


/******************************************************************
*******************************************************************/
void usage(void)
{
/*
  printf("usage:\n"
	  " computeEnergy -b/s sequenceName \n example: \n computeEnergy -b tests/HDV \n"
	  "OR computeEnergy -t -b/s sequenceName (for printing energy trace also) \n example: \n computeEnergy -b tests/HDV \n"
	  "OR \n computeEnergy -d sequence structure \n example: \n computeEnergy -d caaacgaaag \"(...[)...]\" \n"
	  "OR \n computeEnergy -PS -b/s sequenceName \n example: \n computeEnergy -PS tests/HDV \n"
          "OR \n computeEnergy -m model -p paramsFile -b sequenceName (where model= RE, DP (default), CC; and paramsFile= name of parameter file\n");
*/
  printf("usage:\n"
	  " computeEnergy -b/s sequenceName \n example: \n computeEnergy -b tests/HDV \n"
          "OR \n computeEnergy -m model -p paramsFile -s sequence structure (where model= RE, DP (default), CC; and paramsFile= name of parameter file"
	  "\n example: -m \"DP\" -p params/parameters_DP09.txt -s AAACCCUUUGGG \"(((...)))...\"\n"
	  "OR computeEnergy -t -b/s sequenceName (for printing energy trace also) \n example: \n computeEnergy -t -b tests/HDV \n"
	  "OR \n computeEnergy -PS -b/s sequenceName \n example: \n computeEnergy -PS -b/s tests/HDV \n");

  exit(1);
}

