#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

using namespace std;

#include "HotKnotEnergy.h"
#include "initPK.h"
#include "init.h"
#include "params.h"
#include "simfold.h"
#include "constantsPK.h"
#include "externs.h"

int main(int argc, char *argv[])
{

//    char config_file[200] = "params/pairfold.conf";  // CHANGED to lines below
        char config_file[200] = "../bin/params/multirnafold.conf";
        char config_filePK[200] = "../bin/params/pkenergy.conf";
                                  
        int dna_or_rna = RNA;
        double temperature = 37;
                                  
//   init_data (config_file, dna_or_rna, temperature);  // CHANGED to lines below
        init_data(argv[0], config_file, dna_or_rna, temperature);
    fill_data_structures_with_new_parameters ("../bin/params/turner_parameters_fm363_constrdangles.txt");
        init_dataPK(argv[0], config_filePK, dna_or_rna, temperature);

//	cout << "dangle_top[2][1][0] = " << dangle_top[2][1][0] << endl;
//	cout << &dangle_top[2][1][0] << endl;


	printf("Calling Hotknots\n");
	
	char* sequence = new char[100];
//	sequence = "CCCAAAGGGAAA";
	strcpy(sequence, argv[1]);

	int len = strlen(argv[1]);
//	char* sequence = new char[len];
//	for (int i = 0; i < len; i++)
//		sequence[i] = argv[1][i];
//	sequence[len] = '\0';

	printf("len = %d\n", len);

//	char sequence[13] = "CCCAAAGGGAAA";

	printf("%s\n", sequence);

//	int len = strlen(sequence);
	char* structure = new char[len+1];

	char* fileName = new char[100];
	fileName = "pktest1";

	//	double energy =  hotknots_DP_name(sequence, structure, fileName);
//	double energy =  hotknots_CC2006b_name(sequence, structure, fileName);
        double energy =  hotknots_CC2006b(sequence, structure);

	printf("Done call to hotknots: energy = %f\n", energy);

//char * tseq = new char[100];
//tseq = "GCCCCAUGGAGGUGGCUGGGGCCAGCCUCAUGGAGGUGGCUGGGG";
//	char * tstructure = new char[100];
//	tstructure = "(((((.....)).)))...(((((.((((...)))))))))....";

	double free_value = 0;
//	int reset_c = 0;
//	int ignore_dangles = 0;
//	double *c = new double[len];
//	char tseq[100];
//	char tstruct[100];
//	tseq = "";
//	tseq[45] = 

	double fr = 0;
    char cseq[100];
    char cstruct[100];
             
        strcpy(cseq, "GCCCCAUGGAGGUGGCUGGGGCCAGCCUCAUGGAGGUGGCUGGGG");
        cseq[45] = '\0';
        strcpy(cstruct, "(((((.....)).)))...(((((.((((...)))))))))....");
        cstruct[45] = '\0';
             
        printf("SIMFOLD TEST: %f\n", get_feature_counts_restricted (cseq, cstruct, NULL, fr, 1, 0));

	printf("SIMFOLD TEST getenergy: %f\n", free_energy_simfold(cseq, cstruct));

	fr = 0;
	char cseq2[100];
	char cstruct2[100];
	strcpy(cseq2, "GGCCAGCCUCAUGGAGGUGGCU");
	cseq[22] = '\0';
	strcpy(cstruct2, "(((((.((((...)))))))))");
	cstruct2[22] = '\0';

        printf("SIMFOLD TEST getenergy2: %f\n", free_energy_simfold(cseq2, cstruct2));

//	double f = 0;
//	double energy_simfold = get_feature_counts_restricted (tseq, tstructure, c, free_value, reset_c, ignore_dangles);

//	printf("Call to simfold for PtPrp: energy = %f \n", energy_simfold);

	return 0;
}
