#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>

#include "paramsPK.h"

#include "Stack.h"
#include "Loop.h"
#include "Bands.h"
#include "Input.h"
#include "Defines.h"
#include "params.h"		 // for get_num_params()

using namespace std;

///////////////////////// PARAMETER TUNING ///////////////////////////////

void fill_data_structures_with_new_parameters_PK_DP (char *parameters_file){
	// read the simfold parameters
	fill_data_structures_with_new_parameters(parameters_file);

	// read the DP energy model parameters
	FILE *file;
	char buffer[100];
	char  v1[10];
	int line = 0;

	//printf ("FILENAME: %s\n", filename);
	if ((file = fopen (parameters_file, "r")) == NULL){
		giveupPK ( (char *) "Cannot open file", parameters_file);
	}

	// Assumes the parameters are in this order:
	// double Ps;		// pseudoloop initiation energy (kcal/mol)
	// double Psm;		// penalty for introducing pseudoknot inside a multiloop (kcal/mol)
	// double Psp;		// penalty for introducting pseudoknot inside a pseudoloop (kcal/mol)
	// double Pb;		// penalty for band (kcal/mol)
	// double Pup;		// penalty for unpaired base in pseudoloop or band (kcal/mol)
	// double Pps;		// penalty for nested closed region inside a pseudoloop (kcal/mol)
	// double stP;		// multiplicative penalty for stacked pair that spans a band 
	// double intP;		// multiplicative penalty for internal loop that spans a band
	// double a_p;		// penalty for introducing a multiloop that spans a band (kcal/mol)
	// double b_p;		// penalty for multiloop base pair when the multiloop spans a band (kcal/mol)
	// double c_p;		// penalty for unpaired base in multiloop that spans a band (kcal/mol)

	for (int i = 0; i < get_num_params(); i++){
		fgets(buffer, sizeof(buffer), file);  // ignore the simfold parameters
	}

	// PS_penalty: exterior pseudoloop initiation penalty (originally 9.6 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.Ps = ascii_to_doublePK (v1);

	//PSM_penalty: penalty for introducing pseudoknot inside a multiloop (originally 15 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.Psm = ascii_to_doublePK (v1);

	//PSP_penalty: penalty for introducing pseudoknot inside a pseudoloop (originally 15 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.Psp = ascii_to_doublePK (v1);

	//PB_penalty: band penalty (originally 0.2 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.Pb = ascii_to_doublePK (v1);

	//PUP_penalty: penalty for an un-paired base in a pseudoloop or a band (originally 0.1 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.Pup = ascii_to_doublePK (v1);

	//PPS_penalty: penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band(originally 0.1 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.Pps = ascii_to_doublePK (v1);

	// e_stP = 0.83 * e_s
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.stP = ascii_to_doublePK (v1);

	// e_intP = 0.83 * e_int
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.intP = ascii_to_doublePK (v1);

	//a_penalty: penalty for introducing a multiloop (originally 3.4 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.a = ascii_to_doublePK (v1);

	//b_penalty: penalty for base pair in a multiloop (originally 0.4 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.b = ascii_to_doublePK (v1);

	//c_penalty: penalty for un-paired base in a multi-loop (originally 0)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.c = ascii_to_doublePK (v1);

	//ap_penalty: penalty for introducing a multiloop that spans a band (originally 3.4 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.a_p = ascii_to_doublePK (v1);

	//bp_penalty: base pair penalty for a multiloop that spans a band (originally 0.4 Kcal/mol)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.b_p = ascii_to_doublePK (v1);

	//cp_penalty: penalty for unpaired base in a multiloop that spans a band (originally 0)
	fgets (buffer, sizeof(buffer), file);
	line++;
	sscanf (buffer, "%s", v1);
	pkmodelDP.c_p = ascii_to_doublePK (v1);

	fclose (file);
	//printf ("****** we must have 11 lines by now: LINES = %d\n", line);

}


// USED FOR PARAMETER TUNING (DP energy model)
// returns the free energy (in kcal/mol) of sequence folded into structure
double free_energy_PK_DP(char *sequence, char *structure)
// sequence: sequence (input parameter)
// structure: secondary structure in dot-bracket format () and [] (input parameter)
{
	// create an array of shorts which denotes the index of the base pairs:
	// the first base in the sequence has index 1
	// pairseq[0] = 0 always;
	// pairseq[1] = i where i is the index of the base paired with base 1, etc;
	// pairseq[j] = 0 if base j is unpaired;

	int size = strlen (structure);
	short pairseq[size+1];
	detect_original_PKed_pairs_many(structure, pairseq);

	ReadInput* R = new ReadInput(size, sequence, pairseq);

//	printf("DEBUG = After readinput\n");

	Stack * s = new Stack(R);
	Bands * B = new Bands(R, s);

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

	short* secstructure = new short[R->Size+1];
	char *csequence = new char[R->Size+2];
		csequence[R->Size+1] = '\0';
	for (int i = 0; i < R->Size+1; i++) {
		secstructure[i] = (short)(R->Sequence[i]);
		if (secstructure[i] == -1) secstructure[i] = 0;
	  	csequence[i] = R->CSequence[i];
	}
	float totalEnergy = 0;

	// energy returned in kcal

	totalEnergy = -L->Energy(DP)/1000;
	float totalEnergyDang = -L->EnergyDangling()/1000;


//	if (no_pk_dangling_ends == 0)
	return totalEnergy + totalEnergyDang;
//	else
//		return totalEnergy;
}

// USED FOR PARAMETER TUNING (DP energy model)
double get_feature_counts_quadratic_PK_DP (char *sequence, char *structure, double **quadratic_matrix, double *counter, double &free_value)
// sequence: sequence (input parameter)
// structure: secondary structure in dot-bracket format () and [] (input parameter)
// counter: array where counter[i] is the number of times the i-th feature occurs (output parameter)
// quadratic_matrix: TODO
// free_value: TODO
// Note: The counter, free_value, and quadratic_matrix are automatically reset before passed onto other functions
{
	// create an array of shorts which denotes the index of the base pairs:
	// the first base in the sequence has index 1
	// pairseq[0] = 0 always;
	// pairseq[1] = i where i is the index of the base paired with base 1, etc;
	// pairseq[j] = 0 if base j is unpaired;

	int size = strlen (structure);
	short pairseq[size+1];
	detect_original_PKed_pairs_many(structure, pairseq);

	ReadInput* R = new ReadInput(size, sequence, pairseq);

	//printf("DEBUG = After readinput\n");

	Stack * s = new Stack(R);
	Bands * B = new Bands(R, s);

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

	short* secstructure = new short[R->Size+1];
	char *csequence = new char[R->Size+2];
	csequence[R->Size+1] = '\0';
	for (int i = 0; i < R->Size+1; i++) {
		secstructure[i] = (short)(R->Sequence[i]);
		if (secstructure[i] == -1) secstructure[i] = 0;
			csequence[i] = R->CSequence[i];
		//if (DEBUG)
		//	printf("%d %c %d \n", i, csequence[i], secstructure[i]);
	}

	//PlotRna(prefix, &sequence[1], &secstructure[1], outPSFile, L->Energy());
	//printf("Energy Model The total free energy is %f cal/mol\n", L->Energy());

	float totalEnergy = 0;

	//if (DEBUG){
	//	cout << setw(15) << left << "Energy Model" << setw(25) << left << "Free Energy (kcal/mol)" << setw(40) << left << "Free Energy without Dangling (kcal/mol)" << endl;
	//	printf("--------------------------------------------------------------\n");
	//}

	// energy returned in kcal
	
	// PARAMETER TUNING
	// clear the counter, quadratic_matrix, and free_value

	int num_params = get_num_params_PK_DP();
		free_value = 0;
	for (int i=0; i < num_params; i++){
		counter[i] = 0;
		for (int j=i; j < num_params; j++)
			quadratic_matrix[i][j] = 0;
	}

	int reset_c = 0;
	int ignore_dangles = no_pk_dangling_ends;
	int ignore_AU = 0;  // 0 = do include AU penalties

	totalEnergy = -L->Energy(DP, quadratic_matrix, counter, free_value, reset_c, ignore_dangles);
	float totalEnergyDang = -L->EnergyDangling(DP, quadratic_matrix,counter,free_value,reset_c,ignore_dangles,ignore_AU);

	// CHECK VALUES OF COUNTER, ETC
	//int num_params = get_num_params_PK_DP();
	//int num_params_pkfree = get_num_params();
	//printf("Free Value: %f\n", free_value);
	//printf("Counter Values:\n");
	//for (int i = num_params_pkfree; i < num_params; i++)
	//	printf("c[%d]=%f  ", i, counter[i]);

	//printf("Some P_matrix Values:\n");
	//for (int i = 0; i < num_params/10; i++)
	//	printf("P[%d][%d]=%f  ", num_params_pkfree + structure_type_index_PK("stp")-1, i, quadratic_matrix[i][num_params_pkfree + structure_type_index_PK("stp")-1]);


	//if (DEBUG){
	//cout << setw(15) << left << "Dirks&Pierce" << setw(25) << left << totalEnergy + totalEnergyDang << setw(40) << left << totalEnergy << endl;
	//printf("\n");
	//printf("PARAMETER TUNING\n");
	//printf("Ps Psm Psp Pb Pup Pps a b c stP intP a_p b_p c_p\n");
	//cout << g_count_Ps << " " << g_count_Psm << " " << g_count_Psp << " " << g_count_Pb << " " << g_count_Pup << " " <<
	//		g_count_Pps << " " << g_count_a << " " << g_count_b << " " << g_count_c << " " << g_count_stP << " " <<
	//		g_count_intP << " " << g_count_a_p << " " << g_count_b_p << " " << g_count_c_p << endl;
	//}

	// return the free energy:
	//    - deltaG = x' P x + c' x + f
	//        - where x is the vector of parameters
	//        - P is a symmetric matrix of the coefficients for each quadratic term
	//        - c is a vector of counts for each linear term
	//        - c' means c transposed
	//        - f is a constant

	return totalEnergy + totalEnergyDang;

	//int num_params = get_num_params_PK_DP();
	//double * energy_temp1 = new double[num_params];
	//double * energy_temp2 = new double[num_params];
	//double energy = free_value;

	//for (int i = 0; i < num_params; i++)
	//{
	//	energy_temp1[i] = counter[i] * params_all[i];
	//	for (int j = 0; j < num_params; j++)
	//	{
	//		if (i <= j)
	//			energy_temp2[i] += quadratic_matrix[i][j] * params_all[j];
	//		else  // since only the upper triangle of the matrix is filled out
	//			energy_temp2[i] += quadratic_matrix[j][i] * params_all[j];
	//	}
	//}
	//for (int i = 0; i < num_params; i++)
	//{
	//	energy += params_all[i]*energy_temp2[i] + energy_temp1[i];
	//}

	//return energy;
}


// USED FOR TESTING
void get_feature_counts (char *sequence, char *structure, double *counter)
// sequence: sequence (input parameter)
// structure: secondary structure in dot-bracket format () and [] (input parameter)
// counter: array where counter[i] is the number of times the i-th feature occurs (output parameter)
{
	// create an array of shorts which denotes the index of the base pairs:
	// the first base in the sequence has index 1
	// pairseq[0] = 0 always;
	// pairseq[1] = i where i is the index of the base paired with base 1, etc;
	// pairseq[j] = 0 if base j is unpaired;

	int size = strlen (structure);
	short pairseq[size+1];
	detect_original_PKed_pairs_many(structure, pairseq);


	ReadInput* R = new ReadInput(size, sequence, pairseq);

	//printf("DEBUG = After readinput\n");

	Stack * s = new Stack(R);
	Bands * B = new Bands(R, s);

	int printTrace = 0;  // by default, don't print energy trace

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

	L->Print(-1);
	printf("-------------------------------\n");
	/*
	if (DEBUG2){
		for (int i = 1; i <= R->Size; i++){
			if (R->BasePair(i)>= 0){
				s->printPrevStack(i);
			}
		}
		printf("\n");
		if (L != NULL && L->RightChild != NULL)
			printf("L->NumberOfUnpairedInPseudo = %d\n", L->RightChild->NumberOfUnpairedInPseudo);
	}
	*/
	short* secstructure = new short[R->Size+1];
	char *csequence = new char[R->Size+2];
		csequence[R->Size+1] = '\0';
	for (int i = 0; i < R->Size+1; i++) {
		secstructure[i] = (short)(R->Sequence[i]);
		if (secstructure[i] == -1) secstructure[i] = 0;
		csequence[i] = R->CSequence[i];
		printf("%d %c %d \n", i, csequence[i], secstructure[i]);
	}

	//PlotRna(prefix, &sequence[1], &secstructure[1], outPSFile, L->Energy());

	//printf("Energy Model The total free energy is %f cal/mol\n", L->Energy());

	float totalEnergy = 0;

	cout << setw(15) << left << "Energy Model" << setw(25) << left << "Free Energy (kcal/mol)" << setw(40) << left << "Free Energy without Dangling (kcal/mol)" << endl;
	printf("--------------------------------------------------------------\n");
	if (RE_FLAG){
		totalEnergy = -L->Energy(RE);
		if (no_pk_dangling_ends == 0)
			cout << setw(15) << left << "Rivas&Eddy" << setw(25) << left << (totalEnergy - L->EnergyDangling())/1000 << setw(40) << left << totalEnergy/1000 << endl;
		else
			cout << setw(15) << left << "Rivas&Eddy" << setw(25) << left << totalEnergy/1000 << setw(40) << left << totalEnergy/1000 << endl;

		if (printTrace)
			L->printEnergyTrace();
		cout << endl;
	}
	if (DP_FLAG)  // energy returned in kcal
	{
//		double *c = new double[R->Size];
		double ** quadratic_matrix = NULL;
		double f = 0;
		int reset_c = 0;
		int ignore_dangles = no_pk_dangling_ends;
		int ignore_AU = 0;  // 0 = do include AU penalties

		if (DEBUG)
			printf("Before call to Energy for DP model\n");

		totalEnergy = -L->Energy(DP, quadratic_matrix, counter, f, reset_c, ignore_dangles);

		if (DEBUG)		
			printf("After call to Energy for DP model\n");

		cout << setw(15) << left << "Dirks&Pierce" << setw(25) << left << totalEnergy - L->EnergyDangling(DP, quadratic_matrix,counter,f,reset_c,ignore_dangles,ignore_AU) << setw(40) << left << totalEnergy << endl;

		if (printTrace)
			L->printEnergyTrace();
	}

	printf("\n");
	printf("PARAMETER TUNING\n");
	printf("Ps Psm Psp Pb Pup Pps a b c stP intP a_p b_p c_p\n");
	cout << g_count_Ps << " " << g_count_Psm << " " << g_count_Psp << " " << g_count_Pb << " " << g_count_Pup << " " <<
			g_count_Pps << " " << g_count_a << " " << g_count_b << " " << g_count_c << " " << g_count_stP << " " <<
			g_count_intP << " " << g_count_a_p << " " << g_count_b_p << " " << g_count_c_p << endl;


// 	counter[0] = g_count_Ps;
// 	counter[1] = g_count_Psm;
// 	counter[2] = g_count_Psp;
// 	counter[3] = g_count_Pb;
// 	counter[4] = g_count_Pup;
// 	counter[5] = g_count_Pps;
// 	counter[6] = g_count_a;
// 	counter[7] = g_count_b;
// 	counter[8] = g_count_c;
// 	counter[9] = g_count_stP;
// 	counter[10] = g_count_intP;
// 	counter[11] = g_count_a_p;
// 	counter[12] = g_count_b_p;
// 	counter[13] = g_count_c_p;

}

int get_num_params_PK_DP(){
	// Hosna: Feb 15, 2008
	// DONE - Cristina: May 7, 2008
	// this must return the whole parameter counts, i.e. both from simfold and the 11 from the DP model (note: we ignore a, b, c parameters - see Hosna's thesis for what these values are and mean)
	return 14 + get_num_params();
}





////////////////////// CAO & CHEN (CCb) model ////////////////////////

void fill_data_structures_with_new_parameters_PK_CC2006b (char *parameters_file)
// Assumes temperature is 37.0 Celsius
{

/*
	if (DEBUG2)
	{
		printf("BEFORE NEW PARAMS:\n");
		int ii,jj,kk,ll;
		for (ii=0; ii < NUCL; ii++) {
				for (jj=0; jj < NUCL; jj++) {
						for (kk=0; kk < NUCL; kk++) {
								for (ll=0; ll < NUCL; ll++)
								{
										if (coaxstack_f_b[ii][jj][kk][ll] < INF)
										{
											    printf("c_f_b[%d][%d][%d][%d]= %d\n",ii,jj,kk,ll,coaxstack_f_b[ii][jj][kk][ll]);
										}
								}
						}
				}
		}
		for (ii=0; ii < NUCL; ii++) {
				for (jj=0; jj < NUCL; jj++) {
						for (kk=0; kk < NUCL; kk++) {
								for (ll=0; ll < NUCL; ll++)
								{
										if (coaxstack_m1[ii][jj][kk][ll] < INF)
										{
											    printf("c_m1[%d][%d][%d][%d]= %d\n",ii,jj,kk,ll,coaxstack_m1[ii][jj][kk][ll]);
										}
								}
						}
				}
		}
		for (ii=0; ii < NUCL; ii++) {
				for (jj=0; jj < NUCL; jj++) {
						for (kk=0; kk < NUCL; kk++) {
								for (ll=0; ll < NUCL; ll++)
								{
										if (coaxstack_m2[ii][jj][kk][ll] < INF)
										{
											    printf("c_m2[%d][%d][%d][%d]= %d\n",ii,jj,kk,ll,coaxstack_m2[ii][jj][kk][ll]);
										}
								}
						}
				}
		}
		for (ii=0; ii < CC2006_STEMSIZE; ii++) {
			for (jj=0; jj < CC2006_LOOPSIZE; jj++) {
				if (cc2006_s2_l1[ii][jj] < INF)
				{
					printf("c_s2_l1[%d][%d]= %f\n", ii,jj, cc2006_s2_l1[ii][jj]);
				}
			}
		}
		for (ii=0; ii < CC2006_STEMSIZE; ii++) {
			for (jj=0; jj < CC2006_LOOPSIZE; jj++) {
				if (cc2006_s1_l2[ii][jj] < INF)
				{
					printf("c_s1_l2[%d][%d]= %f\n", ii,jj, cc2006_s1_l2[ii][jj]);
				}
			}
		}
		for (ii=0; ii < 4; ii++) {
			for (jj=0; jj < CC2006_STEMSIZE_FORMULA; jj++) {
				if (cc2006_s2_formula[ii][jj] < INF)
				{
					printf("c_s2_f[%d][%d]= %f\n", ii,jj, cc2006_s2_formula[ii][jj]);
				}
			}
		}
		for (ii=0; ii < 4; ii++) {
			for (jj=0; jj < CC2006_STEMSIZE_FORMULA; jj++) {
				if (cc2006_s1_formula[ii][jj] < INF)
				{
					printf("c_s1_f[%d][%d]= %f\n", ii,jj, cc2006_s1_formula[ii][jj]);
				}
			}
		}
	}
*/	

	float temp = 37.0 + 273.15;  // TODO

	// read the simfold and DP parameters

	fill_data_structures_with_new_parameters_PK_DP(parameters_file);

	FILE *file;
	char buffer[100];
	char  v1[10];
	int line = 0;
	double param;

	//printf ("FILENAME: %s\n", filename);
	if ((file = fopen (parameters_file, "r")) == NULL)
	{
		giveupPK ( (char *) "Cannot open file", parameters_file);
	}

	// skip the simfold and DP parameters in the file
	for (int i = 0; i < get_num_params_PK_DP(); i++)
	{
		fgets(buffer, sizeof(buffer), file);
		line++;
	}

	if (DEBUG)
		printf("Skipped %d lines\n", line);

	// read the coaxial stacking parameters
	// assumes parameters are in this order:
	// coaxstack_f_b
	// coaxstack_m1
	// coaxstack_m2

	int i,j,k,l;
	for (i=0; i < NUCL; i++) {
		for (j=0; j < NUCL; j++) {
			for (k=0; k < NUCL; k++) {
				for (l=0; l < NUCL; l++)
				{
					if (coaxstack_f_b[i][j][k][l] < INF)
					{
						// no duplicates here
						fgets (buffer, sizeof(buffer), file);
						line++;
						sscanf (buffer, "%s", v1);
						coaxstack_f_b[i][j][k][l] = ascii_to_intPK (v1);
//						if (DEBUG)
//							printf("read: c_f_b[%d][%d][%d][%d] =%d\n", i,j,k,l,coaxstack_f_b[i][j][k][l]);
					}
				}
			}
		}
	}
	for (i=0; i < NUCL; i++) {
		for (j=0; j < NUCL; j++) {
			for (k=0; k < NUCL; k++) {
				for (l=0; l < NUCL; l++)
				{
					if (coaxstack_m1[i][j][k][l] < INF)
					{
						fgets (buffer, sizeof(buffer), file);
						line++;
						sscanf (buffer, "%s", v1);
						coaxstack_m1[i][j][k][l] = ascii_to_intPK (v1);
//                                              if (DEBUG)
//                                                      printf("read: c_m1[%d][%d][%d][%d] =%d\n", i,j,k,l,coaxstack_m1[i][j][k][l]);
					}
				}
			}
		}
	}
	for (i=0; i < NUCL; i++) {
		for (j=0; j < NUCL; j++) {
			for (k=0; k < NUCL; k++) {
				for (l=0; l < NUCL; l++)
				{
					if (coaxstack_m2[i][j][k][l] < INF)
					{
						fgets (buffer, sizeof(buffer), file);
						line++;
						sscanf (buffer, "%s", v1);
						coaxstack_m2[i][j][k][l] = ascii_to_intPK (v1);
//                                                if (DEBUG)
//                                                        printf("read: c_m2[%d][%d][%d][%d] =%d\n", i,j,k,l,coaxstack_m2[i][j][k][l]);
					}
				}
			}
		}
	}

	if (DEBUG)
		printf("Read %d lines\n", line);

	// read the CC2006b energy model parameters
	// assumes the parameters are in this order:
	// cc2006_s2_l1, cc2006_s1_l2, cc2006_s2_formula, cc2006_s1_formula

	for (i=0; i < CC2006_STEMSIZE; i++) {
		for (j=0; j < CC2006_LOOPSIZE; j++) {
			if (cc2006_s2_l1[i][j] < INF)
			{
				fgets (buffer, sizeof(buffer), file);
				line++;
				sscanf (buffer, "%s", v1);
				cc2006_s2_l1[i][j] = ascii_to_CCdouble_PK(v1,KB*temp);
			}
		}
	}
	for (i=0; i < CC2006_STEMSIZE; i++) {
		for (j=0; j < CC2006_LOOPSIZE; j++) {
			if (cc2006_s1_l2[i][j] < INF)
			{
				fgets (buffer, sizeof(buffer), file);
				line++;
				sscanf (buffer, "%s", v1);
				cc2006_s1_l2[i][j] = ascii_to_CCdouble_PK(v1,KB*temp);
			}
		}
	}
	for (i=1; i < 4; i++) {
		for (j=0; j < CC2006_STEMSIZE_FORMULA; j++) {
			if (cc2006_s2_formula[i][j] < INF)
			{
				fgets (buffer, sizeof(buffer), file);
				line++;
				sscanf (buffer, "%lf", &param);
				cc2006_s2_formula[i][j] = (float) param;
			}
		}
	}
	for (i=1; i < 4; i++) {
		for (j=0; j < CC2006_STEMSIZE_FORMULA; j++) {
			if (cc2006_s1_formula[i][j] < INF)
			{
				fgets (buffer, sizeof(buffer), file);
				line++;
				sscanf (buffer, "%lf", &param);
				cc2006_s1_formula[i][j] = (float) param;
			}
		}
	}

		// Mirela: added dG_assemble as a parameter on Dec 1, 2008
		//fgets (buffer, sizeof(buffer), file);
		//if (!feof (file))
		//{
		//    line++;
		//    sscanf (buffer, "%lf", &param);
		//    pkmodelCC2006.deltaG_assemble = (float) param*100.0;        
		//}


	fclose (file);
	//printf ("****** we must have 11 lines by now: LINES = %d\n", line);

	create_string_params_PK_CC();  // update the string array of parameters

/*
	if (DEBUG2){
		printf("AFTER NEW PARAMS\n");
		for (ii=0; ii < NUCL; ii++) {
				for (jj=0; jj < NUCL; jj++) {
						for (kk=0; kk < NUCL; kk++) {
								for (ll=0; ll < NUCL; ll++){
										if (coaxstack_f_b[ii][jj][kk][ll] < INF){
											    printf("c_f_b[%d][%d][%d][%d]= %d\n",ii,jj,kk,ll,coaxstack_f_b[ii][jj][kk][ll]);    
										}
								}
						}
				}        
		} 
		for (ii=0; ii < NUCL; ii++) {
				for (jj=0; jj < NUCL; jj++) {
						for (kk=0; kk < NUCL; kk++) {
								for (ll=0; ll < NUCL; ll++){
										if (coaxstack_m1[ii][jj][kk][ll] < INF){
											    printf("c_m1[%d][%d][%d][%d]= %d\n",ii,jj,kk,ll,coaxstack_m1[ii][jj][kk][ll]);
										}
								}
						}
				}
		}
		for (ii=0; ii < NUCL; ii++) {
				for (jj=0; jj < NUCL; jj++) {
						for (kk=0; kk < NUCL; kk++) {
								for (ll=0; ll < NUCL; ll++){
										if (coaxstack_m2[ii][jj][kk][ll] < INF){
											    printf("c_m2[%d][%d][%d][%d]= %d\n",ii,jj,kk,ll,coaxstack_m2[ii][jj][kk][ll]);
										}
								}
						}
				}
		}
		for (ii=0; ii < CC2006_STEMSIZE; ii++) {
			for (jj=0; jj < CC2006_LOOPSIZE; jj++) {
				if (cc2006_s2_l1[ii][jj] < INF){
					printf("c_s2_l1[%d][%d]= %f\n", ii,jj, cc2006_s2_l1[ii][jj]);
				}
			}
		}
		for (ii=0; ii < CC2006_STEMSIZE; ii++) {
			for (jj=0; jj < CC2006_LOOPSIZE; jj++) {
				if (cc2006_s1_l2[ii][jj] < INF){
					printf("c_s1_l2[%d][%d]= %f\n", ii,jj, cc2006_s1_l2[ii][jj]);
				}
			}
		}
		for (ii=0; ii < 4; ii++) {
			for (jj=0; jj < CC2006_STEMSIZE_FORMULA; jj++) {
				if (cc2006_s2_formula[ii][jj] < INF){
					printf("c_s2_f[%d][%d]= %f\n", ii,jj, cc2006_s2_formula[ii][jj]);
				}
			}
		}
		for (ii=0; ii < 4; ii++) {
			for (jj=0; jj < CC2006_STEMSIZE_FORMULA; jj++) {
				if (cc2006_s1_formula[ii][jj] < INF){
					printf("c_s1_f[%d][%d]= %f\n", ii,jj, cc2006_s1_formula[ii][jj]);
				}
			}
		}
	}
*/
}

// USED FOR PARAMETER TUNING (CC2006b energy model)
// returns the free energy (in kcal/mol) of sequence folded into structure
double free_energy_PK_CC2006b(char *sequence, char *structure)
// sequence: sequence (input parameter)
// structure: secondary structure in dot-bracket format () and [] (input parameter)
{
	// create an array of shorts which denotes the index of the base pairs:
	// the first base in the sequence has index 1
	// pairseq[0] = 0 always;
	// pairseq[1] = i where i is the index of the base paired with base 1, etc;
	// pairseq[j] = 0 if base j is unpaired;

	int size = strlen (structure);
	short pairseq[size+1];
	detect_original_PKed_pairs_many(structure, pairseq);

	ReadInput* R = new ReadInput(size, sequence, pairseq);

//	printf("DEBUG = After readinput\n");

	Stack * s = new Stack(R);
	Bands * B = new Bands(R, s);

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

	short* secstructure = new short[R->Size+1];
	char *csequence = new char[R->Size+2];
		csequence[R->Size+1] = '\0';
	for (int i = 0; i < R->Size+1; i++) {
		secstructure[i] = (short)(R->Sequence[i]);
		if (secstructure[i] == -1) secstructure[i] = 0;
	  	csequence[i] = R->CSequence[i];
		if (DEBUG)
			 	printf("%d %c %d \n", i, csequence[i], secstructure[i]);
	}

//	PlotRna(prefix, &sequence[1], &secstructure[1], outPSFile, L->Energy());

//	printf("Energy Model The total free energy is %f cal/mol\n", L->Energy());

	float totalEnergy = 0;

	if (DEBUG)
	{
	cout << setw(15) << left << "Energy Model" << setw(25) << left << "Free Energy (kcal/mol)" << setw(40) << left << "Free Energy without Dangling (kcal/mol)" << endl;
	printf("--------------------------------------------------------------\n");
	}

	// energy returned in kcal

	totalEnergy = -L->Energy(CC2006b)/1000;
	float totalEnergyDang = -L->EnergyDangling()/1000;

	if (DEBUG)
		cout << setw(15) << left << "Cao&Chen(b)" << setw(25) << left << totalEnergy + totalEnergyDang << setw(40) << left << totalEnergy << endl;

//	L->printEnergyTrace();

//	if (no_pk_dangling_ends == 0)
		return totalEnergy + totalEnergyDang;
//	else
//		return totalEnergy;
}

// USED FOR PARAMETER TUNING (CC2006b energy model)
double get_feature_counts_quadratic_PK_CC2006b (char *sequence, char *structure, double ** quadratic_matrix, double *counter, double &free_value)
// sequence: sequence (input parameter)
// structure: secondary structure in dot-bracket format () and [] (input parameter)
// counter: array where counter[i] is the number of times the i-th feature occurs (output parameter)
// quadratic_matrix: TODO
// free_value: TODO
// Note: The counter, free_value, and quadratic_matrix are automatically reset before passed onto other functions
{
	// create an array of shorts which denotes the index of the base pairs:
	// the first base in the sequence has index 1
	// pairseq[0] = 0 always;
	// pairseq[1] = i where i is the index of the base paired with base 1, etc;
	// pairseq[j] = 0 if base j is unpaired;

	int size = strlen (structure);
	short pairseq[size+1];
	detect_original_PKed_pairs_many(structure, pairseq);

	ReadInput* R = new ReadInput(size, sequence, pairseq);

	Stack * s = new Stack(R);
	Bands * B = new Bands(R, s);

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


	short* secstructure = new short[R->Size+1];
	char *csequence = new char[R->Size+2];
		csequence[R->Size+1] = '\0';
	for (int i = 0; i < R->Size+1; i++) {
		secstructure[i] = (short)(R->Sequence[i]);
		if (secstructure[i] == -1) secstructure[i] = 0;
	  	csequence[i] = R->CSequence[i];
		if (DEBUG)
			 	printf("%d %c %d \n", i, csequence[i], secstructure[i]);
	}

//	PlotRna(prefix, &sequence[1], &secstructure[1], outPSFile, L->Energy());

//	printf("Energy Model The total free energy is %f cal/mol\n", L->Energy());

	float totalEnergy = 0;

	if (DEBUG)
	{
		cout << setw(15) << left << "Energy Model" << setw(25) << left << "Free Energy (kcal/mol)" << setw(40) << left << "Free Energy without Dangling (kcal/mol)" << endl;
		printf("--------------------------------------------------------------\n");
	}

	// energy returned in kcal
	
	// PARAMETER TUNING
	// clear the counter, quadratic_matrix, and free_value

	int num_params = get_num_params_PK_CC2006b();
		free_value = 0;
	for (int i=0; i < num_params; i++)
	{
		counter[i] = 0;
		for (int j=i; j < num_params; j++)
			quadratic_matrix[i][j] = 0;
	}

	int reset_c = 0;
	int ignore_dangles = no_pk_dangling_ends;
	int ignore_AU = 0;  // 0 = do include AU penalties

	totalEnergy = -L->Energy(CC2006b, quadratic_matrix, counter, free_value, reset_c, ignore_dangles);
	float totalEnergyDang = -L->EnergyDangling(CC2006b, quadratic_matrix, counter,free_value,reset_c,ignore_dangles,ignore_AU);

	// CHECK VALUES OF COUNTER, ETC
	//int num_params = get_num_params_PK_DP();
	//int num_params_pkfree = get_num_params();
	//printf("Free Value: %f\n", free_value);
	//printf("Counter Values:\n");
	//for (int i = num_params_pkfree; i < num_params; i++)
	//	printf("c[%d]=%f  ", i, counter[i]);

	//printf("Some P_matrix Values:\n");
	//for (int i = 0; i < num_params/10; i++)
	//	printf("P[%d][%d]=%f  ", num_params_pkfree + structure_type_index_PK("stp")-1, i, quadratic_matrix[i][num_params_pkfree + structure_type_index_PK("stp")-1]);


	if (DEBUG)
	{
	cout << setw(15) << left << "Cao&Chen(b)" << setw(25) << left << totalEnergy + totalEnergyDang << setw(40) << left << totalEnergy << endl;
	}

	// return the free energy:
	//    - deltaG = x' P x + c' x + f
	//        - where x is the vector of parameters
	//        - P is a symmetric matrix of the coefficients for each quadratic term
	//        - c is a vector of counts for each linear term
	//        - c' means c transposed
	//        - f is a constant

	return totalEnergy + totalEnergyDang;

	//int num_params = get_num_params_PK_DP();
	//double * energy_temp1 = new double[num_params];
	//double * energy_temp2 = new double[num_params];
	//double energy = free_value;

	//for (int i = 0; i < num_params; i++)
	//{
	//	energy_temp1[i] = counter[i] * params_all[i];
	//	for (int j = 0; j < num_params; j++)
	//	{
	//		if (i <= j)
	//			energy_temp2[i] += quadratic_matrix[i][j] * params_all[j];
	//		else  // since only the upper triangle of the matrix is filled out
	//			energy_temp2[i] += quadratic_matrix[j][i] * params_all[j];
	//	}
	//}
	//for (int i = 0; i < num_params; i++)
	//{
	//	energy += params_all[i]*energy_temp2[i] + energy_temp1[i];
	//}

	//return energy;
}


int get_num_params_PK_CC2006b()
{
	// TODO - try to reduce the number of parameters - i.e. are there symmetries in the tables?
	// this must return the whole parameter counts, i.e. both from simfold and the ones in CC model + coaxial stacking
	return 546 + get_num_params_PK_DP();   // + 1;
	// Mirela: added 1 for deltaG_assemble
}
