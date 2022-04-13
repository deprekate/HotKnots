/***************************************************************************
  init.cpp  -  description
  -------------------
begin                : Thu Apr 11 2002
copyright            : (C) 2002 by Mirela Andronescu
email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


// this file contains functions to read the thermodynamic parameters from files
//    into internal data structures

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "initPK.h"
#include "constants.h"
#include "structsPK.h"
#include "globalsPK.h"  // July 16 - necessary because this is the only time globalsPK.h will be included
#include "externs.h"  // July 16 - added instead of globals.h which is included in simfold files
#include "commonPK.h"
//#include "globals.h"  // July 16 - removed
//#include "common.h"  // July 16 - removed - unnecessary; August 5 - needed for giveup function


// from simfold/init.cpp
void read_configuration_filePK (char *filename){
	// PRE:  None
	// POST: read the configuration file, which must be in the standard format
	//      - see the documentation
	char buffer[256], token1[50], token2[5], token3[50];
	FILE *file;
	int i;
	memset(buffer, 0, sizeof(buffer));
	memset(token1, 0, sizeof(token1));
	memset(token2, 0, sizeof(token2));
	memset(token3, 0, sizeof(token3));

	if ((file = fopen (filename,"r")) == NULL){
		giveupPK( (char *) "1aCannot open file", filename);
	}
	fgets(buffer, sizeof(buffer), file);
	while (!feof(file)){
		if (buffer[0] == '#' || buffer[0] == '[' || buffer[0] == '\n'){
			fgets(buffer, sizeof(buffer), file);
			continue;
		}
		sscanf(buffer, "%s%s%s", token1, token2, token3);
		for (i=0; i<nb_params; i++){
			if (strcmp(token1, par_namePK[i]) == 0){
				strcpy(par_valuePK[i], token3);
				break;
			}    
		}
		fgets(buffer, sizeof(buffer), file);
	}
	fclose(file);
}


void read_pkmodelDP_file (char* filename, pkmodelinfoDP &pkmodelDP){
	// initialize the pseudoknotted energy model parameters (Dirks & Pierce)
	char buffer[256];
	FILE *file;
	char  v1[10];
	if ((file = fopen (filename,"r")) == NULL){
		giveupPK ( (char *) "Cannot open file", filename);
	}

	fgets (buffer, sizeof(buffer), file);
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.Ps = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.Psm = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.Psp = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.Pb = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	// efn2 multibranched loops - not in use
	sscanf (buffer, "%s", v1);
	pkmodelDP.Pup = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.Pps = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.stP = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);            
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.intP = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.a = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.b = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.c = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.a_p = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.b_p = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelDP.c_p = ascii_to_doublePK (v1);

}

void read_pkmodelRE_file (char* filename, pkmodelinfoRE &pkmodelRE)
	// initialize the pseudoknotted energy model parameters (Rivas & Eddy)
{
	char buffer[256];
	FILE *file;
	char  v1[10];
	if ((file = fopen (filename,"r")) == NULL)
	{
		giveupPK ( (char *) "Cannot open file", filename);
	}

	fgets (buffer, sizeof(buffer), file);
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.g_interiorPseudo = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.P_tilda = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.P_i = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.Q_tilda = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	// efn2 multibranched loops - not in use
	sscanf (buffer, "%s", v1);
	pkmodelRE.M_tilda = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.Gw = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.Gwh = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);            
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.Gwi = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);            
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.q_unpairedMultiPseudo = ascii_to_doublePK (v1);

	fgets (buffer, sizeof(buffer), file);    
	while (buffer[0] == '#' || buffer[0] == '\n')
		fgets (buffer, sizeof(buffer), file);
	sscanf (buffer, "%s", v1);
	pkmodelRE.p_pairedMultiPseudo = ascii_to_doublePK (v1);
}

void read_pkmodelCC2006_file(char *filename, float temp, float cc2006_s2_l1 [CC2006_STEMSIZE] [CC2006_LOOPSIZE], 
		float cc2006_s1_l2 [CC2006_STEMSIZE] [CC2006_LOOPSIZE],
		float cc2006_s2_formula [4] [CC2006_STEMSIZE_FORMULA],
		float cc2006_s1_formula [4] [CC2006_STEMSIZE_FORMULA])
	// PRE:  filename is the file that contains the pkmodel information for RNA;
	//       temp is in Kelvins
	// POST: the values are stored in the given arrays
{
	char buffer[256];
	FILE *file;
	char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
	char  v9[10], v10[10], v11[10], v12[10];
	int i, j, ii, jj;       // ii = i+1; jj = j-1

	i=0; j=0; ii=0; jj=0;
	int counter;
	if ((file = fopen (filename,"r")) == NULL)
	{
		giveupPK ( (char *) "Cannot open file", filename);
	}
	fgets (buffer, sizeof(buffer), file);
	i = 0;
	counter = 0;

	// loop through the 11 rows of the top half of Table 1
	while (counter < CC2006_STEMSIZE && !feof(file))  
	{
		if (buffer[0] == '#' || buffer[0] == '\n')
		{
			fgets (buffer, sizeof(buffer), file);
			continue;
		}

		// *** Add new energy model code here (if Cao&Chen changes i.e. CC2006_LOOPSIZE changes)
		sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s",		
				v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12);
		ii = 0;
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v1,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 1; array value in 10cal/mol
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v2,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 2
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v3,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 3
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v4,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 4
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v5,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 5
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v6,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 6
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v7,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 7
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v8,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 8
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v9,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 9
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v10,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 10
		cc2006_s2_l1[i][ii++] = ascii_to_CCdouble_PK(v11,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 11
		cc2006_s2_l1[i][ii  ] = ascii_to_CCdouble_PK(v12,KB*temp);  // S2 = 2 initially, i afterwards, L1 = 12
		fgets (buffer, sizeof(buffer), file);
		i++;  // the size of the stem (index of row)
		counter++;

		//printf("ascii_to_intPK(v3) = %d\n", ascii_to_intPK(v3));
	}

	i = 0;
	counter = 0;
	// loop through the 11 rows of the bottom half of Table 1    
	while (counter < CC2006_STEMSIZE && !feof(file))
	{
		if (buffer[0] == '#' || buffer[0] == '\n')
		{
			fgets (buffer, sizeof(buffer), file);
			continue;
		}

		// *** Add new energy model code here (if Cao&Chen changes i.e. CC2006_LOOPSIZE changes)
		sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s",		
				v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12);
		ii = 0;
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v1,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 1
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v2,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 2
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v3,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 3
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v4,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 4
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v5,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 5
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v6,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 6
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v7,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 7
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v8,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 8
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v9,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 9
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v10,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 10
		cc2006_s1_l2[i][ii++] = ascii_to_CCdouble_PK(v11,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 11
		cc2006_s1_l2[i][ii  ] = ascii_to_CCdouble_PK(v12,KB*temp);  // S1 = 2 initially, i afterwards, L2 = 12
		fgets (buffer, sizeof(buffer), file);
		i++;  // the size of the stem (index of row)
		counter++;

		//printf("table2 ascii_to_intPK(v3) = %d\n", ascii_to_intPK(v3));
	}

	//	printf("INITPK: cc2006_s1_l2[1][2] = %f\n", cc2006_s1_l2[1][2]);

	i = 0;
	counter = 0;
	// loop through the 4 rows of the top half of Table 2
	while (counter < 4 && !feof(file))
	{
		if (buffer[0] == '#' || buffer[0] == '\n')
		{
			fgets (buffer, sizeof(buffer), file);
			continue;
		}

		// *** Add new energy model code here (if Cao&Chen changes i.e. CC2006_STEMSIZE_FORMULA changes)
		sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s",		
				v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11);
		ii = 0;
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v1);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v2);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v3);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v4);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v5);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v6);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v7);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v8);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v9);  
		cc2006_s2_formula[i][ii++] = ascii_to_doublePK(v10); 
		cc2006_s2_formula[i][ii  ] = ascii_to_doublePK(v11); 
		fgets (buffer, sizeof(buffer), file);
		i++;  // the size of the stem (index of row)
		counter++;
	}

	i = 0;
	counter = 0;
	// loop through the 4 rows of the bottom half of Table 2
	while (counter < 4 && !feof(file))
	{
		if (buffer[0] == '#' || buffer[0] == '\n')
		{
			fgets (buffer, sizeof(buffer), file);
			continue;
		}

		// *** Add new energy model code here (if Cao&Chen changes i.e. CC2006_STEMSIZE_FORMULA changes)
		sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s",		
				v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11);
		ii = 0;
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v1);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v2);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v3);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v4);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v5);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v6);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v7);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v8);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v9);  
		cc2006_s1_formula[i][ii++] = ascii_to_doublePK(v10); 
		cc2006_s1_formula[i][ii  ] = ascii_to_doublePK(v11); 
		fgets (buffer, sizeof(buffer), file);
		i++;  // the size of the stem (index of row)
		counter++;
	}

	fclose (file);
}

void read_coaxstack_file (char *filename, int stack[][NUCL][NUCL][NUCL])
	// PRE:  filename is the file that contains the stacked pairs energy data for DNA    
	// POST: the values are stored in the 4-D array stack
{
	//	printf("filename = %s\n", filename);

	char buffer[256];
	FILE *file;
	char  v1[10],  v2[10],  v3[10],  v4[10],  v5[10],  v6[10],  v7[10],  v8[10];
	char  v9[10], v10[10], v11[10], v12[10], v13[10], v14[10], v15[10], v16[10];
	int i, j, ii, jj;       // ii = i+1; jj = j-1

	i=0; j=0; ii=0; jj=0;       
	if ((file = fopen (filename,"r")) == NULL)
	{
		giveupPK ( (char *) "Cannot open file", filename);
	}
	fgets (buffer, sizeof(buffer), file);
	while (!feof(file))
	{
		if (buffer[0] == '#' || buffer[0] == '\n')
		{
			fgets (buffer, sizeof(buffer), file);
			continue;
		}
		sscanf (buffer, "%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
				v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16);
		if (ii == NUCL)
		{
			ii = 0; i++;
		}
		// v1-v4 - jj changes
		jj = 0; j = 0;
		stack[i][j][ii][jj++] = ascii_to_intPK(v1);
		stack[i][j][ii][jj++] = ascii_to_intPK(v2);
		stack[i][j][ii][jj++] = ascii_to_intPK(v3);
		stack[i][j][ii][jj  ] = ascii_to_intPK(v4);
		// j grows, jj starts again from 0
		jj = 0; j++;
		stack[i][j][ii][jj++] = ascii_to_intPK(v5);
		stack[i][j][ii][jj++] = ascii_to_intPK(v6);
		stack[i][j][ii][jj++] = ascii_to_intPK(v7);
		stack[i][j][ii][jj  ] = ascii_to_intPK(v8);
		// j grows, jj starts again from 0
		jj = 0; j++;
		stack[i][j][ii][jj++] = ascii_to_intPK(v9);
		stack[i][j][ii][jj++] = ascii_to_intPK(v10);
		stack[i][j][ii][jj++] = ascii_to_intPK(v11);
		stack[i][j][ii][jj  ] = ascii_to_intPK(v12);
		// j grows, jj starts again from 0
		jj = 0; j++;
		stack[i][j][ii][jj++] = ascii_to_intPK(v13);
		stack[i][j][ii][jj++] = ascii_to_intPK(v14);
		stack[i][j][ii][jj++] = ascii_to_intPK(v15);
		stack[i][j][ii][jj++] = ascii_to_intPK(v16);
		// go to the next line: grow ii
		ii++;
		fgets (buffer, sizeof(buffer), file);
	}
	fclose (file);

	//   	printf("stack: coaxstack_m2[2][0][1][2] = %d\n", stack[2][0][1][2]);
	//    printf("end of stack function\n");

}


// similar to simfold/init.cpp
void init_dataPK (char *arg, char *config_file, int what, double temperature)
	//void init_data(char *arg, char *config_file, int what, double temperature)
	// the function that must be called by the main program to read data files
	// PRE:  None
	// POST: Read all data and configuration files
{
	char conf[200];

	char pkmodelDP_filename[200];
	char pkmodelRE_filename[200];
	char pkmodelCC2006_filename[200];
	char coaxstack_f_a_energy37_filename[200];
	char coaxstack_f_b_energy37_filename[200];
	char coaxstack_m1_energy37_filename[200];
	char coaxstack_m2_energy37_filename[200];
	// *** Add new energy model code here

	// configuration file
	// first find the path to params
	int i;
	
	if (what != RNA && what != DNA){
		printf ("Please specify what to fold: RNA or DNA\n");
		exit(1);
	}
	if (temperature < 0 || temperature > 100){
		printf ("Temperature must be between 0 and 100 degrees Celsius\n");
		exit(1);
	}

	int len = strlen(config_file);
	char* dir = (char*) malloc(len+1);
	strcpy(dir, config_file);
	while(len > 0){
		len--;
		if(dir[len] == '\\' || dir[len] == '/'){
			dir[len+1] = '\0';
			break;
		}
	}
	read_configuration_filePK(config_file);
	strcpy(std_dir_par, dir);

	// if temperature is 37, there's no need to read the enthalpies
	if (what == RNA){
		strcpy (pkmodelDP_filename, std_dir_par);
		strcat (pkmodelDP_filename, rna_pkmodelDP_par); 
		strcpy (pkmodelRE_filename, std_dir_par);
		strcat (pkmodelRE_filename, rna_pkmodelRE_par);
		strcpy (pkmodelCC2006_filename, std_dir_par);
		strcat (pkmodelCC2006_filename, rna_pkmodelCC2006_par);
		strcpy (coaxstack_f_a_energy37_filename, std_dir_par);
		strcat (coaxstack_f_a_energy37_filename, rna_coaxstack_f_a_par);
		strcpy (coaxstack_f_b_energy37_filename, std_dir_par);
		strcat (coaxstack_f_b_energy37_filename, rna_coaxstack_f_b_par);
		strcpy (coaxstack_m1_energy37_filename, std_dir_par);
		strcat (coaxstack_m1_energy37_filename, rna_coaxstack_m1_par);
		strcpy (coaxstack_m2_energy37_filename, std_dir_par);
		strcat (coaxstack_m2_energy37_filename, rna_coaxstack_m2_par);
		// *** Add new energy model code here
	}
	else if (what == DNA)
	{
		strcpy (pkmodelDP_filename, std_dir_par);
		strcat (pkmodelDP_filename, dna_pkmodelDP_par);
		strcpy (pkmodelRE_filename, std_dir_par);
		strcat (pkmodelRE_filename, dna_pkmodelRE_par);
		strcpy (pkmodelCC2006_filename, std_dir_par);
		strcat (pkmodelCC2006_filename, dna_pkmodelCC2006_par);
		strcpy (coaxstack_f_a_energy37_filename, std_dir_par);
		strcat (coaxstack_f_a_energy37_filename, dna_coaxstack_f_a_par);
		strcpy (coaxstack_f_b_energy37_filename, std_dir_par);
		strcat (coaxstack_f_b_energy37_filename, dna_coaxstack_f_b_par);
		strcpy (coaxstack_m1_energy37_filename, std_dir_par);
		strcat (coaxstack_m1_energy37_filename, dna_coaxstack_m1_par);
		strcpy (coaxstack_m2_energy37_filename, std_dir_par);
		strcat (coaxstack_m2_energy37_filename, dna_coaxstack_m2_par);
		// *** Add new energy model code here
	}

	//	printf("coaxstack_f_a_energy37_filename = %s\n", coaxstack_f_a_energy37_filename);
	//	printf("coaxstack_f_b_energy37_filename = %s\n", coaxstack_f_b_energy37_filename);
	//	printf("coaxstack_m1_energy37_filename = %s\n", coaxstack_m1_energy37_filename);
	//	printf("coaxstack_m2_energy37_filename = %s\n", coaxstack_m2_energy37_filename);

	read_pkmodelDP_file (pkmodelDP_filename, pkmodelDP);  // initialize the pseudoknotted energy model parameters
	read_pkmodelRE_file (pkmodelRE_filename, pkmodelRE);  // initialize the pseudoknotted energy model parameters
	read_pkmodelCC2006_file (pkmodelCC2006_filename, temperature + 273.15, 
			cc2006_s2_l1, cc2006_s1_l2, cc2006_s2_formula, cc2006_s1_formula);  // initialize the pseudoknotted energy model parameters
	read_coaxstack_file(coaxstack_f_a_energy37_filename, coaxstack_f_a);  // initialize the coaxial stacking parameters (as described by Cao & Chen)
	read_coaxstack_file(coaxstack_f_b_energy37_filename, coaxstack_f_b);  // initialize the coaxial stacking parameters (as described by Mathews)
	read_coaxstack_file(coaxstack_m1_energy37_filename, coaxstack_m1);  // initialize the coaxial stacking parameters (as described by Mathews)
	read_coaxstack_file(coaxstack_m2_energy37_filename, coaxstack_m2);  // initialize the coaxial stacking parameters (as described by Mathews)

	// *** Add new energy model code here

	// for Cao & Chen energy model
	pkmodelCC2006.deltaG_assemble = KB * log (9) * (temperature + 273.15) * 100;  // in 10cal/mol
	pkmodelCC2006.temp = temperature + 273.15;

	//	printf("INITPK: cc2006_s1_l2[1][2] = %f\n", cc2006_s1_l2[1][2]);

	//	printf("INITPK: coaxstack_f_a[0][3][2][1] = %d\n", coaxstack_f_a[0][3][2][1]);
	//	printf("INITPK: coaxstack_m2[2][0][1][2] = %d\n", coaxstack_m2[2][0][1][2]);

	create_string_params_PK_CC();  // update the string array of parameters
}
