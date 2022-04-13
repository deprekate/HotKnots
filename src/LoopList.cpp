/*****************************************************************
HotKnot: A heuristic algorithm for RNA secondary 
structure prediction including pseudoknots
File: LoopList.cpp
Description: 
Contains functions for calculating the free energy
of interior-pseudoknotted and multi-pseudoknotted loops.
Also contains a function (FindChildren) for figuring out the tuples 
associated with a multi-pseudokntted loop.


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

#include "LoopList.h"
#include "Loop.h"
#include "commonPK.h"

//#define NULL 0


LoopList::LoopList(ReadInput* R, int b1, int b2){
	base1 = b1;
	base2 = b2;
	Size = 0;
	Input = R;
}

LoopList::~LoopList(){
}

/*********************************************************************************
FincChildren: figures out the tuples associated with a multi-pseudoknotted loop.
 *********************************************************************************/
void LoopList::FindChildren(){
	Input->loops[base1].num_branches = 0;
	Input->loops[base1].pseudo_num_branches = 0;

	int i = Input->Next[base1];
	int j = base2;

	while (i < j){  // check all bases i from starting base pair of multiloop to ending base pair
		if (Input->ClosedRegions[i] != NULL){  // if base i starts a closed region
			// set 'i'th branch to start of 'i'th closed region in multiloop
			Input->loops[base1].bri[Input->loops[base1].num_branches++] = i;
			if (Input->ClosedRegions[i]->type == pseudo)
				Input->loops[base1].pseudo_num_branches += 2;
			else
				Input->loops[base1].pseudo_num_branches += 1;
			// go to start of next closed region (next base pair after end of current closed region)
			i = Input->Next[Input->ClosedRegions[i]->end];
		}else{
			i = Input->Next[i];  // if i is not the start of a closed region, go to the next paired base
		}
	};

	// CHECK: why is this repeated and what does it do (this one starts from Next(base1))

	i = Input->Next[Input->BasePair(base2)];
	j = Input->BasePair(base1);
	while (i < j){
		if (Input->ClosedRegions[i] != NULL){
			Input->loops[base1].bri[Input->loops[base1].num_branches++] = i;
			if (Input->ClosedRegions[i]->type == pseudo)
				Input->loops[base1].pseudo_num_branches += 2;
			else
				Input->loops[base1].pseudo_num_branches += 1;
			i = Input->Next[Input->ClosedRegions[i]->end];
		}else{
			i = Input->Next[i];
		}
	};
};


///////////////////////////////////////////////////////////////////
// DIRKS & PIERCE ENERGY MODEL
///////////////////////////////////////////////////////////////////


/*********************************************************************************
stackPseudoEnergy: Calculates the free energy of an interior-pseudoknotted loop which
is acutally a stacked pair. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float	LoopList::stackPseudoEnergyDP(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = base1;

	float en = float(LEstacked_pair_energy (i, f[i].pair, sequence));

	/* include if non-canonical base pairs must be stripped from input structure
	   if ((sequence[i] == A && sequence[f[i].pair] == G) ||
	   (sequence[i] == G && sequence[f[i].pair] == A) ||
	   (sequence[i] == A && sequence[f[i].pair] == C) ||
	   (sequence[i] == C && sequence[f[i].pair] == A) ||
	   (sequence[i] == C && sequence[f[i].pair] == U) ||
	   (sequence[i] == U && sequence[f[i].pair] == C) )
	   {
	   printf("+++++");
	   return 0;
	   }
	   */

	g_count_stP++;

	return pkmodelDP.stP * en;
}


/*********************************************************************************
interiorPseudoEnergy: Calculates the free energy of an interior-pseudoknotted loop
which is not a stacked pair. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float	LoopList::interiorPseudoEnergyDP(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;

	if (  (base2 == base1 + 1 ) && (f[base1].pair == f[base2].pair + 1) )
		return stackPseudoEnergyDP();

	int i = base1;

	int ip, jp;
	ip = base2;
	jp = f[ip].pair;

	Input->cannot_add_dangling[ip-1] = 1;
	Input->cannot_add_dangling[jp+1] = 1;

	float en = float(LEinternal_loop_energy (i, f[i].pair, ip, jp, sequence));

	g_count_intP++;

	return pkmodelDP.intP * en;
};

/*********************************************************************************
multiPseudoEnergy:  Calculates the free energy of a multi-pseudoknotted Loop.
Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float	LoopList::multiPseudoEnergyDP(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = base1;
	PARAMTYPE AUpen;
	int h, l;
	PARAMTYPE dang;
	double misc_energy;

	dang = 0;
	misc_energy = 0;
	AUpen = 0;
	int special;
	special = 0;

	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int begin = base1; int end = Input->BasePair(base1);
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++){
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	// add penalty for introducing a multiloop that spans a band
	misc_energy += pkmodelDP.a_p;

	g_count_a_p++;

	//	misc_energy += pkmodelDP.Pps * (f[i].num_branches + 2);
	// +1 for the pseudoknotted base pair, +1 for the closing base pair of the multiloop that spans a band
	misc_energy += pkmodelDP.b_p * (f[i].num_branches + 2);
	//	misc_energy += pkmodelDP.b_p * f[i].pseudo_num_branches;  // use this is pseudoloops count as 2 branches
	//	misc_energy += misc.multi_helix_penalty * f[i].pseudo_num_branches  + 2 * p_pairedMultiPseudo;

	g_count_b_p += f[i].num_branches + 2;

	// count the number of unpaired bases inside the multiloop
	int NumUnpairedInMulti = 0;
	int currentNumBranches = 0;  // the branch we are currently looking at

	// if there are no branches between base1 and the pseudoknotted base pair...
	if ( base2 <= Input->loops[base1].bri[0]){
		for (l = base1 + 1; l < base2; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;

			NumUnpairedInMulti++;
		}
	} else {
		// count the number of unpaired bases from base1 to start of first branch
		for (l=base1+1; l < Input->loops[base1].bri[0]; l++){
			misc_energy += pkmodelDP.c_p;
			g_count_c_p++;
			NumUnpairedInMulti++;
		}
		currentNumBranches++;  // currentNumBranches traversed is now 1
	}

	// DONECHECK: this doesn't count unpaired bases inside a branch of the multiloop !!!!!! Correct

	// DONECHECK: what is the value of num_branches? is it double because of what's happening in FindChildren above?
	// DONE: num_branches does not include the branch that has pseudoknotted base pairs

	h = 0;
	// count number of unpaired bases between end of first branch and beginning of next branch, up the 
	// branch that has the pseudoknotted base pair (starting at base2)
	// exit loop when there are no branches in left-side band, when there are no more branches to
	// look at in left-side band, or when all branches have been dealt with (no more in right-side band)
	while (base2 > Input->loops[base1].bri[0] && base2 > Input->loops[base1].bri[h+1] && currentNumBranches < Input->loops[base1].num_branches){
		for (l = Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1; l < Input->loops[base1].bri[h+1]; l++){
			misc_energy += pkmodelDP.c_p;
			g_count_c_p++;
			NumUnpairedInMulti++;
		}
		h++;
		currentNumBranches++;
	}

	if ( base2 <= Input->loops[base1].bri[0]){
		// if there are no branches between base1 and the pseudoknotted base pair...
		// do nothing; these bases were already counted above
	}else{
		// count the number of unpaired bases between end of last branch considered in loop above and beginning of pseudoknotted branch (base2)
		for (l = Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1; l < base2; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;

			NumUnpairedInMulti++;
		}
	}

	if (currentNumBranches >= Input->loops[base1].num_branches){
		// if no branches in the right-side of the band
		for (l = Input->BasePair(base2) + 1; l < Input->BasePair(base1); l++){
			misc_energy += pkmodelDP.c_p;
			g_count_c_p++;
			NumUnpairedInMulti++;
		}
	}else{
		// count the number of unpaired bases between end of pseudoknotted branch (bp(base2)) and beginning of next branch
		for (l = Input->BasePair(base2) + 1; l < Input->loops[base1].bri[currentNumBranches]; l++){
			misc_energy += pkmodelDP.c_p;
			g_count_c_p++;
			NumUnpairedInMulti++;
		}
	}

	// count the number of unpaired bases between the end of the last branch looked at in loop above, and the beginning of the next branch, until the end of the pseudoloop
	while (currentNumBranches + 1 < Input->loops[base1].num_branches){
		for (l = Input->ClosedRegions[Input->loops[base1].bri[currentNumBranches]]->end + 1; l < Input->loops[base1].bri[currentNumBranches+1]; l++)
		{
			misc_energy += pkmodelDP.c_p;
			g_count_c_p++;
			NumUnpairedInMulti++;
		}
		currentNumBranches++;
	}


	if (currentNumBranches >= Input->loops[base1].num_branches)
	{  // if no branches in the right-side of the band
		// do nothing; these bases were already counted above
	}else{
		// count the number of unpaired bases from the end of the last branch to the end of the multiloop
		for (l = Input->ClosedRegions[Input->loops[base1].bri[Input->loops[base1].num_branches - 1]]->end + 1; l < Input->BasePair(base1); l++){
			misc_energy += pkmodelDP.c_p;
			g_count_c_p++;
			NumUnpairedInMulti++;
		}
	}
	misc_energy *= 100;  // multiply by 100 since the remaining values in this function are in 10*cal/mol

	// add AU_penalties for multi-loop
	AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);

	for (h=0; h < f[i].num_branches; h++){
		AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);
	}

	// add dangling energies for multi-loop
	if (no_pk_dangling_ends == 0){
		dang += dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
		for (l=0; l < f[i].num_branches - 1; l++)
			dang += dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
		dang += dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);
	}
	// add "no-dangling" restriction
	for (l=0; l < f[i].num_branches; l++){
		Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
		Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;

	}

	return misc_energy + float(dang + AUpen);


};

/*********************************************************************************
multiPseudoEnergy: FOR PARAMETER TUNING Calculates the free energy of a multi-pseudoknotted Loop.
Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float	LoopList::multiPseudoEnergyDP(double *c, double &fr, int reset_c, int ignore_dangles){
	int num_params = get_num_params_PK_DP();
	if (reset_c == 1 && c != NULL){
		fr = 0;         
		for (int i=0; i < num_params; i++){
			c[i] = 0;
			//        for (int j=i; j < num_params; j++) 
			//                P_matrix[i][j] = 0;
		}
	}


	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = base1;
	PARAMTYPE AUpen;
	int h, l;
	PARAMTYPE dang;
	double misc_energy;

	dang = 0;
	misc_energy = 0;
	AUpen = 0;
	int special;
	special = 0;

	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int begin = base1; int end = Input->BasePair(base1);
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++){
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';

	int num_params_pkfree = get_num_params();

	// add penalty for introducing a multiloop that spans a band
	misc_energy += pkmodelDP.a_p;

	g_count_a_p++;
	if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "a_p") - 1] += 1;

	// +1 for the pseudoknotted base pair, +1 for the closing base pair of the multiloop that spans a band
	misc_energy += pkmodelDP.b_p * (f[i].num_branches + 2); 
	//	misc_energy += pkmodelDP.b_p * f[i].pseudo_num_branches;  // use this is pseudoloops count as 2 branches
	//	misc_energy += misc.multi_helix_penalty * f[i].pseudo_num_branches  + 2 * p_pairedMultiPseudo;

	g_count_b_p += f[i].num_branches + 2;
	if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "b_p") - 1] += f[i].num_branches + 2;


	// count the number of unpaired bases inside the multiloop
	int NumUnpairedInMulti = 0;
	int currentNumBranches = 0;  // the branch we are currently looking at

	if ( base2 <= Input->loops[base1].bri[0]){  // if there are no branches between base1 and the pseudoknotted base pair...
		for (l = base1 + 1; l < base2; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
	}else{
		// count the number of unpaired bases from base1 to start of first branch
		for (l=base1+1; l < Input->loops[base1].bri[0]; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
		currentNumBranches++;  // currentNumBranches traversed is now 1
	}

	// DONECHECK: this doesn't count unpaired bases inside a branch of the multiloop !!!!!! Correct

	// DONECHECK: what is the value of num_branches? is it double because of what's happening in FindChildren above?
	// DONE: num_branches does not include the branch that has pseudoknotted base pairs

	h = 0;
	// count number of unpaired bases between end of first branch and beginning of next branch, up the
	// branch that has the pseudoknotted base pair (starting at base2)
	// exit loop when there are no branches in left-side band, when there are no more branches to look 
	// at in left-side band, or when all branches have been dealt with (no more in right-side band)
	while (base2 > Input->loops[base1].bri[0] && base2 > Input->loops[base1].bri[h+1] && currentNumBranches < Input->loops[base1].num_branches){
		for (l = Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1; l < Input->loops[base1].bri[h+1]; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
		h++;
		currentNumBranches++;

		//		printf("Input->loops[base1].bri[h+1] = %d", Input->loops[base1].bri[h+1]);
	}

	if ( base2 <= Input->loops[base1].bri[0]){
		// if there are no branches between base1 and the pseudoknotted base pair...
		// do nothing; these bases were already counted above
	}else{
		// count the number of unpaired bases between end of last branch considered in loop above and beginning of pseudoknotted branch (base2)
		for (l = Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1; l < base2; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
	}

	if (currentNumBranches >= Input->loops[base1].num_branches){
		// if no branches in the right-side of the band
		for (l = Input->BasePair(base2) + 1; l < Input->BasePair(base1); l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
	}else{
		// count the number of unpaired bases between end of pseudoknotted branch (bp(base2)) and beginning of next branch
		for (l = Input->BasePair(base2) + 1; l < Input->loops[base1].bri[currentNumBranches]; l++){
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
	}

	// count the number of unpaired bases between the end of the last branch looked at in loop above, and the beginning of the next branch, until the end of the pseudoloop
	while (currentNumBranches + 1 < Input->loops[base1].num_branches){
		for (l = Input->ClosedRegions[Input->loops[base1].bri[currentNumBranches]]->end + 1; l < Input->loops[base1].bri[currentNumBranches+1]; l++)
		{
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
		currentNumBranches++;
	}


	if (currentNumBranches >= Input->loops[base1].num_branches){
		// if no branches in the right-side of the band
		// do nothing; these bases were already counted above
	}else{
		// count the number of unpaired bases from the end of the last branch to the end of the multiloop
		for (l = Input->ClosedRegions[Input->loops[base1].bri[Input->loops[base1].num_branches - 1]]->end + 1; l < Input->BasePair(base1); l++)
		{
			misc_energy += pkmodelDP.c_p;

			g_count_c_p++;
			if (c != NULL)    c[num_params_pkfree + structure_type_index_PK( (char *) "c_p") - 1] += 1;

			NumUnpairedInMulti++;
		}
	}

	misc_energy *= 100;  // multiply by 100 since the remaining values in this function are in 10*cal/mol

	// add AU_penalties for multi-loop
	AUpen += AU_penalty(sequence[i], sequence[f[i].pair]);
	if (c != NULL)    count_AU_penalty(sequence[i], sequence[f[i].pair], c);

	for (h=0; h < f[i].num_branches; h++){
		// AU penalty for pk-free children was already added in their energy calculation
		if (Input->ClosedRegions[f[i].bri[h]]->type == pseudo){
			AUpen += AU_penalty(sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);
			if (c != NULL)    count_AU_penalty(sequence[f[i].bri[h]], sequence[f[f[i].bri[h]].pair], c);
		}
	}

	if (ignore_dangles == 0){
		// add dangling energies for multi-loop
		dang += dangling_energy_left(sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
		if (c != NULL)    count_LEdangling_energy_left (sequence, structure, -1, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair, c);

		for (l=0; l < f[i].num_branches - 1; l++){
			dang += dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
			if (c != NULL)    count_LEdangling_energy (sequence, structure, -1, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair, c);

		}
		dang += dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);
		if (c != NULL)    count_LEdangling_energy_right (sequence, structure, -1, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair, c);

	}

	// add "no-dangling" restriction
	for (l=0; l < f[i].num_branches; l++){
		Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
		Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
	}


	return misc_energy + float(dang + AUpen);

};


///////////////////////////////////////////////////////////////////
// RIVAS & EDDY ENERGY MODEL
///////////////////////////////////////////////////////////////////

float	LoopList::stackPseudoEnergyRE(){
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = base1;

	int en = LEstacked_pair_energy (i, f[i].pair, sequence);

	return pkmodelRE.g_interiorPseudo * en;

}


float	LoopList::interiorPseudoEnergyRE(){

	pk_str_features * f = Input->loops;
	int * sequence = Input->type;

	if (  (base2 == base1 + 1 ) && (f[base1].pair == f[base2].pair + 1) )
		return stackPseudoEnergyRE();

	int i = base1;

	int ip, jp;
	ip = base2;
	jp = f[ip].pair;

	Input->cannot_add_dangling[ip-1] = 1;
	Input->cannot_add_dangling[jp+1] = 1;
	int en = LEinternal_loop_energy (i, f[i].pair, ip, jp, sequence);
	return pkmodelRE.g_interiorPseudo * en;
};


//	Cristina: added structure as a new parameter to dangling_energy functions in commonPK.cpp
float	LoopList::multiPseudoEnergyRE(){

	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = base1;
	int AUpen, h, l;
	int dang;
	float misc_energy;



	dang = 0;
	misc_energy = 0;
	AUpen = 0;
	int special;
	special = 0;


	// create structure: a string of dot-brackets to represent the region
	// NOTE: we don't want this to be the actual structure, involving (,[,<,etc
	//       since we just want to use the old simfold function before < meant
	//       something special
	int begin = base1; int end = Input->BasePair(base1);
	int numbases = end - begin + 1;
	char structure[numbases+1];
	for (int i = 0; i < numbases; i++)
	{
		if (Input->Sequence[begin + i] <= 0)
			structure[i] = '.';
		else if (Input->Sequence[begin + i] > (begin+i))
			structure[i] = '(';
		else
			structure[i] = ')';
	}
	structure[numbases] = '\0';


	misc_energy += 100 * pkmodelRE.M_tilda;

	// CHECK_RE: can change this to P_tilda instead of p_pairedMultiPseudo? Also, why is the extra 2 *p_pairedMultiPseudo required?
	misc_energy += misc.multi_helix_penalty * f[i].pseudo_num_branches  + 2 * 100 * pkmodelRE.p_pairedMultiPseudo;

	// add AU_penalties for multi-loop
	AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);

	for (h=0; h < f[i].num_branches; h++)
		AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);

	// add dangling energies for multi-loop
	dang += dangling_energy_left (sequence, structure, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
	for (l=0; l < f[i].num_branches - 1; l++)
		dang += dangling_energy (sequence, structure, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
	dang += dangling_energy_right (sequence, structure, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);

	// add "no-dangling" restriction
	for (l=0; l < f[i].num_branches; l++)
	{
		Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
		Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
	}
	return misc_energy + dang + AUpen;


};


///////////////////////////////////////////////////////////////////
// CAO & CHEN 2006 ENERGY MODEL
///////////////////////////////////////////////////////////////////


/*********************************************************************************
stackPseudoEnergy: Calculates the free energy of an interior-pseudoknotted loop which
is acutally a stacked pair. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float	LoopList::stackPseudoEnergyCC2006a(){
	if (DEBUG)
		printf("[stackPseudoEnergy CC2006] start\n");
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;
	int i = base1;

	int en = LEstacked_pair_energy (i, f[i].pair, sequence);

	return en;
}


/*********************************************************************************
interiorPseudoEnergy: Calculates the free energy of an interior-pseudoknotted loop
which is not a stacked pair. Energy value is returned in 10*cal/mol.
 *********************************************************************************/
float	LoopList::interiorPseudoEnergyCC2006a(){

	if (DEBUG)
		printf("[interiorPseudoEnergy CC2006] start\n");
	pk_str_features * f = Input->loops;
	int * sequence = Input->type;

	if (  (base2 == base1 + 1 ) && (f[base1].pair == f[base2].pair + 1) )
		return stackPseudoEnergyCC2006a();

	int i = base1;

	int ip, jp;
	ip = base2;
	jp = f[ip].pair;

	Input->cannot_add_dangling[ip-1] = 1;
	Input->cannot_add_dangling[jp+1] = 1;
	int en = LEinternal_loop_energy (i, f[i].pair, ip, jp, sequence);
	return en;
}

// UNECESSARY
/*********************************************************************************
multiPseudoEnergy:  Calculates the free energy of a multi-pseudoknotted Loop.
Energy value is returned in 10*cal/mol.
 *********************************************************************************/
/*
   float	LoopList::multiPseudoEnergyCC2006(){
   if (DEBUG)
   printf("[multiPseudoEnergy] start\n");

   pk_str_features * f = Input->loops;
   int * sequence = Input->type;
   int i = base1;
   int AUpen, h, l;
   int dang;
   double misc_energy;

   dang = 0;
   misc_energy = 0;
   AUpen = 0;
   int special;
   special = 0;

// add penalty for introducing a multiloop that spans a band
misc_energy += pkmodelDP.a_p;

//	misc_energy += pkmodelDP.Pps * (f[i].num_branches + 2);  // +1 for the pseudoknotted base pair, +1 for the closing base pair of the multiloop that spans a band
misc_energy += pkmodelDP.b_p * (f[i].num_branches + 2);  // +1 for the pseudoknotted base pair, +1 for the closing base pair of the multiloop that spans a band
//	misc_energy += pkmodelDP.b_p * f[i].pseudo_num_branches;  // use this is pseudoloops count as 2 branches
//	misc_energy += misc.multi_helix_penalty * f[i].pseudo_num_branches  + 2 * p_pairedMultiPseudo;

if (DEBUG)
printf("[misc_energy]1 is %f kcal/mol\n", misc_energy);

if (DEBUG2)
printf("[multiPseudoEnergyDP]: f[i].num_branches + 2 = %d, f[i].pseudo_num_branches = %d\n", f[i].num_branches + 2, f[i].pseudo_num_branches);

// count the number of unpaired bases inside the multiloop
int NumUnpairedInMulti = 0;
int currentNumBranches = 0;  // the branch we are currently looking at

if ( base2 <= Input->loops[base1].bri[0])
{  // if there are no branches between base1 and the pseudoknotted base pair...
for (l = base1 + 1; l < base2; l++)
{
misc_energy += pkmodelDP.c_p;
NumUnpairedInMulti++;
}
if (DEBUG2)
printf(" (1a) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", base1 + 1, base2, NumUnpairedInMulti);
}
else
{
// count the number of unpaired bases from base1 to start of first branch
for (l=base1+1; l < Input->loops[base1].bri[0]; l++)
{
misc_energy += pkmodelDP.c_p;
NumUnpairedInMulti++;
}
currentNumBranches++;  // currentNumBranches traversed is now 1
if (DEBUG2)
{
printf(" (1b) count unpaired from %d to %d --> ", base1+1, Input->loops[base1].bri[0]);
printf("NumUnpairedInMulti = %d, num_branches = %d, base2 = %d\n", NumUnpairedInMulti, Input->loops[base1].num_branches, base2);
}
}

// DONECHECK: this doesn't count unpaired bases inside a branch of the multiloop !!!!!! Correct

// DONECHECK: what is the value of num_branches? is it double because of what's happening in FindChildren above?
// DONE: num_branches does not include the branch that has pseudoknotted base pairs

h = 0;
// count number of unpaired bases between end of first branch and beginning of next branch, up the branch that has the pseudoknotted base pair (starting at base2)
while (base2 > Input->loops[base1].bri[0] && base2 > Input->loops[base1].bri[h+1] && currentNumBranches < Input->loops[base1].num_branches)
{  // exit loop when there are no branches in left-side band, when there are no more branches to look at in left-side band, or when all branches have been dealt with (no more in right-side band)
for (l = Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1; l < Input->loops[base1].bri[h+1]; l++)
{
	misc_energy += pkmodelDP.c_p;
	NumUnpairedInMulti++;
}
if (DEBUG2)
	printf(" (2) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1, Input->loops[base1].bri[h+1], NumUnpairedInMulti);
	h++;
	currentNumBranches++;

	//		printf("Input->loops[base1].bri[h+1] = %d", Input->loops[base1].bri[h+1]);
	}

if ( base2 <= Input->loops[base1].bri[0])
{  // if there are no branches between base1 and the pseudoknotted base pair...
	// do nothing; these bases were already counted above
	if (DEBUG2)
		printf(" (3a) counted these bases in (1a)\n");
}
else
{
	// count the number of unpaired bases between end of last branch considered in loop above and beginning of pseudoknotted branch (base2)
	for (l = Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1; l < base2; l++)
	{
		misc_energy += pkmodelDP.c_p;
		NumUnpairedInMulti++;
	}
	if (DEBUG2)
		printf(" (3b) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", Input->ClosedRegions[Input->loops[base1].bri[h]]->end + 1, base2, NumUnpairedInMulti);
}

if (currentNumBranches >= Input->loops[base1].num_branches)
{  // if no branches in the right-side of the band
	for (l = Input->BasePair(base2) + 1; l < Input->BasePair(base1); l++)
	{
		misc_energy += pkmodelDP.c_p;
		NumUnpairedInMulti++;
	}
	if (DEBUG2)
		printf(" (4a) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", Input->BasePair(base2) + 1, Input->BasePair(base1), NumUnpairedInMulti);
}
else
{
	// count the number of unpaired bases between end of pseudoknotted branch (bp(base2)) and beginning of next branch
	for (l = Input->BasePair(base2) + 1; l < Input->loops[base1].bri[currentNumBranches]; l++)
	{
		misc_energy += pkmodelDP.c_p;
		NumUnpairedInMulti++;
	}
	if (DEBUG2)
		printf(" (4b) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", Input->BasePair(base2) + 1, Input->loops[base1].bri[currentNumBranches], NumUnpairedInMulti);
}

//	printf("numbranches = %d\n", currentNumBranches);

// count the number of unpaired bases between the end of the last branch looked at in loop above, and the beginning of the next branch, until the end of the pseudoloop
while (currentNumBranches + 1 < Input->loops[base1].num_branches)
{
	for (l = Input->ClosedRegions[Input->loops[base1].bri[currentNumBranches]]->end + 1; l < Input->loops[base1].bri[currentNumBranches+1]; l++)
	{
		misc_energy += pkmodelDP.c_p;
		NumUnpairedInMulti++;
	}
	if (DEBUG2)
		printf(" (5) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", Input->ClosedRegions[Input->loops[base1].bri[currentNumBranches]]->end + 1, Input->loops[base1].bri[currentNumBranches+1], NumUnpairedInMulti);
	currentNumBranches++;
}


if (currentNumBranches >= Input->loops[base1].num_branches)
{  // if no branches in the right-side of the band
	// do nothing; these bases were already counted above
	if (DEBUG2)
		printf(" (6a) counted these bases in (4a)\n");
}
else
{
	// count the number of unpaired bases from the end of the last branch to the end of the multiloop
	for (l = Input->ClosedRegions[Input->loops[base1].bri[Input->loops[base1].num_branches - 1]]->end + 1; l < Input->BasePair(base1); l++)
	{
		misc_energy += pkmodelDP.c_p;
		NumUnpairedInMulti++;
	}
	if (DEBUG2)
		printf(" (6b) count unpaired from %d to %d --> sum(NumUnpairedInMulti) = %d\n", Input->ClosedRegions[Input->loops[base1].bri[Input->loops[base1].num_branches - 1]]->end + 1, Input->BasePair(base1), NumUnpairedInMulti);
}


if (DEBUG)
{
	printf("[misc_energy] with unpaired bases is %f kcal/mol\n", misc_energy);
	printf("NumUnpairedInMulti (multiloop that spans a band is %d\n", NumUnpairedInMulti);
}	

misc_energy *= 100;  // multiply by 100 since the remaining values in this function are in 10*cal/mol

// add AU_penalties for multi-loop
AUpen += AU_penalty (sequence[i], sequence[f[i].pair]);

for (h=0; h < f[i].num_branches; h++)
AUpen += AU_penalty (sequence[f[i].bri[h]],sequence[f[f[i].bri[h]].pair]);

// add dangling energies for multi-loop
dang += dangling_energy_left (sequence, i, f[i].pair, f[i].bri[0], f[f[i].bri[0]].pair);
for (l=0; l < f[i].num_branches - 1; l++)
dang += dangling_energy (sequence, f[i].bri[l], f[f[i].bri[l]].pair, f[i].bri[l+1], f[f[i].bri[l+1]].pair);
dang += dangling_energy_right (sequence, i, f[i].pair, f[i].bri[f[i].num_branches-1], f[f[i].bri[f[i].num_branches-1]].pair);

// add "no-dangling" restriction
for (l=0; l < f[i].num_branches; l++)
{
	Input->cannot_add_dangling [f[i].bri[l] -1] = 1;
	Input->cannot_add_dangling [f[f[i].bri[l]].pair + 1] = 1;
}
return misc_energy + dang + AUpen;


}
*/
