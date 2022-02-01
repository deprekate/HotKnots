/***************************************************************************
                          structsPK.h  -  description
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

#ifndef STRUCTSPK_H
#define STRUCTSPK_H

#include <math.h>
#include "constants.h"
#include "structs.h"
//#include "Loop.h"

// VERSION 1 (below) before June 28 meeting (when b' still existed)
// VERSION 2 (see above) has replaced b' with Pps (nested region in multiloop that spans a band)
// However, VERSION 1 (below) will be used

// info for energy model for pseudoknotted structure
typedef struct
{
    double Ps;		// pseudoloop initiation energy
    double Psm;		// penalty for introducing pseudoknot inside a multiloop
    double Psp;		// penalty for introducting pseudoknot inside a pseudoloop
    double Pb;		// penalty for band
    double Pup;		// penalty for unpaired base in pseudoloop or band
    double Pps;		// penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band
    double stP;		// multiplicative penalty for stacked pair that spans a band 
    double intP;	// multiplicative penalty for internal loop that spans a band
    double a;		// penalty for introducing a multiloop
    double a_p;		// penalty for introducing a multiloop that spans a band
	double b;		// penalty for multiloop base pair 
	double b_p;		// penalty for multiloop base pair when the multiloop spans a band
	double c;		// penalty for unpaired base in multiloop
	double c_p;		// penalty for unpaired base in multiloop that spans a band
} pkmodelinfoDP;


typedef struct
{
	float g_interiorPseudo; // penalty for bulges, stems, and internal loops in a gap matrix (span-band loops)
	float P_tilda; 			// pair in a pseudoknot
	float P_i; 				// for E&R energy model  // CHECK_RE: ? used in pseudoEnergy() (but in R&E this is for pair in a non-nested multiloop)
	float Q_tilda; 			// unpaired bases in a pseudoknot
	float M_tilda;  		// non-nested multiloop (M_tilda)
	float Gw;  				// starting a pseudoknot loop
	float Gwh;  			// (overlapping pseudoknots) each more band region

	float Gwi;  // CHECK_RE: NOT CONSIDERED IN CURRENT CODE! generating a pseudoknot in a multiloop (Gwi)

	float q_unpairedMultiPseudo;  // CHECK_RE: NEVER USED IN CURRENT CODE! 
	float p_pairedMultiPseudo;    // CHECK_RE: where is this used, and convert it to pkmodel parameter

} pkmodelinfoRE;

// Cao & Chen 2006 energy model parameters
typedef struct
{
	float deltaG_assemble;  // pseudoknot penalty in 10cal/mol = KB*T*(ln 9), where T is multiplied after init_dataPK() is 
							//  executed from initPK.cpp
	float temp;  // temperature in Kelvin; used in Loop.cpp pseudoEnergy() function; set in initPK.cpp
} pkmodelinfoCC2006;

const int CC2006_STEMSIZE = 11;                         // max size of stem in Table 1 (starting from 1)
const int CC2006_LOOPSIZE = 12;                         // max size of loops in Table 1 (starting from 1)
const int CC2006_STEMSIZE_FORMULA = 11;         // max size of stem in Table 2 (starting from 1)
const float KB = 0.001987;  // Gas constant/Boltzmann constant [units: cal/(10 mol K)] used to multiply input file parameters from Table 1
							// use these units because entropy parameters are multiplied by 100 in initPK.cpp

// *** Add new energy model code here

typedef struct pk_str_features
{
    short int pair;
    char type;                   // type can be 'H', 'S', 'I', 'M'
    short int num_branches;      // number of branches in loop (for multipseudoknotted loops, the pseudoknotted base pair doesn't count)
	short int pseudo_num_branches;
    int bri[MAX_BRANCHES];       // the i of each branch
    
    pk_str_features()
    {
        pair = -1;
        type = NONE;
        num_branches = 0;
    }
} pk_str_features;

// THE STRUCT BELOW IS LOCAL TO LOOP.CPP  -- defined in Loop.h
//typedef struct pk_coax_features
//{
//	int index0;  // index of first loop in coaxial stacks
//	int index1;  // index of last loop in coaxial stacks; if there is more than one pair, this might not refer to loop that stacks with index0
//	Loop * loop0;  // pointer to first loop in the coaxial stacking; only set if this stacking involves one pair of loops
//	Loop * loop1;  // pointer to first loop in the coaxial stacking; only set if this stacking involves one pair of loops
//	int lastPairAdded;  // indexes the array pk_coax_features and points to the index of the last pair added in the coaxial stacking (between 0 and NumberOfChildren)
//	int lastArrayIndex; // indexes the array pk_coax_features and points to the index of the combination of pairs without lastPairAdded
//	float coax_energy;  // energy of the coaxial stacking of all combinations of pairs in this struct
//	int valid;  // represents whether this struct is a valid coaxial stacking 
//		// (can remove this if space is allocated dynamically in coaxial stacking function of Loop.cpp)
//		// can remove this even now, since the minimum is taken two at a time, not once at the end by looping through whole array
//} pk_coax_features;

#endif

