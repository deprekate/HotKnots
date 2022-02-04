/***************************************************************************
                          externs.h  -  description
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

#ifndef EXTERNSPK_H
#define EXTERNSPK_H

#include "constants.h"
#include "structsPK.h"

//#include "structs.h"  // July 16 - removed - unnecessary

// energies information
extern pkmodelinfoDP pkmodelDP;
extern pkmodelinfoRE pkmodelRE;  
extern pkmodelinfoCC2006 pkmodelCC2006;
extern float cc2006_s2_l1 [CC2006_STEMSIZE] [CC2006_LOOPSIZE];  // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
extern float cc2006_s1_l2 [CC2006_STEMSIZE] [CC2006_LOOPSIZE];  // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
extern float cc2006_s2_formula [4] [CC2006_STEMSIZE_FORMULA];    // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
extern float cc2006_s1_formula [4] [CC2006_STEMSIZE_FORMULA];    // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
extern int coaxstack_f_a [NUCL] [NUCL] [NUCL] [NUCL];  // coaxial stacking energy parameters (as described by Cao & Chen)
extern int coaxstack_f_b [NUCL] [NUCL] [NUCL] [NUCL];  // flush coaxial stacking energy parameters (as described by Mathews)
extern int coaxstack_m1 [NUCL] [NUCL] [NUCL] [NUCL];  // mismatch coaxial stacking for continuous backbone (as described by Mathews)
extern int coaxstack_m2 [NUCL] [NUCL] [NUCL] [NUCL];  // mismatch coaxial stacking for discontinuous backbone (as described by Mathews)
// *** Add new energy model code here (extern the variables in globalsPK.h)

extern char string_params_PK_CC[5000][MAXPNAME];    // for playing with the parameters

// PARAMETER TUNING

extern int g_count_Ps;		// pseudoloop initiation energy
extern int g_count_Psm;		// penalty for introducing pseudoknot inside a multiloop
extern int g_count_Psp;		// penalty for introducting pseudoknot inside a pseudoloop
extern int g_count_Pb;		// penalty for band
extern int g_count_Pup;		// penalty for unpaired base in pseudoloop or band
extern int g_count_Pps;		// penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band
extern int g_count_stP;		// multiplicative penalty for stacked pair that spans a band 
extern int g_count_intP;	// multiplicative penalty for internal loop that spans a band
extern int g_count_a;		// penalty for introducing a multiloop
extern int g_count_a_p;		// penalty for introducing a multiloop that spans a band
extern int g_count_b;		// penalty for multiloop base pair 
extern int g_count_b_p;		// penalty for multiloop base pair when the multiloop spans a band
extern int g_count_c;		// penalty for unpaired base in multiloop
extern int g_count_c_p;	

// parameters from the configuration file
extern char par_namePK [14] [100];   	// *** Add new energy model code here -- change first index to +2
extern char par_valuePK [14] [100];   // *** Add new energy model code here  -- change first index to +2

extern char * rna_pkmodelDP_par;
extern char * rna_pkmodelRE_par;
extern char * rna_pkmodelCC2006_par;
extern char * rna_coaxstack_f_a_par;
extern char * rna_coaxstack_f_b_par;
extern char * rna_coaxstack_m1_par;
extern char * rna_coaxstack_m2_par;
// *** Add new energy model code here

extern char * dna_pkmodelDP_par;
extern char * dna_pkmodelRE_par;
extern char * dna_pkmodelCC2006_par;
extern char * dna_coaxstack_f_a_par;
extern char * dna_coaxstack_f_b_par;
extern char * dna_coaxstack_m1_par;
extern char * dna_coaxstack_m2_par;
// *** Add new energy model code here


#endif
