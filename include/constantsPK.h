#ifndef CONSTANTSPK_H_
#define CONSTANTSPK_H_

/////////////////////////////////////////////////////////////////////////////
//		 				USER DEFINED PARAMETERS
/////////////////////////////////////////////////////////////////////////////

// *** Change these flags to 1 to allow debugging statements to print out
const int DEBUG = 0;   // comphrehensive debugging statements for all energy functions
const int DEBUG2 = 0;  // further debugging information (simfold functions + pre-energy calculations)
const int DEBUGH = 0;  // print statements from the original HotKnots

// *** Specify the location of configuration file for input parameter filenames
//char config_file[200] = "/cs/beta/People/Cpop/2007/PKEnergy/Version1.0/params/pairfold.conf";

// *** Change these flags to 1 to output energy for one or more energy models; 0 to disable energy computation
const int RE_FLAG = 0;  // compute energy using Rivas & Eddy energy model
const int DP_FLAG = 1;  // compute energy using Dirks & Pierce energy model
const int CC2006a_FLAG = 0;  // compute energy using Cao & Chen 2006 energy model
const int CC2006b_FLAG = 0;  // compute energy using Cao & Chen 2006 energy model with variation
const int CC2006c_FLAG = 0;  // compute energy using Cao & Chen 2006 energy model with variation

// *** Add other energy models here (use next consecutive integer for new energy models)
const int RE = 0;  // id for Rivas&Eddy energy computation in Loop.cpp
const int DP = 1;  // id for Dirks&Pierce energy computation in Loop.cpp
const int CC2006a = 2;  // id for Cao&Chen 2006 energy compution in Loop.cpp
const int CC2006b = 3;  // id for Cao&Chen 2006 energy compution with variation in Loop.cpp
const int CC2006c = 4;  // id for Cao&Chen 2006 energy compution with variation in Loop.cpp

// *** For mismatch coaxial stacking; change to 1 if want to restrict only the dangling end involved in the stacking;
//									change to 0 if want to restrict both dangling ends close to the stacking helices
const int RES_STACK_DANGLE = 0;

// *** Maximum number of bases in the RNA strand
const int MaxN = 1600;

const int MAXCOAXSTACK = 100;  // maximum number of coaxial stacking combinations between children of a loop

// *** Change to 0 to include dangling ends for pk regions; 1 otherwise
//     NOTE: currently the 0 case is not handled properly - March 6, 2008
const int no_pk_dangling_ends = 0;

// *** Change to 0 to include coaxial stacking energies for pkfree regions; 1 only for pseudoknots (CC model)
//     NOTE: currently the 0 case is not handled properly - June 11, 2008
const int no_coax_in_pkfree = 1;

/////////////////////////////////////////////////////////////////////////////
// Do not change


// April 29, 2008 - Cristina
// the following matches the declaration at the top of Simfold/include/constants.h
// *** Add Simfold changes here
#ifdef DOUBLEPARAMS
#define PRINT_CHAR 0    // double
#define PRINT_MULTIPLIER 100

#elif LDOUBLEPARAMS
#define PRINT_CHAR 1    // long double
#define PRINT_MULTIPLIER 100

#else
#define PRINT_CHAR 2    // integer
#define PRINT_MULTIPLIER 100
#endif


const int LOWINF = 16000;  // infinity used by HotKnots (100 times smaller than that used by Simfold)

const int COAX_MULTI = 0;
const int COAX_PSEUDO = 1;
const int COAX_OTHER = 2;
// *** Add new energy model code here -- only applies to that dealing with coaxial stacking

/*
// Rivas & Eddy energy parameters -- implemented in structsPK.h
const float g_interiorPseudo = 0.83;
const int multi_OffsetPseudo = 843;
const int q_unpairedMultiPseudo = 0;
const int p_pairedMultiPseudo = 100;
*/



#endif /*CONSTANTSPK_H_*/
