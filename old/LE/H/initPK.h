#ifndef INITPK_
#define INITPK_

/*
#include "constants.h"
#include "structsPK.h"
#include "globalsPK.h"
#include "externs.h"  // July 16 - added instead of globals.h
#include "commonPK.h"
//#include "globals.h"  // July 16 - remove
//#include "common.h"
*/  // July 16 - remove and place in initPK.cpp

/*
double ascii_to_doublePK (char *string);
void read_configuration_filePK (char *filename);
void giveupPK (char *string1, char *string2);
// to add: variable nb of parameters, as in scanf, printf

void read_pkmodelDP_file (char* filename, pkmodelinfoDP &pkmodelDP);
void read_pkmodelRE_file (char* filename, pkmodelinfoRE &pkmodelRE);
*/ // July 16 - remove since internal, helper functions for initPK.cpp


void init_dataPK(char* arg, char *config_file, int what, double temperature);
// the function that must be called by the main program to read data files
// PRE:  None
// POST: Read all data and configuration files

#endif /*INITPK_*/
