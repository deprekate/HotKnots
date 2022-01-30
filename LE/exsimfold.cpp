/* example of how to use the library libpairfold.dll */

// differences comparing to mfold: UAAUUGAUGACUGGGCCCAACUCCCUUGCUUGCUUCCCGC - simfold and RNAfold gives -4, mfold gives -2.8



#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// include the simfold library
#include "simfold.h"

#define MAXSLEN 3000

void create_random_sequence (int length, char *sequence)
// function to create uniformly random sequences - for demonstration purposes
{
    int rnumber, i;
    char base;
    for (i=0; i< length; i++)
    {
        rnumber = rand() % 100;
        if (rnumber >=0 && rnumber < 25)
            base = 'A';
        else if (rnumber >= 25 && rnumber < 50 )
            base = 'C';
        else if (rnumber >= 50 && rnumber < 75)
            base = 'G';
        else
            base = 'U';
        sequence[i] = base;
    }
    sequence[i] = '\0';
}




int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    double energy;
    int i;

    // configuration file
    char config_file[200] = "params/pairfold.conf";

    // what to fold: RNA or DNA
    int dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37;

    // seed initialization - needed by the random sequence generator
    srand((unsigned int)time((time_t *)(NULL)));

    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again (see below)
    init_data (config_file, dna_or_rna, temperature);

    strcpy (structure, "xxxxx........xxxxxxxxxx...........................xxxxx.............................................");
    printf ("%s\n------\n",structure);

    for (i=0; i < 10; i++)
      {
	create_random_sequence (100, sequence);
	strcpy (structure, "xxxxx........xxxxxxxxxx...........................xxxxx.............................................");
	energy = simfold (sequence, structure);
	printf ("%s\n%s\n%.2lf\n-----\n",sequence, structure, energy);
      }
    return 0;
}

