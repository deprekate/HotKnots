
// a simple driver for the simfold

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "constants.h"
#include "params.h"

char sequence[MAXSLEN];
char restricted[MAXSLEN] = "";
int num_subopt = 1;
char parameter_filename[500] = "";

void usage (const char *prog)
// PRE:  None
// POST: Prints the usage message and exits
{
    printf ("Usage: %s -s <sequence> [options]\n\n", prog);
    printf ("  -s <sequence>\n");
    printf ("\tThe RNA sequence to fold, max length is %d\n\n", MAXSLEN);
    printf ("Options:\n");
    printf ("  -r \"<restricted_structure>\"\n");
    printf ("\tThe restricted structure, if any. Restricted structure symbols are:\n");
    printf ("\t\t() restricted base pair\n");
    printf ("\t\t. base restricted to be unpaired\n\t\t_ no restriction\n\n");
    printf ("  -n <number_suboptimals>\n");
    printf ("\tThe number of suboptimal structures, including the optimal one. \n");
    printf ("\tIgnore if only the MFE (optimal) structure is desired.\n");
    printf ("\tMaximum allowed number of suboptimal structures is %d.\n\n", MAXSUBSTR);
    printf ("  -p <parameter_file>\n");
    printf ("\tIf none specified, it is using the Turner99 parameters.\n\n");
    printf ("  -d\n\tTurn debug on, by default it is off\n");
    printf ("  -h\n\tPrint this help message\n\n");
    printf ("Examples:\n");
    printf ("\t%s -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" -n 10 \n", prog);
    printf ("\t%s -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" -n 10 -r \"(____________________________)\"\n", prog);
    printf ("\t%s -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" -p params/CG_best_parameters_ISMB2007.txt\n", prog);
    exit (0);
}


void get_arguments (int argc, char *argv[])
// PRE:  None
// POST: Get the parameters from the command line
{
    int c;
    extern char *optarg;
    extern int optind;
    int errflag = 0;

    while ((c = getopt (argc, argv, "s:r:n:p:dh?")) != -1)
    {
        switch (c)
        {
            case 's':
                strcpy (sequence, optarg);
                break;
            case 'r':
                strcpy (restricted, optarg);
                break;
            case 'n':
                num_subopt = atoi (optarg);
                break;
            case 'p':
                strcpy (parameter_filename, optarg);
                break;
            case 'd':
                debug = 1;
                break;
            case 'h':
            case '?':
            default:
                errflag = 1;
        }
    }
    if (errflag || optind != argc || argc < 3)
        usage (argv[0]);

}



int main (int argc, char *argv[])
{
    char structure[MAXSLEN];    
    double energy;
    char structures[MAXSUBSTR][MAXSLEN];
    double energies[MAXSUBSTR];    
    int actual_num_str;
    //debug = 1;

    get_arguments (argc, argv);
    
    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory
    
    // configuration file, the path should be relative to the location of this executable
    char config_file[200];
    strcpy (config_file, "params/multirnafold.conf");

    // what to fold: RNA or DNA
    int dna_or_rna;
    dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37.0;
    
    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again   
    init_data (argv[0], config_file, dna_or_rna, temperature);         

    if (dna_or_rna == RNA && temperature == 37.0)
    {
        if (strlen (parameter_filename) > 0)
            fill_data_structures_with_new_parameters (parameter_filename);
        else
            fill_data_structures_with_new_parameters ("params/turner_parameters_fm363_constrdangles.txt");
    }

    printf ("Seq: %s\n", sequence);

    // if this structure must be restricted
    if (strlen (restricted) == strlen (sequence))
    {
        // only the MFE structure
        if (num_subopt == 1)
        {
            energy = simfold_restricted (sequence, restricted, structure);
            printf ("RES: %s  %.2lf\n", structure, energy);
        }
        else    // also suboptimal structures
        {
            // Just compute the first num_subopt restricted suboptimal structures.
            //  They may not be exactly the true ones for this model, because the model
            //      used here is a bit different from the one used for MFE prediction.
            actual_num_str = simfold_restricted_unordered_suboptimals (sequence, restricted, num_subopt, structures, energies);
            for (int i=0; i < actual_num_str; i++)
            {
                printf ("S %d: %s  %.2lf\n", i, structures[i], energies[i]);
            }        
        }        
    }
    else
    {
        // only the MFE structure
        if (num_subopt == 1)
        {
            energy = simfold (sequence, structure);
            printf ("MFE: %s  %.2lf\n", structure, energy);
        }
        else    // also suboptimal structures
        {
            // Here we actually compute 2*num_subopt suboptimal structures, compute free energy and reorder.
            //  That's because the energy model for the suboptimal structures is not exactly the same as the other one.
            actual_num_str = simfold_ordered_suboptimals (sequence, num_subopt, structures, energies);
            for (int i=0; i < actual_num_str; i++)
            {
                printf ("S %d: %s  %.2lf\n", i, structures[i], energies[i]);
            }
        }
    }
    
    return 0;
}
    
  
