#ifndef PARAMS_PK_H
#define PARAMS_PK_H

#include <stdio.h>
#include <stdlib.h>

void fill_data_structures_with_new_parameters_PK_DP (char *parameters_file);
void fill_data_structures_with_new_parameters_PK_CC2006b (char *parameters_file);
double free_energy_PK_DP(char *sequence, char *structure);
double free_energy_PK_CC2006b(char *sequence, char *structure);
double get_feature_counts_quadratic_PK_DP (char *sequence, char *structure, double **quadratic_matrix, double *counter, double &free_value);
double get_feature_counts_quadratic_PK_CC2006b (char *sequence, char *structure, double **quadratic_matrix, double *counter, double &free_value);
void get_feature_counts (char *sequence, char *structure, double *counter);
int get_num_params_PK_DP();
int get_num_params_PK_CC2006b();


#endif
