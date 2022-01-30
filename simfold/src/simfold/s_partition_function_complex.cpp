/***************************************************************************
                          s_partition_function_complex.cpp  -  description
                             -------------------
    begin                : Wed Mar 15 2006
    copyright            : (C) 2006 by Mirela Andronescu
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>
#include <math.h>
#include <ctype.h>

#include "constants.h"
#include "structs.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "s_energy_matrix.h"
#include "s_hairpin_loop.h"
#include "s_stacked_pair.h"
#include "s_partition_function_complex.h"
#include "s_min_folding.h"
#include "s_sub_folding.h"
#include "params.h"

#define MERR 1.0e-14
// machine error


#ifdef DOUBLEPARAMS

Complex s_partition_function_complex::simple_function (int which, Complex addition)
{
    printf ("stack1=%d\n",  stack[1][2][1][2]);
    printf ("stack2=%d\n",  stack[0][3][0][3]);
    if (which == 0)
    {        
        stack_complex[1][2][1][2] = (Complex)(stack[1][2][1][2]) + addition;
        stack_complex[0][3][0][3] = (Complex)(stack[0][3][0][3]);
    }
    else
    {
        stack_complex[1][2][1][2] = (Complex)(stack[1][2][1][2]);
        stack_complex[0][3][0][3] = (Complex)(stack[0][3][0][3]) + addition;    
    }
    return (stack_complex[1][2][1][2]*stack_complex[1][2][1][2] + stack_complex[0][3][0][3]*stack_complex[0][3][0][3]);
}

// test constructor
void s_partition_function_complex::test_simple_function ()
// The constructor
{
    const Complex imaginary(0,1);   //sqrt(-1)
    PARAMTYPE hc = 1e-20;
    PARAMTYPE derivative1 = simple_function(0, hc*imaginary).imag()/hc;
    PARAMTYPE derivative2 = simple_function(1, hc*imaginary).imag()/hc;
    printf ("complex derivative1: %g\n", derivative1);
    printf ("complex derivative2: %g\n", derivative2);
}






s_partition_function_complex::s_partition_function_complex (char *cseq, int ignore, int compute_g_dangles, int pdouble)
// same as the constructor of s_partition_function, but with complex variables
// The constructor
{
    int i;
    seqlen = strlen (cseq);
    ignore_dangles = ignore;    
    compute_gradient_dangles = compute_g_dangles;
    params_are_double = pdouble;
    
    // for the exhaustive calculations
    IFD
        no_dangling_ends = 1;
    else
        simple_dangling_ends = 1;
            
    csequence = cseq;
    //strcpy (csequence, cseq);    // just refer it from where it is in memory
    
    for (i=0; i < seqlen; i++)
    {
        toupper(csequence[i]);
    }          
    sequence = new int[seqlen];
    if (sequence == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    for (i=0; i < seqlen; i++) 
    {
        sequence[i] = nuc_to_int(csequence[i]);
    }
    
    oneoverRT = -10.0/(1.98717*310.15);
    //oneoverRT = 1000.0/(1.98717*310.15);
    RT = -1.98717*310.15/10.0;
    
    index = new int [seqlen];
    // an array with indexes, such that we don't work with a 2D array, but with a 1D array of length (n*(n+1))/2
    int total_length = (seqlen *(seqlen+1))/2;
    index[0] = 0;
    for (int i=1; i < seqlen; i++)
        index[i] = index[i-1]+seqlen-i+1;

    up = new Complex [total_length];
    if (up == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");

    upm = new Complex [total_length];
    if (upm == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");   

    p = new Complex [total_length];
    if (p == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");

    IFD
    {
        u = new Complex[total_length];
        if (u == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        u1 = new Complex[total_length];
        if (u1 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        s1 = new Complex[total_length];
        if (s1 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        s2 = new Complex[total_length];
        if (s2 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        s3 = new Complex[total_length];
        if (s3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        pm = new Complex[total_length];
        if (pm == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        pm1 = new Complex[total_length];
        if (pm1 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        pm2 = new Complex[total_length];
        if (pm2 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    }
    else
    {
        u_ip_jp = new Complex [total_length];
        if (u_ip_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    
        u_ip_ju = new Complex [total_length];
        if (u_ip_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
        u_iu_jp = new Complex [total_length];
        if (u_iu_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
        u_iu_ju = new Complex [total_length];
        if (u_iu_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
        
        s1_jp = new Complex [total_length];
        if (s1_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");    
    
        s1_ju = new Complex [total_length];
        if (s1_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");        
                                    
        u1_ip_jp = new Complex [total_length];
        if (u1_ip_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    
        u1_ip_ju_jm1p = new Complex [total_length];
        if (u1_ip_ju_jm1p == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
        u1_ip_ju = new Complex [total_length];
        if (u1_ip_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
        u1_iu_jp = new Complex [total_length];
        if (u1_iu_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    
        u1_iu_ju_jm1p = new Complex [total_length];
        if (u1_iu_ju_jm1p == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");    
            
        u1_iu_ju = new Complex [total_length];
        if (u1_iu_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
    
        s2_jp = new Complex [total_length];
        if (s2_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
        s2_ju = new Complex [total_length];
        if (s2_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");    
        
        s3_jp = new Complex [total_length];
        if (s3_jp == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
        s3_ju_jm1p = new Complex [total_length];
        if (s3_ju_jm1p == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
            
        s3_ju = new Complex [total_length];
        if (s3_ju == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
                               
        pmnod3_needmidd3 = new Complex[total_length];
        if (pmnod3_needmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    
        pmnod3_noneedmidd3 = new Complex[total_length];
        if (pmnod3_noneedmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");    
            
        pmd3_needmidd3 = new Complex[total_length];
        if (pmd3_needmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    
        pmd3_noneedmidd3 = new Complex[total_length];
        if (pmd3_noneedmidd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");    
            
        pm1nod3_needendd3 = new Complex[total_length];
        if (pm1nod3_needendd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");        
        
        pm1d3_needendd3 = new Complex[total_length];
        if (pm1d3_needendd3 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
                
        pm2d5_needmidd5 = new Complex[total_length];
        if (pm2d5_needmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
    
        pm2d5_noneedmidd5 = new Complex[total_length];
        if (pm2d5_noneedmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
    
        pm2nod5_needmidd5 = new Complex[total_length];
        if (pm2nod5_needmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
    
        pm2nod5_noneedmidd5 = new Complex[total_length];
        if (pm2nod5_noneedmidd5 == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");            
    }                        
    
    
    // TODO: to erase in the final version
    num_params = create_string_params();
    
    // added this, which wasn't in s_partition_function
    GlogZ_numerical = new PFTYPE[num_params];
    if (GlogZ_numerical == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    for (i=0; i < num_params; i++)
        GlogZ_numerical[i] = 0.0;
    
    GlogZ_finite_differences = new PFTYPE[num_params];
    if (GlogZ_finite_differences == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    for (i=0; i < num_params; i++)
        GlogZ_finite_differences[i] = 0.0;
    
    //printf ("num_params=%d\n", num_params);
    GlogZ = new Complex[num_params];
    if (GlogZ == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");                        
    for (i=0; i < num_params; i++)
        GlogZ[i] = 0.0;

    parameters_complex = new Complex*[num_params];
    if (parameters_complex == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
        
    HlogZ_numerical = new PFTYPE*[num_params];
    for (i=0; i < num_params; i++)
    {
        HlogZ_numerical[i] = new PFTYPE[num_params];
    }
                   
    HlogZ_finite_differences = new PFTYPE*[num_params];
    for (i=0; i < num_params; i++)
    {
        HlogZ_finite_differences[i] = new PFTYPE[num_params];
    }

            
    initialize_arrays ();
    if (params_are_double)
        fill_data_structures_complex_double();
    else    
        fill_data_structures_complex();
    initialize_parameters_complex();

    // faster to compute eAU
    // TODO: it looks like EXPA, EXPB1 etc used in this form doesn't work, it gives the gradient of multi_helix_penalty 0.
    // TODO: to change it for the case ignore_dangles = 0

    // apparently expl doesn't exist in the complex class
    eAU = exp (AU_penalty_complex(A,U)* oneoverRT);
    EXPA = exp (misc_complex.multi_offset * oneoverRT);
    // EXPB is always called as EXPB1 or EXPB2
    Complex scaled_helix_penalty = misc_complex.multi_helix_penalty * oneoverRT; //
    EXPB1 = (Complex) exp (scaled_helix_penalty);
    //(Complex) exp (2*misc_complex.multi_helix_penalty * oneoverRT)
    EXPB2 = (Complex) exp (2 * scaled_helix_penalty.real());
    EXPC = new Complex [seqlen];
    if (EXPC == NULL) giveup ("Cannot allocate memory", "s_partition_function_complex");
    Complex scaled_free_base_penalty = misc_complex.multi_free_base_penalty * oneoverRT; //
    for (int i=0; i < seqlen; i++)
    {
        EXPC[i] = (Complex) exp (i * scaled_free_base_penalty.real());
    }    
    //(Complex) exp (i.0 * misc_complex.multi_free_base_penalty * oneoverRT)
    
    if (!ignore_dangles)
    {
        // now fill the edangle3 and edangle5 arrays
        // it's actually a looot faster (about 5 times on length 200)
        for (int i=0; i < NUCL; i++)
        {
            for (int j=0; j < NUCL; j++)
            {
                if (can_pair (i,j))
                {
                    for (int k=0; k < NUCL; k++)
                    {                
                        edangle5[i][j][k] = exp (dangle_bot_complex[i][j][k]*oneoverRT);
                        edangle3[i][j][k] = exp (dangle_top_complex[i][j][k]*oneoverRT);
                    }
                }
            }
        }
    }    
}


void s_partition_function_complex::fill_data_structures_complex ()
// Mirela: Mar 12, 2007
// fill the complex parameter data structures from the global ones
{
    int index;
    int i, j, k, l, m, n, o, p;
    index = 0;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                dangle_top_complex[i][j][k] = (Complex) dangle_top[i][j][k];
                dangle_bot_complex[i][j][k] = (Complex) dangle_bot[i][j][k];
                for (l=0; l < NUCL; l++)
                {
                    stack_complex[i][j][k][l] = (Complex) stack[i][j][k][l];
                    tstackh_complex[i][j][k][l] = (Complex) tstackh[i][j][k][l];
                    tstacki_complex[i][j][k][l] = (Complex) tstacki[i][j][k][l];
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            int11_complex[i][j][k][l][m][n] = (Complex) int11[i][j][k][l][m][n];
                            for (o=0; o < NUCL; o++)
                            {
                                int21_complex[i][j][k][l][m][n][o] = (Complex) int21[i][j][k][l][m][n][o];
                                for (p=0; p < NUCL; p++)
                                    int22_complex[i][j][k][l][m][n][o][p] = (Complex) int22[i][j][k][l][m][n][o][p];
                            }
                        }
                }
            }


    for (i=1; i <= MAXLOOP_I; i++)
    {
        internal_penalty_by_size_complex[i] = (Complex) internal_penalty_by_size[i];
    }
    for (i=1; i <= MAXLOOP_B; i++)
    {
        bulge_penalty_by_size_complex[i] = (Complex) bulge_penalty_by_size[i];
    }
    for (i=1; i <= MAXLOOP_H; i++)
    {
        hairpin_penalty_by_size_complex[i] = (Complex) hairpin_penalty_by_size[i];
    }

    #if (MODEL == SIMPLE)
    for(i=0; i < nb_triloops; i++)
    {
        triloop_complex[i].energy = (Complex) triloop[i].energy;
        //triloop_complex[i].seq = triloop[i].seq;
    }
    for(i=0; i < nb_tloops; i++)
    {
        tloop_complex[i].energy = (Complex) tloop[i].energy;
        //tloop_complex[i].seq = tloop[i].seq;
    }
    #endif
                    
    misc_complex.param_greater30 = misc.param_greater30;    // double, not complex yet
    misc_complex.terminal_AU_penalty = (Complex) misc.terminal_AU_penalty;
    misc_complex.hairpin_GGG = (Complex) misc.hairpin_GGG;
    misc_complex.hairpin_c1 = (Complex) misc.hairpin_c1;
    misc_complex.hairpin_c2 = (Complex) misc.hairpin_c2;
    misc_complex.hairpin_c3 = (Complex) misc.hairpin_c3;

    misc_complex.asymmetry_penalty_max_correction = (Complex) misc.asymmetry_penalty_max_correction;
    for (i=0; i < 4; i++)
        misc_complex.asymmetry_penalty_array[i] = (Complex) misc.asymmetry_penalty_array[i];
    misc_complex.gail_rule = misc.gail_rule;    // this is int
  
    misc_complex.multi_offset = (Complex) misc.multi_offset;
    misc_complex.multi_helix_penalty = (Complex) misc.multi_helix_penalty;
    misc_complex.multi_free_base_penalty = (Complex) misc.multi_free_base_penalty;
    misc_complex.intermolecular_initiation = (Complex) misc.intermolecular_initiation;

    // 3 params instead of the whole tstacki table
    misc_complex.internal_AU_closure = (Complex) misc.internal_AU_closure;
    misc_complex.internal_AG_mismatch = (Complex) misc.internal_AG_mismatch;
    misc_complex.internal_UU_mismatch = (Complex) misc.internal_UU_mismatch;

   
    #if (MODEL == SIMPLE)
    // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
    misc_complex.internal11_basic_mismatch = (Complex) misc.internal11_basic_mismatch;
    misc_complex.internal11_GG_mismatch = (Complex) misc.internal11_GG_mismatch;    
    misc_complex.internal21_match = (Complex) misc.internal21_match;
    misc_complex.internal21_AU_closure = (Complex) misc.internal21_AU_closure;    
    #endif

    misc_complex.internal22_delta_same_size = (Complex) misc.internal22_delta_same_size;
    misc_complex.internal22_delta_different_size = (Complex) misc.internal22_delta_different_size;
    misc_complex.internal22_delta_1stable_1unstable = (Complex) misc.internal22_delta_1stable_1unstable;
    misc_complex.internal22_delta_AC = (Complex) misc.internal22_delta_AC;
    misc_complex.internal22_match = (Complex) misc.internal22_match;
}


void s_partition_function_complex::fill_data_structures_complex_double ()
// DEPRICATED
// Mirela: Sep 11, 2007
// fill the complex parameter data structures from the global ones, but with the double parameters
{
    int index;
    int i, j, k, l, m, n, o, p;
    index = 0;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                dangle_top_complex[i][j][k] = (Complex) dangle_top_double[i][j][k];
                dangle_bot_complex[i][j][k] = (Complex) dangle_bot_double[i][j][k];
                for (l=0; l < NUCL; l++)
                {
                    stack_complex[i][j][k][l] = (Complex) stack_double[i][j][k][l];
                    tstackh_complex[i][j][k][l] = (Complex) tstackh_double[i][j][k][l];
                    tstacki_complex[i][j][k][l] = (Complex) tstacki_double[i][j][k][l];
                    for (m=0; m < NUCL; m++)
                        for (n=0; n < NUCL; n++)
                        {
                            int11_complex[i][j][k][l][m][n] = (Complex) int11_double[i][j][k][l][m][n];
                            for (o=0; o < NUCL; o++)
                            {
                                int21_complex[i][j][k][l][m][n][o] = (Complex) int21_double[i][j][k][l][m][n][o];
                                for (p=0; p < NUCL; p++)
                                    int22_complex[i][j][k][l][m][n][o][p] = (Complex) int22_double[i][j][k][l][m][n][o][p];
                            }
                        }
                }
            }


    for (i=1; i <= MAXLOOP_I; i++)
    {
        internal_penalty_by_size_complex[i] = (Complex) internal_penalty_by_size_double[i];
    }
    for (i=1; i <= MAXLOOP_B; i++)
    {
        bulge_penalty_by_size_complex[i] = (Complex) bulge_penalty_by_size_double[i];
    }
    for (i=1; i <= MAXLOOP_H; i++)
    {
        hairpin_penalty_by_size_complex[i] = (Complex) hairpin_penalty_by_size_double[i];
    }

    #if (MODEL == SIMPLE)
    for(i=0; i < nb_triloops; i++)
    {
        triloop_complex[i].energy = (Complex) triloop_double[i].energy;
        //triloop_complex[i].seq = triloop[i].seq;
    }
    for(i=0; i < nb_tloops; i++)
    {
        tloop_complex[i].energy = (Complex) tloop_double[i].energy;
        //tloop_complex[i].seq = tloop[i].seq;
    }
    #endif    
                    
    misc_complex.param_greater30 = misc_double.param_greater30;    // double, not complex yet
    misc_complex.terminal_AU_penalty = (Complex) misc_double.terminal_AU_penalty;
    misc_complex.hairpin_GGG = (Complex) misc_double.hairpin_GGG;
    misc_complex.hairpin_c1 = (Complex) misc_double.hairpin_c1;
    misc_complex.hairpin_c2 = (Complex) misc_double.hairpin_c2;
    misc_complex.hairpin_c3 = (Complex) misc_double.hairpin_c3;

    misc_complex.asymmetry_penalty_max_correction = (Complex) misc_double.asymmetry_penalty_max_correction;
    for (i=0; i < 4; i++)
        misc_complex.asymmetry_penalty_array[i] = (Complex) misc_double.asymmetry_penalty_array[i];
    misc_complex.gail_rule = misc_double.gail_rule;    // this is int
  
    misc_complex.multi_offset = (Complex) misc_double.multi_offset;
    misc_complex.multi_helix_penalty = (Complex) misc_double.multi_helix_penalty;
    misc_complex.multi_free_base_penalty = (Complex) misc_double.multi_free_base_penalty;
    misc_complex.intermolecular_initiation = (Complex) misc_double.intermolecular_initiation;

    // 3 params instead of the whole tstacki table
    misc_complex.internal_AU_closure = (Complex) misc_double.internal_AU_closure;
    misc_complex.internal_AG_mismatch = (Complex) misc_double.internal_AG_mismatch;
    misc_complex.internal_UU_mismatch = (Complex) misc_double.internal_UU_mismatch;

    // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
    misc_complex.internal11_basic_mismatch = (Complex) misc_double.internal11_basic_mismatch;
    misc_complex.internal11_GG_mismatch = (Complex) misc_double.internal11_GG_mismatch;
    
    misc_complex.internal21_match = (Complex) misc_double.internal21_match;
    misc_complex.internal21_AU_closure = (Complex) misc_double.internal21_AU_closure;

    misc_complex.internal22_delta_same_size = (Complex) misc_double.internal22_delta_same_size;
    misc_complex.internal22_delta_different_size = (Complex) misc_double.internal22_delta_different_size;
    misc_complex.internal22_delta_1stable_1unstable = (Complex) misc_double.internal22_delta_1stable_1unstable;
    misc_complex.internal22_delta_AC = (Complex) misc_double.internal22_delta_AC;
    misc_complex.internal22_match = (Complex) misc_double.internal22_match;
    //printf ("stack complex = %lf\n", stack_complex[1][2][1][2].real());
}



void s_partition_function_complex::compute_logZ_gradient_numerical ()
{
    int index;

    // variables for the complex trick
    PARAMTYPE hc = 1e-20;
    const Complex imaginary(0,1);   // sqrt(-1)
    Complex saved;

    Complex logZ;
    // TODO: should go to num_params
    for (index = 0; index < num_params; index++)
    {      
        saved = *parameters_complex[index];
        *parameters_complex[index] += hc*imaginary;

        compute_partition_function();
        logZ = log(Z);
        GlogZ_numerical[index] = RT.real() * logZ.imag() / hc;
        *parameters_complex[index] = saved;
    }
}


void s_partition_function_complex::compute_logZ_gradient_numerical (int index)
// only compute the partial derivative wrt index
{
    // variables for the complex trick
    PARAMTYPE hc = 1e-20;
    const Complex imaginary(0,1);   // sqrt(-1)
    Complex saved;

    Complex logZ;
    saved = *parameters_complex[index];
    *parameters_complex[index] += hc*imaginary;

    compute_partition_function();
    logZ = log(Z);
    GlogZ_numerical[index] = RT.real() * logZ.imag() / hc;
    *parameters_complex[index] = saved;
}


void s_partition_function_complex::compute_logZ_gradient_finite_differences ()
{
    int index;

    // variables for the finite differences method
    PARAMTYPE mu = 1e-7;
    Complex saved;

    PFTYPE regular_logZ = log(Z.real());
    PFTYPE new_logZ;
    // TODO: should go to num_params
    for (index = 0; index < num_params; index++)
    {
        saved = *parameters_complex[index];
        *parameters_complex[index] += mu;

        compute_partition_function();
        new_logZ = log(Z.real());
        GlogZ_finite_differences[index] = RT.real() * (new_logZ - regular_logZ) / mu;
        //printf ("regular=%e, new=%e, diff=%e, G[%d] finite diff = %g\n", regular_logZ, Z.real(), new_logZ-regular_logZ, index, GlogZ_finite_differences[index]);
        *parameters_complex[index] = saved;
    }
}


void s_partition_function_complex::compute_logZ_second_derivatives_numerical ()
// We can evalute this by also computing the finite differences of it
// compute the logZ 
{
    int index1, index2;

    // variables for the complex trick
    PARAMTYPE hc = 1e-20;
    const Complex imaginary(0,1);   // sqrt(-1)
    Complex saved;

    for (index1 = 0; index1 < num_params; index1++)
    {    
        IFD
        {
            if (index1 >= 260 && index1 <= 307)
                continue;
        }    
        saved = *parameters_complex[index1];
        *parameters_complex[index1] += hc*imaginary;

        compute_partition_function();
        compute_base_pair_probabilities();
        compute_logZ_gradient();
        // it's a symmetric matrix, so it's enough to compute the upper triangle
        for (index2=index1; index2 < num_params; index2++)
        {
            HlogZ_numerical[index1][index2] = RT.real() * GlogZ[index2].imag() / hc;
//             if (HlogZ_numerical[index1][index2] != 0)
//             {
//                 printf ("H[%d][%d] = %g\t%d=%s, %d=%s\n", index1, index2, HlogZ_numerical[index1][index2], index1, string_params[index1], index2, string_params[index2]);
//             }
            if (index1 != index2)
                HlogZ_numerical[index2][index1] = HlogZ_numerical[index1][index2];
        }
        *parameters_complex[index1] = saved;
    }
}


void s_partition_function_complex::compute_logZ_second_derivatives_finite_differences ()
{
    
    // variables for the finite differences method
    PARAMTYPE mu = 1e-7;
    Complex saved;

    /*
    PFTYPE regular_logZ = log(Z.real());
    PFTYPE new_logZ;
    // TODO: should go to num_params
    for (index = 0; index < num_params; index++)
    {
        saved = *parameters_complex[index];
        *parameters_complex[index] += mu;

        compute_partition_function();
        new_logZ = log(Z.real());
        GlogZ_finite_differences[index] = RT.real() * (new_logZ - regular_logZ) / mu;
        //printf ("regular=%e, new=%e, diff=%e, G[%d] finite diff = %g\n", regular_logZ, Z.real(), new_logZ-regular_logZ, index, GlogZ_finite_differences[index]);
        *parameters_complex[index] = saved;
    }


    int index1, index2;

    // variables for the finite differences method
    PARAMTYPE mu = 1e-7;
    Complex saved;

    // TODO: IFD, don't do the dangling end derivatives
    for (index1 = 0; index1 < num_params; index1++)
    {    
        IFD
        {
            if (index1 >= 260 && index1 <= 307)
                continue;
        }    
        saved = *parameters_complex[index1];
        *parameters_complex[index1] += mu;

        compute_partition_function();
        compute_base_pair_probabilities();
        compute_logZ_gradient();
        // it's a symmetric matrix, so it's enough to compute the upper triangle
        for (index2=index1; index2 < num_params; index2++)
        {
            HlogZ_numerical[index1][index2] = RT.real() * GlogZ[index2].imag() / hc;
//             if (HlogZ_numerical[index1][index2] != 0)
//             {
//                 printf ("H[%d][%d] = %g\t%d=%s, %d=%s\n", index1, index2, HlogZ_numerical[index1][index2], index1, string_params[index1], index2, string_params[index2]);
//             }
            if (index1 != index2)
                HlogZ_numerical[index2][index1] = HlogZ_numerical[index1][index2];
        }
        *parameters_complex[index1] = saved;
    }
    */
}

        
void s_partition_function_complex::compute_logZ_second_derivatives_numerical (int index1)
// compute the second derivatives of logZ wrt theta_index1
// this way it is easier to parallelize
// unfortunately we can't really evaluate this
{
    int index2;

    // variables for the complex trick
    PARAMTYPE hc = 1e-20;
    const Complex imaginary(0,1);   // sqrt(-1)
    Complex saved;
    
    saved = *parameters_complex[index1];
    *parameters_complex[index1] += hc*imaginary;

    compute_partition_function();
    compute_base_pair_probabilities();
    compute_logZ_gradient();
    // it's a symmetric matrix, so it's enough to compute the upper triangle
    for (index2=index1; index2 < num_params; index2++)
    {
        HlogZ_numerical[index1][index2] = RT.real() * GlogZ[index2].imag() / hc;
        // in this case, I only care about the current row, so I don't need what's below
        //if (index1 != index2)
        //    HlogZ_numerical[index2][index1] = HlogZ_numerical[index1][index2];
    }
    *parameters_complex[index1] = saved;

}


void s_partition_function_complex::print_Hessian ()
// print the second derivative matrix to the screen
{
    int index1, index2;
    for (index1 = 0; index1 < num_params; index1++)
    {
        for (index2=0; index2 < num_params; index2++)
        {
            printf ("%g\t", HlogZ_numerical[index1][index2]);
        }
        printf ("\n");
    }
}


void s_partition_function_complex::print_Hessian (char *filename)
// print the second derivative matrix into the file filename
{
    FILE *file;
    if ((file = fopen (filename, "w")) == NULL)
    {
        printf ("Cannot open file %s\n", filename);
        exit (0);
    }     

    int index1, index2;
    for (index1 = 0; index1 < num_params; index1++)
    {
        for (index2=0; index2 < num_params; index2++)
        {
            fprintf (file, "%g\t", HlogZ_numerical[index1][index2]);
        }
        fprintf (file, "\n");
    }
    fclose (file);
}


void s_partition_function_complex::print_Hessian (int index1, char *filename)
// print the second derivative matrix into the file filename
// only print the row index1
{
    FILE *file;
    if ((file = fopen (filename, "w")) == NULL)
    {
        printf ("Cannot open file %s\n", filename);
        exit (0);
    }     

    int index2;
    // in this case, I don't know what's before index1, so I don't start from 0, but from index1
    for (index2=index1; index2 < num_params; index2++)
    {
        fprintf (file, "%g\t", HlogZ_numerical[index1][index2]);
    }
    fprintf (file, "\n");
    fclose (file);
}


void s_partition_function_complex::initialize_parameters_complex ()
// make the parameters_complex array to point to the right parameters;
{
    int index;
    int i,j,k,l,m,n,o,p;

    index = 0;
    
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (stack[i][j][k][l] < INF)
                    {
                        // exclude duplicates
                        // stack[i][j][k][l] is the same as stack[l][k][j][i]
                        if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                        {
                            parameters_complex[index++] = &stack_complex[i][j][k][l];
                        }
                    }
                }
                
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (tstackh[i][j][k][l] < INF)
                    {
                        // no duplicates here
                        //fprintf (file, "%.2lf\n", tstackh[i][j][k][l]/100.0);
                        //sprintf (string_params[index++], "tstackh[%d][%d][%d][%d]", i, j, k, l);
                        parameters_complex[index++] = &tstackh_complex[i][j][k][l];
                    }
                }
                
    
    // we only use 3 parameters here, instead of the whole tstacki table
    //fprintf (file, "%.2lf\n", misc.internal_AU_closure/100.0);
    parameters_complex[index++] = &misc_complex.internal_AU_closure;
    
    //fprintf (file, "%.2lf\n", misc.internal_AG_mismatch/100.0);
    parameters_complex[index++] = &misc_complex.internal_AG_mismatch;
    
    //fprintf (file, "%.2lf\n", misc.internal_UU_mismatch/100.0);
    parameters_complex[index++] = &misc_complex.internal_UU_mismatch;


    if (!simple_internal_energy)
    {
        // do the few int 11 params: only those enclosed by CG and CG (any order), + 2 more params
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {
                                if ( (((i==C && j==G) || (i==G && j==C)) &&
                                     ((m==C && n==G) || (m==G && n==C)) &&
                                     !can_pair(k,l)) ||
                                     (watson_crick(i,j) && watson_crick(m,n) && k==U && l==U))
                                {
                                    if (int11[i][j][k][l][m][n] < INF)
                                    {
                                        // exclude duplicates
                                        // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                        if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                        {
                                            //fprintf (file, "%.2lf\n", int11[i][j][k][l][m][n]/100.0);
                                            parameters_complex[index++] = &int11_complex[i][j][k][l][m][n];
                                        }
                                    }
                                }
                            }
                            
        //fprintf (file, "%.2lf\n", misc.internal11_basic_mismatch/100.0);
        parameters_complex[index++] = &misc_complex.internal11_basic_mismatch;
        
        //fprintf (file, "%.2lf\n", misc.internal11_GG_mismatch/100.0);
        parameters_complex[index++] = &misc_complex.internal11_GG_mismatch;
    
        // go with few parameters, as in Mathews et al 1999
        // closed by CG
        i=C; j=G; m=C; n=G;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                    if (!can_pair(k,l) && !can_pair(k,o))
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            //fprintf (file, "%.2lf\n", int21[i][j][k][l][m][n][o]/100.0);
                            parameters_complex[index++] = &int21_complex[i][j][k][l][m][n][o];
                        }
                        
        // closed by GC                        
        i=G; j=C; m=G; n=C;
        for (k=0; k < NUCL; k++)
            for (l=0; l < NUCL; l++)
                for (o=0; o < NUCL; o++)
                    if (!can_pair(k,l) && !can_pair(k,o))
                        if (int21[i][j][k][l][m][n][o] < INF)
                        {
                            // no duplicates here
                            //fprintf (file, "%.2lf\n", int21[i][j][k][l][m][n][o]/100.0);
                            parameters_complex[index++] = &int21_complex[i][j][k][l][m][n][o];
                        }

                        
        //fprintf (file, "%.2lf\n", misc.internal21_match/100.0);
        parameters_complex[index++] = &misc_complex.internal21_match;
        
        //fprintf (file, "%.2lf\n", misc.internal21_AU_closure/100.0);
        parameters_complex[index++] = &misc_complex.internal21_AU_closure;

        // go with the 53 params, like in Mathews et al 1999
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                    {
                        n = i;
                        m = j;
                        p = k;
                        o = l;                    
                        if (watson_crick(i,j) && !watson_crick(k,l))
                        {
                            // exclude duplicates
                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + p <=
                                n*10000000 + m*1000000 + p*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                            {
                                //fprintf (file, "%.2lf\n", int22[i][j][k][l][m][n][o][p]/100.0);
                                //sprintf (string_params[index++], "int22[%d][%d][%d][%d][%d][%d][%d][%d]", i, j, k, l, m, n, o, p);
                                parameters_complex[index++] = &int22_complex[i][j][k][l][m][n][o][p];
                            }
                        }
                    }

        
        // add the 4 deltas
        //fprintf (file, "%.2lf\n", misc.internal22_delta_same_size/100.0);
        parameters_complex[index++] = &misc_complex.internal22_delta_same_size;
        
        //fprintf (file, "%.2lf\n", misc.internal22_delta_different_size/100.0);
        parameters_complex[index++] = &misc_complex.internal22_delta_different_size;
        
        //fprintf (file, "%.2lf\n", misc.internal22_delta_1stable_1unstable/100.0);
        parameters_complex[index++] = &misc_complex.internal22_delta_1stable_1unstable;
        
        //fprintf (file, "%.2lf\n", misc.internal22_delta_AC/100.0);
        parameters_complex[index++] = &misc_complex.internal22_delta_AC;
        
        //fprintf (file, "%.2lf\n", misc.internal22_match/100.0);
        parameters_complex[index++] = &misc_complex.internal22_match;
        
    }    // end if (!simple_internal_energy)

    
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                if (dangle_top[i][j][k] < INF)
                {
                    parameters_complex[index++] = &dangle_top_complex[i][j][k];
                }
            }
            
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                if (dangle_bot[i][j][k] < INF)
                {
                    parameters_complex[index++] = &dangle_bot_complex[i][j][k];
                }
            }
            
    int start;        
    if (!simple_internal_energy)
        start = 4;
    else
        start = 1;                
    for (i=start; i <= MAXLOOP_I; i++)
    {
        if (internal_penalty_by_size[i] < INF)
        {
            // no duplicates here
            //fprintf (file, "%.2lf\n", internal_penalty_by_size[i]/100.0);
            //sprintf (string_params[index++], "internal_penalty_by_size[%d]", i);
            parameters_complex[index++] = &internal_penalty_by_size_complex[i];
        }
    }
    
    for (i=1; i <= MAXLOOP_B; i++)
    {
        if (bulge_penalty_by_size[i] < INF)
        {
            // no duplicates here
            //fprintf (file, "%.2lf\n", bulge_penalty_by_size[i]/100.0);
            //sprintf (string_params[index++], "bulge_penalty_by_size[%d]", i);
            parameters_complex[index++] = &bulge_penalty_by_size_complex[i];
        }
    }
    
    for (i=1; i <= MAXLOOP_H; i++)
    {
        if (hairpin_penalty_by_size[i] < INF)
        {
            // no duplicates here
            //fprintf (file, "%.2lf\n", hairpin_penalty_by_size[i]/100.0);
            //sprintf (string_params[index++], "hairpin_penalty_by_size[%d]", i);
            parameters_complex[index++] = &hairpin_penalty_by_size_complex[i];
        }
    }
    
    //fprintf (file, "%.2lf\n", misc.param_greater30/100.0);
    //sprintf (string_params[index++], "misc.param_greater30");
    
    //fprintf (file, "%.2lf\n", misc.terminal_AU_penalty/100.0);  
    //sprintf (string_params[index++], "misc.terminal_AU_penalty");
    parameters_complex[index++] = &misc_complex.terminal_AU_penalty;
    
    //fprintf (file, "%.2lf\n", misc.hairpin_GGG/100.0);
    //sprintf (string_params[index++], "misc.hairpin_GGG");
    parameters_complex[index++] = &misc_complex.hairpin_GGG;
    
    //fprintf (file, "%.2lf\n", misc.hairpin_c1/100.0);
    //sprintf (string_params[index++], "misc.hairpin_c1");
    parameters_complex[index++] = &misc_complex.hairpin_c1;
    
    //fprintf (file, "%.2lf\n", misc.hairpin_c2/100.0);    
    //sprintf (string_params[index++], "misc.hairpin_c2");
    parameters_complex[index++] = &misc_complex.hairpin_c2;
    
    //fprintf (file, "%.2lf\n", misc.hairpin_c3/100.0);
    //sprintf (string_params[index++], "misc.hairpin_c3");
    parameters_complex[index++] = &misc_complex.hairpin_c3;

    
    // TODO
    //fprintf (file, "%d\n", misc.asymmetry_penalty_max_correction);
    //fprintf (file, "%d\n", misc.asymmetry_penalty_array[0]);
    //fprintf (file, "%d\n", misc.asymmetry_penalty_array[1]);
    
    //fprintf (file, "0\n");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[1]");
    //fprintf (file, "0\n");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[2]");
    //fprintf (file, "0\n");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[3]");
    //fprintf (file, "0\n");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[4]");
    //fprintf (file, "0\n");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[5]");
    //fprintf (file, "0\n");
    //sprintf (string_params[index++], "misc.asymmetry_penalty[6]");  
    
    
    //fprintf (file, "%.2lf\n", misc.gail_rule);
    //sprintf (string_params[index++], "misc.gail_rule");
    
    //fprintf (file, "%.2lf\n", misc.multi_offset/100.0);
    //sprintf (string_params[index++], "misc.multi_offset");
    parameters_complex[index++] = &misc_complex.multi_offset;
    
    //fprintf (file, "%.2lf\n", misc.multi_helix_penalty/100.0);
    //sprintf (string_params[index++], "misc.multi_helix_penalty");
    parameters_complex[index++] = &misc_complex.multi_helix_penalty;
    
    //fprintf (file, "%.2lf\n", misc.multi_free_base_penalty/100.0);
    //sprintf (string_params[index++], "misc.multi_free_base_penalty");
    parameters_complex[index++] = &misc_complex.multi_free_base_penalty;
    
    //fprintf (file, "%.2lf\n", misc.intermolecular_initiation/100.0);
    //sprintf (string_params[index++], "misc.intermolecular_initiation");
    parameters_complex[index++] = &misc_complex.intermolecular_initiation;
    
    #if (MODEL == SIMPLE)
    for(i=0; i < nb_triloops; i++)
    {
        //fprintf (file, "%.2lf\n", triloop[i].energy/100.0);
        //sprintf (string_params[index++], "triloop[%d].energy", i);
        parameters_complex[index++] = &triloop_complex[i].energy;
    }
    
    for(i=0; i < nb_tloops; i++)
    {
        //fprintf (file, "%.2lf\n", tloop[i].energy/100.0);
        //sprintf (string_params[index++], "tloop[%d].energy", i);
        parameters_complex[index++] = &tloop_complex[i].energy;
    }
    #endif
    
    if (index != num_params)
    {
        printf ("index=%d should be num_params=%d, ABORT\n", index, num_params);
        exit(1);
    }

}



int s_partition_function_complex::validate_gradient_numerical (s_partition_function *pf)
// return 1 if correct, 0 otherwise
{
    int i;
    int all_correct = 1;
    // TODO: to go all the way to num_params
    for (i=0; i < num_params; i++)
    {
        if (fabs (GlogZ_numerical[i] - pf->GlogZ[i]) > MERR)
        {
            #ifdef DOUBLEPARAMS
            printf ("Gr[%d] fd=%g\tc=%g\tr=%g,\tdiff=%g,\ttype=%s\n", i, GlogZ_finite_differences[i], GlogZ_numerical[i], pf->GlogZ[i], GlogZ_numerical[i]-pf->GlogZ[i], string_params[i]);
            #else
            printf ("Gr[%d] fd=%Lg\tc=%Lg\tr=%Lg,\tdiff=%Lg,\ttype=%s\n", i, GlogZ_finite_differences[i], GlogZ_numerical[i], pf->GlogZ[i], GlogZ_numerical[i]-pf->GlogZ[i], string_params[i]);
            #endif                        
            all_correct = 0;
        }
        //else if (i==11)
        //    printf ("Gr[%d] c=%g\tr=%g, type=%s\n", i, GlogZ_numerical[i], pf->GlogZ[i], string_params[i]);
    }
    if (all_correct)
    {
        printf ("Numerical gradient CORRECT!\n");
        return 1;
    }
    return 0;
}


int s_partition_function_complex::validate_gradient_analytical_vs_numerical ()
// return 1 if correct, 0 otherwise
{
    int i;
    int all_correct = 1;
    // TODO: to go all the way to num_params
    for (i=0; i < num_params; i++)
    {
        if (fabs (GlogZ_numerical[i] - GlogZ[i].real()) > MERR)
        {
            #ifdef DOUBLEPARAMS
            printf ("Gr[%d] c=%g\tr=%g,\tdiff=%g,\ttype=%s\n", i, GlogZ_numerical[i], GlogZ[i].real(), GlogZ_numerical[i]-GlogZ[i].real(), string_params[i]);
            #else
            printf ("Gr[%d] c=%Lg\tr=%Lg,\tdiff=%Lg,\ttype=%s\n", i, GlogZ_numerical[i], GlogZ[i].real(), GlogZ_numerical[i]-GlogZ[i].real(), string_params[i]);
            #endif                        
            all_correct = 0;
        }
        //else if (i==11)
        //    printf ("Gr[%d] c=%g\tr=%g, type=%s\n", i, GlogZ_numerical[i], pf->GlogZ[i], string_params[i]);
    }
    if (all_correct)
    {
        printf ("Numerical and analytical gradients are the SAME!\n");
        return 1;
    }
    return 0;
}


int s_partition_function_complex::validate_gradient_dp (s_partition_function *pf)
// return 1 if correct, 0 otherwise
{
    int i;
    int all_correct = 1;
    for (i=0; i < num_params; i++)
    {
        if (fabs (GlogZ[i].real() - pf->GlogZ[i]) > MERR)
        {
            printf ("Gr[%d] c=%g\tr=%g,\tdiff=%g,\ttype=%s\n", i, GlogZ[i].real(), pf->GlogZ[i], GlogZ[i].real()-pf->GlogZ[i], string_params[i]);
            all_correct = 0;
        }
        //else if (i==11)
        //    printf ("Gr[%d] c=%g\tr=%g, type=%s\n", i, GlogZ_numerical[i], pf->GlogZ[i], string_params[i]);
    }
    if (all_correct)
    {
        printf ("Dp gradient CORRECT!\n");
        return 1;
    }
    return 0;
}



s_partition_function_complex::~s_partition_function_complex ()
// The destructor
{
    delete [] sequence;
    delete [] index;     
    
    delete [] up;
    delete [] upm;
    delete [] p;
        
    IFD
    {
        delete [] u;
        delete [] u1;
        delete [] s1;
        delete [] s2;
        delete [] s3;
        delete [] pm;
        delete [] pm1;
        delete [] pm2;
    }
    else
    {        
        delete [] u_ip_jp;
        delete [] u_iu_jp;
        delete [] u_ip_ju;
        delete [] u_iu_ju;        
        delete [] s1_jp;
        delete [] s1_ju;    
        delete [] u1_ip_jp;
        delete [] u1_ip_ju_jm1p;
        delete [] u1_iu_jp;
        delete [] u1_iu_ju_jm1p;
        delete [] u1_ip_ju;
        delete [] u1_iu_ju;
        delete [] s2_jp;
        delete [] s2_ju;    
        delete [] s3_jp;
        delete [] s3_ju_jm1p;
        delete [] s3_ju;
        
        delete [] pmnod3_needmidd3;
        delete [] pmnod3_noneedmidd3;
        delete [] pmd3_needmidd3;
        delete [] pmd3_noneedmidd3;
        delete [] pm1nod3_needendd3;
        delete [] pm1d3_needendd3;

        delete [] pm2d5_needmidd5;
        delete [] pm2d5_noneedmidd5;
        delete [] pm2nod5_needmidd5;
        delete [] pm2nod5_noneedmidd5;        
    }
        
    delete [] EXPC;
    delete [] GlogZ;    
    delete [] GlogZ_numerical;
    delete [] GlogZ_finite_differences;
    delete [] parameters_complex;
    
    for (int i=0; i < num_params; i++)
    {
        delete [] HlogZ_numerical[i];
        delete [] HlogZ_finite_differences[i];
    }
    delete [] HlogZ_numerical;
    delete [] HlogZ_finite_differences;
}


void s_partition_function_complex::copy_gradient_numerical (PFTYPE *grad)
// write the gradient values into grad
{
    int i;
    for (i=0; i < num_params; i++)
    {
        grad[i] = GlogZ_numerical[i];
    }
}


void s_partition_function_complex::copy_gradient (PFTYPE *grad)
// write the gradient values into grad
{
    int i;
    for (i=0; i < num_params; i++)
    {
        grad[i] = GlogZ[i].real();
    }
}


void s_partition_function_complex::initialize_arrays ()
// U
{

    int i, j, ij;
    for (i = 0; i < seqlen; i++)
    {
        for (j = i; j < seqlen; j++)
        //for (j = i; j <= MIN(i + TURN, seqlen-1); j++)
        {
            ij = index[i] + j - i;
            up[ij] = 0;
            upm[ij] = 0;
            p[ij] = 0;
            
            IFD
            {
                u[ij] = 1;
                u1[ij] = 0;
                s1[ij] = 0;
                s2[ij] = 0;
                s3[ij] = 0;
                pm[ij] = 0;
                pm1[ij] = 0;
                pm2[ij] = 0;
            }
            else
            {
                u_ip_jp[ij] = 0;
                u_ip_ju[ij] = 0;
                u_iu_jp[ij] = 0;
                u_iu_ju[ij] = 1;
                s1_jp[ij] = 0;
                s1_ju[ij] = 0;
                            
                u1_ip_jp[ij] = 0;
                u1_ip_ju_jm1p[ij] = 0;
                u1_ip_ju[ij] = 0;
                u1_iu_jp[ij] = 0;
                u1_iu_ju_jm1p[ij] = 0;
                u1_iu_ju[ij] = 0;
                
                s2_jp[ij] = 0;
                s2_ju[ij] = 0;
                
                s3_jp[ij] = 0;
                s3_ju_jm1p[ij] = 0;
                s3_ju[ij] = 0;
                
                pmnod3_needmidd3[ij] = 0;
                pmnod3_noneedmidd3[ij] = 0;
                pmd3_needmidd3[ij] = 0;
                pmd3_noneedmidd3[ij] = 0;
                pm1nod3_needendd3[ij] = 0;
                pm1d3_needendd3[ij] = 0;
    
                pm2d5_needmidd5[ij] = 0;
                pm2d5_noneedmidd5[ij] = 0;
                pm2nod5_needmidd5[ij] = 0;
                pm2nod5_noneedmidd5[ij] = 0;        
            }            
        }
    }
}


Complex s_partition_function_complex::compute_partition_function ()
// the recursions are taken from Ding and Lawrence, "A statistical sampling algorithm for RNA secondary structure prediction", NAR 2003
{
    int i, j, ij;
    for (j=TURN+1; j < seqlen; j++)
    {
        for (i=j-TURN-1; i>=0; i--)
        {
            ij = index[i] + j - i;

            IFD
            {
                compute_upm_nodangles (i, j);    // doesn't matter where it is, all dependencies have been computed at previous steps
                compute_up (i, j);    // must be after upm            
                compute_s1 (i,j);   
                compute_u (i,j);    // must be after s1
                compute_s3 (i,j);   // Mirela: moved this before u1 on Aug 17, 2007
                compute_u1 (i,j);   // must be after s3!! 
                compute_s2 (i,j);
            }
            else
            {            
                compute_upm(i, j);    // doesn't matter where it is, all dependencies have been computed at previous steps
                compute_up (i, j);    // must be after upm                
                compute_u_ip_jp (i, j);    // must be after up
                compute_u_ip_ju (i, j);    // must be after up
                compute_u_iu_jp (i, j);    // must be after up
                compute_u_iu_ju (i, j);    // must be after up
                compute_s1_jp (i, j);    // must be after up
                compute_s1_ju (i, j);    // must be after up
                            
                compute_u1_ip_jp (i, j);    // must be after up
                compute_u1_ip_ju_jm1p (i, j);    // must be after up
                compute_u1_ip_ju (i, j);    // must be after up
                compute_u1_iu_jp (i, j);    // must be after up
                compute_u1_iu_ju_jm1p (i, j);    // must be after up
                compute_u1_iu_ju (i, j);    // must be after up
                
                compute_s2_jp (i, j);    // must be after up
                compute_s2_ju (i, j);    // must be after up
                
                compute_s3_jp (i, j);    // must be after up
                compute_s3_ju_jm1p (i, j);    // must be after up
                compute_s3_ju (i, j);    // must be after up
                
                //printf ("u[%d,%d] = %g\n", i, j, u_ip_jp[ij] + u_ip_ju[ij] + u_iu_jp[ij] + u_iu_ju[ij]);
    //                 printf ("u_ip_jp[%d,%d] = %g\n", i, j, u_ip_jp[ij]);
    //                 printf ("u_ip_ju[%d,%d] = %g\n", i, j, u_ip_ju[ij]);
    //                 printf ("u_iu_jp[%d,%d] = %g\n", i, j, u_iu_jp[ij]);
    //                 printf ("u_iu_ju[%d,%d] = %g\n", i, j, u_iu_ju[ij]);
            }            
        }
    }    

    int firstlast = index[0]+seqlen-1;                    
//     Z = u[firstlast];
//     return u[firstlast];
    IFD
    {
        Z = u[firstlast];
    }
    else
    {
        Z = u_ip_jp[firstlast] + u_ip_ju[firstlast] + u_iu_jp[firstlast] + u_iu_ju[firstlast];
    }
    
    //printf ("parition in complex = %lf\n", Z.real());
    return Z;
}


void s_partition_function_complex::validate_partition_function (s_partition_function *pf)
{
    int i, j, ij;
    for (j=TURN+1; j < seqlen; j++)
    {
        for (i=j-TURN-1; i>=0; i--)
        {
            ij = index[i]+j-i;
            if (upm[ij].real() != 0)
            {
                if ((upm[ij].real() - pf->upm[ij])/upm[ij].real() > MERR || (- upm[ij].real() + pf->upm[ij])/upm[ij].real() > MERR)
                    printf ("upm[%d,%d] c = %g\tr = %g\tdiff = %g\tratio = %g\n", i, j, upm[ij].real(), pf->upm[ij], upm[ij].real()-pf->upm[ij], (upm[ij].real() - pf->upm[ij])/upm[ij].real());
            }
            if ((up[ij].real() - pf->up[ij])/up[ij].real() > MERR || (-up[ij].real()+pf->up[ij])/up[ij].real() > MERR)
                printf ("up[%d,%d] c = %g\tr = %g\tdiff = %g\tratio = %g\n", i, j, up[ij].real(), pf->up[ij], up[ij].real()-pf->up[ij], (up[ij].real() - pf->up[ij])/up[ij].real());
            
            IFD
            {
            }
            else
            {
                if (fabs (u_ip_jp[ij].real() - pf->u_ip_jp[ij]) > MERR)
                    printf ("u_ip_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u_ip_jp[ij].real(), pf->u_ip_jp[ij], u_ip_jp[ij].real()-pf->u_ip_jp[ij]);
                if (fabs (u_ip_ju[ij].real() - pf->u_ip_ju[ij]) > MERR)
                    printf ("u_ip_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u_ip_ju[ij].real(), pf->u_ip_ju[ij], u_ip_ju[ij].real()-pf->u_ip_ju[ij]);
                if (fabs (u_iu_jp[ij].real() - pf->u_iu_jp[ij]) > MERR)
                    printf ("u_iu_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u_iu_jp[ij].real(), pf->u_iu_jp[ij], u_iu_jp[ij].real()-pf->u_iu_jp[ij]);
                if (fabs (u_iu_ju[ij].real() - pf->u_iu_ju[ij]) > MERR)
                    printf ("u_iu_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u_iu_ju[ij].real(), pf->u_iu_ju[ij], u_iu_ju[ij].real()-pf->u_iu_ju[ij]);
                if (fabs (s1_jp[ij].real() - pf->s1_jp[ij]) > MERR)
                    printf ("s1_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s1_jp[ij].real(), pf->s1_jp[ij], s1_jp[ij].real()-pf->s1_jp[ij]);
                if (fabs (s1_ju[ij].real() - pf->s1_ju[ij]) > MERR)
                    printf ("s1_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s1_ju[ij].real(), pf->s1_ju[ij], s1_ju[ij].real()-pf->s1_ju[ij]);
    
                if (fabs (u1_ip_jp[ij].real() - pf->u1_ip_jp[ij]) > MERR)
                    printf ("u1_ip_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u1_ip_jp[ij].real(), pf->u1_ip_jp[ij], u1_ip_jp[ij].real()-pf->u1_ip_jp[ij]);
                if (fabs (u1_ip_ju_jm1p[ij].real() - pf->u1_ip_ju_jm1p[ij]) > MERR)
                    printf ("u1_ip_ju_jm1p[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u1_ip_ju_jm1p[ij].real(), pf->u1_ip_ju_jm1p[ij], u1_ip_ju_jm1p[ij].real()-pf->u1_ip_ju_jm1p[ij]);
                if (fabs (u1_ip_ju[ij].real() - pf->u1_ip_ju[ij]) > MERR)
                    printf ("u1_ip_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u1_ip_ju[ij].real(), pf->u1_ip_ju[ij], u1_ip_ju[ij].real()-pf->u1_ip_ju[ij]);
                if (fabs (u1_iu_jp[ij].real() - pf->u1_iu_jp[ij]) > MERR)
                    printf ("u1_iu_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u1_iu_jp[ij].real(), pf->u1_iu_jp[ij], u1_iu_jp[ij].real()-pf->u1_iu_jp[ij]);
                if (fabs (u1_iu_ju_jm1p[ij].real() - pf->u1_iu_ju_jm1p[ij]) > MERR)
                    printf ("u1_iu_ju_jm1p[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u1_iu_ju_jm1p[ij].real(), pf->u1_iu_ju_jm1p[ij], u1_iu_ju_jm1p[ij].real()-pf->u1_iu_ju_jm1p[ij]);
                if (fabs (u1_iu_ju[ij].real() - pf->u1_iu_ju[ij]) > MERR)
                    printf ("u1_iu_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, u1_iu_ju[ij].real(), pf->u1_iu_ju[ij], u1_iu_ju[ij].real()-pf->u1_iu_ju[ij]);
    
                if (fabs (s2_jp[ij].real() - pf->s2_jp[ij]) > MERR)
                    printf ("s2_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s2_jp[ij].real(), pf->s2_jp[ij], s2_jp[ij].real()-pf->s2_jp[ij]);
                if (fabs (s2_ju[ij].real() - pf->s2_ju[ij]) > MERR)
                    printf ("s2_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s2_ju[ij].real(), pf->s2_ju[ij], s2_ju[ij].real()-pf->s2_ju[ij]);
    
                if (fabs (s3_jp[ij].real() - pf->s3_jp[ij]) > MERR)
                    printf ("s3_jp[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s3_jp[ij].real(), pf->s3_jp[ij], s3_jp[ij].real()-pf->s3_jp[ij]);
                if (fabs (s3_ju_jm1p[ij].real() - pf->s3_ju_jm1p[ij]) > MERR)
                    printf ("s3_ju_jm1p[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s3_ju_jm1p[ij].real(), pf->s3_ju_jm1p[ij], s3_ju_jm1p[ij].real()-pf->s3_ju_jm1p[ij]);
                if (fabs (s3_ju[ij].real() - pf->s3_ju[ij]) > MERR)
                    printf ("s3_ju[%d,%d] c = %g\tr = %g\tdiff = %g\n", i, j, s3_ju[ij].real(), pf->s3_ju[ij], s3_ju[ij].real()-pf->s3_ju[ij]);
            }
        }
    }
}

void s_partition_function_complex::print_partition_function()
// PRE: Z was computed. i.e. the function compute_partition_function has been called
{
    #ifdef DOUBLEPARAMS
    printf ("PF complex: %g\n", Z.real());
    #else
    printf ("PF complex: %Lg\n", Z.real());
    #endif
}


void s_partition_function_complex::print_u()
// print the u array
{
    int i, j, ij;
    for (j=TURN+1; j < seqlen; j++)
    {
        for (i=j-TURN-1; i>=0; i--)
        {
            ij = index[i] + j - i;
            printf ("u(%d,%d) = %g\n", i, j, u[ij].real());
        }
    }
}



Complex s_partition_function_complex::asymmetry_penalty_complex (int size1, int size2)
{
    Complex penalty;
    Complex m1 = misc_complex.asymmetry_penalty_max_correction;
    Complex m2 = (Complex)(abs(size1-size2)) * misc_complex.asymmetry_penalty_array [MIN (2, MIN ((size1), (size2)))-1];
    if (m1.real() < m2.real())  penalty = m1;
    else penalty = m2;
    //printf ("Asym penalty complex: %g\n", penalty.real());
    return penalty;
}


Complex s_partition_function_complex::exp_AUpenalty (int i, int j)
{
    Complex AUpen = AU_penalty_complex (sequence[i], sequence[j]);
    return (Complex)(exp (AUpen * oneoverRT));
}


// this doesn't work, for some reason. The gradient wrt AUpenalty gives 0 in this case
// Complex s_partition_function_complex::exp_AUpenalty (int i, int j)
// {
//     //Complex AUpen = AU_penalty_complex (sequence[i], sequence[j]);
//     //return (Complex)(exp (AUpen * oneoverRT));
//     // Note: doing exp(sum log) is much slower than doing exp(product)
// 
//     // this way is faster than computing exp(AUpen*oneoverRT) every time
//     if (has_AU_penalty(sequence[i], sequence[j]))
//     {
//         //printf ("Has AU penalty, eAU=%g\n", eAU.real());
//         return (Complex) eAU;
//     }
//     return (Complex) 1.0;
// }



/*
Complex s_partition_function_complex::exp_dangle5 (int i, int j, int k)
// dangle_bot
{    
    // added the next lines on Oct 7, 2006
    // make sure we are not out of bounds
    if (i < 0 || i >= seqlen || j < 0 || j >= seqlen || k < 0 || k >= seqlen)
        return (Complex) 1.0;    
    Complex dang = dangle_bot_complex[sequence[i]][sequence[j]][sequence[k]];
    return (Complex) (exp (dang * oneoverRT));
}
*/

Complex s_partition_function_complex::exp_dangle5 (int i, int j, int k)
// dangle_bot
{    
    // added the next lines on Oct 7, 2006
    // make sure we are not out of bounds
    if (i < 0 || i >= seqlen || j < 0 || j >= seqlen || k < 0 || k >= seqlen)
        return (Complex) 1.0;
    // make sure this dangling end makes sense
    if (dangle_bot_complex[sequence[i]][sequence[j]][sequence[k]].real() == INF)
        return (Complex) 1.0;
    // try calling the arrays directly, it's actually a looot faster (about 5 times on length 200)
    return edangle5[sequence[i]][sequence[j]][sequence[k]];
}



Complex s_partition_function_complex::penalty_by_size_complex (int size, char type)
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop
{
    //return (Complex)500.0;
    Complex penalty30, penalty;
    double logval;
    
    // the penalties for size <= MAXLOOP _H, _B, _I should be read from the file "loop"
    //if (size <= MAXLOOP)
    if (type == 'H' && size <= MAXLOOP_H)
    {
        //printf ("complex: size=%d, penalty=%g\n", size, hairpin_penalty_by_size_complex[size].real());
        //return (Complex)50.0;
        return hairpin_penalty_by_size_complex[size];
        //return hairpin_penalty_by_size_complex[size]/100.0;
    }
    if (type == 'I' && size <= MAXLOOP_I)
    {
        //return (Complex)50.0;
        return internal_penalty_by_size_complex[size];
    }
    if (type == 'B' && size <= MAXLOOP_B)
    {
        //return (Complex)50.0;        
        return bulge_penalty_by_size_complex[size];
    }

    //return (Complex)50.0;
    // size > MAXLOOP _H, _B, _I
    if (type == 'H')
    {
        penalty30 = hairpin_penalty_by_size_complex[MAXLOOP_H];
        logval = log (1.0*size/MAXLOOP_H);
    }        
    else if (type == 'I')
    {
        penalty30 = internal_penalty_by_size_complex[MAXLOOP_I];
        logval = log (1.0*size/MAXLOOP_I);
    }
    else
    {
        penalty30 = bulge_penalty_by_size_complex[MAXLOOP_B];
        logval = log (1.0*size/MAXLOOP_B);
    }

    // TODO: not sure if we need the (int)
    penalty = penalty30 + (Complex)(100.0*misc_complex.param_greater30 * logval);
    
    //penalty = penalty30 + (Complex)((int)(round(100.0*misc_complex.param_greater30 * logval)));
    //if (type == 'H')
    //    printf ("complex: size=%d, penalty=%g\n", size, penalty.real());

    return penalty;
}


/*
Complex s_partition_function_complex::exp_dangle3 (int i, int j, int k)
// dangle_top
{
    if (i < 0 || i >= seqlen || j < 0 || j >= seqlen || k < 0 || k >= seqlen)
        return 1.0;        
    Complex dang = dangle_top_complex[sequence[i]][sequence[j]][sequence[k]];
    return (Complex) (exp (dang * oneoverRT));
}
*/

Complex s_partition_function_complex::exp_dangle3 (int i, int j, int k)
// dangle_top
{
    if (i < 0 || i >= seqlen || j < 0 || j >= seqlen || k < 0 || k >= seqlen)
        return (Complex) 1.0;
    // make sure this dangling end makes sense
    if (dangle_top_complex[sequence[i]][sequence[j]][sequence[k]].real() == INF)
        return (Complex) 1.0;       
    // try calling the arrays directly, it's actually a looot faster (about 5 times on length 200)
    return edangle3[sequence[i]][sequence[j]][sequence[k]];    
}


void s_partition_function_complex::compute_u (int i, int j)
// called only when no dangling ends
// U
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;

    // must be initialized - took me a while to fix this bug
    u[ij] = 1.0;

    // we don't need separate cases to compute u, because we don't care about dangling ends

    for (h = i; h < j; h++)                 // case ...(...)...---
    {
        hj = index[h]+j-h;
        u[ij] += s1[hj];
    }
}


void s_partition_function_complex::compute_s1 (int h, int j)
// called only when no dangling ends
// U
{
    // if (h < 1)  return;  // we don't want this here, because we don't care about dangling ends
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s1[hj] = 0;

    // I split it into 2 for loops to avoid the if inside. Maybe this is faster.
    for (l = h+1; l < j-2; l++)        // (...)---
    {
        hl = index[h] + l - h;
        lp1j = index[l+1]+j-l-1;
        s1[hj] += up[hl] * exp_AUpenalty (h, l) * u[lp1j];
    }
    
    for (l = j-2; l <= j; l++)        // (...)---
    {
        hl = index[h] + l - h;
        s1[hj] += up[hl] * exp_AUpenalty (h, l);
    }


//     for (l = h+1; l < j; l++)        // .(...)---
//     {
//         hl = index[h] + l - h;
//         if (l+2 < j)
//         {
//             lp1j = index[l+1]+j-l-1;
//             s1[hj] += up[hl] * exp_AUpenalty (h, l) * u[lp1j];
//         }
//         else
//         {
//             s1[hj] += up[hl] * exp_AUpenalty (h, l);
//         }
//     }
    
}


void s_partition_function_complex::compute_u1 (int i, int j)
// called only when no dangling ends
// contains at least one branch of a multi-loop
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1[ij] = 0;

    for (h=i; h <= j-1; h++)    // ...(...)---
    {
        hj = index[h]+j-h;
        u1[ij] += (Complex) exp ((PARAMTYPE)(h-i) * misc_complex.multi_free_base_penalty * oneoverRT) *  //EXPC[h-i] *
            s3[hj];
    }
    u1[ij] *= (Complex) exp (misc_complex.multi_helix_penalty * oneoverRT);     //EXPB1;
}

void s_partition_function_complex::compute_s3 (int h, int j)
// called only when no dangling ends
// must contain at least one branch
// s3 doesn't contain the helix penalty
// U
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3[hj] = 0;   
    
    for (l = h+1; l <= j; l++)
    {
        hl = index[h] + l - h;        
        
        if (l+2 < j)                  
        {
            lp1j = index[l+1]+j-l-1;
            s3[hj] += up[hl] * exp_AUpenalty (h, l) * (u1[lp1j] +
                (Complex) exp ((PARAMTYPE)(j-l) * misc_complex.multi_free_base_penalty * oneoverRT));    //EXPC[j-l]);
        }
        else
        {
            s3[hj] += up[hl] * exp_AUpenalty (h, l) *
               (Complex) exp ((PARAMTYPE)(j-l) * misc_complex.multi_free_base_penalty * oneoverRT);    // EXPC[j-l];
        }
    }
}


void s_partition_function_complex::compute_upm_nodangles (int i, int j)
// i and j close a multi-loop
{
    int ij = index[i]+j-i;
    int l, h;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    Complex upm_temp;
    
    if (!can_pair(sequence[i], sequence[j]))
        return;
    upm[ij] = 0;
    // took this out for speed
    Complex upm_common;
    upm_common = exp_AUpenalty (i,j) * exp (misc_complex.multi_offset * oneoverRT) *    //EXPA *
            (Complex) exp ((PARAMTYPE)2.0*misc_complex.multi_helix_penalty * oneoverRT);   //EXPB2;

    /*
    // I don't think we need all these cases, since we don't care about dangling ends
    // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1: l+1+TURN <= j-1
    for (l=i+2; l < j-TURN-2; l++)    // case ((...)--(--)-)
    {
        ip1l = index[i+1]+l-i-1;
        lp1jm1 = index[l+1]+j-1-l-1;

        upm[ij] += up[ip1l] * exp_AUpenalty (i+1,l) * u1[lp1jm1];
    }
    for (l=i+3; l < j-TURN-2; l++)    // case [.(...)--(--)-]
    {
        ip2l = index[i+2]+l-i-2;
        lp1jm1 = index[l+1]+j-1-l-1;
        upm[ij] += up[ip2l] * EXPC[1] * exp_AUpenalty (i+2, l) * u1[lp1jm1];
    }
    */
    
    for (h=i+1; h < j-TURN-2; h++)
    //for (h=i+3; h < j-TURN-2; h++)    // case (....(...)--(--)-)
    {
        hj = index[h]+j-h;
        int hjm1 = index[h] + j-1 -h;
        upm[ij] += (Complex) exp ((PARAMTYPE)(h-i-1) * misc_complex.multi_free_base_penalty * oneoverRT) *   //EXPC[h-i-1] *
            s2[hjm1];
    }
    upm[ij] *= upm_common;
    //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
}

void s_partition_function_complex::compute_s2 (int h, int j)
// called only when no dangling ends
// must contain at least 2 branches
// helper for computing upm, which closes a multi-loop
// modified from Ding and Lawrence, to be able to add the right most d5 dangling end: (--(...)..)
// has at least 2 branches, and h is paired with some l in between h and j
// U
{
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s2[hj] = 0;
    
    for (l = h+1; l < j-3; l++)    // .(...)
    {
        hl = index[h] + l - h;
        lp1j = index[l+1]+j-l-1;
        s2[hj] += up[hl] * exp_AUpenalty (h, l) * u1[lp1j];
    }
}




void s_partition_function_complex::compute_u_ip_jp (int i, int j)
// i paired
// j unpaired or j-1 unpaired
// U
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;    
    int ijm1 = index[i]+j-1-i;

    u_ip_jp[ij] = up[ij] * exp_AUpenalty (i,j);    // case (...)
    u_ip_jp[ij] += up[ijm1] * exp_AUpenalty (i, j-1) * exp_dangle3 (j-1, i, j);    // case (...).
    
    for (l = i+1; l < j-2; l++)                    // case (...)-(--)
    {
        il = index[i]+l-i;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        u_ip_jp[ij] += up[il] * exp_AUpenalty (i, l) * (u_ip_jp[lp1j] + exp_dangle3 (l, i, l+1) * (u_ip_jp[lp2j] + u_iu_jp[lp2j]));
        //if (u_ip_jp[lp2j] != 1)            
        //    u_ip_jp[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1) * u_ip_jp[lp2j];
        //if (u_iu_jp[lp2j] != 1)
        //    u_ip_jp[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1) * u_iu_jp[lp2j];
    }                                   
}


void s_partition_function_complex::compute_u_ip_ju (int i, int j)
// i paired
// j must be unpaired and j-1 must be unpaired
// U
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;    
    
    u_ip_ju[ij] = 0;
    
    
    l=j-2;
    il = index[i]+l-i;
    u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1);
    
    for (l = i+1; l < j-2; l++)                    // case (...)...---
    {
        il = index[i]+l-i;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
    
        //u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1);    
        
        //if (u_ip_ju[lp1j] != 1)
            u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) * u_ip_ju[lp1j];
        //if (u_ip_ju[lp2j] != 1)
            u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1) * u_ip_ju[lp2j];
        //if (u_iu_ju[lp2j] != 1)
            u_ip_ju[ij] += up[il] * exp_AUpenalty (i, l) * exp_dangle3 (l, i, l+1) * u_iu_ju[lp2j];
        
    }                                   
    
}


void s_partition_function_complex::compute_u_iu_jp (int i, int j)
// i unpaired
// j paired or j-1 paired
// U
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j, hjm1;
    int h, l;    
    
    u_iu_jp[ij] = 0;
    for (h = i+1; h < j-1; h++)                    // case ...(...) or ...(...).
    {
        hj = index[h]+j-h;
        hjm1 = index[h]+j-1-h;
        //u_iu_jp[ij] += up[hj] * exp_dangle5 (j, h, h-1) * exp_AUpenalty (h, j);
        //u_iu_jp[ij] += up[hjm1] * exp_dangle5 (j-1, h, h-1) * exp_AUpenalty (h, j-1) * exp_dangle3(j-1,h,j);
    }                    
    for (h = i+1; h < j-1; h++)                 // case ...(...)...---
    {
        hj = index[h]+j-h;
        u_iu_jp[ij] += s1_jp[hj];
    }
}

void s_partition_function_complex::compute_u_iu_ju (int i, int j)
// i unpaired
// j unpaired and j-1 unpaired
// U
{
    int ij = index[i]+j-i;
    int hj, il, lp2j, lp1j;
    int h, l;    
    
    u_iu_ju[ij] = 1;
    for (h = i+1; h < j-1; h++)                 // case ...(...)...---
    {
        hj = index[h]+j-h;
        u_iu_ju[ij] += s1_ju[hj];
    }
}


void s_partition_function_complex::compute_s1_jp (int h, int j)
// j paired or j-1 paired
// U
{
    if (h < 1)  return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    int hjm1 = index[h]+j-1-h;
    s1_jp[hj] = 0;
    
    // case .(...)
    s1_jp[hj] += up[hj] * exp_dangle5 (j, h, h-1) * exp_AUpenalty (h, j);        
    // case .(...).
    s1_jp[hj] += up[hjm1] * exp_dangle5 (j-1, h, h-1) * exp_AUpenalty (h, j-1) * exp_dangle3 (j-1, h, j);        
    
    for (l = h+1; l < j-2; l++)        // .(...)---
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;        
        //if (u_ip_jp[lp1j] != 1)
            s1_jp[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * u_ip_jp[lp1j];
        //if (u_ip_jp[lp2j] != 1)
            s1_jp[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * exp_dangle3(l,h,l+1) * u_ip_jp[lp2j];
        //if (u_iu_jp[lp2j] != 1)
            s1_jp[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * exp_dangle3(l,h,l+1) * u_iu_jp[lp2j];
            
    }
}


void s_partition_function_complex::compute_s1_ju (int h, int j)
// j unpaired and j-1 unpaired
// U
{
    if (h < 1)  return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s1_ju[hj] = 0;
    
    l=j-2;
    hl = index[h] + l - h;
    s1_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * exp_dangle3 (l, h, l+1);
    
    for (l = h+1; l < j-2; l++)        // .(...)---
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;        
        //s1_ju[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * exp_dangle3(l,h,l+1);
        //if (u_ip_ju[lp1j] != 1)
            s1_ju[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * u_ip_ju[lp1j];
        //if (u_ip_ju[lp2j] != 1)
            s1_ju[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * exp_dangle3(l,h,l+1) * u_ip_ju[lp2j];
        //if (u_iu_ju[lp2j] != 1)
            s1_ju[hj] += up[hl] * exp_dangle5(l,h,h-1) * exp_AUpenalty(h,l) * exp_dangle3(l,h,l+1) * u_iu_ju[lp2j];
    }
}



void s_partition_function_complex::compute_up (int i, int j)
// need to use the complex functions here
{
    int ij = index[i]+j-i;
    int ip1jm1 = index[i+1]+j-1-i-1;
    int ip, jp, minq;
    Complex en_hairpin, en_stack, en_internal;

    up[ij] = 0;
    
    if (! can_pair (sequence[i], sequence[j]))
    {
        return;
    }    
    
    // hairpin loop
    en_hairpin = get_hairpin_energy (i, j, sequence, csequence);
    //printf ("COMPLE en_hairpin=%g, size=%d\n", en_hairpin.real(), j-i-1);
    

    // weirdly, if this has the GGG hairpin, it's added here, although it shouldn't
    // so try to remove it in case it's added
    //    that would only give me equality between the exhaustive and the true u's and up's
//     if (j-i-1 == 3 && i > 1)
//     {
//         if (sequence[i-2] == G && sequence[i-1] == G &&
//             sequence[i] == G && sequence[j] == U)
//             en_hairpin -= misc.hairpin_GGG;
//     }

    
    if (en_hairpin.real() >= INF/2)
    {
        //printf ("** Infinite hairpin (%d, %d) !\n", i, j);
    }
    else
        up[ij] += exp (en_hairpin * oneoverRT);

//     if (i==0 && j==8)    printf ("c1 up[0,8] = %g\n", up[ij].real());

    // stack pair
    if (can_pair (sequence[i+1], sequence[j-1]))    
    {
        en_stack = get_stacked_energy (i, j, sequence);
        //printf ("COMPLEX i=%d, j=%d, en_stack=%g\n", i, j, en_stack.real());
        if (en_stack.real() >= INF/2)
        {
            //printf ("** Infinite stack   (%d, %d) !\n", i, j);
        }
        else            
            up[ij] += exp (en_stack * oneoverRT) * up [ip1jm1];
    }             
//     if (i==0 && j==8)    printf ("c2 up[0,8] = %g\n", up[ij].real());            

    // TODO: remove return;
    if (!ignore_internal)
    {
        // internal loop/bulge
        for (ip = i+1; ip <= MIN(j-2-TURN,i+MAXLOOP+1) ; ip++)  // j-2-TURN
        {
            minq = MAX (j-i+ip-MAXLOOP-2, ip+1+TURN);    // ip+1+TURN);
            for (jp = minq; jp < j; jp++)
            {        
                if (sequence[ip]+sequence[jp] == 3 ||
                    sequence[ip]+sequence[jp] == 5)        
                {
                    if (ip == i+1 && jp == j-1) continue;    // we don't want stacked pairs here
                    int ipjp = index[ip] + jp - ip;
                    en_internal = get_internal_energy (i, j, ip, jp, sequence);
                    //en_internal = (Complex) s_internal_loop::get_energy (i, j, ip, jp, sequence);
                    if (en_internal.real() >= INF/2)
                    {
                        //printf ("** Infinite internal(%d, %d) !\n", i, j);
                    }
                    else
                    {
                        up[ij] += exp (en_internal * oneoverRT) * up[ipjp];
    //                     if (i==0 && j==8)    printf ("c3 ip=%d, jp=%d, up[0,8] = %g, en_internal=%g\n", ip, jp, up[ij].real(), en_internal.real());
                    }
                }
            }        
        }
    }

    if (ignore_multi)
        return;
        
    // multi-loop
    up[ij] += upm[ij];
//     if (i==0 && j==8)    printf ("c4 up[0,8] = %g\n", up[ij].real());
}


void s_partition_function_complex::compute_upm (int i, int j)
// i and j close a multi-loop
{
    int ij = index[i]+j-i;
    int l, h;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    Complex upm_temp;
    
    if (!can_pair(sequence[i], sequence[j]))
        return;
    upm[ij] = 0;
    // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1: l+1+TURN <= j-1
    for (l=i+2; l < j-TURN-2; l++)    // case ((...)--(--)-)
    {
        ip1l = index[i+1]+l-i-1;
        lp2jm1 = index[l+2]+j-1-l-2;
        lp1jm1 = index[l+1]+j-1-l-1;

//         if (ignore_dangles)
//             upm[ij] += exp_AUpenalty (i,j) * up[ip1l] * EXPA * EXPB2 *
//                         exp_AUpenalty (i+1,l) *
//                         u1[lp1jm1];
//         else        
            upm[ij] += exp_AUpenalty (i,j) * up[ip1l] * EXPA * EXPB2 *
                        exp_AUpenalty (i+1,l) *
                    (  u1_ip_jp[lp1jm1]      // [(...)(.-..)] or [(...)(.-..).]                    
                       + exp_dangle3 (l, i+1, l+1) * EXPC[1] *
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1] +     // [(...).-(...)] or [(...).-(...).] 
                                exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])) +    // [(...).-(...)-..]
                       + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);    // [(...)(...)-..]
    }
    for (l=i+3; l < j-TURN-2; l++)    // case [.(...)--(--)-]
    {
        ip2l = index[i+2]+l-i-2;
        lp2jm1 = index[l+2]+j-1-l-2;
        lp1jm1 = index[l+1]+j-1-l-1;
        lp2jm2 = index[l+2]+j-2-l-2;
        lp1jm2 = index[l+1]+j-2-l-1;        
        lp2jm3 = index[l+2]+j-3-l-2;
        lp1jm3 = index[l+1]+j-3-l-1;      
//         if (ignore_dangles)
//             upm[ij] += exp_AUpenalty (i,j) * up[ip2l] * EXPA * EXPB2 * EXPC[1] *
//                     exp_AUpenalty (i+2, l) * 
//                     u1[lp1jm1];                    
//         else
            upm[ij] += exp_AUpenalty (i,j) * up[ip2l] * EXPA * EXPB2 * EXPC[1] *
                    exp_dangle3 (i, j, i+1) * exp_AUpenalty (i+2, l) * 
                    (  u1_ip_jp[lp1jm1]      // [.(...)(.-..)] or [.(...)(.-..).]                    
                       + exp_dangle3 (l, i+2, l+1) * EXPC[1] *
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1] +     // [.(...).-(...)] or [.(...).-(...).] 
                                exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])) +    // [.(...).-(...)-..]
                       + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);    // [.(...)(...)-..] 
    }
    upm_temp = 0;
    for (h=i+3; h < j-TURN-2; h++)    // case (....(...)--(--)-)
    {
        hj = index[h]+j-h;
        int hjm1 = index[h] + j-1 -h;
        int hjm2 = index[h] + j-2 -h;
        int hjm3 = index[h] + j-3 -h;
//         if (ignore_dangles)
//             upm_temp += EXPA * EXPB2 * EXPC[h-i-1] * s2[hjm1];
//         else
            upm_temp += EXPA * EXPB2 * EXPC[h-i-1] *
                ( s2_jp[hjm1]    // --(...)) or --(...).)
                  + s2_ju[hjm1] * exp_dangle5(i,j,j-1));    // --(...)..)
    }   
//     if (ignore_dangles) 
//         upm[ij] += exp_AUpenalty (i,j) * upm_temp;
//     else
        upm[ij] += exp_AUpenalty (i,j) * exp_dangle3 (i, j, i+1) * upm_temp;
    //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
}


void s_partition_function_complex::compute_s2_jp (int h, int j)
// helper for computing upm, which closes a multi-loop
// modified from Ding and Lawrence, to be able to add the right most d5 dangling end: (--(...)..)
// has at least 2 branches, and h is paired with some l in between h and j
// U
{
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s2_jp[hj] = 0;
    
    for (l = h+1; l < j-3; l++)    // .(...)
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        s2_jp[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                ( u1_ip_jp[lp1j] + exp_dangle3 (l, h, l+1) * EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
                //u1[lp2j] + u1nod5[lp1j] - EXPC[1]*u1nod5[lp2j] );
    }
}



void s_partition_function_complex::compute_s2_ju (int h, int j)
// helper for computing upm, which closes a multi-loop
// modified from Ding and Lawrence, to be able to add the right most d5 dangling end: (--(...)..)
// has at least 2 branches, and h is paired with some l in between h and j
// U
{
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s2_ju[hj] = 0;
    
    for (l = h+1; l < j-3; l++)    // .(...)
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        s2_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                ( u1_ip_ju[lp1j] + exp_dangle3 (l, h, l+1) * EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
                //u1[lp2j] + u1nod5[lp1j] - EXPC[1]*u1nod5[lp2j] );
    }
}


void s_partition_function_complex::compute_u1_ip_jp (int i, int j)
// contains at least one branch of a multi-loop
// i must be paired
// j must be paired, or j-1 must be paired
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_ip_jp[ij] = 0;


            
    for (l=j-1; l <= j; l++)    // (...) or (...).
    {
        il = index[i]+l-i;        
        u1_ip_jp[ij] += up[il] * EXPB1 *
            exp_AUpenalty(i,l) * ( fd3(j+1,i,l) * EXPC[j-l] );
    }
                
    for (l=i+1; l < j-2; l++)    // (...)-(---) or (...)-(---).
    {
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;                      
        il = index[i]+l-i;        
        u1_ip_jp[ij] += up[il] * EXPB1 *
                    exp_AUpenalty(i,l) *
                    (u1_ip_jp[lp1j] + exp_dangle3 (l,i,l+1)*EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }
}


void s_partition_function_complex::compute_u1_ip_ju_jm1p (int i, int j)
// contains at least one branch of a multi-loop
// i must be paired
// j must be paired, or j-1 must be paired
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_ip_ju_jm1p[ij] = 0;
    int ijm1 = index[i]+j-1-i;       
    u1_ip_ju_jm1p[ij] += up[ijm1] * EXPB1 *
            exp_AUpenalty(i,j-1)* fd3(j+1,i,j-1) * EXPC[1];
                
    for (l=i+1; l < j-2; l++)    // (...)-(---).
    {
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;                      
        il = index[i]+l-i;        
        u1_ip_ju_jm1p[ij] += up[il] * EXPB1 *
                    exp_AUpenalty(i,l) *
                    (u1_ip_ju_jm1p[lp1j] + exp_dangle3 (l,i,l+1)*EXPC[1]*(u1_ip_ju_jm1p[lp2j] + u1_iu_ju_jm1p[lp2j]));
    }
}


void s_partition_function_complex::compute_u1_ip_ju (int i, int j)
// contains at least one branch of a multi-loop
// i must be paired
// j must be unpaired, and j-1 must be unpaired
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_ip_ju[ij] = 0;

    for (l=i+1; l <= j-2; l++)    // (...)----
    {
        il = index[i]+l-i;        // (...)....
        u1_ip_ju[ij] += up[il] * EXPB1 *
                exp_AUpenalty(i,l) *
                ( fd3(j+1,i,l) * EXPC[j-l] );
                  
        if (l+2 < j)            // (...)-(--)-
        {          
            lp1j = index[l+1]+j-l-1;
            lp2j = index[l+2]+j-l-2;                  
            u1_ip_ju[ij] += up[il] * EXPB1 *
                        exp_AUpenalty(i,l) *
                        (u1_ip_ju[lp1j] + exp_dangle3 (l,i,l+1)*EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
        }                        
    }
}

    
void s_partition_function_complex::compute_u1_iu_jp (int i, int j)
// contains at least one branch of a multi-loop
// i must be unpaired
// j must be paired, or j-1 must be paired
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_iu_jp[ij] = 0;

    for (l=j-1; l <= j; l++)    // .(...) or .(...).
    {
        ip1l = index[i+1]+l-i-1;        
        u1_iu_jp[ij] += up[ip1l] * EXPB1 *
                exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                ( fd3(j+1,i+1,l) * EXPC[j-l] );
    }
    
    for (l=i+2; l < j-2; l++)    // .(...)-(---) or .(...)-(---).
    {
        ip1l = index[i+1]+l-i-1;
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;  
        u1_iu_jp[ij] += up[ip1l] * EXPB1 *
                    exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                    (u1_ip_jp[lp1j] + exp_dangle3 (l,i+1,l+1)*EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }

    for (h=i+2; h <= j-1; h++)    // ..-(...)-(--) or ..-(...)-(--).
    {
        hj = index[h]+j-h;
        // d5 is included in s3, but helix penalty is not 
        u1_iu_jp[ij] += EXPB1 *
            EXPC[h-i] * s3_jp[hj];
    }
}


void s_partition_function_complex::compute_u1_iu_ju_jm1p (int i, int j)
// contains at least one branch of a multi-loop
// i must be unpaired
// j must be paired, or j-1 must be paired
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_iu_ju_jm1p[ij] = 0;

    int ip1jm1 = index[i+1]+j-1-i-1;        // .(...).
    u1_iu_ju_jm1p[ij] += up[ip1jm1] * EXPB1 *
        exp_AUpenalty(i+1,j-1) * exp_dangle5 (j-1, i+1, i) * fd3(j+1,i+1,j-1) * EXPC[2];
    
    for (l=i+2; l < j-2; l++)    // .(...)-(---) or .(...)-(---).
    {
        ip1l = index[i+1]+l-i-1;
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;  
        u1_iu_ju_jm1p[ij] += up[ip1l] * EXPB1 * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                    (u1_ip_ju_jm1p[lp1j] + exp_dangle3 (l,i+1,l+1)*EXPC[1]*(u1_ip_ju_jm1p[lp2j] + u1_iu_ju_jm1p[lp2j]));
    }

    for (h=i+2; h <= j-1; h++)    // ..-(...)-(--) or ..-(...)-(--).
    {
        hj = index[h]+j-h;
        // d5 is included in s3, but helix penalty is not 
        u1_iu_ju_jm1p[ij] += EXPB1 * EXPC[h-i] * s3_ju_jm1p[hj];
    }
}



void s_partition_function_complex::compute_u1_iu_ju (int i, int j)
// contains at least one branch of a multi-loop
// i must be unpaired
// j must be unpaired, and j-1 must be unpaired
// U
{
    int ij, il, l, h, lp1j, lp2j, ip1l, hj;
    ij = index[i]+j-i;
    u1_iu_ju[ij] = 0;

    for (l=i+2; l < j-1; l++)    // .(...)-..
    {
        ip1l = index[i+1]+l-i-1;        
        u1_iu_ju[ij] += up[ip1l] * EXPB1 * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                ( fd3(j+1,i+1,l) * EXPC[j-l] );
    }
        
    for (l=i+2; l < j-2; l++)    // .(...)-(--)-..
    {
        ip1l = index[i+1]+l-i-1;
        lp1j = index[l+1]+j-l-1;
        lp2j = index[l+2]+j-l-2;                 
        u1_iu_ju[ij] += up[ip1l] * EXPB1 * exp_AUpenalty(i+1,l) * exp_dangle5 (l, i+1, i) * EXPC[1] *
                    (u1_ip_ju[lp1j] + exp_dangle3 (l,i+1,l+1)*EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
    }

    for (h=i+2; h <= j-1; h++)    // ..-(...)-(--)-..
    {
        hj = index[h]+j-h;
        // d5 is included in s3, but helix penalty is not 
        u1_iu_ju[ij] += EXPB1 * EXPC[h-i] * s3_ju[hj];
    }
}



Complex s_partition_function_complex::fd3 (int jplus1, int h, int l)
{
    if (l > jplus1-1)    giveup ("Error, l > j", "f function, partition_function");
    if (l == jplus1-1)
        return 1.0;
    return exp_dangle3 (l, h, l+1);    
}



void s_partition_function_complex::compute_s3_jp (int h, int j)
// s3 doesn't contain the helix penalty
// j paired or j-1 paired
// U
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3_jp[hj] = 0;   

    for (l = j-1; l <= j; l++)
    {
        hl = index[h] + l - h;
        s3_jp[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                    ( fd3 (j+1, h, l) * EXPC[j-l] );
    }
    for (l = h+1; l < j-2; l++)
    {
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        hl = index[h] + l - h;
        s3_jp[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) *
                (u1_ip_jp[lp1j] + exp_dangle3 (l,h,l+1)*EXPC[1]*(u1_ip_jp[lp2j] + u1_iu_jp[lp2j]));
    }
}


void s_partition_function_complex::compute_s3_ju_jm1p (int h, int j)
// s3 doesn't contain the helix penalty
// j paired or j-1 paired
// U
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3_ju_jm1p[hj] = 0;   

    int hjm1 = index[h]+j-1-h;
    s3_ju_jm1p[hj] += up[hjm1] * exp_dangle5 (j-1, h, h-1) * exp_AUpenalty (h, j-1) * fd3 (j+1, h, j-1) * EXPC[1];

    for (l = h+1; l < j-2; l++)
    {
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        hl = index[h] + l - h;
        s3_ju_jm1p[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) *
                (u1_ip_ju_jm1p[lp1j] + exp_dangle3 (l,h,l+1)*EXPC[1]*(u1_ip_ju_jm1p[lp2j] + u1_iu_ju_jm1p[lp2j]));
    }
}


void s_partition_function_complex::compute_s3_ju (int h, int j)
// s3 doesn't contain the helix penalty
// j unpaired and j-1 unpaired
// U
{
    if (h < 1) return;
    int hj, hl, l, lp2j, lp1j;
    hj = index[h]+j-h;
    s3_ju[hj] = 0;   
    
    for (l = h+1; l < j-1; l++)
    {
        hl = index[h] + l - h;
        s3_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) * 
                    ( fd3 (j+1, h, l) * EXPC[j-l] );
    }
    
    for (l = h+1; l < j-2; l++)
    {
        hl = index[h] + l - h;
        lp2j = index[l+2]+j-l-2;
        lp1j = index[l+1]+j-l-1;
        s3_ju[hj] += up[hl] * exp_dangle5 (l, h, h-1) * exp_AUpenalty (h, l) *
                (u1_ip_ju[lp1j] + exp_dangle3 (l,h,l+1)*EXPC[1]*(u1_ip_ju[lp2j] + u1_iu_ju[lp2j]));
    }
}




void s_partition_function_complex::compute_base_pair_probabilities ()
// Nov 9, 2006. Computes base pair probabilities
// PRE: the partition function arrays have been filled
// U
{
    int h, l;
    for (h=0; h < seqlen; h++)   
    {
        for (l=seqlen-1; l > h+TURN; l--)
        {
            if (can_pair (sequence[h], sequence[l]))    
            {
                compute_p (h, l);
            }
            IFD
            {
                compute_pm (h, l);
                compute_pm1 (h, l);
                // moved to the gradient
                //compute_pm2 (h, l);
            }
            else
             {
                compute_pmnod3_needmidd3 (h, l);
                compute_pmnod3_noneedmidd3 (h, l);
                compute_pmd3_needmidd3 (h, l);
                compute_pmd3_noneedmidd3 (h, l);
                compute_pm1nod3_needendd3 (h, l);
                compute_pm1d3_needendd3 (h, l);            
                
                // the next 4 are only needed to compute the partial derivative wrt multi_free_base_penalty
                // moved to the gradient                
//                 compute_pm2d5_needmidd5 (h, l);
//                 compute_pm2d5_noneedmidd5 (h, l);
//                 compute_pm2nod5_needmidd5 (h, l);
//                 compute_pm2nod5_noneedmidd5 (h, l);
            }
        }
    }
}


void s_partition_function_complex::print_base_pair_probabilities (PFTYPE threshold)
//prints all the base pair probabilities above the given threshold
{
    int i, j;
    printf ("Base pair probabilities above threshold %.2lf\nIndeces start from 0.\n", threshold);
    printf ("i\tj\tprobability\n");
    for (i=0; i < seqlen-TURN-1; i++)
    {
        for (j=i+TURN+1; j < seqlen; j++)
        {
            int ij = index[i]+j-i;
            if (p[ij].real() >= threshold)
            {
                printf ("%d\t%d\t%g\n", i, j, p[ij].real());
            }
        }
    }
}



void s_partition_function_complex::compute_p (int h, int l)
{
    int hl, zeronminus1, hm1lp1;
    Complex en_stack;
    int i,j;
    Complex en_internal;
    Complex term1, term2;
    
    hl = index[h]+l-h;
    
    // initialize, just in case
    p[hl] = 0;

    //exterior base pair
    // first, 5' end
    term1 = 0.0;
    
    IFD
    {
        if (h > 0)      term1 = u[h-1];
        else            term1 = 1;
    }
    else
    {
        if (h > TURN)
        {
            term1 +=  (u_ip_jp[h-1] + u_iu_jp[h-1] + exp_dangle5(l,h,h-1)*(u_ip_ju[h-1] + u_iu_ju[h-1]));
        }
        else if (h > 0)       
        {
            term1 += exp_dangle5(l,h,h-1);  //-.[...]
        }
        else                  term1 = 1;  
    }
    
    // then, 3' end
    term2 = 0.0;
    IFD
    {
        if (l < seqlen-1)   term2 = u[index[l+1]+seqlen-1-(l+1)];
        else                term2 = 1;
    }
    else
    {
        if (l < seqlen-3)
        {
            int lp1n = index[l+1]+seqlen-1-(l+1);
            int lp2n = index[l+2]+seqlen-1-(l+2);            
            term2 += (u_ip_jp[lp1n] + u_ip_ju[lp1n] + 
                        exp_dangle3(l,h,l+1)*(u_ip_jp[lp2n] + u_ip_ju[lp2n] + u_iu_jp[lp2n] + u_iu_ju[lp2n]));
        }
        else if (l < seqlen-1) 
        {    
            term2 += exp_dangle3(l,h,l+1);
            //printf ("Add4 d3(%d,%d,%d)\n", l, h, l+1);
        }
        else                  term2 = 1;
    }
    
    zeronminus1 = index[0] + seqlen-1;
    IFD
    {
        p[hl] = term1 * up[hl] * term2 * exp_AUpenalty (h, l) / u[zeronminus1];
    }
    else
    {
        p[hl] = term1 * up[hl] * term2 * exp_AUpenalty (h, l) /
            (u_ip_jp[zeronminus1] + u_ip_ju[zeronminus1] + u_iu_jp[zeronminus1] + u_iu_ju[zeronminus1]);
    }
    
    // case stacking energies
    if (h > 0 && l < seqlen-1)
    {
        hm1lp1 = index[h-1] + l+1 - h+1;
        if (can_pair (sequence[h-1], sequence[l+1]))    
        {
            en_stack = get_stacked_energy (h-1, l+1, sequence);
            if (en_stack.real() >= INF/2)
            {
                //printf ("** Infinite stack   (%d, %d) !\n", h-1, l+1);
            }
            else            
                p[hl] += p[hm1lp1] * up[hl] / up[hm1lp1] * exp (en_stack * oneoverRT);
        }
    }
    
    // case internal loop
    for (i = MAX(0, h-MAXLOOP-1); i < h; i++)
    {
        for (j = l+1; j <= MIN(MAXLOOP+l-h+i-2, seqlen-1); j++)
        {        
            if (sequence[i]+sequence[j] == 3 ||
                sequence[i]+sequence[j] == 5)        
            {
                if (i == h-1 && j == l+1) continue;    // we don't want stacked pairs here
                int ij = index[i] + j - i;
                en_internal = get_internal_energy (i, j, h, l, sequence);
                if (en_internal.real() >= INF/2)
                {
                    //printf ("** Infinite internal(%d, %d) !\n", i, j);
                }
                else
                    p[hl] += p[ij] * up[hl] /up[ij] * exp (en_internal * oneoverRT);
            }
        }        
    }

    // case multi-loop
    //     u1 has multi_helix_penalty                    
        

    Complex pml = 0;
    for (i=0; i < h; i++)
    {
        int il = index[i] + l - i;
                
        // the case when h-l is the first branch of the multi-loop    //  (-[...] ..      
        IFD
            pml += EXPC[h-i-1] * pm[il];
        else
        {
            if (i < h-1)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
            {
                 pml += EXPC[h-i-1] * (i<h-2?exp_dangle5(l,h,h-1):1) *
                    (pmd3_noneedmidd3[il] + pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1]);
                    // (.-[...](---)-)    i.-h...l(---)-j 
                    // (.-[...].-(---)-)    i.-h...l.-(---)-j
            }
            else    // case ((... , no dangling end
            {
                pml += pmnod3_noneedmidd3[il] + pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1];
                    // ([...](---)-)    ih...l(---)-j 
                    // ([...].-(---)-)    ih...l.-(---)-j
            }
        }
                
        if (i < h - TURN -2)    // only now we can have a branch to the left
        {      
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            IFD
            {
                pml += (pm1[il] + pm[il]) * u1[ip1hm1];
                // no branch to the right of h-l
                // branch to the left and to the right of h-l
            }
            else
            {   
                int ip2hm1 = index[i+2] + h-1 - (i+2);
                int ilp1 = index[i] + l+1 - i;
                Complex term1 = 0.0;
                if (up[ilp1].real() != 0)    // j is l+1
                {
                    term1 = p[ilp1] / up[ilp1] * exp_AUpenalty (i,l+1) *
                                (  // first, the case ((..-)-[...])
                                  ( u1_ip_jp[ip1hm1]    // ((..-)[...]) or ((..-).[...])
                                    + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1] )   // ((..-).-[...])  
                                        // no need to add EXPC[1], it's in u1_ip_ju
                                        
                                    // next, the case (.-(...)-[...])      i..-(...)h...lj
                                  + exp_dangle3(i,l+1,i+1) * EXPC[1] *
                                    ( u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1] 
                                            + exp_dangle5(l,h,h-1)*(u1_ip_ju[ip2hm1]+u1_iu_ju[ip2hm1]) ) 
                                        // (..-(..-).-[...])     
                                );
                }
                
                pml +=
                
                // first, when h-l is the last branch to the right
                ( term1 + 
                
                // case ((..-)-[...].-)  
                pm1nod3_needendd3[il] * exp_dangle3(l,h,l+1) *    // EXPC added in pm1nod3
                    (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                            // c is added in pm1nod3_needendd3
                            // ((...).-[...].-)         i(...)h...l.-j
                // case (.(..-)-[...].-)
                + pm1d3_needendd3[il] * exp_dangle3(l,h,l+1) *     // don't need EXPC[1]
                    (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)* u1_ip_ju[ip2hm1])
                            // (.(...).-[...].-)        i.(...)h...l.-j
                            
                // case (..-(..-)-[...].-)            
                + pm1d3_needendd3[il]  * exp_dangle3(l,h,l+1) *     // don't need EXPC[1]
                    (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1])
                            // (..-(...)..-[...].-)      i..-(...)h...l.-j

                            
                // we have branches to the left and to the right of h-l
                // let's do same as above, but with pm instead of pm1
                                // case ((..-)-[...].-(--)-) or  ((..-)-[...](--)-)
                + (pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1] + pmnod3_noneedmidd3[il]) *
                    (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                            // c is added in pm1nod3_needendd3
                            // ((...).-[...].-(--)-)         i(...)h...l.-j
                // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                + EXPC[1] * (pmd3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il]) *
                        (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip2hm1])
                            // (.(...).-[...].-)        i.(...)h...l.-j
                            
                // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                + EXPC[1] * (pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il])*
                    (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1]) ); 
                            
            }    // end if-else ignore_dangles
        }      // end if (i < h - TURN -2)  
    }
    p[hl] +=  up[hl] * EXPA * EXPB2 * exp_AUpenalty (h,l) * pml;
}


void s_partition_function_complex::compute_pm (int i, int l)
// called only when no dangling ends
//  region l+1  - j-1 has at least one branch
//  i.-h    l-(---)-j
//  (.-(    )-(---)-)
// U
{
    int j;
    int ij, lp1jm1, il;
    il = index[i]+l-i;
    // intialize, just in case
    pm[il] = 0;
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp1jm1 = index[l+1] + j-1 - l-1;
            pm[il] += exp_AUpenalty (i,j) * p[ij] / up[ij] * u1[lp1jm1];            
        }
    }    
}


void s_partition_function_complex::compute_pmd3_needmidd3 (int i, int l)
// add the dangling end d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  i.-h    l-(---)-j
//  (.-(    )-(---)-)
// U
{
    int j;
    int ij, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;
    // intialize, just in case
    pmd3_needmidd3[il] = 0;
    
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;            
            pmd3_needmidd3[il] += exp_AUpenalty (i,j) * exp_dangle3(i,j,i+1) * p[ij] / up[ij] * // don't need to add EXPC[1] here, it's added in the calling function
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1]      // ..))  or ..).)
                                + exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1])  );    //)-..)            
        }
    }
}





void s_partition_function_complex::compute_pmd3_noneedmidd3 (int i, int l)
// add the dangling end d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  i.-h    l-(---)-j
//  (.-(    )-(---)-)
// U
{
    int j;
    int ij, lp1jm1, lp1jm2, lp1jm3, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;
    // intialize, just in case
    pmd3_noneedmidd3[il] = 0;
    
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp1jm1 = index[l+1] + j-1 - l-1;
            lp1jm2 = index[l+1] + j-2 - l-1;
            lp1jm3 = index[l+1] + j-3 - l-1;            
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;

            pmd3_noneedmidd3[il] += exp_AUpenalty (i,j) * exp_dangle3(i,j,i+1) * p[ij] / up[ij] * // shouldn't add EXPC[1] here because it's in the calling function
                            (u1_ip_jp[lp1jm1] + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);  
                                // ..))  or ..).),  and then )-..)            
        }
    }
}

void s_partition_function_complex::compute_pmnod3_needmidd3 (int i, int l)
// the pm from McCaskill, which does not have d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  ih    l-(---)-j
//  ((    )-(---)-)
// U
{
    int j;
    int ij, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;
    // intialize, just in case
    pmnod3_needmidd3[il] = 0;
    
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;
            pmnod3_needmidd3[il] += exp_AUpenalty (i,j) * p[ij] / up[ij] * // no need to add EXPC[1] here, it's added in the calling function
                            ( u1_ip_jp[lp2jm1] + u1_iu_jp[lp2jm1]      // ..))  or ..).)
                                + exp_dangle5(i,j,j-1) * (u1_ip_ju[lp2jm1] + u1_iu_ju[lp2jm1]));    //)-..)
        }
    }
}


void s_partition_function_complex::compute_pmnod3_noneedmidd3 (int i, int l)
// the pm from McCaskill, which does not have d3(i,j,i+1)
//  region l+1  - j-1 has at least one branch
//  ih    l(---)-j
//  ((    )(---)-)    // right branch follows right after, no free base in between
// U
{
    int j;
    int ij, lp1jm1, lp1jm2, lp1jm3, lp2jm1, lp2jm2, lp2jm3, il;
    il = index[i]+l-i;
    // intialize, just in case
    pmnod3_noneedmidd3[il] = 0;
    
    for (j = l+TURN+3; j < seqlen; j++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            lp1jm1 = index[l+1] + j-1 - l-1;
            lp1jm2 = index[l+1] + j-2 - l-1;
            lp1jm3 = index[l+1] + j-3 - l-1;
            lp2jm1 = index[l+2] + j-1 - l-2;
            lp2jm2 = index[l+2] + j-2 - l-2;
            lp2jm3 = index[l+2] + j-3 - l-2;            

            pmnod3_noneedmidd3[il] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                            (u1_ip_jp[lp1jm1] + exp_dangle5(i,j,j-1) * u1_ip_ju[lp1jm1]);                               
        }
    }
}


void s_partition_function_complex::compute_pm1 (int i, int l)
// called only when no dangling ends
// the pm1 from McCaskill, no dangling ends
// region l+1 ... j-1 is unpaired
//  i h    l.-j   
//  ( (    ).-)
// U
{
    int j;
    int ij, il;
    il = index[i]+l-i;
    // intialize, just in case
    pm1[il] = 0;
    
    // the case j=l+1 is dealt with separately, directly in compute_p
    for (j = l+1; j < seqlen; j++)    
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;           
            pm1[il] += p[ij] / up[ij] * exp_AUpenalty (i,j) * EXPC[j-l-1];
        }
    }
}


void s_partition_function_complex::compute_pm1nod3_needendd3 (int i, int l)
// all free bases between l and j
// the pm1 from McCaskill, which does not have d3(i,j,i+1)
// region l+1 ... j-1 is unpaired
//  i h    l.-j    - we need to add d3(l,h,l+1)
//  ( (    ).-)
// U
{
    int j;
    int ij, il;
    il = index[i]+l-i;
    // intialize, just in case
    pm1nod3_needendd3[il] = 0;
        
    // the case j=l+1 is dealt with separately, directly in compute_p
    for (j = l+2; j < seqlen; j++)    
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;           
            pm1nod3_needendd3[il] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                                     EXPC[j-l-1] * (l+2<j?exp_dangle5(i,j,j-1):1);
        }
    }
}


void s_partition_function_complex::compute_pm1d3_needendd3 (int i, int l)
// all free bases between l and j
// add the dangling end d3(i,j,i+1)
// region l+1 ... j-1 is unpaired
//  i h    l........j
//  ( (    )........)
// U
{
    int j;
    int ij, il;
    il = index[i]+l-i;
    // intialize, just in case
    pm1d3_needendd3[il] = 0;
    
    // the case j=l+1 is dealt with separately, directly in compute_p
    for (j = l+2; j < seqlen; j++)    
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;           
            pm1d3_needendd3[il] += exp_dangle3(i,j,i+1) * EXPC[1] * p[ij] / up[ij] * exp_AUpenalty (i,j) *
                                   EXPC[j-l-1] * (l+2<j?exp_dangle5(i,j,j-1):1);
        }
    }
}


void s_partition_function_complex::compute_pm2 (int h, int j)
// called only when no dangling ends
//  region i+1  - h-1 has at least one branch
// the 2 free bases are not aplicable here because we don't care about dangling ends
// i---------h...l.-.j
// (-(---)-..[...].-.)
// U
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    
    // intialize, just in case
    pm2[hj] = 0;
    
    for (i=0; i < h-TURN-2; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            pm2[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * u1[ip1hm1];
        }
    }
}


void s_partition_function_complex::compute_pm2d5_needmidd5 (int h, int j)
//  region i+1  - h-3 has at least one branch, and there must at least 2 free bases before h
// also add the dangling end d5(i,j,j-1)
// i---------h...l.-.j
// (-(---)-..[...].-.)
// U
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    // intialize, just in case
    pm2d5_needmidd5[hj] = 0;
    
    for (i=0; i < h-TURN-3; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2d5_needmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * exp_dangle5(i,j,j-1) *
                // the following is the same in compute_pm2nod5_needmidd5
                // we have to make sure the middle d3 is included
                // first, ((---)-..[...].-.)
                ( u1_ip_ju[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_ju[ip2hm1] + u1_iu_ju[ip2hm1]));
        }
    }
}


void s_partition_function_complex::compute_pm2d5_noneedmidd5 (int h, int j)
//  region i+1  - h-1 has at least one branch, and there must be 0 or 1 free bases left of h
// also add the dangling end d5(i,j,j-1)
// i------.h...l.-.j
// (-(---).[...].-.) or (-(---)[...].-.)
// U
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    // intialize, just in case
    pm2d5_noneedmidd5[hj] = 0;
    
    for (i=0; i < h-TURN-2; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2d5_noneedmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * exp_dangle5(i,j,j-1) * 
                // first, ((---)[...].-.) or ((---).[...].-.)
                ( u1_ip_jp[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1]));
        }
    }
}


void s_partition_function_complex::compute_pm2nod5_needmidd5 (int h, int j)
//  region i+1  - h-3 has at least one branch, and there must at least 2 free bases before h
// DO NOT add the dangling end d5(i,j,j-1)
// i---------h...l.j
// (-(---)-..[...]) or (-(---)-..[...].)
// exactly the same as compute_pm2d5_needmidd5, but withour exp_dangle5
// U
{
    int i, ij;
    int hj, ip1;
    hj = index[h]+j-h;   
    // intialize, just in case
    pm2nod5_needmidd5[hj] = 0;
    
    for (i=0; i < h-TURN-2; i++)
    {
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2nod5_needmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                // we have to make sure the middle d3 is included
                // first, ((---)-..[...]) or ((---)-..[...].)
                ( u1_ip_ju[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_ju[ip2hm1] + u1_iu_ju[ip2hm1]));
        }
    }
}


void s_partition_function_complex::compute_pm2nod5_noneedmidd5 (int h, int j)
//  region i+1  - h-1 has at least one branch, and there must be 0 or 1 free bases left of h
// also add the dangling end d5(i,j,j-1)
// i------.h...l.-.j
// (-(---).[...]) or (-(---)[...])         or   (-(---).[...].) or (-(---)[...].)
// same as compute_pm2d5_noneedmidd5, just remove the exp_dangle5
// U
{
    int i;
    int hj, ip1, ij;
    hj = index[h]+j-h;   
    // intialize, just in case
    pm2nod5_noneedmidd5[hj] = 0;
    
    for (i=0; i < h-TURN-2; i++)
    {
    
        if (can_pair (sequence[i], sequence[j]))
        {
            ij = index[i]+j-i;   
            int ip1hm1 = index[i+1] + h-1 - (i+1);
            int ip2hm1 = index[i+2] + h-1 - (i+2);
            pm2nod5_noneedmidd5[hj] += p[ij] / up[ij] * exp_AUpenalty (i,j) * 
                // first, ((---)[...].-.) or ((---).[...].-.)
                ( u1_ip_jp[ip1hm1] + exp_dangle3(i,j,i+1) * EXPC[1] * (u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1]));
        }
    }
}


int s_partition_function_complex::correct_gradient_nan ()
// return 1 if no derivative is nan
{
    int i;
    int correct = 1;
    for (i=0; i < num_params; i++)
    {
        if (isnan(GlogZ[i].real()))
        {
            printf ("Glog[%d]=%g\n", i, GlogZ[i].real());
            correct = 0;
        }
    }
    return correct;
}



void s_partition_function_complex::print_gradient ()
// print the gradient
{
    int i;
    for (i=0; i < num_params; i++)
    {
        printf ("Glog[%d]=%g\n", i, GlogZ[i].real());
    }
}



void s_partition_function_complex::compute_logZ_gradient ()
// u, up and p arrays are filled
{

    for (int h=0; h < seqlen; h++)
    {
        for (int l=seqlen-1; l > h+TURN; l--)
        {
            IFD
            {
                compute_pm2 (h, l);
            }
            else
             {
                // the next 4 are only needed to compute the partial derivative wrt multi_free_base_penalty
                compute_pm2d5_needmidd5 (h, l);
                compute_pm2d5_noneedmidd5 (h, l);
                compute_pm2nod5_needmidd5 (h, l);
                compute_pm2nod5_noneedmidd5 (h, l);
            }
        }
    }



    int index_param;
    int i, j, k, l, m, n, o, r, h;
    int ii, jj, iip, jjp, iijj, iipjjp, iip1jjm1;
    Complex en_hairpin, en_stack, en_internal;
    Complex AUpen, dangle, dangle5, dangle3;
    char type[100];
    int tindex;    // type index
    int il_i, il_j, il_ip1, il_jm1;
    int index_should_be;
    
    num_params = create_string_params();

    int AUpen_index = structure_type_index ("misc.terminal_AU_penalty");
   
    for (i=0; i < num_params; i++)    GlogZ[i] = 0;
    
    // stacking energies
    index_param = 0;
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (stack[i][j][k][l] < INF)
                    {
                        // exclude duplicates
                        // stack[i][j][k][l] is the same as stack[l][k][j][i]
                        if (i*1000 + j*100 + k*10 + l <= l*1000 + k*100 + j*10 + i)
                        {
                            // traverse the whole sequence and look for this building block
                            // first look from the left
                            for (ii = 0; ii < seqlen; ii++)
                            {
                                for (jj = ii+TURN+3; jj < seqlen; jj++)
                                {
                                    if (sequence[ii] == i && sequence[jj] == j)
                                    {
                                        iijj = index[ii] + jj -ii;
                                        // stacking energies are also involved in bulges of size 1
                                        for (iip = ii + 1; iip <= ii+2; iip++)
                                        {
                                            for (jjp = jj-1; jjp >= jj-2; jjp--)
                                            {
                                                if (iip == ii+2 && jjp == jj-2)     // this is not stacking energy, is internal loop 1x1
                                                    continue;                                                
                                                if (sequence[iip] == k && sequence[jjp] == l) 
                                                {                                                        
                                                    iipjjp = index[iip]+jjp-iip;
                                                    if (iip == ii+1 && jjp == jj-1)    // stack pair
                                                        en_stack = get_stacked_energy (ii, jj, sequence);
                                                    else    // bulge
                                                        en_stack = get_internal_energy (ii, jj, iip, jjp,     sequence);
                                                    GlogZ[index_param] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_stack * oneoverRT);
                                                }
                                            }                                    
                                        }
                                    }
                                }
                            }
                            // check if it's symmetric. If so, it was already considered when we looked from the left
                            if (i!=l || j!=k)    // then it's not symmetric
                            {
                                // second, look from the right
                                for (ii = seqlen-1; ii > TURN; ii--)
                                {
                                    for (jj = ii-TURN-1; jj > 0; jj--)
                                    {
                                        if (sequence[ii] == i && sequence[jj] == j)
                                        {
                                            iijj = index[jj] + ii -jj;
                                            
                                            // stacking energies are also involved in bulges of size 1
                                            for (iip = ii + 1; iip <= MIN(ii+2,seqlen-1); iip++)
                                            {
                                                for (jjp = jj-1; jjp >= MAX(jj-2,0); jjp--)
                                                {
                                                    if (iip == ii+2 && jjp == jj-2)     // this is not stacking energy, is internal loop 1x1
                                                        continue;                                            
                                                     
                                                    if(sequence[iip] == k && sequence[jjp] == l)       
                                                    {                                                        
                                                        iipjjp = index[jjp]+iip-jjp;
                                                        if (iip == ii+1 && jjp == jj-1)    // stack pair
                                                            en_stack = get_stacked_energy (jjp, iip, sequence);
                                                        else    // bulge
                                                            en_stack = get_internal_energy (jjp, iip, jj, ii,     sequence);
                                                        GlogZ[index_param] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_stack * oneoverRT);
                                                    }
                                                }
                                            
                                            }       
                                        }
                                    }
                                }                                
                            }    
                            index_param++;                                                        
                        }                        
                    }                    
                }
                 
    index_should_be = structure_type_index("tstackh[0][3][0][0]");                
    if (index_param != index_should_be)
    {
        printf ("Index param after stack = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }                
    
    // tstackh and hairpin energies
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (tstackh[i][j][k][l] < INF)
                    {
                        // no duplicates here
                        for (ii = 0; ii < seqlen; ii++)
                        {                            
                            for (jj = ii+TURN+1; jj < seqlen; jj++)
                            {
                                if (sequence[ii] == i && sequence[jj] == j && 
                                    sequence[ii+1] == k && sequence[jj-1] == l)
                                {
                                
                                    // compute the derivative for hairpin penalty by size
                                    char s[100];                                    
                                    sprintf (s, "hairpin_penalty_by_size[%d]", jj-ii-1);
                                    // this is very slow, it's linear in the number of parameters
                                    int sizeindex = structure_type_index (s);
                                    iijj = index[ii] + jj -ii;
                                    en_hairpin = get_hairpin_energy (ii, jj, sequence, csequence);
                                    Complex grad = p[iijj] / up[iijj] * exp (en_hairpin * oneoverRT);
                                    GlogZ[sizeindex] += grad;  
                                    
                                    // add the AU_penalty
                                    
                                    if (jj-ii-1 == 3)
                                    {
                                        if (has_AU_penalty (sequence[ii], sequence[jj]))
                                        {
                                            GlogZ[AUpen_index] += grad;
                                        }                                        
                                    }
                                    
                                    // check to see if it's a special triloop or tetraloop
                                        // check if it is a triloop
                                    #if (MODEL == SIMPLE)
                                    if (jj-ii-1 == 3)
                                    {
                                        char seq[10];
                                        substr (csequence, ii, jj, seq);
                                        for (int kk=0; kk < nb_triloops; kk++)
                                        {
                                            if (strcmp (seq, triloop[kk].seq) == 0)
                                            {
                                                sprintf (s, "triloop[%d].energy", kk);
                                                sizeindex = structure_type_index (s);
                                                GlogZ[sizeindex] += grad;  
                                            }
                                        }
                                    }
                                    
                                    // check to see it is a tetraloop in tloop
                                    else if (jj-ii-1 == 4)
                                    {
                                        char seq[10];
                                        substr (csequence, ii, jj, seq);                                        
                                        for (int kk=0; kk < nb_tloops; kk++)
                                        {
                                            if (strcmp (seq, tloop[kk].seq) == 0)
                                            {
                                                sprintf (s, "tloop[%d].energy", kk);
                                                sizeindex = structure_type_index (s);
                                                GlogZ[sizeindex] += grad;  
                                            }
                                        }
                                    }
                                    #endif

                                    // compute the derivatives for misc.hairpin_GGG
                                    if (ii > 1)
                                    {
                                        if (sequence[ii]==G && sequence[ii-1]==G && sequence[ii-2]==G and sequence[jj]==U)
                                        {
                                            sizeindex = structure_type_index ("misc.hairpin_GGG");
                                            GlogZ[sizeindex] += grad;  
                                        }
                                    }
                                    
                                    // check for the special case of "poly-C" hairpin loop
                                    int is_poly_C = 1;
                                    for (int kk=ii+1; kk<jj; kk++)
                                    {
                                        if (sequence[kk] != C)
                                        {
                                            is_poly_C = 0;
                                            break;
                                        }
                                    }
                                    if (is_poly_C)
                                    {
                                        if (jj-ii-1 == 3)
                                        {
                                            sizeindex = structure_type_index ("misc.hairpin_c3");
                                            GlogZ[sizeindex] += grad;  
                                        }
                                        else
                                        {
                                            sizeindex = structure_type_index ("misc.hairpin_c1");
                                            // we have to multiply by size because this param is multiplied by size in the hairpin loop calculation
                                            GlogZ[sizeindex] += ((Complex)(jj-ii-1)) * grad;
                                            sizeindex = structure_type_index ("misc.hairpin_c2");
                                            GlogZ[sizeindex] += grad;
                                        }
                                    }
                                    
                                    // compute the derivatives of tstackh
                                    // tstackh is not added to hairpin loops of size 3
                                    if (jj-ii-1 > 3)
                                        GlogZ[index_param] += grad;
                                }
                            }
                        }
                        index_param++;
                    }
                }

    index_should_be = structure_type_index("misc.internal_AU_closure");                
    if (index_param != index_should_be)
    {
        printf ("Index param after tstackh = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }                
                                
    // tstacki energies                
    // in fact, we only have 3 parameters                            
    
    char s[100];
    //strcpy (s, "misc.terminal_AU_penalty");
    // this is very slow, it's linear in the number of parameters
    //int index_terminal_AU_penalty = structure_type_index (s);                        
    strcpy (s, "misc.internal_AU_closure");
    int index_internal_AU_closure = structure_type_index (s);                        
    strcpy (s, "misc.internal_AG_mismatch");
    int index_internal_AG_mismatch = structure_type_index (s);
    strcpy (s, "misc.internal_UU_mismatch");
    int index_internal_UU_mismatch = structure_type_index (s);    
                                
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
                for (l=0; l < NUCL; l++)
                {
                    if (tstacki[i][j][k][l] < INF)
                    {
                        // no duplicates here
                        for (ii = 0; ii < seqlen; ii++)
                        {
                            for (jj = ii+TURN+4; jj < seqlen; jj++)
                            {
                                // first check the weird case of k==0 and l==0
                                if (sequence[ii] == i && sequence[jj] == j && k==0 && l==0)
                                {
                                    iijj = index[ii] + jj -ii;
                                    
                                    for (iip = ii+1; iip <= MIN(jj-2-TURN,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, iip+1+TURN);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                // compute the derivative for internal penalty by size
                                                    
                                                iipjjp = index[iip] + jjp - iip;  
                                                if (((branch1 == 1 && branch2 > 2) || (branch1 > 2 && branch2 == 1)) && misc.gail_rule)
                                                {
                                                    //char s[100];
                                                    sprintf (s, "internal_penalty_by_size[%d]", branch1+branch2);
                                                    // this is very slow, it's linear in the number of parameters
                                                    int sizeindex = structure_type_index (s);
                                                    //printf ("sizeindex=%d\n", sizeindex);
                                                
                                                    en_internal = get_internal_energy (ii, jj, iip, jjp, sequence);
                                                    Complex grad = p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                    GlogZ[sizeindex] += grad; 
                                                    
                                                    // now check if any of the 3 internal params is involved
                                                    if (((i == A || i == G) && j == U) ||
                                                        ((j == A || j == G) && i == U))
                                                    {
                                                        //GlogZ[index_terminal_AU_penalty] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                        GlogZ[index_internal_AU_closure] += grad;
                                                    }
                                                    // AG and UU mismatch can't happen here, because k is A and l is A
                                                }
                                            }
                                        }
                                    }                                
                                }

                                // left tstacki, the right one to be added
                                if (sequence[ii] == i && sequence[jj] == j && 
                                    sequence[ii+1] == k && sequence[jj-1] == l)
                                {
                                    iijj = index[ii] + jj -ii;
                                    
                                    for (iip = ii+1; iip <= MIN(jj-2-TURN,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, iip+1+TURN);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                if (branch1 < 0 || branch2 < 0)    // not valid
                                                    continue;
                                                if (branch1 == 0 && branch2 == 0)    // stack pair, not good here
                                                    continue;
                                                if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)
                                                    continue;                                                    
                                                if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 1 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
                                                    continue;                                                

                                                if (branch1 == 0 || branch2 == 0)    // bulge
                                                {
                                                    // compute the derivative for bulge penalty by size
                                                    //char s[100];
                                                    sprintf (s, "bulge_penalty_by_size[%d]", branch1+branch2);
                                                    // this is very slow, it's linear in the number of parameters
                                                    int sizeindex = structure_type_index (s);
                                                    en_internal = get_internal_energy (ii, jj, iip, jjp, sequence);
                                                    iipjjp = index[iip] + jjp - iip;
                                                    Complex grad = p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                    GlogZ[sizeindex] += grad; 
                                                    
                                                    if (branch1 + branch2 > 1)    // add AU_penalty
                                                    {                                                        
                                                        if (has_AU_penalty (sequence[ii], sequence[jj]))
                                                        {
                                                            GlogZ[AUpen_index] += grad;
                                                        }
                                                        if (has_AU_penalty (sequence[iip], sequence[jjp]))
                                                        {
                                                            GlogZ[AUpen_index] += grad;
                                                        }
                                                    }
                                                    continue;
                                                }
                                                    
                                                // compute the derivative for internal penalty by size
                                                //char s[100];
                                                sprintf (s, "internal_penalty_by_size[%d]", branch1+branch2);
                                                // this is very slow, it's linear in the number of parameters
                                                int sizeindex = structure_type_index (s);
                                                //printf ("sizeindex=%d\n", sizeindex);
                                                    
                                                iipjjp = index[iip] + jjp - iip;  
                                                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                                                    continue;
                                                    
                                                en_internal = get_internal_energy (ii, jj, iip, jjp, sequence);
                                                Complex grad = p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                GlogZ[sizeindex] += grad; 
                                                
                                                // now check if any of the 3 internal params is involved
                                                if (((i == A || i == G) && j == U) ||
                                                    ((j == A || j == G) && i == U))
                                                {
                                                    //GlogZ[index_terminal_AU_penalty] += p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT); 
                                                    GlogZ[index_internal_AU_closure] += grad;
                                                }
                                                if ((k == A && l == G) ||
                                                    (l == A && k == G))
                                                {
                                                    GlogZ[index_internal_AG_mismatch] += grad; 
                                                }
                                                if (k == U && l == U)
                                                {
                                                    GlogZ[index_internal_UU_mismatch] += grad; 
                                                }  
                                            }
                                        }
                                    }                                    
                                }
                            }
                        }           
                        // the upside down orientation
                        // no duplicates here
                        for (ii = seqlen-2; ii > TURN; ii--)
                        {
                            for (jj = ii-TURN-1; jj > 0; jj--)
                            {
                                // first check the weird case of k==0 si l==0
                                if (sequence[ii] == i && sequence[jj] == j && k==0 && l==0)
                                {
                                    // jj < ii
                                    iijj = index[jj] + ii -jj;
                                    
                                    for (iip = ii+1; iip <= MIN(seqlen-1,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, 0);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                // compute the derivative for internal penalty by size
                                                    
                                                iipjjp = index[jjp] + iip - jjp; 
                                                if (((branch1 == 1 && branch2 > 2) || (branch1 > 2 && branch2 == 1)) && misc.gail_rule)
                                                {                                                
                                                    en_internal = get_internal_energy (jjp, iip, jj, ii, sequence);
                                                    // now check if any of the 3 internal params is involved
                                                    if (((i == A || i == G) && j == U) ||
                                                        ((j == A || j == G) && i == U))
                                                    {
                                                        //GlogZ[index_terminal_AU_penalty] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                        GlogZ[index_internal_AU_closure] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                    }
                                                    // AG and UU mismatch can't happen here, because k is A and l is A                                               
                                                }
                                            }
                                        }
                                    }                                
                                }                            
                            
                                // the general case
                                // left tstacki, the right one to be added
                                if (sequence[ii] == i && sequence[jj] == j && 
                                    sequence[ii+1] == k && sequence[jj-1] == l)
                                {
                                    // jj < ii
                                    iijj = index[jj] + ii -jj;
                                    for (iip = ii+1; iip <= MIN(seqlen-1,ii+MAXLOOP+1) ; iip++)  // j-2-TURN
                                    {
                                        int minq = MAX (jj-ii+iip-MAXLOOP-2, 0);    // ip+1+TURN);
                                        for (jjp = minq; jjp < jj; jjp++)
                                        {        
                                            if (sequence[iip]+sequence[jjp] == 3 ||
                                                sequence[iip]+sequence[jjp] == 5)        
                                            {
                                                // same as before
                                                int branch1 = iip-ii-1;
                                                int branch2 = jj-jjp-1;
                                                if (branch1 < 0 || branch2 < 0)    // not valid
                                                    continue;
                                                if (branch1 == 0 && branch2 == 0)    // stack pair, not good here
                                                    continue;
                                                if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 1 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
                                                    continue;
                                                if (branch1 == 0 || branch2 == 0)    // bulge
                                                    continue;
                                                    
                                                //jjp < iip
                                                iipjjp = index[jjp] + iip - jjp;  
                                                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                                                    continue;
                                                en_internal = get_internal_energy (jjp, iip, jj, ii, sequence);
                                                Complex grad = p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                
                                                // now check if any of the 3 internal params is involved
                                                if (((i == A || i == G) && j == U) ||
                                                    ((j == A || j == G) && i == U))
                                                {
                                                    //GlogZ[index_terminal_AU_penalty] += p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT); 
                                                    GlogZ[index_internal_AU_closure] += grad; 
                                                }
                                                if ((k == A && l == G) ||
                                                    (l == A && k == G))
                                                {
                                                    GlogZ[index_internal_AG_mismatch] += grad; 
                                                }
                                                if (k == U && l == U)
                                                {
                                                    GlogZ[index_internal_UU_mismatch] += grad; 
                                                }                                                 
                                            }
                                        }
                                    }                                    
                                }
                            }
                        }
                        //index_param++;   // there are only 3 tstacki parameters                              
                    }    // end if (tstacki[i][j][k][l] < INF)                                            
                }
                
    index_param += 3;           
    index_should_be = structure_type_index("int11[0][3][3][3][0][3]");     
    if (index_param != index_should_be)
    {
        printf ("Index param after tstacki = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }                
                                      
    // internal loops 1x1                        
    if (!simple_internal_energy)          
    {
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                            {                            
                                if (int11[i][j][k][l][m][n] < INF)
                                {
                                    // exclude duplicates
                                    // int11[i][j][k][l][m][n] is the same as int11[n][m][l][k][j][i]
                                    if (i*100000 + j*10000 + k*1000 + l*100 + m*10 + n <= n*100000 + m*10000 + l*1000+ k*100 + j*10 + i)
                                    {             
                                       
                                        // look from the left
                                        for (ii = 0; ii < seqlen-TURN-5; ii++)
                                        {
                                            if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                            {
                                                for (jj = ii + TURN+5; jj < seqlen; jj++)
                                                {
                                                    if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == n)
                                                    {                      
                                                        iijj = index[ii] + jj -ii;    
                                                        iipjjp = index[ii+2] + jj-2 - (ii+2); 
                                                        en_internal = get_internal_energy (ii, jj, ii+2, jj-2, sequence);
                                                        Complex grad = p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                              
                                                        int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc;
                                                        iinuc = sequence[ii]; 
                                                        jjnuc = sequence[jj];
                                                        kknuc = sequence[ii+1];
                                                        llnuc = sequence[jj-1];
                                                        mmnuc = sequence[ii+2];
                                                        nnnuc = sequence[jj-2];
                                        
                                                        if ( ((iinuc==C && jjnuc==G) || (iinuc==G && jjnuc==C)) && ((mmnuc==C && nnnuc==G) || (mmnuc==G && nnnuc==C))) 
                                                        {
                                                            if (!can_pair(kknuc,llnuc))
                                                                sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                            else
                                                                sprintf (type, "misc.internal11_basic_mismatch");
                                                            tindex = structure_type_index (type);                                                            
                                                            GlogZ[tindex] += grad;
                                                        }        
                                                        else if (watson_crick(iinuc,jjnuc) && watson_crick(mmnuc,nnnuc) && kknuc==U && llnuc==U)
                                                        {
                                                            sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += grad;
                                                        }
                                                        else
                                                        {
                                                            if (kknuc==G && llnuc==G)
                                                                sprintf (type, "misc.internal11_GG_mismatch");
                                                            else
                                                                sprintf (type, "misc.internal11_basic_mismatch");
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += grad;
                                                                                        
                                                            if (has_AU_penalty(iinuc,jjnuc))
                                                            {
                                                                sprintf (type, "misc.internal_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;         
                                                            }
                                                            if (has_AU_penalty(mmnuc,nnnuc))
                                                            {
                                                                sprintf (type, "misc.internal_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }
                                                        }                                                                                                       
                                                    }    
                                                }
                                            }                   
                                        }
                        
                                        // look from the right, if it's not symmetric
                                        if (!(i==n && k==l && m==j))
                                        {
                                            for (ii = seqlen-3; ii >= TURN+3; ii--)
                                            {
                                                if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                                {
                                                    for (jj = ii-TURN-1; jj >= 2; jj--)
                                                    {
                                                        if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == n)
                                                        {
                                                            //jj < ii
                                                            iijj = index[jj] + ii -jj;    
                                                            iipjjp = index[jj-2] + ii+2 - (jj-2); 
                                                            en_internal = get_internal_energy (jj-2, ii+2, jj, ii, sequence);
                                                            Complex grad = p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                            
                                                            int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc;
                                                            iinuc = sequence[ii]; 
                                                            jjnuc = sequence[jj];
                                                            kknuc = sequence[ii+1];
                                                            llnuc = sequence[jj-1];
                                                            mmnuc = sequence[ii+2];
                                                            nnnuc = sequence[jj-2];
                                            
                                                            if ( ((iinuc==C && jjnuc==G) || (iinuc==G && jjnuc==C)) && ((mmnuc==C && nnnuc==G) || (mmnuc==G && nnnuc==C))) 
                                                            {
                                                                if (!can_pair(kknuc,llnuc))
                                                                    sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                                else
                                                                    sprintf (type, "misc.internal11_basic_mismatch");
                                                                tindex = structure_type_index (type);                                                            
                                                                GlogZ[tindex] += grad;
                                                            }        
                                                            else if (watson_crick(iinuc,jjnuc) && watson_crick(mmnuc,nnnuc) && kknuc==U && llnuc==U)
                                                            {
                                                                sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }
                                                            else
                                                            {
                                                                if (kknuc==G && llnuc==G)
                                                                    sprintf (type, "misc.internal11_GG_mismatch");
                                                                else
                                                                    sprintf (type, "misc.internal11_basic_mismatch");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                                                            
                                                                if (has_AU_penalty(iinuc,jjnuc))
                                                                {
                                                                    sprintf (type, "misc.internal_AU_closure");
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += grad;  
                                                                }
                                                                if (has_AU_penalty(mmnuc,nnnuc))
                                                                {
                                                                    sprintf (type, "misc.internal_AU_closure");
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += grad;
                                                                }
                                                            }                                                                                                       
                                                                
                                                        }    
                                                    }
                                                }                   
                                            }
                                        }                         
                                    }                                            
                                } // end if (int11[i][j][k][l][m][n] < INF)
                            }    // end int11
        index_param += 33;
        index_should_be = structure_type_index("int21[1][2][0][0][1][2][0]");
        if (index_param != index_should_be)
        {
            printf ("Index param after int11 = %d, should be %d\n", index_param, index_should_be);
            exit(1);
        }                 
        
        // int12 energies
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for (o=0; o < NUCL; o++)
                                {
                                    if (int21[i][j][k][l][m][n][o] < INF)
                                    {
                                        // no duplicates here
                                        // look from the left
                                        for (ii = 0; ii < seqlen-TURN-6; ii++)
                                        {
                                            if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                            {
                                                for (jj = ii + TURN+6; jj < seqlen; jj++)
                                                {
                                                    if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == o && sequence[jj-3] == n)
                                                    {
                                                        iijj = index[ii] + jj -ii;    
                                                        iipjjp = index[ii+2] + jj-3 - (ii+2); 
                                                        en_internal = get_internal_energy (ii, jj, ii+2, jj-3, sequence);
                                                        Complex grad = p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                        
                                                        int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc;
                                                        iinuc = sequence[ii];
                                                        jjnuc = sequence[jj];
                                                        kknuc = sequence[ii+1];
                                                        llnuc = sequence[jj-1];
                                                        mmnuc = sequence[ii+2];
                                                        nnnuc = sequence[jj-3];
                                                        oonuc = sequence[jj-2];            
                                            
                                                        if ((iinuc==C && jjnuc==G && mmnuc==C && nnnuc==G) ||  // these are already filled above, except what can pair inside
                                                            (iinuc==G && jjnuc==C && mmnuc==G && nnnuc==C))
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                                sprintf (type, "misc.internal21_match");
                                                            else
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc);
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += grad;
                                                        }
                                                        else
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                            {
                                                                sprintf (type, "misc.internal21_match");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }
                                                            else
                                                            {
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kknuc, llnuc, C, G, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                                
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kknuc, llnuc, G, C, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                            }    
                                                            if (has_AU_penalty(iinuc,jjnuc))
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }    
                                                            if (has_AU_penalty(mmnuc,nnnuc))    
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }    
                                                        }                                                                                                                                                    
                                                    }    
                                                }
                                            }                   
                                        }
                        
                                        // int21 cannot be symmetric
                                        for (ii = seqlen-3; ii >= TURN+4; ii--)
                                        {
                                            if (sequence[ii] == i && sequence[ii+1]==k && sequence[ii+2] == m)
                                            {
                                                for (jj = ii-TURN-1; jj >= 3; jj--)
                                                {
                                                    if (sequence[jj] == j && sequence[jj-1] == l && sequence[jj-2] == o && sequence[jj-3] == n)
                                                    {
                                                        //jj < ii
                                                        iijj = index[jj] + ii -jj;    
                                                        iipjjp = index[jj-3] + ii+2 - (jj-3); 
                                                        en_internal = get_internal_energy (jj-3, ii+2, jj, ii, sequence);
                                                        Complex grad = p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                        
                                                        int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc;
                                                        iinuc = sequence[ii];
                                                        jjnuc = sequence[jj];
                                                        kknuc = sequence[ii+1];
                                                        llnuc = sequence[jj-1];
                                                        mmnuc = sequence[ii+2];
                                                        nnnuc = sequence[jj-3];
                                                        oonuc = sequence[jj-2];            
                                            
                                                        if ((iinuc==C && jjnuc==G && mmnuc==C && nnnuc==G) ||  // these are already filled above, except what can pair inside
                                                            (iinuc==G && jjnuc==C && mmnuc==G && nnnuc==C))
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                                sprintf (type, "misc.internal21_match");
                                                            else
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc);
                                                            tindex = structure_type_index (type);
                                                            GlogZ[tindex] += grad;
                                                        }
                                                        else
                                                        {
                                                            if (can_pair(kknuc,llnuc) || can_pair(kknuc,oonuc))
                                                            {
                                                                sprintf (type, "misc.internal21_match");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }
                                                            else
                                                            {
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kknuc, llnuc, C, G, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                                
                                                                sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kknuc, llnuc, G, C, oonuc);
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                            }
                                                            if (has_AU_penalty(iinuc,jjnuc))
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }    
                                                            if (has_AU_penalty(mmnuc,nnnuc))    
                                                            {
                                                                sprintf (type, "misc.internal21_AU_closure");
                                                                tindex = structure_type_index (type);
                                                                GlogZ[tindex] += grad;
                                                            }    
                                                        }                
                                                    }    
                                                }
                                            }                   
                                        }
                                    //index_param++;                            
                                    }                                            
                                } // end int21
                                
        index_param += 54;                                
        index_should_be = structure_type_index("int22[0][3][0][0][3][0][0][0]");                                
        if (index_param != index_should_be)
        {
            printf ("Index param after int21 = %d, should be %d\n", index_param, index_should_be);
            exit(1);
        }                                     
        // energies int22                                
        for (i=0; i < NUCL; i++)
            for (j=0; j < NUCL; j++)
                for (k=0; k < NUCL; k++)
                    for (l=0; l < NUCL; l++)
                        for (m=0; m < NUCL; m++)
                            for (n=0; n < NUCL; n++)
                                for(o=0; o < NUCL; o++)
                                    for (r=0; r < NUCL; r++)
                                    {
                                        if (int22[i][j][k][l][m][n][o][r] < INF)
                                        {
                                            // exclude duplicates
                                            // int22[i][j][k][l][m][n][o][p] is the same as int22[n][m][p][o][j][i][l][k]
                                            if (i*10000000 + j*1000000 + k*100000 + l*10000 + m*1000 + n*100 + o*10 + r <= n*10000000 + m*1000000 + r*100000 + o*10000 + j*1000 + i*100 + l*10 + k)
                                            {
                                                // look from the left
                                                for (ii = 0; ii < seqlen-TURN-7; ii++)
                                                {
                                                    if (sequence[ii]==i && sequence[ii+1]==k && sequence[ii+2]==o && sequence[ii+3]==m)
                                                    {
                                                        for (jj = ii + TURN+7; jj < seqlen; jj++)
                                                        {
                                                            if (sequence[jj]==j && sequence[jj-1]==l && sequence[jj-2]==r && sequence[jj-3]==n)
                                                            {
                                                                iijj = index[ii] + jj -ii;    
                                                                iipjjp = index[ii+3] + jj-3 - (ii+3); 
                                                                en_internal = get_internal_energy (ii, jj, ii+3, jj-3, sequence);
                                                                Complex grad = p[iijj] * up[iipjjp] / up[iijj] * exp (en_internal * oneoverRT);
                                                                
                                                                int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc;
                                                                iinuc = sequence[ii];
                                                                jjnuc = sequence[jj];
                                                                kknuc = sequence[ii+1];
                                                                llnuc = sequence[jj-1];
                                                                mmnuc = sequence[ii+3];
                                                                nnnuc = sequence[jj-3];
                                                                oonuc = sequence[ii+2];            
                                                                ppnuc = sequence[jj-2];
                                                                
                                                                if (nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc & watson_crick(iinuc,jjnuc) && !watson_crick(kknuc,llnuc))
                                                                {
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += grad;
                                                                }
                                                            
                                                                int iinuc2, jjnuc2, mmnuc2, nnnuc2;                
                                                                if (iinuc==G && jjnuc==U)   iinuc2 = A;     else iinuc2 = iinuc;
                                                                if (iinuc==U && jjnuc==G)   jjnuc2 = A;     else jjnuc2 = jjnuc;
                                                                if (mmnuc==G && nnnuc==U)   mmnuc2 = A;     else mmnuc2 = mmnuc;
                                                                if (mmnuc==U && nnnuc==G)   nnnuc2 = A;     else nnnuc2 = nnnuc;
                                                                
                                                                if (watson_crick(kknuc,llnuc) || watson_crick(oonuc,ppnuc))
                                                                {
                                                                    sprintf (type, "misc.internal22_match");
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += grad;
                                                                }
                                                                else if ( ((iinuc==G && jjnuc==U) || (iinuc==U && jjnuc==G) || (mmnuc==G && nnnuc==U) || (mmnuc==U && nnnuc==G)) &&  
                                                                        (nnnuc2==iinuc2 && mmnuc2==jjnuc2 && ppnuc==kknuc && oonuc==llnuc))  // the UG closing pairs are the same as UA
                                                                {
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += grad;
                                                                }                
                                                                else if (!(nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc))   // was already filled above
                                                                {
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, jjnuc2, iinuc2, llnuc, kknuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                                    sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nnnuc2, mmnuc2, ppnuc, oonuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                                    
                                                                    int result = check_stability_and_size (kknuc, llnuc, oonuc, ppnuc);                    
                                                                    switch (result)
                                                                    {
                                                                        case 1: 
                                                                            sprintf (type, "misc.internal22_delta_same_size");
                                                                            break;
                                                                        case 2: 
                                                                            sprintf (type, "misc.internal22_delta_different_size");
                                                                            break;
                                                                        case 3: 
                                                                            sprintf (type, "misc.internal22_delta_1stable_1unstable");
                                                                            break;
                                                                        case 4: 
                                                                            sprintf (type, "misc.internal22_delta_AC");
                                                                            break;
                                                                        default: 
                                                                            printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kknuc, llnuc, oonuc, ppnuc);
                                                                            exit(1);
                                                                    }
                                                                    tindex = structure_type_index (type);
                                                                    GlogZ[tindex] += grad;
                                                                }                                                                
                                                            }    
                                                        }
                                                    }                   
                                                }
                                
                                                // look from the right, if it's not symmetric
                                                if (!(i==n && k==r && o==l && m==j))
                                                {
                                                    for (ii = seqlen-4; ii >= TURN+4; ii--)
                                                    {
                                                        if (sequence[ii]==i && sequence[ii+1]==k && sequence[ii+2]==o && sequence[ii+3]==m)
                                                        {
                                                            for (jj = ii-TURN-1; jj >= 3; jj--)
                                                            {
                                                                if (sequence[jj]==j && sequence[jj-1]==l && sequence[jj-2]==r && sequence[jj-3]==n)
                                                                {
                                                                    //jj < ii
                                                                    iijj = index[jj] + ii -jj;    
                                                                    iipjjp = index[jj-3] + ii+3 - (jj-3); 
                                                                    en_internal = get_internal_energy (jj-3, ii+3, jj, ii, sequence);
                                                                    Complex grad = p[iipjjp] * up[iijj] / up[iipjjp] * exp (en_internal * oneoverRT);
                                                                    
                                                                    int iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc;
                                                                    iinuc = sequence[ii];
                                                                    jjnuc = sequence[jj];
                                                                    kknuc = sequence[ii+1];
                                                                    llnuc = sequence[jj-1];
                                                                    mmnuc = sequence[ii+3];
                                                                    nnnuc = sequence[jj-3];
                                                                    oonuc = sequence[ii+2];            
                                                                    ppnuc = sequence[jj-2];
                                                                    
                                                                    if (nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc & watson_crick(iinuc,jjnuc) && !watson_crick(kknuc,llnuc))
                                                                    {
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc, jjnuc, kknuc, llnuc, mmnuc, nnnuc, oonuc, ppnuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += grad;
                                                                    }
                                                                
                                                                    int iinuc2, jjnuc2, mmnuc2, nnnuc2;                
                                                                    if (iinuc==G && jjnuc==U)   iinuc2 = A;     else iinuc2 = iinuc;
                                                                    if (iinuc==U && jjnuc==G)   jjnuc2 = A;     else jjnuc2 = jjnuc;
                                                                    if (mmnuc==G && nnnuc==U)   mmnuc2 = A;     else mmnuc2 = mmnuc;
                                                                    if (mmnuc==U && nnnuc==G)   nnnuc2 = A;     else nnnuc2 = nnnuc;
                                                                    
                                                                    if (watson_crick(kknuc,llnuc) || watson_crick(oonuc,ppnuc))
                                                                    {
                                                                        sprintf (type, "misc.internal22_match");
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += grad;
                                                                    }
                                                                    else if ( ((iinuc==G && jjnuc==U) || (iinuc==U && jjnuc==G) || (mmnuc==G && nnnuc==U) || (mmnuc==U && nnnuc==G)) &&  
                                                                            (nnnuc2==iinuc2 && mmnuc2==jjnuc2 && ppnuc==kknuc && oonuc==llnuc))  // the UG closing pairs are the same as UA
                                                                    {
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += grad;
                                                                    }                
                                                                    else if (!(nnnuc==iinuc && mmnuc==jjnuc && ppnuc==kknuc && oonuc==llnuc))   // was already filled above
                                                                    {
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", iinuc2, jjnuc2, kknuc, llnuc, jjnuc2, iinuc2, llnuc, kknuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                                        sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nnnuc2, mmnuc2, ppnuc, oonuc, mmnuc2, nnnuc2, oonuc, ppnuc);
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += (PARAMTYPE)0.5 * grad;
                                                                        
                                                                        int result = check_stability_and_size (kknuc, llnuc, oonuc, ppnuc);                    
                                                                        switch (result)
                                                                        {
                                                                            case 1: 
                                                                                sprintf (type, "misc.internal22_delta_same_size");
                                                                                break;
                                                                            case 2: 
                                                                                sprintf (type, "misc.internal22_delta_different_size");
                                                                                break;
                                                                            case 3: 
                                                                                sprintf (type, "misc.internal22_delta_1stable_1unstable");
                                                                                break;
                                                                            case 4: 
                                                                                sprintf (type, "misc.internal22_delta_AC");
                                                                                break;
                                                                            default: 
                                                                                printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kknuc, llnuc, oonuc, ppnuc);
                                                                                exit(1);
                                                                        }
                                                                        tindex = structure_type_index (type);
                                                                        GlogZ[tindex] += grad;
                                                                    }
                                                                }    
                                                            }
                                                        }                   
                                                    }
                                                }    
                                            //index_param++;                            
                                            }                                            
                                        } // end if (int11[i][j][k][l][m][n] < INF)
                                    }    // end int11
                                        
        index_param += 53;                                    
        index_should_be = structure_type_index("dangle_top[0][3][0]");                                                
        if (index_param != index_should_be)
        {
            printf ("Index param after int22 = %d, should be %d\n", index_param, index_should_be);
            exit(1);
        }              
                                    
    }    // end int11, int12, int22,     // end if (!simple_internal_energy)     
    

    

    if (compute_gradient_dangles)    // we don't need to compute dangling ends
    {
    
    // dangling 3', or dangle_top    
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                if (dangle_top[i][j][k] < INF)
                {
                    // no duplicates here                    
                    // traverse the whole sequence and look for this building block                    

                    // it never gets here but anyway
                    IFD
                    {
                        GlogZ[index_param] = 0.0;
                        index_param++;
                        continue;
                    }
                    
                    // first add the dangling energies that participate to the external loop                     
                    for (ii = 0; ii < seqlen; ii++)
                    {
                        for (jj = ii+TURN+1; jj < seqlen-1; jj++)
                        {
                            if (sequence[jj] == i && sequence[ii] == j && sequence[jj+1] == k)
                            {            
                                Complex term1 = 0.0;
                                Complex term2 = 0.0;
                                                            
                                iijj = index[ii]+jj-ii;
                                
                                if (ii > TURN)        
                                    term1 += u_ip_jp[ii-1] + u_iu_jp[ii-1] + exp_dangle5(jj,ii,ii-1)*(u_ip_ju[ii-1] + u_iu_ju[ii-1]);                                  
                                else if (ii > 0)      term1 += exp_dangle5(jj,ii,ii-1);
                                else                  term1 = 1;  
    
                                if (jj < seqlen-3)
                                {
                                    int jjp2n = index[jj+2]+seqlen-1-(jj+2);
                                    term2 += 
                                            exp_dangle3(jj,ii,jj+1) * (u_ip_jp[jjp2n] + u_ip_ju[jjp2n] + u_iu_jp[jjp2n] + u_iu_ju[jjp2n]);
                                }
                                else if (jj < seqlen-1) term2 += exp_dangle3(jj,ii,jj+1);
                                // it will never get here
                                //else                  term2 = 1;
                                
                                GlogZ[index_param] += term1*up[iijj]*term2* exp_AUpenalty (ii,jj);
                            }
                        }
                    }
                    // moved the denominator outside, it's the same everywhere
                    GlogZ[index_param] /= (u_ip_jp[seqlen-1] + u_ip_ju[seqlen-1] + u_iu_jp[seqlen-1] + u_iu_ju[seqlen-1]) ;

                    
                    // contribution from multi-loops
                    // the first free bases
                    // taken from compute_upm  
                                      
                    Complex temp = 0;
                    for (jj=TURN+1; jj < seqlen; jj++)
                    {
                        for (ii=jj-TURN-1; ii>=0; ii--)
                        {                     
                            iijj = index[ii] + jj -ii;  
                            if (sequence[ii] == i && sequence[jj] == j && sequence[ii+1] == k)
                            {                 
                                temp = 0;
                                if (upm[iijj].real() == 0)
                                    continue;                                                              
                                for (l=ii+3; l < jj-TURN-2; l++)    // case (.(...)--(--)-)
                                {
                                    int iip2l = index[ii+2]+l-ii-2;
                                    int lp2jjm1 = index[l+2]+jj-1-l-2;
                                    int lp1jjm1 = index[l+1]+jj-1-l-1;
                                    int lp2jjm2 = index[l+2]+jj-2-l-2;
                                    int lp1jjm2 = index[l+1]+jj-2-l-1;        
                                    int lp2jjm3 = index[l+2]+jj-3-l-2;
                                    int lp1jjm3 = index[l+1]+jj-3-l-1;                
                                    temp += up[iip2l] *
                                            exp_AUpenalty (ii+2, l) *
                                            ( u1_ip_jp[lp1jjm1]      // [.(...)(.-..)] or [.(...)(.-..).]                    
                                                + exp_dangle3 (l, ii+2, l+1) * EXPC[1] *
                                                    ( u1_ip_jp[lp2jjm1] + u1_iu_jp[lp2jjm1] +     // [.(...).-(...)] or [.(...).-(...).] 
                                                       exp_dangle5(ii,jj,jj-1) * (u1_ip_ju[lp2jjm1] + u1_iu_ju[lp2jjm1])) +    // [.(...).-(...)-..]
                                                + exp_dangle5(ii,jj,jj-1) * u1_ip_ju[lp1jjm1]);    // [.(...)(...)-..]                                             
                                }
                                for (h=ii+3; h < jj-TURN-2; h++)    // case (....(...)--(--)-)
                                {
                                    int hjj = index[h]+jj-h;
                                    int hjjm1 = index[h] + jj-1 -h;
                                    int hjjm2 = index[h] + jj-2 -h;
                                    int hjjm3 = index[h] + jj-3 -h;
                                    temp += EXPC[h-ii-2] *
                                        ( s2_jp[hjjm1]    // --(...)) or --(...).) 
                                          + s2_ju[hjjm1]* exp_dangle5(ii,jj,jj-1));    // --(...)..)
                                    //s2[hj];
                                }    
                                //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);

                                GlogZ[index_param] += temp * EXPA * EXPB2 * EXPC[1] *exp_AUpenalty (ii,jj) * exp_dangle3 (ii, jj, ii+1) * p[iijj] / up[iijj];
                            }
                        }
                    }
                    

                    // the multi-loop 3' dangling end which can appear to the right of branches
                    // taken from compute_p
                    temp = 0;
                    for (h=0; h < seqlen; h++)   
                    {
                        for (l=seqlen-2; l > h+TURN; l--)    // it's seqlen-2 because l+1 must be valid
                        {
                            //printf ("h=%d, l=%d, seqlen=%d\n", h, l, seqlen);
                            if (can_pair (sequence[h], sequence[l]) && 
                                sequence[l] == i && sequence[h] == j && sequence[l+1] == k)    
                            {
                                Complex temp = 0;
                                int hl = index[h] + l - h;
                                
                                for (ii=0; ii < h; ii++)
                                {
                                    int iil = index[ii] + l - ii;                                    
                                    // the case when h-l is the first branch of the multi-loop
                                    if (ii < h-1)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
                                    {
                                        // removed one line, no dangling end of interest
                                        temp += EXPC[h-ii] * pmd3_needmidd3[iil]  *   // EXPC[1] * // added it to the left
                                                    (ii<h-2?exp_dangle5(l,h,h-1):1);    
                                            // (.-[...].-(---)-)    i.-h...l.-(---)-j                 
                                    }
                                    else    // case ((... , no dangling end
                                    {
                                        // removed one line, no dangling end of interest
                                        temp += pmnod3_needmidd3[iil] * EXPC[h-ii]; // * EXPC[1];
                                            // ([...].-(---)-)    ih...l.-(---)-j                 
                                    }
                                    
                                    if (ii < h - TURN -2)    // only now we can have a branch to the left
                                    {        
                                        int iip1hm1 = index[ii+1] + h-1 - (ii+1);
                                        int iip1hm2 = index[ii+1] + h-2 - (ii+1);
                                        int iip1hm3 = index[ii+1] + h-3 - (ii+1);
                                        int iip2hm1 = index[ii+2] + h-1 - (ii+2);
                                        int iip2hm2 = index[ii+2] + h-2 - (ii+2);
                                        int iip2hm3 = index[ii+2] + h-3 - (ii+2);
                                        int iip3hm1 = index[ii+3] + h-1 - (ii+3);
                                        int iip3hm2 = index[ii+3] + h-2 - (ii+3);
                                        int iip3hm3 = index[ii+3] + h-3 - (ii+3);
                                        int iilp1 = index[ii] + l+1 - ii;
                                        // removed term1, no dangling end of interest
                                        
                                        temp +=
                                        
                                        // first, when h-l is the last branch to the right
                                                                                            
                                        // case ((..-)-[...].-)  
                                        pm1nod3_needendd3[iil] *    // EXPC added in pm1nod3
                                            (u1_ip_jp[iip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1])
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-)
                                        + pm1d3_needendd3[iil] * //  don't need EXPC[1]
                                            (u1_ip_jp[iip2hm1] + u1_iu_jp[iip2hm1] +
                                                exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1] + u1_iu_ju[iip2hm1]))
                                                    // (.(...).-[...].-)        i.(...)h...l.-j                                                    
                                        // case (..-(..-)-[...].-)            
                                                    // (..-(...)[...].-)      i..-(...)h...l.-j
                                                    
                                        // we have branches to the left and to the right of h-l
                                        // if this is changed in compute_p, it should be changed here too
                                        //+ pmnod3_noneedmidd3[iil]  * u1[iip1hm1]     
                                                    // ((...)[...]-(---)...)    i(...)h...l-(---)...j 
                                        //+ pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * u1[iip1hm1]);        
                                                    // ((...).-[...]-(---)...)    i(...).-h...l-(---)...j 
                                                    
                                        + pmnod3_needmidd3[iil] * EXPC[1]*    //+ pmnod3_noneedmidd3[il]) *
                                            (u1_ip_jp[iip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[iip1hm1])
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-(--)-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                                        + EXPC[2]* pmd3_needmidd3[iil] * //EXPC[1]*     //+ pmd3_noneedmidd3[il]) *
                                                (u1_ip_jp[iip2hm1] + u1_iu_jp[iip2hm1] +
                                                    exp_dangle5(l,h,h-1)*(u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1]));
                                                    // (.(...).-[...].-)        i.(...)h...l.-j
                                                    
                                        // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                                    }
                                }
                                GlogZ[index_param] += up[hl] * EXPA * EXPB2 * exp_AUpenalty (h,l) *
                                        exp_dangle3(l,h,l+1) * temp;                                
                            }
                        }
                    }
                    
                    index_param++;
                }
            }
    index_should_be = structure_type_index("dangle_bot[0][3][0]");            
    if (index_param != index_should_be)
    {
        printf ("Index param after dangle_top = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }      

    // dangling 5', or dangle_bot
    for (i=0; i < NUCL; i++)
        for (j=0; j < NUCL; j++)
            for (k=0; k < NUCL; k++)
            {
                if (dangle_bot[i][j][k] < INF)
                {
                    // it never gets here
                    IFD
                    {
                        GlogZ[index_param] = 0.0;
                        index_param++;
                        continue;
                    }
                    // no duplicates here
                    // first add the dangling energies that participate to the external loop
                    // traverse the whole sequence and look for this building block
                    
                    
                    for (ii = 1; ii < seqlen; ii++)
                    {
                        for (jj = ii+TURN+1; jj < seqlen; jj++)
                        {
                            if (sequence[jj] == i && sequence[ii] == j && sequence[ii-1] == k)
                            {                            
                                Complex term1 = 0.0;
                                Complex term2 = 0.0;
                                iijj = index[ii]+jj-ii;                                
                                if (ii > TURN)        
                                        // case i-1 paired
                                    term1 += exp_dangle5(jj,ii,ii-1)*(u_ip_ju[ii-1] + u_iu_ju[ii-1]);
                                else if (ii > 0)      term1 += exp_dangle5(jj,ii,ii-1);
                                // this case doesn't occur anyway
                                //else                  term1 = 1;  
                                if (jj < seqlen-3)
                                {
                                    int jjp1n = index[jj+1]+seqlen-1-(jj+1);
                                    int jjp2n = index[jj+2]+seqlen-1-(jj+2);
                                    term2 += (u_ip_jp[jjp1n] + u_ip_ju[jjp1n] + 
                                        exp_dangle3(jj,ii,jj+1)*(u_ip_jp[jjp2n] + u_ip_ju[jjp2n] + u_iu_jp[jjp2n] + u_iu_ju[jjp2n]));
                                }    
                                else if (jj < seqlen-1) term2 += exp_dangle3(jj,ii,jj+1);
                                else                  term2 = 1;
                                
                                GlogZ[index_param] += term1 * up[iijj] * term2 * exp_AUpenalty (ii,jj);
                                //printf ("GlogZ[%d]=%g\n", index_param, GlogZ[index_param]);
                            }
                        }
                    }
                    GlogZ[index_param] /= (u_ip_jp[seqlen-1] + u_ip_ju[seqlen-1] + u_iu_jp[seqlen-1] + u_iu_ju[seqlen-1]);                   
                    
                    
                    // contribution from multi-loops
                    // the last dangling 5'
                    // taken from compute_upm                    
                    Complex temp = 0;
                    
                    for (jj=TURN+1; jj < seqlen; jj++)
                    {
                        for (ii=jj-TURN-1; ii>=0; ii--)
                        {                     
                            iijj = index[ii] + jj -ii;  
                            if (upm[iijj].real() == 0)
                                continue;
                            if (sequence[ii] == i && sequence[jj] == j && sequence[jj-1] == k)
                            {       
                                temp = 0;
                                // there must be at least one more branch (i.e. 5 nucleotides) from l+1 to j-1, plus at least 2 free bases: l+1+TURN <= j-3
                                // removed the cases where there isn't a 5' dangling end on the last free base
                                for (l=ii+2; l < jj-TURN-4; l++)    // case ((...)--(--)-..)
                                {
                                    int iip1l = index[ii+1]+l-ii-1;
                                    int lp1jjm1 = index[l+1]+jj-1-l-1;        
                                    int lp2jjm1 = index[l+2]+jj-1-l-2;
                                    temp += up[iip1l] * exp_AUpenalty (ii+1,l) *
                                            (  //u1_ip_jp[lp1jjm1] +     // no dangling end of interest [(...)(.-..)] or [(...)(.-..).]                    
                                             exp_dangle3 (l, ii+1, l+1) * EXPC[1] *
                                                    ( //u1_ip_jp[lp2jjm1] + u1_iu_jp[lp2jjm1] +     // [(...).-(...)] or [(...).-(...).] 
                                                        (u1_ip_ju[lp2jjm1] + u1_iu_ju[lp2jjm1])) +    // [(...).-(...)-..]
                                            + u1_ip_ju[lp1jjm1]);    // [(...)(...)-..]
                                            // TODO shouldn't there be an EXPC[1] here in the last line?
                                }   
                                                             
                                Complex upm_temp = 0;
                                for (l=ii+3; l < jj-TURN-4; l++)    // case (.(...)--(--)-..)
                                {
                                    int iip2l = index[ii+2]+l-ii-2;
                                    int lp1jjm1 = index[l+1]+jj-1-l-1;        
                                    int lp2jjm1 = index[l+2]+jj-1-l-2;
                                    upm_temp += up[iip2l] * exp_AUpenalty (ii+2, l) *
                                            (  //u1_ip_jp[lp1jjm1] +     // no dangling end of interest [.(...)(.-..)] or [.(...)(.-..).]                    
                                              exp_dangle3 (l, ii+2, l+1) * EXPC[1] *
                                                    ( //u1_ip_jp[lp2jjm1] + u1_iu_jp[lp2jjm1] +     // no dangling end [.(...).-(...)] or [.(...).-(...).] 
                                                        (u1_ip_ju[lp2jjm1] + u1_iu_ju[lp2jjm1]))   // [.(...).-(...)-..]
                                              + u1_ip_ju[lp1jjm1]);    // [.(...)(...)-..]
                                }
                                                                
                                for (h=ii+3; h < jj-TURN-4; h++)    // case (....(...)--(--)-..)
                                {
                                    int hjjm1 = index[h] + jj-1 -h;
                                    upm_temp +=  EXPC[h-ii-2] * s2_ju[hjjm1];
                                         //s2_jp[hjjm1] +   // no dangling end --(...)) or --(...).)
                                         // --(...)..)
                                }    
                                temp += exp_dangle3 (ii, jj, ii+1) * EXPC[1] * upm_temp;                                
                                //printf ("i=%d, j=%d, upm=%g\n", i, j, upm[ij]);
                            
                                GlogZ[index_param] += temp * EXPA * EXPB2 * exp_dangle5(ii,jj,jj-1) *
                                    exp_AUpenalty (ii,jj) * p[iijj] / up[iijj];
                            }
                        }
                    }                                                                   

                    
                    // the multi-loop 5' dangling end which can appear to the left of branches
                    // taken from compute_p
                    temp = 0;
                    
                    for (h=1; h < seqlen; h++)   // starts from 1 because h-1 must be valid
                    {
                        for (l=seqlen-1; l > h+TURN; l--)
                        {
                            if (can_pair (sequence[h], sequence[l]) && 
                                sequence[l] == i && sequence[h] == j && sequence[h-1] == k)    
                            {                    
                                // from now on, taken from compute_p      
                                int hl = index[h] + l - h;
                                
                                for (ii=0; ii < h; ii++)
                                {
                                    int iil = index[ii] + l - ii;
                                    // the case when h-l is the first branch of the multi-loop
                                    if (ii < h-2)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
                                    {
                                        temp += EXPC[h-ii-1] * 
                                            (pmd3_noneedmidd3[iil] + pmd3_needmidd3[iil]  * exp_dangle3(l,h,l+1) * EXPC[1]);
                                            // (.-[...](---)-)    i.-h...l(---)-j 
                                            // (.-[...].-(---)-)    i.-h...l.-(---)-j
                                    }
                                    
                                    if (ii < h - TURN -2)    // only now we can have a branch to the left
                                    {        
                                        int iip1hm1 = index[ii+1] + h-1 - (ii+1);
                                        int iip1hm2 = index[ii+1] + h-2 - (ii+1);
                                        int iip1hm3 = index[ii+1] + h-3 - (ii+1);
                                        int iip2hm1 = index[ii+2] + h-1 - (ii+2);
                                        int iip2hm2 = index[ii+2] + h-2 - (ii+2);
                                        int iip2hm3 = index[ii+2] + h-3 - (ii+2);
                                        int iip3hm1 = index[ii+3] + h-1 - (ii+3);
                                        int iip3hm2 = index[ii+3] + h-2 - (ii+3);
                                        int iip3hm3 = index[ii+3] + h-3 - (ii+3);
                                        int iilp1 = index[ii] + l+1 - ii;
                                        Complex term1 = 0.0;
                                        
                                        if (up[iilp1].real() != 0)
                                        {
                                            //printf ("EXPC=%g\n", exp ((misc.multi_free_base_penalty) * oneoverRT) );
                                            term1 = p[iilp1] / up[iilp1] *  exp_AUpenalty (ii,l+1) *
                                            
                                                        (  // first, the case ((..-)-[...])
                                                             u1_ip_ju[iip1hm1]    // ((..-).-[...])
                                                                
                                                            // next, the case (.-(...)-[...])      i..-(...)h...lj
                                                        + exp_dangle3(ii,l+1,ii+1) * EXPC[1] *
                                                            (u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1])
                                                                // (..-(..-).-[...])     
                                                        );
                                        }
                                        
                                        temp += 
                                        
                                        // first, when h-l is the last branch to the right
                                        ( term1 + 
                                        
                                        
                                        // case ((..-)-[...].-)  
                                        pm1nod3_needendd3[iil] * exp_dangle3(l,h,l+1) *   // EXPC added in pm1nod3
                                            u1_ip_ju[iip1hm1]
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-)
                                        + pm1d3_needendd3[iil] * exp_dangle3(l,h,l+1) *     // don't need EXPC[1]
                                            (u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1])
                                                    // (.(...).-[...].-)        i.(...)h...l.-j                                                    
                                        // case (..-(..-)-[...].-)            
                                                    // (..-(...)..-[...].-)      i..-(...)h...l.-j
                                                    
                                        // we have branches to the left and to the right of h-l
                                        //+ pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * u1[iip1hm1]);
                                                    // ((...).-[...]-(---)...)    i(...).-h...l-(---)...j 
                                                                                                        
                                        + (pmnod3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1]+ pmnod3_noneedmidd3[iil]) *
                                            u1_ip_ju[iip1hm1]
                                                    // c is added in pm1nod3_needendd3
                                                    // ((...).-[...].-(--)-)         i(...)h...l.-j
                                        // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                                        + EXPC[1]*(pmd3_needmidd3[iil] * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[iil]) *
                                            (u1_ip_ju[iip2hm1]+u1_iu_ju[iip2hm1]));
                                                    // (.(...).-[...].-)        i.(...)h...l.-j                                                    
                                        // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                                                    // (..-(...)..-[...].-)      i..-(...)h...l.-j
                                                    
                                    }        
                                }
                                GlogZ[index_param] += up[hl] * EXPA * EXPB2 * exp_dangle5(l,h,h-1) * exp_AUpenalty (h,l) * temp;
                            }
                        }
                    }             
                           
                    index_param++;
                }
            }
    index_should_be = structure_type_index("internal_penalty_by_size[4]");            
    if (index_param != index_should_be)
    {
        printf ("Index param after dangle_bot = %d, should be %d\n", index_param, index_should_be);
        exit(1);
    }
          


    }       // end if (compute_gradient_dangles)
    
    // compute derivatives for AU_penalties
    // the contribution from hairpin loop of size 3 was computed above
    // now consider the constributions from the exterior loop

    Complex term1 = 0.0;
    Complex term2 = 0.0;
        
    
    for (ii = 0; ii < seqlen; ii++)
    {
        for (jj = ii+TURN+1; jj < seqlen; jj++)
        {
            if (can_pair (sequence[ii], sequence[jj]) && has_AU_penalty (sequence[ii], sequence[jj]))
            {   
                iijj = index[ii]+jj-ii;
                AUpen = AU_penalty (sequence[ii], sequence[jj]);
    
                
                // first, exterior loop - same as in compute_p            
                
                term1 = 0.0;
                term2 = 0.0;
                    
                IFD
                {
                    if (ii > 0)     term1 = u[ii-1];
                    else            term1 = 1.0;
                }
                else
                {
                    if (ii > TURN)        
                            // case i-1 paired
                        term1 += (u_ip_jp[ii-1] + u_iu_jp[ii-1] + exp_dangle5(jj,ii,ii-1)*(u_ip_ju[ii-1] + u_iu_ju[ii-1]));
                    else if (ii > 0)      term1 += exp_dangle5(jj,ii,ii-1);
                    else                  term1 = 1;  
                }
                                
                IFD
                {
                    if (jj < seqlen-1)    term2 = u[index[jj+1]+seqlen-1-(jj+1)];
                    else                  term2 = 1.0;
                }
                else
                {
                    int jjp1n = index[jj+1] + seqlen-1 - jj-1;
                    int jjp2n = index[jj+2] + seqlen-1 - jj-2;
                    if (jj < seqlen-3)    
                        term2 += u_ip_jp[jjp1n] + u_ip_ju[jjp1n] + 
                                exp_dangle3(jj,ii,jj+1)*(u_ip_jp[jjp2n] + u_ip_ju[jjp2n] + u_iu_jp[jjp2n] + u_iu_ju[jjp2n]);
                    else if (jj < seqlen-1) term2 += exp_dangle3(jj,ii,jj+1);
                    else                  term2 = 1;
                }

                IFD
                {
                    GlogZ[AUpen_index] += term1 *up[iijj] * term2 * exp_AUpenalty(ii,jj) / u[seqlen-1];                        
                }
                else
                {
                    GlogZ[AUpen_index] += term1 *up[iijj] * term2 * exp_AUpenalty(ii,jj) /
                        (u_ip_jp[seqlen-1] + u_ip_ju[seqlen-1] + u_iu_jp[seqlen-1] + u_iu_ju[seqlen-1]);
                }
                // if it closes a multi loop, or closes a helix of a multi-loop - it's considered at helix below
            }                                        
        }
    }               
    
    int sizeindex = structure_type_index ("misc.multi_offset");
    for (ii = 0; ii < seqlen; ii++)
    {
        for (jj = ii+TURN+1; jj < seqlen; jj++)
        {
            iijj = index[ii]+jj-ii;
            // TODO: not so sure test upm[iijj].real() > 0 is right, how about imag()?
            if (upm[iijj].real() > 0)
            {   
                GlogZ[sizeindex] += upm[iijj] * p[iijj] / up[iijj];
            }
        }
    }
    
    sizeindex = structure_type_index ("misc.multi_helix_penalty");
    // took the multi-loop contribution of p, and upm, like above
    // taken from compute_upm and compute_p
    Complex temp;
    for (h = 0; h < seqlen; h++)
    {
        for (l = h+TURN+1; l < seqlen; l++)
        {
            if (!can_pair (sequence[h], sequence[l]))
                continue;
            int hl = index[h]+l-h;
            // the exterior multi-loop pair      
            temp = upm[hl] * p[hl] / up[hl];      
            GlogZ[sizeindex] += temp;
            if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;

            temp = 0;
            for (i=0; i < h; i++)
            {
                int il = index[i] + l - i;  
                IFD
                {
                    temp += EXPC[h-i-1] * pm[il];
                }
                else
                {                          
                    if (i < h-1)  // we must add d3 (which is added in pmd3); also add d5 here, if i < h-2
                    {
                        temp += EXPC[h-i-1] * (i<h-2?exp_dangle5(l,h,h-1):1) *
                            (pmd3_noneedmidd3[il] + pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1]);
                            // (.-[...](---)-)    i.-h...l(---)-j                                 
                            // (.-[...].-(---)-)    i.-h...l.-(---)-j  
                    }
                    else    // case ((... , no dangling end
                    {
                        temp +=  pmnod3_noneedmidd3[il];
                            // ([...](---)-)    ih...l(---)-j
                            
                        temp +=  pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1];
                            // ([...].-(---)-)    ih...l.-(---)-j                 
                    }
                }                    
                if (i < h - TURN -2)    // only now we can have a branch to the left
                {        
                    int ip1hm1 = index[i+1] + h-1 - (i+1);                
                    
                    IFD
                    {
                        temp += u1[ip1hm1] * (pm1[il] + pm[il]);
                    }
                    else
                    {
                        int ip1hm2 = index[i+1] + h-2 - (i+1);
                        int ip1hm3 = index[i+1] + h-3 - (i+1);
                        int ip2hm1 = index[i+2] + h-1 - (i+2);
                        int ip2hm2 = index[i+2] + h-2 - (i+2);
                        int ip2hm3 = index[i+2] + h-3 - (i+2);
                        int ip3hm1 = index[i+3] + h-1 - (i+3);
                        int ip3hm2 = index[i+3] + h-2 - (i+3);
                        int ip3hm3 = index[i+3] + h-3 - (i+3);
                        int ilp1 = index[i] + l+1 - i;
                        Complex term1 = 0.0;
                        if (up[ilp1].real() != 0)
                        {
                            //printf ("EXPC=%g\n", exp ((misc.multi_free_base_penalty) * oneoverRT) );
                            term1 = p[ilp1] / up[ilp1] *  exp_AUpenalty (i,l+1) *
                                        (  // first, the case ((..-)-[...])
                                          ( u1_ip_jp[ip1hm1]    // ((..-)[...]) or ((..-).[...])
                                            + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1] )   // ((..-).-[...])  
                                                
                                            // next, the case (.-(...)-[...])      i..-(...)h...lj
                                          + exp_dangle3(i,l+1,i+1) * EXPC[1] *
                                            ( u1_ip_jp[ip2hm1] + u1_iu_jp[ip2hm1] 
                                                    + exp_dangle5(l,h,h-1)*(u1_ip_ju[ip2hm1]+u1_iu_ju[ip2hm1]) ) 
                                                // (..-(..-).-[...])     
                                        ); 
                        }
                        
                        temp += 
                        
                        // first, when h-l is the last branch to the right
                        ( term1 + 
                        
                        // case ((..-)-[...].-)  
                        pm1nod3_needendd3[il] * exp_dangle3(l,h,l+1) * // we don't need EXPC[1] here, because it's in pm1nod3_needendd3
                            (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                                    // c is added in pm1nod3_needendd3
                                    // ((...).-[...].-)         i(...)h...l.-j
                        // case (.(..-)-[...].-)
                        + pm1d3_needendd3[il] * exp_dangle3(l,h,l+1) * // we don't need EXPC[1] here, because it's in pm1d3_needendd3
                            (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip2hm1])
                                    // (.(...).-[...].-)        i.(...)h...l.-j
                                    
                        // case (..-(..-)-[...].-)            
                        + pm1d3_needendd3[il]  * exp_dangle3(l,h,l+1) * // don't need EXPC[1]
                            (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1])
                                    // (..-(...)..-[...].-)      i..-(...)h...l.-j
                                    
                        // we have branches to the left and to the right of h-l
                        //+ pmnod3_noneedmidd3[il]  * u1[ip1hm1] +    
                                    // ((...)[...]-(---)...)    i(...)h...l-(---)...j 
                        //pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * u1[ip1hm1]);        
                                    // ((...).-[...]-(---)...)    i(...).-h...l-(---)...j 
                          
                        + (pmnod3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1]+ pmnod3_noneedmidd3[il]) *
                            (u1_ip_jp[ip1hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip1hm1])
                                    // c is added in pm1nod3_needendd3
                                    // ((...).-[...].-(--)-)         i(...)h...l.-j
                        // case (.(..-)-[...].-(--)-) or (.(..-)-[...](--)-)
                        + EXPC[1]*(pmd3_needmidd3[il] * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il]) *
                                (u1_ip_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_ip_ju[ip2hm1])
                                    
                        // case (..-(..-)-[...].-(--)-) or (..-(..-)-[...](--)-)
                        + EXPC[1]*(pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) * EXPC[1] + pmd3_noneedmidd3[il])*
                            (u1_iu_jp[ip2hm1] + exp_dangle5(l,h,h-1)*u1_iu_ju[ip2hm1]) ); 
                                    
                    }                
                }        
            }
            temp *= up[hl] * EXPA * EXPB2 * exp_AUpenalty (h,l);
            GlogZ[sizeindex] += temp;
            if (has_AU_penalty (sequence[h],sequence[l]))    GlogZ[AUpen_index] += temp;
        }
    }
                
    sizeindex = structure_type_index ("misc.multi_free_base_penalty");
    // traverse the sequence from left to right         
    int ij, hl;
    int ip1l, ip2l, lp2jm1, lp1jm1, lp2jm2, lp1jm2, lp2jm3, lp1jm3, hj;
    Complex upm_temp;

    // the first free bases
    for (h=0; h < seqlen; h++)   
    {
        for (l=seqlen-1; l > h+TURN; l--)
        {    
            if (!can_pair(sequence[h], sequence[l]))
                continue;
            hl = index[h]+l-h;        
                
            Complex temp = 0;
            for (i=0; i < h-1; i++)
            {
                int il = index[i] + l - i;                
                IFD
                {
                    temp += ((Complex)(h-i-1))* EXPC[h-i-1] * pm[il];
                }
                else
                {
                    // the case when h-l is the first branch of the multi-loop
                    temp += ((Complex)(h-i-1)) * EXPC[h-i-1] * (i<h-2?exp_dangle5(l,h,h-1):1) *
                        (pmd3_noneedmidd3[il] + pmd3_needmidd3[il]  * exp_dangle3(l,h,l+1) *  EXPC[1]);
                        // (.-[...](---)-)    i.-h...l(---)-j
                        // (.-[...].-(---)-)    i.-h...l.-(---)-j
                }
            }
            GlogZ[sizeindex]  += up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB2 * temp;
        }
    }
    
    //printf ("GlogZ[%d] = %g\n", sizeindex, GlogZ[sizeindex]);
       
    // second, there is at least a branch to the left of k, and a branch to the right
    // Now, k is the free base inside a multi-loop, closed be i,j     
    for (i=0; i < seqlen-2*TURN-5; i++)
    {
        for (j=i+2*TURN+5; j < seqlen; j++)
        {
            ij = index[i]+j-i;
            if (upm[ij].real() > 0)
            {
                Complex temp = 0;
                for (k=i+TURN+3; k < j-TURN-2; k++)
                {
                    int ip1km1 = index[i+1] + k-1 - (i+1);
                    int kp1jm1 = index[k+1] + j-1 - (k+1);
                    
                    IFD
                    {
                        temp += EXPC[1] * u1[ip1km1] * u1[kp1jm1];
                    }
                    else
                    {
                        int ip1k = index[i+1] + k - (i+1);                    
                        int ip1km2 = index[i+1] + k-2 - (i+1);
                        int ip2k = index[i+2] + k - (i+2);
                        int ip2km1 = index[i+2] + k-1 - (i+2);
                        int ip2km2 = index[i+2] + k-2 - (i+2);
                        int kjm1 = index[k] + j-1 - (k);
                        int kjm2 = index[k] + j-2 - (k);
                        int kjm3 = index[k] + j-3 - (k);                    
                        int kp1jm2 = index[k+1] + j-2 - (k+1);
                        int kp1jm3 = index[k+1] + j-3 - (k+1);
                        // first situation (---)-.k---
                                          
                        temp +=
                            
                            // the branch(es) to the left of k
                            
                            ( exp_dangle3 (i, j, i+1) * EXPC[1]* (u1_ip_ju[ip2k] + u1_iu_ju[ip2k])    // (.-(---)-.k(
                              + u1_ip_ju[ip1k] ) *   // ((---).k(
                            // the branch(es) to the right of k
                            ( u1_iu_jp[kjm1] / EXPC[1] + u1_iu_ju[kjm1]/EXPC[1] *exp_dangle5(i,j,j-1));
                                           
                        // second situation (---)k---
                        temp +=
                            // the branch(es) to the left of k
                            (exp_dangle3 (i, j, i+1) * EXPC[1]*(u1_ip_ju_jm1p[ip2k] + u1_iu_ju_jm1p[ip2k]) + u1_ip_ju_jm1p[ip1k]) *    // only j should be paired, not j-1
                            // the branch(es) to the right of k
                            (u1_ip_jp[kp1jm1] + u1_iu_jp[kp1jm1] + (u1_ip_ju[kp1jm1] + u1_iu_ju[kp1jm1])*exp_dangle5(i,j,j-1));                                                        
                    }
                }
                GlogZ[sizeindex] += p[ij] / up[ij] * EXPA * EXPB1 * exp_AUpenalty(i,j) * temp;
            }
        }
    }



    // last, no branch to the right of k
    // similar to the first free bases 
    for (h=0; h < seqlen; h++)   
    {
        int hj;
        for (l=seqlen-3; l > h+TURN; l--)
        {    
            if (!can_pair(sequence[h], sequence[l]))
                continue;
            hl = index[h]+l-h;        
                
            Complex temp = 0;
            IFD
            {
                //GlogZ[sizeindex] += (j-l-1)*up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB2 * EXPC[j-l-1] * pm2[hj];
                for (j=l+2; j < seqlen; j++)
                {
                    hj = index[h] + j - h;
                    temp += ((Complex)(j-l-1)) * EXPC[j-l-1] * pm2[hj];
                }
                GlogZ[sizeindex] += up[hl] * exp_AUpenalty (h,l) * EXPA * EXPB2 * temp;
            }
            else
            {
                j = l+2;
                hj = index[h] + j - h;

                temp += ((Complex)(j-l-1))* EXPC[j-l-1] *
                    (pm2nod5_noneedmidd5[hj] + pm2nod5_needmidd5[hj] * exp_dangle5(l,h,h-1));
                    // (-(---)[...].) or (-(---).[...].)
                    // (-(---)-..[...].)

                for (j=l+3; j < seqlen; j++)
                {
                    // needs the last 5' dangling end
                    hj = index[h] + j - h;
                    temp += ((Complex)(j-l-1))* EXPC[j-l-1] *
                                (pm2d5_noneedmidd5[hj] + pm2d5_needmidd5[hj] * exp_dangle5(l,h,h-1));
                                // (-(---)[...].-.) or (-(---).[...].-.)
                                // (-(---)-..[...].-.)
                }
                GlogZ[sizeindex] += up[hl] * exp_dangle3(l,h,l+1) * exp_AUpenalty (h,l) * EXPA * EXPB2 * temp;
            }            
        }
    }    
    
}



Complex s_partition_function_complex::get_hairpin_energy (int i, int j, int* sequence, char *csequence)
// returns the COMPLEX free energy of the hairpin loop closed at (i,j)
{
    Complex energy=0;

    Complex terminal_mismatch_energy = 0, bonus = 0, special_bonus = 0, AU_pen = 0;
    int k, is_poly_C;
    int size;
    char seq[6] = "";


    size = j-i-1;

    // TODO 
//     if (size < 3)
//         return INF;
//     return (Complex)500;


    if (size < 3)
        return INF;            
    else if (size == 3)
    {
        terminal_mismatch_energy = 0;
        AU_pen = AU_penalty_complex (sequence[i], sequence[j]);
    }
    else
    {
        terminal_mismatch_energy =
             tstackh_complex[sequence[i]]
                    [sequence[j]]
                    [sequence[i+1]]
                    [sequence[j-1]];
    }

    #if (MODEL == SIMPLE)
    // check if it is a triloop
    if (size == 3)
    {
        substr (csequence, i, j, seq);
        for (k=0; k < nb_triloops; k++)
        {
            // don't need triloop_complex[].seq
            if (strcmp (seq, triloop[k].seq) == 0)
            //if (strcmp (seq, triloop_complex[k].seq) == 0)            
                bonus = triloop_complex[k].energy;
        }
    }

    // check to see it is a tetraloop in tloop
    else if (size == 4)
    {
        substr (csequence, i, j, seq);
        for (k=0; k < nb_tloops; k++)   // nb_tloops should be the same as for energy
        {
            if (strcmp (seq, tloop[k].seq) == 0)
            //if (strcmp (seq, tloop_complex[k].seq) == 0)
                bonus = tloop_complex[k].energy;
        }
    }
    #endif

    // special_bonus from miscloop file
    // check if we have to apply "GGG" loop special bonus
    // Vienna package doesn't have it

    if (i > 1)
    {
        if (sequence[i-2] == G && sequence[i-1] == G &&
            sequence[i] == G && sequence[j] == U)
            special_bonus += misc_complex.hairpin_GGG;
    }


    // check for the special case of "poly-C" hairpin loop
    is_poly_C = 1;
    for (k=i+1; k<j; k++)
    {
        if (sequence[k] != C)
        {
            is_poly_C = 0;
            break;
        }
    }
    if (is_poly_C)
    {
        if (size == 3)
            special_bonus += misc_complex.hairpin_c3;
        else
            special_bonus += misc_complex.hairpin_c2 + misc_complex.hairpin_c1 * (Complex) size;
    }

    energy = penalty_by_size_complex (size, 'H') +
             terminal_mismatch_energy + bonus + special_bonus + AU_pen;

    return (Complex) energy;
}



Complex s_partition_function_complex::get_stacked_energy (int i, int j, int *sequence)
// returns the free energy of the stacked pair closed at (i,j)
{
    int k,l;
    k=i+1;
    l=j-1;
    // TODO

/*    if (sequence[i] == 0 && sequence[j] == 3 && sequence[k] == 0 && sequence[l] == 3)
        return stack_complex[0][3][0][3];
    else if (sequence[l] == 0 && sequence[k] == 3 && sequence[j] == 0 && sequence[i] == 3)
        return stack_complex[0][3][0][3];
    else
        return (Complex)-100;  */      
        
    if (sequence[i]*1000 + sequence[j]*100 + sequence[k]*10 + sequence[l] <=
        sequence[l]*1000 + sequence[k]*100 + sequence[j]*10 + sequence[i])
        return stack_complex [sequence[i]] [sequence[j]] [sequence[k]] [sequence[l]];
    return stack_complex [sequence[l]] [sequence[k]] [sequence[j]] [sequence[i]];
                            
}


Complex s_partition_function_complex::get_internal_energy (int i, int j, int ip, int jp, int *sequence)
// returns the free energy of the internal loop closed at (i,j,ip,jp)
// very similar to function s_internal_loop::count_get_energy
{
    // TODO
    //return (Complex)100;

    Complex energy = 0;
    Complex penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
    int branch1, branch2, len;

    int k,l,m,n,o,p;


    branch1 = ip-i-1;
    branch2 = j-jp-1;

    if (branch1 != 0 || branch2 != 0)
    {
        // check if it is a bulge loop of size 1
        // check if it is int11 or int21 or int22

        if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
        {                        
            if (sequence[i]*100000 + sequence[j]*10000 + sequence[i+1]*1000 + sequence[j-1]*100 + sequence[ip]*10 + sequence[jp] >
                sequence[jp]*100000 + sequence[ip]*10000 + sequence[j-1]*1000+ sequence[i+1]*100 + sequence[j]*10 + sequence[i])
            {
                int ii, jj, kk, ll, mm, nn;
                ii = sequence[jp];
                jj = sequence[ip];
                kk = sequence[j-1];
                ll = sequence[i+1];
                mm = sequence[j];
                nn = sequence[i];
                if ( ((ii==C && jj==G) || (ii==G && jj==C)) && ((mm==C && nn==G) || (mm==G && nn==C)))
                {
                    if (!can_pair(kk,ll))
                        energy += int11_complex[ii][jj][kk][ll][mm][nn];
                    else
                        energy += misc_complex.internal11_basic_mismatch;
                }
                else if (watson_crick(ii,jj) && watson_crick(mm,nn) && kk==U && ll==U)
                {
                    energy += int11_complex[ii][jj][kk][ll][mm][nn];
                }
                else
                {
                    if (kk==G && ll==G)
                        energy += misc_complex.internal11_GG_mismatch;
                    else
                        energy += misc_complex.internal11_basic_mismatch;

                    if (has_AU_penalty(ii,jj))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                    if (has_AU_penalty(mm,nn))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                }
            }
            else
            {
                int ii, jj, kk, ll, mm, nn;
                ii = sequence[i];
                jj = sequence[j];
                kk = sequence[i+1];
                ll = sequence[j-1];
                mm = sequence[ip];
                nn = sequence[jp];

                if ( ((ii==C && jj==G) || (ii==G && jj==C)) && ((mm==C && nn==G) || (mm==G && nn==C)))
                {
                    if (!can_pair(kk,ll))
                        energy += int11_complex[ii][jj][kk][ll][mm][nn];
                    else
                        energy += misc_complex.internal11_basic_mismatch;
                }
                else if (watson_crick(ii,jj) && watson_crick(mm,nn) && kk==U && ll==U)
                {
                    energy += int11_complex[ii][jj][kk][ll][mm][nn];
                }
                else
                {
                    if (kk==G && ll==G)
                        energy += misc_complex.internal11_GG_mismatch;
                    else
                        energy += misc_complex.internal11_basic_mismatch;

                    if (has_AU_penalty(ii,jj))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                    if (has_AU_penalty(mm,nn))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                }
            }
        }
        else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
        {
            int ii, jj, kk, ll, mm, nn, oo;
            ii = sequence[i];
            jj = sequence[j];
            kk = sequence[i+1];
            ll = sequence[j-1];
            mm = sequence[ip];
            nn = sequence[jp];
            oo = sequence[jp+1];            

            if ((ii==C && jj==G && mm==C && nn==G) ||  // these are already filled above, except what can pair inside
                (ii==G && jj==C && mm==G && nn==C))
            {
                if (can_pair(kk,ll) || can_pair(kk,oo))
                    energy += misc_complex.internal21_match;
                else
                    energy += int21_complex[ii][jj][kk][ll][mm][nn][oo];
            }
            else
            {
                if (can_pair(kk,ll) || can_pair(kk,oo))
                {
                    energy += misc_complex.internal21_match;
                }
                else
                {
                    if (params_are_double)
                    {
                        // for double parameters, we don't want conversion to int
                        energy += Complex(0.5*int21_complex[C][G][kk][ll][C][G][oo].real(), 0.5*int21_complex[C][G][kk][ll][C][G][oo].imag());
                        energy += Complex(0.5*int21_complex[G][C][kk][ll][G][C][oo].real(), 0.5*int21_complex[G][C][kk][ll][G][C][oo].imag());
                    }
                    else
                    {
                        // we want the real part to be converted to int, so that it's the same as s_partition_function.
                        energy += Complex((int)(0.5*int21_complex[C][G][kk][ll][C][G][oo].real()), 0.5*int21_complex[C][G][kk][ll][C][G][oo].imag());
                        energy += Complex((int)(0.5*int21_complex[G][C][kk][ll][G][C][oo].real()), 0.5*int21_complex[G][C][kk][ll][G][C][oo].imag());
                    }
                }        
                if (has_AU_penalty(ii,jj))
                {
                    energy += misc_complex.internal21_AU_closure;
                }    
                if (has_AU_penalty(mm,nn))    
                {
                    energy += misc_complex.internal21_AU_closure;
                }    
            }                                    
        }
        else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
        {
            // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]

            int ii, jj, kk, ll, mm, nn, oo;
            ii = sequence[jp];
            jj = sequence[ip];
            kk = sequence[j-1];
            ll = sequence[ip-1];
            mm = sequence[j];
            nn = sequence[i];
            oo = sequence[i+1];            
            if ((ii==C && jj==G && mm==C && nn==G) ||  // these are already filled above, except what can pair inside
                (ii==G && jj==C && mm==G && nn==C))
            {
                if (can_pair(kk,ll) || can_pair(kk,oo))
                    energy += misc_complex.internal21_match;
                else
                    energy += int21_complex[ii][jj][kk][ll][mm][nn][oo];
            }
            else
            {
                if (can_pair(kk,ll) || can_pair(kk,oo))
                {
                    energy += misc_complex.internal21_match;
                }
                else
                {
                    if (params_are_double)
                    {
                        energy += Complex(0.5*int21_complex[C][G][kk][ll][C][G][oo].real(), 0.5*int21_complex[C][G][kk][ll][C][G][oo].imag());
                        energy += Complex(0.5*int21_complex[G][C][kk][ll][G][C][oo].real(), 0.5*int21_complex[G][C][kk][ll][G][C][oo].imag());
                    }
                    else
                    {
                        energy += Complex((int)(0.5*int21_complex[C][G][kk][ll][C][G][oo].real()), 0.5*int21_complex[C][G][kk][ll][C][G][oo].imag());
                        energy += Complex((int)(0.5*int21_complex[G][C][kk][ll][G][C][oo].real()), 0.5*int21_complex[G][C][kk][ll][G][C][oo].imag());
                    }
                }        
                if (has_AU_penalty(ii,jj))
                {
                    energy += misc_complex.internal21_AU_closure;
                }    
                if (has_AU_penalty(mm,nn))    
                {
                    energy += misc_complex.internal21_AU_closure;
                }    
            }
        }
        else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
        {
            // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
            if (sequence[i]*10000000 + sequence[j]*1000000 + sequence[i+1]*100000 + sequence[j-1]*10000 +
                sequence[ip]*1000 + sequence[jp]*100 + sequence[ip-1]*10 + sequence[jp+1] >
                sequence[jp]*10000000 + sequence[ip]*1000000 + sequence[jp+1]*100000 + sequence[ip-1]*10000 +
                sequence[j]*1000 + sequence[i]*100 + sequence[j-1]*10 + sequence[i+1])
            {
                int ii, jj, kk, ll, mm, nn, oo, pp;
                ii = sequence[jp];
                jj = sequence[ip];
                kk = sequence[jp+1];
                ll = sequence[ip-1];
                mm = sequence[j];
                nn = sequence[i];
                oo = sequence[j-1];            
                pp = sequence[i+1];
                
                if (nn==ii && mm==jj && pp==kk && oo==ll & watson_crick(ii,jj) && !watson_crick(kk,ll))
                {
                    energy += int22_complex[ii][jj][kk][ll][mm][nn][oo][pp];
                }
             
                int ii2, jj2, mm2, nn2;                
                if (ii==G && jj==U)   ii2 = A;     else ii2 = ii;
                if (ii==U && jj==G)   jj2 = A;     else jj2 = jj;
                if (mm==G && nn==U)   mm2 = A;     else mm2 = mm;
                if (mm==U && nn==G)   nn2 = A;     else nn2 = nn;
                
                if (watson_crick(kk,ll) || watson_crick(oo,pp))
                {
                    energy += misc_complex.internal22_match;
                }
                else if ( ((ii==G && jj==U) || (ii==U && jj==G) || (mm==G && nn==U) || (mm==U && nn==G)) &&  
                          (nn2==ii2 && mm2==jj2 && pp==kk && oo==ll))  // the UG closing pairs are the same as UA
                {
                    energy += int22_complex[ii2][jj2][kk][ll][mm2][nn2][oo][pp];
                }                
                else if (!(nn==ii && mm==jj && pp==kk && oo==ll))   // was already filled above
                {
                    if (params_are_double)
                    {
                        energy += Complex(0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].real(), 0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].imag());
                        energy += Complex(0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].real(), 0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].imag());
                    }
                    else
                    {
                        energy += Complex((int)(0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].real()), 0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].imag());
                        energy += Complex((int)(0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].real()), 0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].imag());
                    }
                    
                    int result = check_stability_and_size (kk, ll, oo, pp);                    
                    switch (result)
                    {
                        case 1: 
                            energy += misc_complex.internal22_delta_same_size;
                            break;
                        case 2: 
                            energy += misc_complex.internal22_delta_different_size;
                            break;
                        case 3: 
                            energy += misc_complex.internal22_delta_1stable_1unstable;
                            break;
                        case 4: 
                            energy += misc_complex.internal22_delta_AC;
                            break;
                        default: 
                            printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kk, ll, oo, pp);
                            exit(1);
                    }
                }
                         
                         
            }
            else
            {
                int ii, jj, kk, ll, mm, nn, oo, pp;
                ii = sequence[i];
                jj = sequence[j];
                kk = sequence[i+1];
                ll = sequence[j-1];
                mm = sequence[ip];
                nn = sequence[jp];
                oo = sequence[ip-1];            
                pp = sequence[jp+1];
                
                if (nn==ii && mm==jj && pp==kk && oo==ll & watson_crick(ii,jj) && !watson_crick(kk,ll))
                {
                    energy += int22_complex[ii][jj][kk][ll][mm][nn][oo][pp];
                }
             
                int ii2, jj2, mm2, nn2;                
                if (ii==G && jj==U)   ii2 = A;     else ii2 = ii;
                if (ii==U && jj==G)   jj2 = A;     else jj2 = jj;
                if (mm==G && nn==U)   mm2 = A;     else mm2 = mm;
                if (mm==U && nn==G)   nn2 = A;     else nn2 = nn;
                
                if (watson_crick(kk,ll) || watson_crick(oo,pp))
                {
                    energy += misc_complex.internal22_match;
                }
                else if ( ((ii==G && jj==U) || (ii==U && jj==G) || (mm==G && nn==U) || (mm==U && nn==G)) &&  
                          (nn2==ii2 && mm2==jj2 && pp==kk && oo==ll))  // the UG closing pairs are the same as UA
                {
                    energy += int22_complex[ii2][jj2][kk][ll][mm2][nn2][oo][pp];
                }
                else if (!(nn==ii && mm==jj && pp==kk && oo==ll))   // was already filled above
                {
                    if (params_are_double)
                    {
                        energy += Complex(0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].real(), 0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].imag());
                        energy += Complex(0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].real(), 0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].imag());
                    }
                    else
                    {
                        energy += Complex((int)(0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].real()), 0.5*int22_complex[ii2][jj2][kk][ll][jj2][ii2][ll][kk].imag());
                        energy += Complex((int)(0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].real()), 0.5*int22_complex[nn2][mm2][pp][oo][mm2][nn2][oo][pp].imag());
                    }
                    
                    int result = check_stability_and_size (kk, ll, oo, pp);                    
                    switch (result)
                    {
                        case 1: 
                            energy += misc_complex.internal22_delta_same_size;
                            break;
                        case 2: 
                            energy += misc_complex.internal22_delta_different_size;
                            break;
                        case 3: 
                            energy += misc_complex.internal22_delta_1stable_1unstable;
                            break;
                        case 4: 
                            energy += misc_complex.internal22_delta_AC;
                            break;
                        default: 
                            printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kk, ll, oo, pp);
                            exit(1);
                    }
                }
            }
        }

        else
        {
            // this case is not int11, int21, int22

            // check if it is a bulge
            if (branch1 == 0 || branch2 == 0)
            {
                len = branch1+branch2;
                penalty_size = penalty_by_size_complex (len, 'B');
                if (len == 1)
                {
                    // bulge of size 1
                    // stack[i][j][i+1][j-1]
                    if (sequence[i]*1000 + sequence[j]*100 + sequence[ip]*10 + sequence[jp] <=
                        sequence[jp]*1000 + sequence[ip]*100 + sequence[j]*10 + sequence[i])
                        energy = stack_complex [sequence[i]][sequence[j]][sequence[ip]][sequence[jp]] +
                               penalty_size;
                    else
                        energy = stack_complex [sequence[jp]][sequence[ip]][sequence[j]][sequence[i]] +
                               penalty_size;
                    
                }
                else
                {
                    // bulge of size bigger than 1
                    // check if (i,j) and (ip,jp) can pair
                    energy = penalty_size +
                               AU_penalty_complex (sequence[i],sequence[j]) +
                               AU_penalty_complex (sequence[ip], sequence[jp]);
                }
            }
            // it is an internal loop (not a bulge)
            else
            {
                len = branch1+branch2;
                energy += penalty_by_size_complex (len, 'I');
                energy += asymmetry_penalty_complex (branch1, branch2);

                if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
                // If gail_rule is set to 1 in miscloop file,
                // i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                {
                    // actually, just use 3 parameters instead of the tstacki table
                    if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                        ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                    // actually, just use 3 parameters instead of the tstacki table
                    if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                        ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                }
                else
                {
                    // actually, just use 3 parameters instead of the tstacki table
                    if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
                        ((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
                    {
                        energy += misc_complex.internal_AU_closure;
                    }
                    if ((sequence[i+1] == A && sequence[j-1] == G) ||
                        (sequence[j-1] == A && sequence[i+1] == G))
                    {
                        energy += misc_complex.internal_AG_mismatch;
                    }
                    if (sequence[i+1] == U && sequence[j-1] == U)
                    {
                        energy += misc_complex.internal_UU_mismatch;
                    }                                          
                                          
                    // actually, just use 3 parameters instead of the tstacki table
                    if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
                        ((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
                    {  
                        energy += misc_complex.internal_AU_closure;
                    }
                    if ((sequence[ip-1] == A && sequence[jp+1] == G) ||
                        (sequence[jp+1] == A && sequence[ip-1] == G))
                    {
                        energy += misc_complex.internal_AG_mismatch;
                    }
                    if (sequence[ip-1] == U && sequence[jp+1] == U)
                    {
                        energy += misc_complex.internal_UU_mismatch;
                    }                                                                      
                }
            }
        }
    }
    //printf ("COMPLEX: i=%d, j=%d, ip=%d, jp=%d, energy=%g\n", i, j, ip, jp, energy.real());
    return energy;
}


#endif

