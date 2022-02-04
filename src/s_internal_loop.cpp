/***************************************************************************
                          s_internal_loop.cpp  -  description
                             -------------------
    begin                : Fri Apr 12 2002
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

// a class for internal loop related functions
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "s_internal_loop.h"
#include "externs.h"
#include "common.h"
#include "simfold.h"
#include "params.h"


s_internal_loop::s_internal_loop (int *seq, int length)
// The constructor
{
    seqlen = length;
    sequence = seq;
    this->V = NULL;
}


s_internal_loop::~s_internal_loop ()
// The destructor
{
}

PARAMTYPE s_internal_loop::compute_energy (int i, int j)
// computes the MFE of the structure closed by an internal loop closed at (i,j)
{
    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
    int branch1, branch2, l;
    mmin = INF;

    for (ip = i+1; ip <= MIN(j-2-TURN,i+MAXLOOP+1) ; ip++){  // j-2-TURN
        minq = MAX (j-i+ip-MAXLOOP-2, ip+1+TURN);    // ip+1+TURN);
        for (jp = minq; jp < j; jp++){        
            // could replace the code from now to the end of the function, by get_energy_str
            // it's very redundant (duplicated code) as it is now, but it's faster.
            if (sequence[ip]+sequence[jp] == 3 || sequence[ip]+sequence[jp] == 5){
                branch1 = ip-i-1;
                branch2 = j-jp-1;
                if (branch1 == 0 && branch2 == 0)
                    continue;

                if (branch1 != 0 || branch2 != 0){
                    // check if it is a bulge loop of size 1
                    // check if it is int11 or int21 or int22
                    // 1x1 internal loop
                    if (branch1 == 1 && branch2 == 1){ // && !simple_internal_energy){     // it is int11
                        // int11[i][j][i+1][j-1][ip][jp]
                        en = int11 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin){
                            mmin = ttmp;
                        }
                        continue;
                    // 1x2 internal loop
					}else if (branch1 == 1 && branch2 == 2){ // && !simple_internal_energy){
                        // int21[i][j][i+1][j-1][ip][jp][jp+1]
                        en = int21 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]]
                                   [sequence[jp+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin){
                            mmin = ttmp;
                        }
                        continue;
                    // 1x2 internal loop
					}else if(branch1 == 2 && branch2 == 1){ // && !simple_internal_energy){
                        // after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
                        en = int21 [sequence[jp]][sequence[ip]]
                                   [sequence[j-1]][sequence[ip-1]]
                                   [sequence[j]][sequence[i]]
                                   [sequence[i+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin){
                            mmin = ttmp;
                        }
                        continue;
                    // 2x2 internal loop
					}else if (branch1 == 2 && branch2 == 2){ // && !simple_internal_energy){
                        // int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
                        en = int22 [sequence[i]][sequence[j]]
                                   [sequence[i+1]][sequence[j-1]]
                                   [sequence[ip]][sequence[jp]]
                                   [sequence[ip-1]][sequence[jp+1]];
                        ttmp = en + V->get_energy (ip, jp);
                        if (ttmp < mmin){
                            mmin = ttmp;
                        }
                        continue;
                    }else{
                        // this case is not int11, int21, int22
                        
                        // check if it is a bulge
                        if (branch1 == 0 || branch2 == 0){
                            l = branch1+branch2;
                            penalty_size = penalty_by_size (l, 'B');
                            if (l == 1){
                                // bulge of size 1
                                // stack[i][j][i+1][j-1]
                                en = stack [sequence[i]][sequence[j]]
                                           [sequence[ip]][sequence[jp]];
                                ttmp = en + penalty_size + V->get_energy (ip, jp);
                                if (ttmp < mmin){
                                    mmin = ttmp;
                                }
                                continue;
                            }else{
                                // bulge of size bigger than 1
                                // check if (i,j) and (ip,jp) can pair
                                ttmp = penalty_size + AU_penalty(sequence[i],sequence[j]) + AU_penalty (sequence[ip],sequence[jp]) + V->get_energy(ip, jp);
                                if (ttmp < mmin){
                                    mmin = ttmp;
                                }
                                continue;
                            }
                        // it is a regular internal loop (not a bulge)
						}else{
							//printf("HERE\n");
                            l = branch1 + branch2;
                            penalty_size = penalty_by_size (l, 'I');
                            asym_penalty = asymmetry_penalty(branch1, branch2);

                            if ((branch1 == 1 || branch2 == 1) && misc.gail_rule){
                            	// If gail_rule is set to 1 in miscloop file,
                            	// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
                                i_j_energy  =  tstacki[sequence[i]][sequence[j]] [0][0];
                                ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [0][0];
                            }else{
                                i_j_energy   = tstacki[sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]];
                                ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [sequence[jp+1]][sequence[ip-1]];
                            }
                            ttmp = i_j_energy + ip_jp_energy + penalty_size + asym_penalty + V->get_energy(ip, jp);
                            if (ttmp < mmin){
                                mmin = ttmp;
                            }
                            continue;
                        }
                    }
                }
            }
        }
    }
    return mmin;
}


PARAMTYPE s_internal_loop::compute_energy_restricted (int i, int j, str_features *fres){
// computes the MFE of the structure closed by a restricted internal loop closed by (i,j)

    int ip, jp, minq;
    PARAMTYPE mmin, ttmp;
    mmin = INF;

    for (ip = i+1; ip <= MIN(j-2,i+MAXLOOP+1) ; ip++){  // the -TURN shouldn't be there
        minq = MAX (j-i+ip-MAXLOOP-2, ip+1);    // without TURN
        for (jp = minq; jp < j; jp++){
            if (exists_restricted (i,ip,fres) || exists_restricted (jp,j,fres))
                continue;
            ttmp = get_energy_str (i, j, ip, jp);
            if (ttmp < mmin){
                mmin = ttmp;
            }
        }
    }
    return mmin;
}

PARAMTYPE s_internal_loop::get_energy_str (int i, int j, int ip, int jp){
	// returns the free energy of the structure closed by the internal loop (i,j,ip,jp)
	// This function is just most of what is inside the second for loop of compute_energy
	// TODO
	//return 0;
	//printf ("\n2\n");

	PARAMTYPE mmin, ttmp;
	PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy, en;
	int branch1, branch2, l;
	mmin = INF;

	if (sequence[ip]+sequence[jp] == 3 || sequence[ip]+sequence[jp] == 5){
		branch1 = ip-i-1;
		branch2 = j-jp-1;
		if (branch1 != 0 || branch2 != 0){
			// check if it is a bulge loop of size 1
			// check if it is int11 or int21 or int22
			if (branch1 == 1 && branch2 == 1 && !simple_internal_energy){     // it is int11
				// int11[i][j][i+1][j-1][ip][jp]
				en = int11 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]];
				ttmp = en + V->get_energy(ip, jp);
				if (ttmp < mmin){
					mmin = ttmp;
				}
			}else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy){
				// int21[i][j][i+1][j-1][ip][jp][jp+1]
				en = int21 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]] [sequence[jp+1]];
				ttmp = en + V->get_energy(ip, jp);
				if (ttmp < mmin){
					mmin = ttmp;
				}
			}else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy){
				// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
				en = int21 [sequence[jp]][sequence[ip]] [sequence[j-1]][sequence[ip-1]] [sequence[j]][sequence[i]] [sequence[i+1]];
				ttmp = en + V->get_energy(ip, jp);
				if (ttmp < mmin){
					mmin = ttmp;
				}
			}else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy){
				// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
				en = int22 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]] [sequence[ip-1]][sequence[jp+1]];
				ttmp = en + V->get_energy (ip, jp);
				if (ttmp < mmin){
					mmin = ttmp;
				}
			}else{
				// this case is not int11, int21, int22

				// check if it is a bulge
				if (branch1 == 0 || branch2 == 0){
					l = branch1+branch2;
					penalty_size = penalty_by_size (l, 'B');
					if (l == 1){
						// bulge of size 1
						// stack[i][j][i+1][j-1]
						en = stack [sequence[i]][sequence[j]] [sequence[ip]][sequence[jp]];
						ttmp = en + penalty_size + V->get_energy (ip, jp);
						if (ttmp < mmin){
							mmin = ttmp;
						}
					}else{
						// bulge of size bigger than 1
						// check if (i,j) and (ip,jp) can pair
						ttmp = penalty_size + AU_penalty(sequence[i],sequence[j]) + AU_penalty(sequence[ip],sequence[jp]) + V->get_energy(ip, jp);
						if (ttmp < mmin){
							mmin = ttmp;
						}
					}
				// it is an internal loop (not a bulge)
				}else{
					l = branch1+branch2;
					penalty_size = penalty_by_size (l, 'I');
					asym_penalty = asymmetry_penalty (branch1, branch2);

					if ((branch1 == 1 || branch2 == 1) && misc.gail_rule){
						// If gail_rule is set to 1 in miscloop file,
						// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
						i_j_energy  =  tstacki[sequence[i]][sequence[j]] [0][0];
						ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [0][0];
					}else{
						i_j_energy   = tstacki[sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]];
						ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [sequence[jp+1]][sequence[ip-1]];
					}
					ttmp = i_j_energy + ip_jp_energy + penalty_size + asym_penalty + V->get_energy(ip, jp);
					if (ttmp < mmin){
						mmin = ttmp;
					}
				}
			}
		}
	}    
	if (mmin < INF/2)
		return mmin;
	return INF;            
}


// Not used. Tried to see if the gradient in s_partition_function gets computed faster
//  if you make a different function for each situation.
// It doesn't seem to be the case.
PARAMTYPE s_internal_loop::get_energy_00 (int i, int j, int ip, int jp, int *sequence)
{
	PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
	PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
	int branch1, branch2, l;

	branch1 = ip-i-1;
	branch2 = j-jp-1;
	l = branch1+branch2;
	penalty_size = penalty_by_size (l, 'I');
	asym_penalty = asymmetry_penalty (branch1, branch2);

	i_j_energy  =  tstacki[sequence[i]][sequence[j]]
		[0][0];
	ip_jp_energy = tstacki[sequence[jp]][sequence[ip]]
		[0][0];
	energy = i_j_energy + ip_jp_energy + penalty_size +
		asym_penalty;
	return energy;
}                            


PARAMTYPE s_internal_loop::get_energy (int i, int j, int ip, int jp, int *sequence, int *ptable){
	// returns the free energy of the internal loop closed at (i,j,ip,jp)

	PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
	PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
	int branch1, branch2, l;
	if (exists_restricted_ptable (i,ip,ptable) || exists_restricted_ptable (jp,j,ptable))
		return INF;

	branch1 = ip-i-1;
	branch2 = j-jp-1;

	if (branch1 != 0 || branch2 != 0){
		// check if it is a bulge loop of size 1
		// check if it is int11 or int21 or int22
		if (branch1 == 1 && branch2 == 1 && !simple_internal_energy){     // it is int11
			// int11[i][j][i+1][j-1][ip][jp]
			energy = int11 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]];
		}else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy){
			// int21[i][j][i+1][j-1][ip][jp][jp+1]
			energy = int21 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]] [sequence[jp+1]];
		}else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy){
			// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
			energy = int21 [sequence[jp]][sequence[ip]] [sequence[j-1]][sequence[ip-1]] [sequence[j]][sequence[i]] [sequence[i+1]];
		}else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy){
			// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
			energy = int22 [sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]] [sequence[ip]][sequence[jp]] [sequence[ip-1]][sequence[jp+1]];
		}else{
			// this case is not int11, int21, int22

			// check if it is a bulge
			if (branch1 == 0 || branch2 == 0){
				l = branch1+branch2;
				penalty_size = penalty_by_size (l, 'B');
				if (l == 1){
					// bulge of size 1
					// stack[i][j][i+1][j-1]
					energy = stack [sequence[i]][sequence[j]] [sequence[ip]][sequence[jp]] + penalty_size;
				}else{
					// bulge of size bigger than 1
					// check if (i,j) and (ip,jp) can pair
					energy = penalty_size + AU_penalty(sequence[i],sequence[j]) + AU_penalty(sequence[ip], sequence[jp]);
				}
			// it is an internal loop (not a bulge)
			}else{
				l = branch1 + branch2;
				penalty_size = penalty_by_size(l, 'I');
				asym_penalty = asymmetry_penalty(branch1, branch2);

				if ((branch1 == 1 || branch2 == 1) && misc.gail_rule){
					// If gail_rule is set to 1 in miscloop file,
					// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
					i_j_energy  =  tstacki[sequence[i]][sequence[j]] [0][0];
					ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [0][0];
				}else{
					i_j_energy   = tstacki[sequence[i]][sequence[j]] [sequence[i+1]][sequence[j-1]];
					ip_jp_energy = tstacki[sequence[jp]][sequence[ip]] [sequence[jp+1]][sequence[ip-1]];
				}
				energy = i_j_energy + ip_jp_energy + penalty_size + asym_penalty;
			}
		}
	}
	//printf ("REAL:   i=%d, j=%d, ip=%d, jp=%d, energy=%g\n", i, j, ip, jp, energy);
	return energy;
}

PARAMTYPE s_internal_loop::get_enthalpy (int i, int j, int ip, int jp, int *sequence){
	// returns the enthalpy of the internal loop closed at (i,j,ip,jp)
	PARAMTYPE energy = INF;   // just in case i,j,ip,jp are not closing an internal loop
	PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
	int branch1, branch2, l;

	branch1 = ip-i-1;
	branch2 = j-jp-1;

	printf("HERE\n");
	if (branch1 != 0 || branch2 != 0){
		// check if it is a bulge loop of size 1
		// check if it is int11 or int21 or int22
		if (branch1 == 1 && branch2 == 1 && !simple_internal_energy){     // it is int11
			// int11[i][j][i+1][j-1][ip][jp]
			energy = enthalpy_int11 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]];
		}else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy){
			// int21[i][j][i+1][j-1][ip][jp][jp+1]
			energy = enthalpy_int21 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]]
					[sequence[jp+1]];
		}else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy){
			// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
			energy = enthalpy_int21 [sequence[jp]][sequence[ip]]
				[sequence[j-1]][sequence[ip-1]]
				[sequence[j]][sequence[i]]
					[sequence[i+1]];
		}else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy){
			// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
			energy = enthalpy_int22 [sequence[i]][sequence[j]]
				[sequence[i+1]][sequence[j-1]]
				[sequence[ip]][sequence[jp]]
					[sequence[ip-1]][sequence[jp+1]];
		}else{
			// this case is not int11, int21, int22

			// check if it is a bulge
			if (branch1 == 0 || branch2 == 0){
				l = branch1+branch2;
				penalty_size = penalty_by_size_enthalpy (l, 'B');
				if (l == 1){
					// bulge of size 1
					// stack[i][j][i+1][j-1]
					energy = enthalpy_stack [sequence[i]][sequence[j]]
						[sequence[ip]][sequence[jp]] +
						penalty_size;
				}else{
					// bulge of size bigger than 1
					// check if (i,j) and (ip,jp) can pair
					energy = penalty_size +
						AU_penalty_enthalpy (sequence[i],sequence[j]) +
						AU_penalty_enthalpy (sequence[ip], sequence[jp]);
				}
			// it is an internal loop (not a bulge)
			}else{
				l = branch1+branch2;
				penalty_size = penalty_by_size_enthalpy (l, 'I');
				asym_penalty = asymmetry_penalty_enthalpy (branch1, branch2);

				if ((branch1 == 1 || branch2 == 1) && enthalpy_misc.gail_rule){
					// If gail_rule is set to 1 in miscloop file,
					// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
					i_j_energy  =  enthalpy_tstacki[sequence[i]][sequence[j]]
						[0][0];
					ip_jp_energy = enthalpy_tstacki[sequence[jp]][sequence[ip]]
						[0][0];
				}else{
					i_j_energy   = enthalpy_tstacki[sequence[i]][sequence[j]]
						[sequence[i+1]][sequence[j-1]];
					ip_jp_energy = enthalpy_tstacki[sequence[jp]][sequence[ip]]
						[sequence[jp+1]][sequence[ip-1]];
				}
				energy = i_j_energy + ip_jp_energy + penalty_size + asym_penalty;
			}
		}
	}
	return energy;
}


void count_int21_MODEL_EXTENDED (double *counter, int ii, int jj, int kk, int ll, int mm, int nn, int oo){
	// do the counts for int21, when MODEL is EXTENDED
#if (MODEL == EXTENDED)
	char type[100];
	int index;
	// exactly the same as above
	// try to follow the model proposed by Badhwar_Znosko_2007
	// the initiation appears in all cases
	sprintf (type, "misc.internal21_initiation");
	index = structure_type_index (type);
	counter[index]++;
	// look for AU closure
	if ((ii==A && jj==U) || (ii==U && jj==A))
	{
		sprintf (type, "misc.internal21_AU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	if ((mm==A && nn==U) || (mm==U && nn==A))
	{
		sprintf (type, "misc.internal21_AU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	// look for GU closure
	if ((ii==G && jj==U) || (ii==U && jj==G))
	{
		sprintf (type, "misc.internal21_GU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	if ((mm==G && nn==U) || (mm==U && nn==G))
	{
		sprintf (type, "misc.internal21_GU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	// look for AG mismatch - but not applied to 5'RA/3'YG loops
	if ((kk==A && ll==G &&   ii!=A && ii!=G   &&   jj!=U && jj!=C) ||
			(kk==G && ll==A) ||
			(kk==G && oo==A &&   mm!=U && mm!=C   &&   nn!=A && nn!=G) ||
			(kk==A && oo==G))
	{
		sprintf (type, "misc.internal21_AG_mismatch");
		index = structure_type_index (type);
		counter[index]++;
	}        
	// look for GG mismatch
	if (kk==G && (ll==G || oo==G))
	{
		sprintf (type, "misc.internal21_GG_mismatch");
		index = structure_type_index (type);
		counter[index]++;
	}
	// look for UU mismatch
	if (kk==U && (ll==U || oo==U))
	{
		sprintf (type, "misc.internal21_UU_mismatch");
		index = structure_type_index (type);
		counter[index]++;
	}
	// IN THIS CASE, int21 is just the additional value, on top of the above values, in order to make them be what the experiments say, and not an approximation
	// HMM - MAYBE THIS IS NOT SUCH A GOOD IDEA ACTUALLY, BUT THAT'S WHAT'S RECOMMENDED BY THE OPTICAL MELTING PAPERS
	// TODO: separate into experimental_addition and not

	// the following is counted only if it's part of experimental_addition, otherwise it's not
	//if (int21_experimental_addition[ii][jj][kk][ll][mm][nn][oo] < INF)
	{
		sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
		//sprintf (type, "int21_experimental_addition[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
		index = structure_type_index (type);
		counter[index]++;
	}
#endif
}


void count_int11_MODEL_EXTENDED (double *counter, int ii, int jj, int kk, int ll, int mm, int nn)
	// do the counts for int21, when MODEL is EXTENDED
{
#if (MODEL == EXTENDED)
	char type[100];
	int index;
	// try to follow the model proposed by Davis_Znosko_2007
	// look for AU closure
	if ((ii==A && jj==U) || (ii==U && jj==A))
	{
		sprintf (type, "misc.internal11_AU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	if ((mm==A && nn==U) || (mm==U && nn==A))
	{
		sprintf (type, "misc.internal11_AU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	// look for GU closure
	if ((ii==G && jj==U) || (ii==U && jj==G))
	{
		sprintf (type, "misc.internal11_GU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	if ((mm==G && nn==U) || (mm==U && nn==G))
	{
		sprintf (type, "misc.internal11_GU_closure");
		index = structure_type_index (type);
		counter[index]++;
	}
	// look for AG mismatch
	if ((kk==A && ll==G) || (kk==G && ll==A))
	{
		sprintf (type, "misc.internal11_AG_mismatch");
		index = structure_type_index (type);
		counter[index]++;
	}    
	// look for GG mismatch
	if (kk==G && ll==G)
	{
		sprintf (type, "misc.internal11_GG_mismatch");
		index = structure_type_index (type);
		counter[index]++;
	}
	// look for UU mismatch
	if (kk==U && ll==U)
	{
		sprintf (type, "misc.internal11_UU_mismatch");
		index = structure_type_index (type);
		counter[index]++;
	}
	// check if it is internal11_5YRR_5YRR

	if (isY(ii) && isR(jj) && isR(kk) && isR(ll) && isR(mm) && isY(nn))
	{
		sprintf (type, "misc.internal11_5YRR_5YRR");
		index = structure_type_index (type);
		counter[index]++;
	}    
	if ( isR(ii) && isY(jj) && isY(kk) && isY(ll) && isY(mm) && isR(nn) )
	{
		sprintf (type, "misc.internal11_5RYY_5RYY");
		index = structure_type_index (type);
		counter[index]++;
	}    
	if ( isY(ii) && isR(jj) && isY(kk) && isY(ll) && isR(mm) && isY(nn) )
	{
		sprintf (type, "misc.internal11_5YYR_5YYR");
		index = structure_type_index (type);
		counter[index]++;
	}    
	if ( isY(ii) && isR(jj) && isR(kk) && isY(ll) && isY(mm) && isR(nn) )
	{
		sprintf (type, "misc.internal11_5YRY_5RYR");
		index = structure_type_index (type);
		counter[index]++;
	}    
	if ( isR(ii) && isY(jj) && isR(kk) && isY(ll) && isY(mm) && isR(nn) )
	{
		sprintf (type, "misc.internal11_5RRY_5RYY");
		index = structure_type_index (type);
		counter[index]++;
	}    

	// IN THIS CASE, int21 is just the additional value, on top of the above values, in order to make them be what the experiments say, and not an approximation
	// HMM - MAYBE THIS IS NOT SUCH A GOOD IDEA ACTUALLY, BUT THAT'S WHAT'S RECOMMENDED BY THE OPTICAL MELTING PAPERS
	// TODO: separate into experimental_addition and not

	// the following is counted only if it's part of experimental_addition, otherwise it's not
	//if (int11_experimental_addition[ii][jj][kk][ll][mm][nn] < INF)
	{
		//sprintf (type, "int11_experimental_addition[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
		sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
		index = structure_type_index (type);
		counter[index]++;
	}
#endif
}

void s_internal_loop::count_get_energy (int i, int j, int ip, int jp, int *sequence, double *counter){
	// this function is needed for parameter learning, not for folding
	// fill the counter vectro accordingly
	// Mirela: Nov 23, 2003
	PARAMTYPE energy;
	PARAMTYPE penalty_size, asym_penalty, ip_jp_energy, i_j_energy;
	int branch1, branch2, l;
	char type[100];
	int index;
	branch1 = ip-i-1;
	branch2 = j-jp-1;    

	if (branch1 != 0 || branch2 != 0)
	{
		// check if it is a bulge loop of size 1
		// check if it is int11 or int21 or int22
		if (branch1 == 1 && branch2 == 1 && !simple_internal_energy)     // it is int11
		{
			// int11[i][j][i+1][j-1][ip][jp]
			energy = IGINF(int11 [sequence[i]][sequence[j]]
					[sequence[i+1]][sequence[j-1]]
					[sequence[ip]][sequence[jp]]);
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
#if (MODEL == SIMPLE)
				if ( ((ii==C && jj==G) || (ii==G && jj==C)) && ((mm==C && nn==G) || (mm==G && nn==C)))
				{
					if (!can_pair(kk,ll))
						sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
					else
						sprintf (type, "misc.internal11_basic_mismatch");
					index = structure_type_index (type);
					counter[index]++;                        
				}        
				else if (watson_crick(ii,jj) && watson_crick(mm,nn) && kk==U && ll==U)
				{
					sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
					index = structure_type_index (type);
					counter[index]++;                                                    
				}
				else
				{
					if (kk==G && ll==G)
						sprintf (type, "misc.internal11_GG_mismatch");
					else
						sprintf (type, "misc.internal11_basic_mismatch");
					index = structure_type_index (type);
					counter[index]++;                        

					if (has_AU_penalty(ii,jj))
					{
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                                                    
					}
					if (has_AU_penalty(mm,nn))
					{
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                        
					}
				}
#endif
				count_int11_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn);
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

#if (MODEL == SIMPLE)
				if ( ((ii==C && jj==G) || (ii==G && jj==C)) && ((mm==C && nn==G) || (mm==G && nn==C))) 
				{
					if (!can_pair(kk,ll))
						sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
					else
						sprintf (type, "misc.internal11_basic_mismatch");
					index = structure_type_index (type);
					counter[index]++;                        
				}        
				else if (watson_crick(ii,jj) && watson_crick(mm,nn) && kk==U && ll==U)
				{
					sprintf (type, "int11[%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn);
					index = structure_type_index (type);
					counter[index]++;                                                    
				}
				else
				{
					if (kk==G && ll==G)
						sprintf (type, "misc.internal11_GG_mismatch");
					else
						sprintf (type, "misc.internal11_basic_mismatch");
					index = structure_type_index (type);
					counter[index]++;                        

					if (has_AU_penalty(ii,jj))
					{
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                                                    
					}
					if (has_AU_penalty(mm,nn))
					{
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                        
					}
				}
#endif
				count_int11_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn);
			}
		}
		else if (branch1 == 1 && branch2 == 2 && !simple_internal_energy)
		{
			// int21[i][j][i+1][j-1][ip][jp][jp+1]
			energy = IGINF(int21 [sequence[i]][sequence[j]]
					[sequence[i+1]][sequence[j-1]]
					[sequence[ip]][sequence[jp]]
					[sequence[jp+1]]);

			int ii, jj, kk, ll, mm, nn, oo;
			ii = sequence[i];
			jj = sequence[j];
			kk = sequence[i+1];
			ll = sequence[j-1];
			mm = sequence[ip];
			nn = sequence[jp];
			oo = sequence[jp+1];            

#if (MODEL == SIMPLE)
			if ((ii==C && jj==G && mm==C && nn==G) ||  // these are already filled above, except what can pair inside
					(ii==G && jj==C && mm==G && nn==C))
			{
				if (can_pair(kk,ll) || can_pair(kk,oo))
					sprintf (type, "misc.internal21_match");
				else
					sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
				index = structure_type_index (type);
				counter[index]++;
			}
			else
			{
				if (can_pair(kk,ll) || can_pair(kk,oo))
				{
					sprintf (type, "misc.internal21_match");
					index = structure_type_index (type);
					counter[index] ++;                
				}
				else
				{
					sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kk, ll, C, G, oo);
					index = structure_type_index (type);
					counter[index] += 0.5;                

					sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kk, ll, G, C, oo);
					index = structure_type_index (type);
					counter[index] += 0.5;  
				}        
				if (has_AU_penalty(ii,jj))
				{
					sprintf (type, "misc.internal21_AU_closure");
					index = structure_type_index (type);
					counter[index]++;                        
				}    
				if (has_AU_penalty(mm,nn))    
				{
					sprintf (type, "misc.internal21_AU_closure");
					index = structure_type_index (type);
					counter[index]++;                        
				}    
			}
#endif
			count_int21_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn, oo);
		}
		else if(branch1 == 2 && branch2 == 1 && !simple_internal_energy)
		{
			// after rotation: int21[jp][ip][j-1][ip-1][j][i][i+1]
			energy = IGINF(int21 [sequence[jp]][sequence[ip]]
					[sequence[j-1]][sequence[ip-1]]
					[sequence[j]][sequence[i]]
					[sequence[i+1]]);
			int ii, jj, kk, ll, mm, nn, oo;
			ii = sequence[jp];
			jj = sequence[ip];
			kk = sequence[j-1];
			ll = sequence[ip-1];
			mm = sequence[j];
			nn = sequence[i];
			oo = sequence[i+1];
#if (MODEL == SIMPLE)
			if ((ii==C && jj==G && mm==C && nn==G) ||  // these are already filled above, except what can pair inside
					(ii==G && jj==C && mm==G && nn==C))
			{
				if (can_pair(kk,ll) || can_pair(kk,oo))
					sprintf (type, "misc.internal21_match");
				else
					sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo);
				index = structure_type_index (type);
				counter[index]++;
			}
			else
			{
				if (can_pair(kk,ll) || can_pair(kk,oo))
				{
					sprintf (type, "misc.internal21_match");
					index = structure_type_index (type);
					counter[index]++;
				}
				else
				{
					sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", C, G, kk, ll, C, G, oo);
					index = structure_type_index (type);
					counter[index] += 0.5;                

					sprintf (type, "int21[%d][%d][%d][%d][%d][%d][%d]", G, C, kk, ll, G, C, oo);
					index = structure_type_index (type);
					counter[index] += 0.5;  
				}        
				if (has_AU_penalty(ii,jj))
				{
					sprintf (type, "misc.internal21_AU_closure");
					index = structure_type_index (type);
					counter[index]++;                        
				}    
				if (has_AU_penalty(mm,nn))    
				{
					sprintf (type, "misc.internal21_AU_closure");
					index = structure_type_index (type);
					counter[index]++;                        
				}    
			}
#endif
			count_int21_MODEL_EXTENDED (counter, ii, jj, kk, ll, mm, nn, oo);
		}
		else if (branch1 == 2 && branch2 == 2 && !simple_internal_energy)
		{
			// int22[i][j][i+1][j-1][ip][jp][ip-1][jp+1]
			energy = IGINF(int22 [sequence[i]][sequence[j]]
					[sequence[i+1]][sequence[j-1]]
					[sequence[ip]][sequence[jp]]
					[sequence[ip-1]][sequence[jp+1]]);
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

#if (MODEL == SIMPLE)
				if ( (nn==ii) && (mm==jj) && (pp==kk) && (oo==ll) & watson_crick(ii,jj) && !watson_crick(kk,ll))
				{
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
					index = structure_type_index (type);
					counter[index]++;
				}

				int ii2, jj2, mm2, nn2;                
				if (ii==G && jj==U)   ii2 = A;     else ii2 = ii;
				if (ii==U && jj==G)   jj2 = A;     else jj2 = jj;
				if (mm==G && nn==U)   mm2 = A;     else mm2 = mm;
				if (mm==U && nn==G)   nn2 = A;     else nn2 = nn;

				if (watson_crick(kk,ll) || watson_crick(oo,pp))
				{
					sprintf (type, "misc.internal22_match");
					index = structure_type_index (type);
					counter[index]++;
				}
				else if ( ((ii==G && jj==U) || (ii==U && jj==G) || (mm==G && nn==U) || (mm==U && nn==G)) &&  
						(nn2==ii2 && mm2==jj2 && pp==kk && oo==ll))  // the UG closing pairs are the same as UA
				{
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii2, jj2, kk, ll, mm2, nn2, oo, pp);
					index = structure_type_index (type);
					counter[index]++;
				}                
				else if (!(nn==ii && mm==jj && pp==kk && oo==ll))   // was already filled above
				{
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii2, jj2, kk, ll, jj2, ii2, ll, kk);
					index = structure_type_index (type);
					counter[index] += 0.5;
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nn2, mm2, pp, oo, mm2, nn2, oo, pp);
					index = structure_type_index (type);
					counter[index] += 0.5;

					int result = check_stability_and_size (kk, ll, oo, pp);                    
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
							printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kk, ll, oo, pp);
							exit(1);
					}
					index = structure_type_index (type);
					counter[index]++;
				}
#elif (MODEL == EXTENDED)
				// in the EXTENDED model, all the int22 features are listed
				sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
				index = structure_type_index (type);
				counter[index]++;     
#endif                                                  
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

#if (MODEL == SIMPLE)
				if ( (nn==ii) && (mm==jj) && (pp==kk) && (oo==ll) & watson_crick(ii,jj) && !watson_crick(kk,ll))
				{
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
					index = structure_type_index (type);
					counter[index]++;
				}

				int ii2, jj2, mm2, nn2;                
				if (ii==G && jj==U)   ii2 = A;     else ii2 = ii;
				if (ii==U && jj==G)   jj2 = A;     else jj2 = jj;
				if (mm==G && nn==U)   mm2 = A;     else mm2 = mm;
				if (mm==U && nn==G)   nn2 = A;     else nn2 = nn;

				if (watson_crick(kk,ll) || watson_crick(oo,pp))
				{
					sprintf (type, "misc.internal22_match");
					index = structure_type_index (type);
					counter[index]++;
				}
				else if ( ((ii==G && jj==U) || (ii==U && jj==G) || (mm==G && nn==U) || (mm==U && nn==G)) &&  
						(nn2==ii2 && mm2==jj2 && pp==kk && oo==ll))  // the UG closing pairs are the same as UA
				{
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii2, jj2, kk, ll, mm2, nn2, oo, pp);
					index = structure_type_index (type);
					counter[index]++;
				}
				else if (!(nn==ii && mm==jj && pp==kk && oo==ll))   // was already filled above
				{
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii2, jj2, kk, ll, jj2, ii2, ll, kk);
					index = structure_type_index (type);
					counter[index] += 0.5;
					sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", nn2, mm2, pp, oo, mm2, nn2, oo, pp);
					index = structure_type_index (type);
					counter[index] += 0.5;

					int result = check_stability_and_size (kk, ll, oo, pp);                    
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
							printf ("ERROR: result %d for k=%d, l=%d, o=%d, p=%d, ABORT!\n", result, kk, ll, oo, pp);
							exit(1);
					}
					index = structure_type_index (type);
					counter[index]++;
				}
#elif (MODEL == EXTENDED)
				// in the EXTENDED model, all the int22 features are listed
				sprintf (type, "int22[%d][%d][%d][%d][%d][%d][%d][%d]", ii, jj, kk, ll, mm, nn, oo, pp);
				index = structure_type_index (type);
				counter[index]++;  
#endif
			}

		}       // end if 2x2
		else
		{
			// this case is not int11, int21, int22                                                                   
			// check if it is a bulge
			if (branch1 == 0 || branch2 == 0)
			{
				l = branch1+branch2;                            
				if (l == 1)
				{
#if (MODEL == SIMPLE)
					// bulge of size 1
					// stack[i][j][i+1][j-1]
					penalty_size = penalty_by_size (l, 'B');    
					count_penalty_by_size (l, 'B', counter);
					energy = IGINF(stack [sequence[i]][sequence[j]]
							[sequence[ip]][sequence[jp]]) +
						penalty_size;
					if (1000*sequence[i] + 100*sequence[j] + 10*sequence[ip] + sequence[jp] >
							1000*sequence[jp] + 100*sequence[ip] + 10*sequence[j] + sequence[i])
						sprintf (type, "stack[%d][%d][%d][%d]", sequence[jp], sequence[ip], sequence[j], sequence[i]);
					else
						sprintf (type, "stack[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[ip], sequence[jp]);
					index = structure_type_index (type);
					counter[index]++;
#elif (MODEL == EXTENDED)
					int i2, j2, k2, ip2, jp2;       //  bulge1[i2][j2][k2][ip2][jp2], the bulged nucleotide on top
					if (branch1 == 1)
					{
						i2 = sequence[i];
						j2 = sequence[j];
						k2 = sequence[i+1];
						ip2 = sequence[ip];
						jp2 = sequence[jp];
					}
					else        // it's upside down
					{
						i2 = sequence[jp];
						j2 = sequence[ip];
						k2 = sequence[j-1];
						ip2 = sequence[j];
						jp2 = sequence[i];
					}
					sprintf (type, "bulge1[%d][%d][%d][%d][%d]", i2, j2, k2, ip2, jp2);
					index = structure_type_index (type);
					counter[index]++;                    
#endif
				}
				else
				{
					// bulge of size bigger than 1
					// check if (i,j) and (ip,jp) can pair
					penalty_size = penalty_by_size (l, 'B');    
					count_penalty_by_size (l, 'B', counter);
					energy = penalty_size +
						AU_penalty (sequence[i],sequence[j]) +
						AU_penalty (sequence[ip], sequence[jp]);
					count_AU_penalty (sequence[i],sequence[j], counter);
					count_AU_penalty (sequence[ip], sequence[jp], counter);
				}
			}
			// it is an internal loop (not a bulge)
			else
			{
				l = branch1+branch2;
#if (MODEL == SIMPLE)
				penalty_size = penalty_by_size (l, 'I');
				count_penalty_by_size (l, 'I', counter);
				// TODO: to add asym penalty
				asym_penalty = asymmetry_penalty (branch1, branch2);
				count_asymmetry_penalty (branch1, branch2, counter);
#elif (MODEL == EXTENDED)
				penalty_size = IL_penalty_by_size_2D (branch1, branch2);
				IL_count_penalty_by_size_2D (branch1, branch2, counter);
#endif

				if ((branch1 == 1 || branch2 == 1) && misc.gail_rule)
					// If gail_rule is set to 1 in miscloop file,
					// i_j_energy and ip_jp_energy will be calculated as if it was a loop of As
				{
					i_j_energy  =  IGINF(tstacki[sequence[i]][sequence[j]]
							[0][0]);

					// actually, just use 3 parameters instead of the tstacki table
					if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
							((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
					{
						//internal_AU_closure includes terminal_AU_penalty
						//sprintf (type, "misc.terminal_AU_penalty");
						//index = structure_type_index (type);
						//counter[index]++;                        
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                                                
					}
					/*  // the full tstacki table
						sprintf (type, "tstacki[%d][%d][0][0]", sequence[i], sequence[j]);
						index = structure_type_index (type);
						counter[index]++;
						*/
					ip_jp_energy = IGINF(tstacki[sequence[jp]][sequence[ip]]
							[0][0]);
					// actually, just use 3 parameters instead of the tstacki table
					if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
							((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
					{
						//internal_AU_closure includes terminal_AU_penalty
						//sprintf (type, "misc.terminal_AU_penalty");
						//index = structure_type_index (type);
						//counter[index]++;                        
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                                                
					}
					/*  // the full tstacki table                      
						sprintf (type, "tstacki[%d][%d][0][0]", sequence[jp], sequence[ip]);
						index = structure_type_index (type);
						counter[index]++;
						*/
				}
				else
				{
					i_j_energy   = IGINF(tstacki[sequence[i]][sequence[j]]
							[sequence[i+1]][sequence[j-1]]);

#if (MODEL == SIMPLE)
					// actually, just use 3 parameters instead of the tstacki table
					if (((sequence[i] == A || sequence[i] == G) && sequence[j] == U) ||
							((sequence[j] == A || sequence[j] == G) && sequence[i] == U))
					{
						//internal_AU_closure includes terminal_AU_penalty
						//sprintf (type, "misc.terminal_AU_penalty");
						//index = structure_type_index (type);
						//counter[index]++;                        
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                                                
					}
					if ((sequence[i+1] == A && sequence[j-1] == G) ||
							(sequence[j-1] == A && sequence[i+1] == G))
					{
						sprintf (type, "misc.internal_AG_mismatch");
						index = structure_type_index (type);
						counter[index]++;                                            
					}
					if (sequence[i+1] == U && sequence[j-1] == U)
					{
						sprintf (type, "misc.internal_UU_mismatch");
						index = structure_type_index (type);
						counter[index]++;                                            
					}          
#elif (MODEL == EXTENDED)                                    

					// the full tstacki table                                                                
					sprintf (type, "tstacki[%d][%d][%d][%d]", sequence[i], sequence[j], sequence[i+1], sequence[j-1]);
					index = structure_type_index (type);
					counter[index]++;
#endif                    

					ip_jp_energy = IGINF(tstacki[sequence[jp]][sequence[ip]]
							[sequence[jp+1]][sequence[ip-1]]);

#if (MODEL == SIMPLE)
					// actually, just use 3 parameters instead of the tstacki table
					if (((sequence[ip] == A || sequence[ip] == G) && sequence[jp] == U) ||
							((sequence[jp] == A || sequence[jp] == G) && sequence[ip] == U))
					{
						//internal_AU_closure includes terminal_AU_penalty
						//sprintf (type, "misc.terminal_AU_penalty");
						//index = structure_type_index (type);
						//counter[index]++;                        
						sprintf (type, "misc.internal_AU_closure");
						index = structure_type_index (type);
						counter[index]++;                                                
					}
					if ((sequence[ip-1] == A && sequence[jp+1] == G) ||
							(sequence[jp+1] == A && sequence[ip-1] == G))
					{
						sprintf (type, "misc.internal_AG_mismatch");
						index = structure_type_index (type);
						counter[index]++;                                            
					}
					if (sequence[ip-1] == U && sequence[jp+1] == U)
					{
						sprintf (type, "misc.internal_UU_mismatch");
						index = structure_type_index (type);
						counter[index]++;                                            
					}   
#elif (MODEL == EXTENDED)                                                                   
					// the full tstacki table                                                                
					sprintf (type, "tstacki[%d][%d][%d][%d]", sequence[jp], sequence[ip], sequence[jp+1], sequence[ip-1]);
					index = structure_type_index (type);
					counter[index]++;
#endif
				}
				energy = i_j_energy + ip_jp_energy + penalty_size +
					asym_penalty;
			}
		}
	}
	//    return energy;
}


