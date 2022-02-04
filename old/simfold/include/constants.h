/***************************************************************************
                          constants.h  -  description
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

// define the model
#define MODEL SIMPLE
#define SIMPLE 0        // the model with 363 parameters that I used in the ISMB 2007 paper
#define EXTENDED 1


// define the parameters type
#ifdef DOUBLEPARAMS
#define PARAMTYPE double
#define PFTYPE double
#define EXP exp

#elif LDOUBLEPARAMS
#define PARAMTYPE long double
#define PFTYPE long double
#define EXP expl

#else
#define PARAMTYPE int
#define PFTYPE double
#define EXP exp
#endif


#ifdef DOUBLEPARAMS
#define INCLUDE_FUNCTIONS
#elif LDOUBLEPARAMS
#define INCLUDE_FUNCTIONS
#endif

// the type of the partition function arrays

#define PFFORMAT g

#define RNA   0
#define DNA   1

// change these if necessary
#define MAXSLEN         1000  // maximum of total sequence length
#define MAXENERGY       0     // the maximum energy that this program returns (usually 0)
#define MAX_BRANCHES    MAXSLEN/2    // maximum # of branches in multiloops   // should be no more than 50 or so, but if the parameters go nuts, I may have many branches
#define MAXSUBSTR       2000  // maximum number of suboptimal structures

// No need to change the following

#define MAXNUMPARAMS    15000  // max # of thermodynamic parameters
#define MAXPNAME        500    // max length of a parameter name

// for MultiFold
#define MAXNUMSEQ      50

#define NUM_DANG       48        // there are 48 dangling ends 

#define A               0
#define C               1
#define G               2
#define U               3
#define T               3

#define INF             1600000      // a very big value (infinity)
#define NUCL            4            // number of nucleotides: 4: A, C, G, T
#define MAXLOOP         30           // max size for internal loops

// max loop size for which tabulated values are experimentally determined 
//    for len > MAX, the following formula is applied: lenpen (len) = lenpen(MAX) + 1.75 RT*log(1.0*len/MAX));
//    RT = 1.75 * 1.98717 / 1000 * 310.15;

#if (MODEL == SIMPLE)
#define MAXLOOP_H     9    // hairpin
#elif (MODEL == EXTENDED)
#define MAXLOOP_H     30    // hairpin
#endif

#if (MODEL == SIMPLE)
#define MAXLOOP_I     6    // internal
#elif (MODEL == EXTENDED)
#define MAXLOOP_I     30    // internal
#endif

#if (MODEL == SIMPLE)
#define MAXLOOP_B     6    // bulge
#elif (MODEL == EXTENDED)
#define MAXLOOP_B     30    // bulge
#endif

#define MAXTRILOOPNO    100          // max number of hairpin loops of size 3
#define MAXTLOOPNO      400          // max number of hairpin loops of size 4
#define MAX_SPECIAL_LOOP_NO 400      // max number of special hairpin loops for MODEL EXTENDED
#define EPSILON         0.0001       // a very small value
#define TURN            3            // the minimum number of free bases in a hairpin loop
// for partition function for restricted sequences, TURN should be 0, otherwise I can't get internal loops etc next to ()
// However, this screws up the partition function and gradient for the unrestricted case
#define MIN_HAIRPIN     3
// TURN used to be 3, but for the restricted cases it can be 0. So I changed it to 0.


#define NONE            'N'         // no structure
#define HAIRP           'H'         // closes a hairpin loop
#define STACK           'S'         // closes a stacked pair
#define INTER           'I'         // closes an internal loop
#define MULTI           'M'         // closes a regular multi-loop
#define MULTI_LINK      'O'         // closes a special multi-loop

#define M_WM            'B'         // closes a regular partial multi-loop
#define M_WM_LINK       'C'         // closes a special partial multi-loop

#define FREE            'W'         // this base is free to be paired or unpaired to other base
#define LOOP            'V'         // closes a loop

#define M_FM            'F'         // regular FM from the multiloop decomposition of Wutchy et al 
#define M_FM1           'A'         // regular FM1 from the multiloop decomposition of Wutchy et al

#define M_FM_LINK       'P'         // special FM adapted from FM = the multiloop decomposition of Wutchy et al
#define M_FM1_LINK      'Q'         // special FM1 adapted from FM1 = the multiloop decomposition of Wutchy et al


#endif

