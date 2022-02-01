/*****************************************************************
         HotKnot: A heuristic algorithm for RNA secondary 
            structure prediction including pseudoknots 
    Date        : Oct. 16, 2004
    copyright   : (C) 2004 by Jihong Ren, Baharak Rastegari  
    email       : jihong@cs.ubc.ca, baharak@cs.ubc.ca       
******************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include "goodStem.h"
#include "simfold.h"
#include "hotspot.h"
#include "utils.h"
#include "initPK.h"  // for init_dataPK()
#include "params.h"  // for init_dataPK()
#include "init.h"  // for init_dataPK()
#include "constantsPK.h"
//#include "externs.h"

struct Hotspot* hotspots;
int noPS=0;
int noGU = 0;
int TRACE=0;    // Mirela: to change back to 0
int *sp;
float THR=400.0;  // Mirela: was 400.0. //threshold for sub-sequence
float T_RATIO = 0.8; //Mirela: was 0.8
/*Consider only subsequent 2ndary structures that have energy lower than 
  T_RATIO*energy of the best non-pseudoknotted structure */
int MaxSubOpt = 20;
struct Node* listOfNodes[50];   // for Rivas&Eddy
struct Node* listOfNodes2[50];  // for Dirks&Pierce
struct Node* listOfNodes3[50];  // for Cao&Chen (a)
struct Node* listOfNodes4[50];  // for Cao&Chen (b)
struct Node* listOfNodes5[50];  // for Cao&Chen (c)
// *** Add new energy model code here (add another listOfNodes)

int count;//number of nodes
int numRnaStruct;  //total number of different Rna structures

extern void PlotRna(char* seqName, char *sequence, short *structure, char *filename, float score);


#define PRIVATE static

static char  scale1[] = "....,....1....,....2....,....3....,....4";
static char  scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);

double hotknots_DP_name(char *seq, char *structure, char* fileName)
// fileName: input - does not include path, just the name of the file without extension (sequence name)
// sequence: input
// structure: output
// By Mirela: THE NAME iS NOT USED ANY MORE, I wrote a function called bpseq2dp in utils.c
//            TO CLEAN UP THE MESS IN THIS FILE.
{
    if (DEBUG)
        printf("Input parameter seq = %s, strlen = %d\n", seq, (int) strlen(seq));

    int length = strlen(seq);
    char* sequence = new char[length+1];
    for (int i = 0; i < length; i++)
        sequence[i] = toupper(seq[i]);
    sequence[length] = '\0';

    if (DEBUG)
        printf("\nInput to hotknots_energy: len=%d, seq = %s, fileName = %s\n", length, sequence, fileName);

    char *string;
    char line[5000];
    char *cstructure=NULL, *cstruc=NULL;  // constraint structure
    char  fname[13], ffname[20], gfname[20];

    char  *ns_bases=NULL, *c;
    int   i, l, sym, r;
    double energy, min_en;
    double kT, sfact=1.07;
    int   pf=0, istty;
    int noconv=0;
    int InputFile = 1;
    int StrCon = 0;    
    int endFlag;
    char seqName[1000];  // includes path, used for PlotRNA
    char inFile[1000];
//  char ctFile[100], outPSFile[100], bpseqFile[100];
    char ctFile2[1000], outPSFile2[1000], bpseqFile2[1000];
//  char ctFile3[100], outPSFile3[100], bpseqFile3[100];
//  char ctFile4[100], outPSFile4[100], bpseqFile4[100];
//  char ctFile5[100], outPSFile5[100], bpseqFile5[100];
//  // *** Add new energy model code here

//  struct Node *rootNode;
    struct Node *rootNode2;
//      struct Node *rootNode3;
//      struct Node *rootNode4;
//      struct Node *rootNode5;
//  // *** Add new energy model code here

    int len = 0;
    int index = -1;  
    char separator;
//    char fileName[200];  // filename of input file sequence, without the path
    
    // *** Add new energy model code here (also change output path as desired)
    
    //char outpath[30] = "../HotKnots/bin/TestSeq/";  // starting part of path for output files assuming directory of the parameter tuning files are in ../tools
//    char outpath1[20] = "RivasEddy/";  // directory to put files generated using Rivas&Eddy (preceeded by outpath)
    //char outpath2[20] = "DirksPierce/";  // see above
    //char outpath3[20] = "CaoChen_a/";  // see above
    //char outpath4[20] = "CaoChen_b/";  // see above
    //char outpath5[20] = "CaoChen_c/";  // see above

    // Mirela: added the following two lines: I want the path to be only what fileName tells me to be
    // Mirela: Actually I removed those variables altogether
    //outpath[0] = '\0';
    //outpath2[0] = '\0';
    //sprintf(seqName, "%s%s%s", outpath, outpath2, fileName);
    sprintf(seqName, "%s", fileName);

    char tempFile[10]="hello";
    FILE* input_file;
    int MaxHotspots = 200;

//  string=NULL;
//      length = (int) strlen(sequence);

    cstructure = (char *) malloc((length+1)*sizeof(char));

    for (i = 0; i < length; i++) 
        cstructure[i] = '.';
    cstructure[length] = 0;

    int first = 1;  // don't print suboptimals
    int noPS = 1;   // don't do the PS file
    int noGU = 0;   // include GU base pairs

//  char bpseqFile_temp[] = "RNASeq_temp.bpseq";

//        cout << "after init, COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;


    /////////////////////////// RIVAS and EDDY ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode=(struct Node *)malloc(sizeof(struct Node));
    rootNode->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode->constraint[length] = 0;
    rootNode->numChild = 0;
    rootNode->length = length;
    rootNode->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode->secStructure[i]=0;
        rootNode->fixedPairs[i]=0;
        //rootNode->constraint[i]='.';
    }
    strcpy(rootNode->constraint, cstructure);
    //===end of initialization    
 
    if (DEBUGH)
    {
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", string);    
    }
    
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);
    GenerateStemList(length, string, cstructure);
    
    endFlag=secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, RE);  // ADDED: for RE energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;

    if (DEBUGH)
    {
    printf("In total, %d nodes created. \n", count);
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    }

    sprintf(ctFile, "%s%s%s%s", outpath, outpath1, fileName, "_RE.ct");  // ADDED: for RE energy model
    FILE *cfile = fopen(ctFile, "w");
    if (cfile == NULL) {
        printf("can't open %s \n", ctFile);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.ps");  // ADDED: for RE energy model
            sprintf(bpseqFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.bpseq");  // ADDED: for RE energy model
        
        if (DEBUGH)
        {
            printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
            printRnaStruct(listOfNodes[i]->secStructure,length);
            printf("\n");
        }

        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile, listOfNodes[i]->score);
            FILE *bfile = fopen(bpseqFile, "w");
            if (bfile == NULL) {
                printf("can't open %s \n", bpseqFile);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);            
            }          
        }
        fprintf(cfile, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/

    /////////////////////////// DIRKS and PIERCE ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode2=(struct Node *)malloc(sizeof(struct Node));
    rootNode2->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode2->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode2->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode2->constraint[length] = 0;
    rootNode2->numChild = 0;
    rootNode2->length = length;
    rootNode2->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode2->secStructure[i]=0;
        rootNode2->fixedPairs[i]=0;
        //rootNode2->constraint[i]='.';
    }
    strcpy(rootNode2->constraint, cstructure);
    //===end of initialization    

    if (DEBUGH)
    {
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", sequence); 
    }
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section

//  printf("Input to Generate Stem List: len = %d, seq = %s, scruct = %s", length, sequence, cstructure);

    GenerateStemList(length, sequence, cstructure);  // depends only on hairpin, stacked, internal loop parameters

//  printf("After generate stem list\n");
//              cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//      cout << &dangle_top[2][1][0] << endl;

    endFlag=secondaryStruct(sequence,length,rootNode2,rootNode2,MaxHotspots, DP);  // ADDED: for DP energy model

//  printf("After secondaryStruct\n");

    ClearHotspots(MaxHotspots);

//                cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//      cout << &dangle_top[2][1][0] << endl;

//  printf("After ClearHotspots\n");

    min_en=endFlag;

    if (DEBUGH)
    {
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    }

    //sprintf(ctFile2, "%s%s%s%s", outpath, outpath2, fileName, "_DP.ct");
    // Mirela: don't write the ct file
/*    sprintf(ctFile2, "%s%s", fileName, "_DP.ct");
    FILE *cfile2 = fopen(ctFile2, "w");
    if (cfile2 == NULL) {
        printf("can't open %s \n", ctFile2);
        exit(1);
    }*/
    

    // Mirela: commented out this for now
    /*
    for (i=0; i < numRnaStruct; i++) {
        //sprintf(outPSFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i,"_DP.ps");  // ADDED: for DP energy model
        sprintf(outPSFile2, "%s%d%s", fileName, i,"_DP.ps");  // ADDED: for DP energy model
        //sprintf(bpseqFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i, "_DP.bpseq");  // ADDED: for DP energy model
        sprintf(bpseqFile2, "%s%d%s", fileName, i, "_DP.bpseq");  // ADDED: for DP energy model

        if (DEBUGH)
        {
            printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
            printRnaStruct(listOfNodes[i]->secStructure,length);
            printf("\n");
        }

        if ((first == 1 && i==0) || (first == 0))
        {
        
            // Mirela: that's where the bpseq file is written
            //printf ("--- INSIDE if first = 1... i=%d\n", i);
            if (noPS == 0) PlotRna(seqName, sequence, listOfNodes[i]->secStructure, outPSFile2, listOfNodes[i]->score);
            FILE *bfile2 = fopen(bpseqFile2, "w");

        // Used later:
//      FILE *bfile_temp = fopen(bpseqFile_temp, "w");  // erase previous contents

            if (bfile2 == NULL) {
                    printf("can't open %s \n", bpseqFile2);
                    exit(1);
            }

//          if (bfile_temp == NULL) {
//              printf("can't open %s \n", bpseqFile_temp);
//                  exit(1);
//          }

            for (int b = 0; b < length; b++) {
                    fprintf(bfile2, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
//              fprintf(bfile_temp, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
            }          
//          fclose(bfile_temp);
            fclose(bfile2);
        }

//      fprintf(cfile2, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
//          for (int b = 0; b < length; b++) {
//              fprintf(cfile2, "%5d %c    %4d %4d %4d %4d\n", b+1, sequence[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
//          }
            
        }
    //fclose(cfile2);
 */    

    /////////////////////////// CAO and CHEN (a) ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode3=(struct Node *)malloc(sizeof(struct Node));
    rootNode3->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode3->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode3->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode3->constraint[length] = 0;
    rootNode3->numChild = 0;
    rootNode3->length = length;
    rootNode3->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode3->secStructure[i]=0;
        rootNode3->fixedPairs[i]=0;
        //rootNode3->constraint[i]='.';
    }
    strcpy(rootNode3->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode3,rootNode3, MaxHotspots, CC2006a);  // ADDED: for CC2006 energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile3, "%s%s%s%s", outpath, outpath3, fileName, "_CCa.ct");  // ADDED: for CCa energy model (filename)
    FILE *cfile3 = fopen(ctFile3, "w");
    if (cfile3 == NULL) {
        printf("can't open %s \n", ctFile3);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i,"_CCa.ps");  // ADDED: for CCa energy model (filename)
            sprintf(bpseqFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i, "_CCa.bpseq");  // ADDED: for CCa energy model (filename)
        printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
        printRnaStruct(listOfNodes[i]->secStructure,length);
        printf("\n");
        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile3, listOfNodes[i]->score);
            FILE *bfile3 = fopen(bpseqFile3, "w");
            if (bfile3 == NULL) {
                printf("can't open %s \n", bpseqFile3);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile3, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
            }
        }
        fprintf(cfile3, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile3, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }

    /////////////////////////// CAO and CHEN (b) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode4=(struct Node *)malloc(sizeof(struct Node));
    rootNode4->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode4->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode4->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode4->constraint[length] = 0;
    rootNode4->numChild = 0;
    rootNode4->length = length;
    rootNode4->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode4->secStructure[i]=0;
        rootNode4->fixedPairs[i]=0;
        //rootNode4->constraint[i]='.';
    }
    strcpy(rootNode4->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode4,rootNode4, MaxHotspots, CC2006b);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile4, "%s%s%s%s", outpath, outpath4, fileName, "_CCb.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile4 = fopen(ctFile4, "w");
    if (cfile4 == NULL) {
        printf("can't open %s \n", ctFile4);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i,"_CCb.ps");  // ADDED: for CC2006b energy model (filename)
            sprintf(bpseqFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i, "_CCb.bpseq");  // ADDED: for CCb energy model (filename)
        printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
        printRnaStruct(listOfNodes[i]->secStructure,length);
        printf("\n");
        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile4, listOfNodes[i]->score);
            FILE *bfile4 = fopen(bpseqFile4, "w");
            if (bfile4 == NULL) {
                printf("can't open %s \n", bpseqFile4);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile4, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
            }
        }
        fprintf(cfile4, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile4, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }

    /////////////////////////// CAO and CHEN (c) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode5=(struct Node *)malloc(sizeof(struct Node));
    rootNode5->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode5->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode5->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode5->constraint[length] = 0;
    rootNode5->numChild = 0;
    rootNode5->length = length;
    rootNode5->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode5->secStructure[i]=0;
        rootNode5->fixedPairs[i]=0;
        //rootNode4->constraint[i]='.';
    }
    strcpy(rootNode5->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode5,rootNode5, MaxHotspots, CC2006c);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile5, "%s%s%s%s", outpath, outpath5, fileName, "_CCc.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile5 = fopen(ctFile5, "w");
    if (cfile5 == NULL) {
        printf("can't open %s \n", ctFile5);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i,"_CCc.ps");  // ADDED: for CC2006b energy model (filename)
            sprintf(bpseqFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i, "_CCc.bpseq");  // ADDED: for CCb energy model (filename)
        printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
        printRnaStruct(listOfNodes[i]->secStructure,length);
        printf("\n");
        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile5, listOfNodes[i]->score);
            FILE *bfile5 = fopen(bpseqFile5, "w");
            if (bfile5 == NULL) {
                printf("can't open %s \n", bpseqFile5);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile5, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
            }
        }
        fprintf(cfile5, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile5, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/
    /////////////// *** Add new energy model code here (similar to above section) /////////////


    //(void) fflush(stdout);
    //free(string);
    free(cstructure);

    // TODO
    // fill out structure with the predicted structure
    
    /*    
    //  Mirela: commented out this part
    char dpFile_temp[1000];
    sprintf(dpFile_temp, "%s%s", seqName, "0_DP.dp");

    // make sure the bpseq file is closed
 
    char bpseqFileTemp[1000];
    sprintf(bpseqFileTemp, "%s%s", seqName, "0_DP.bpseq");
    
//  FILE* ftemp = fopen (bpseqFileTemp, "r");
//  if (ftemp == NULL)
//  {
//      printf("Before system call, cannot open file %s\n", bpseqFileTemp);
//      exit(1);
//  }
//  else
//      fclose(ftemp);
  

    char command[200];
    // assumes bpseq2dp.pl and Analyser (needed for pl script) are in same folder as parameter tuning executables (ie. tools)
    sprintf(command, "/cs/local/bin/perl ./bpseq2dp.pl %s%s %s%s", seqName, "0_DP.bpseq", seqName, "0_DP.dp");

        if (DEBUG)
                printf("System call: %s\n", command);

    // Assumes bpseq2dp.pl is in the same folder as the executable Hotknots
    system(command);

    if (DEBUG)
        printf("Done system call to bpseq2dp.pl\n");

    FILE *dpfile_temp;
        char buffer[MaxN];
        char  v1[MaxN];
//  char dpFile_temp[] = "RNASeq_temp.dp";

    //printf ("FILENAME: %s\n", filename);
    if ((dpfile_temp = fopen (dpFile_temp, "r")) == NULL)
    {
        printf ("Cannot open file %s\n", dpFile_temp);
        exit(1);
    }

    fgets (buffer, sizeof(buffer), dpfile_temp);
    sscanf (buffer, "%s", v1);
    strcpy(structure, v1);
    fclose(dpfile_temp);

    if (DEBUGH)
        printf("Structure returned: %s\n", structure);

    // Mirela: added the following 3 lines to clean up the mess
    remove(dpFile_temp);
    remove(bpseqFileTemp);
    //remove (ctFile2);
*/
    bpseq2dp (length, listOfNodes[0]->secStructure, structure);
    if (DEBUGH)
        printf("Structure returned: %s\n", structure);

    return -listOfNodes[0]->score /1000.0;   //return best suboptimal energy
}


// Prototype: double hotknots_DP (char *sequence, char *structure)
//      - sequence is the given RNA sequence to be folded (input parameter)
//      - structure is the predicted secondary structure (output parameter)
//      - the function may return the energy
// Uses the Dirks&Pierce model
// Defaults: no PS is output, GU is included (if want to change this, set it manually in ? TODO )
double hotknots_DP(char *sequence, char *structure)
{       
        return hotknots_DP_name (sequence, structure, (char *) "Default");
}


// Prototype: hotknots_DP_suboptimals(char *seq, char **structures, double* energies)
double hotknots_DP_suboptimals(char *seq, char **structures, double* energies, int& num_suboptimals)
// fileName: input - does not include path, just the name of the file without extension (sequence name)
// sequence: input
// structure: output
// By Mirela: THE NAME iS NOT USED ANY MORE, I wrote a function called bpseq2dp in utils.c
//            TO CLEAN UP THE MESS IN THIS FILE.
{
    if (DEBUG)
        printf("Input parameter seq = %s, strlen = %d\n", seq, (int) strlen(seq));

    int length = strlen(seq);
    char* sequence = new char[length+1];
    for (int i = 0; i < length; i++)
        sequence[i] = toupper(seq[i]);
    sequence[length] = '\0';

    if (DEBUG)
        printf("\nInput to hotknots_energy: len=%d, seq = %s\n", length, sequence);

    char *string;
    char line[5000];
    char *cstructure=NULL, *cstruc=NULL;  // constraint structure
    char  fname[13], ffname[20], gfname[20];

    char  *ns_bases=NULL, *c;
    int   i, l, sym, r;
    double energy, min_en;
    double kT, sfact=1.07;
    int   pf=0, istty;
    int noconv=0;
    int InputFile = 1;
    int StrCon = 0;    
    int endFlag;
    char seqName[1000];  // includes path, used for PlotRNA
    char inFile[1000];
//  char ctFile[100], outPSFile[100], bpseqFile[100];
    char ctFile2[1000], outPSFile2[1000], bpseqFile2[1000];
//  char ctFile3[100], outPSFile3[100], bpseqFile3[100];
//  char ctFile4[100], outPSFile4[100], bpseqFile4[100];
//  char ctFile5[100], outPSFile5[100], bpseqFile5[100];
//  // *** Add new energy model code here

//  struct Node *rootNode;
    struct Node *rootNode2;
//      struct Node *rootNode3;
//      struct Node *rootNode4;
//      struct Node *rootNode5;
//  // *** Add new energy model code here

    int len = 0;
    int index = -1;  
    char separator;
//    char fileName[200];  // filename of input file sequence, without the path
    
    // *** Add new energy model code here (also change output path as desired)
    
    //char outpath[30] = "../HotKnots/bin/TestSeq/";  // starting part of path for output files assuming directory of the parameter tuning files are in ../tools
//    char outpath1[20] = "RivasEddy/";  // directory to put files generated using Rivas&Eddy (preceeded by outpath)
    //char outpath2[20] = "DirksPierce/";  // see above
    //char outpath3[20] = "CaoChen_a/";  // see above
    //char outpath4[20] = "CaoChen_b/";  // see above
    //char outpath5[20] = "CaoChen_c/";  // see above

    // Mirela: added the following two lines: I want the path to be only what fileName tells me to be
    // Mirela: Actually I removed those variables altogether
    //outpath[0] = '\0';
    //outpath2[0] = '\0';
    //sprintf(seqName, "%s%s%s", outpath, outpath2, fileName);
//    sprintf(seqName, "%s", fileName);
    sprintf(seqName, "%s", "default");

    char tempFile[10]="hello";
    FILE* input_file;
    int MaxHotspots = 200;

//  string=NULL;
//      length = (int) strlen(sequence);

    cstructure = (char *) malloc((length+1)*sizeof(char));

    for (i = 0; i < length; i++) 
        cstructure[i] = '.';
    cstructure[length] = 0;

    int first = 1;  // don't print suboptimals
    int noPS = 1;   // don't do the PS file
    int noGU = 0;   // include GU base pairs

//  char bpseqFile_temp[] = "RNASeq_temp.bpseq";

//        cout << "after init, COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;


    /////////////////////////// RIVAS and EDDY ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode=(struct Node *)malloc(sizeof(struct Node));
    rootNode->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode->constraint[length] = 0;
    rootNode->numChild = 0;
    rootNode->length = length;
    rootNode->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode->secStructure[i]=0;
        rootNode->fixedPairs[i]=0;
        //rootNode->constraint[i]='.';
    }
    strcpy(rootNode->constraint, cstructure);
    //===end of initialization    
 
    if (DEBUGH)
    {
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", string);    
    }
    
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);
    GenerateStemList(length, string, cstructure);
    
    endFlag=secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, RE);  // ADDED: for RE energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;

    if (DEBUGH)
    {
    printf("In total, %d nodes created. \n", count);
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    }

    sprintf(ctFile, "%s%s%s%s", outpath, outpath1, fileName, "_RE.ct");  // ADDED: for RE energy model
    FILE *cfile = fopen(ctFile, "w");
    if (cfile == NULL) {
        printf("can't open %s \n", ctFile);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.ps");  // ADDED: for RE energy model
            sprintf(bpseqFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.bpseq");  // ADDED: for RE energy model
        
        if (DEBUGH)
        {
            printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
            printRnaStruct(listOfNodes[i]->secStructure,length);
            printf("\n");
        }

        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile, listOfNodes[i]->score);
            FILE *bfile = fopen(bpseqFile, "w");
            if (bfile == NULL) {
                printf("can't open %s \n", bpseqFile);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);            
            }          
        }
        fprintf(cfile, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/

    /////////////////////////// DIRKS and PIERCE ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode2=(struct Node *)malloc(sizeof(struct Node));
    rootNode2->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode2->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode2->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode2->constraint[length] = 0;
    rootNode2->numChild = 0;
    rootNode2->length = length;
    rootNode2->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode2->secStructure[i]=0;
        rootNode2->fixedPairs[i]=0;
        //rootNode2->constraint[i]='.';
    }
    strcpy(rootNode2->constraint, cstructure);
    //===end of initialization    

    if (DEBUGH)
    {
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", sequence); 
    }
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section

//  printf("Input to Generate Stem List: len = %d, seq = %s, scruct = %s", length, sequence, cstructure);

    GenerateStemList(length, sequence, cstructure);  // depends only on hairpin, stacked, internal loop parameters

//  printf("After generate stem list\n");
//              cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//      cout << &dangle_top[2][1][0] << endl;

    endFlag=secondaryStruct(sequence,length,rootNode2,rootNode2,MaxHotspots, DP);  // ADDED: for DP energy model

//  printf("After secondaryStruct\n");

    ClearHotspots(MaxHotspots);

//                cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//      cout << &dangle_top[2][1][0] << endl;

//  printf("After ClearHotspots\n");

    min_en=endFlag;

    if (DEBUGH)
    {
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    }

    //sprintf(ctFile2, "%s%s%s%s", outpath, outpath2, fileName, "_DP.ct");
    // Mirela: don't write the ct file
/*    sprintf(ctFile2, "%s%s", fileName, "_DP.ct");
    FILE *cfile2 = fopen(ctFile2, "w");
    if (cfile2 == NULL) {
        printf("can't open %s \n", ctFile2);
        exit(1);
    }*/
    

    // Mirela: commented out this for now
    /*
    for (i=0; i < numRnaStruct; i++) {
        //sprintf(outPSFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i,"_DP.ps");  // ADDED: for DP energy model
        sprintf(outPSFile2, "%s%d%s", fileName, i,"_DP.ps");  // ADDED: for DP energy model
        //sprintf(bpseqFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i, "_DP.bpseq");  // ADDED: for DP energy model
        sprintf(bpseqFile2, "%s%d%s", fileName, i, "_DP.bpseq");  // ADDED: for DP energy model

        if (DEBUGH)
        {
            printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
            printRnaStruct(listOfNodes[i]->secStructure,length);
            printf("\n");
        }

        if ((first == 1 && i==0) || (first == 0))
        {
        
            // Mirela: that's where the bpseq file is written
            //printf ("--- INSIDE if first = 1... i=%d\n", i);
            if (noPS == 0) PlotRna(seqName, sequence, listOfNodes[i]->secStructure, outPSFile2, listOfNodes[i]->score);
            FILE *bfile2 = fopen(bpseqFile2, "w");

        // Used later:
//      FILE *bfile_temp = fopen(bpseqFile_temp, "w");  // erase previous contents

            if (bfile2 == NULL) {
                    printf("can't open %s \n", bpseqFile2);
                    exit(1);
            }

//          if (bfile_temp == NULL) {
//              printf("can't open %s \n", bpseqFile_temp);
//                  exit(1);
//          }

            for (int b = 0; b < length; b++) {
                    fprintf(bfile2, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
//              fprintf(bfile_temp, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
            }          
//          fclose(bfile_temp);
            fclose(bfile2);
        }

//      fprintf(cfile2, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
//          for (int b = 0; b < length; b++) {
//              fprintf(cfile2, "%5d %c    %4d %4d %4d %4d\n", b+1, sequence[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
//          }
            
        }
    //fclose(cfile2);
 */    

    /////////////////////////// CAO and CHEN (a) ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode3=(struct Node *)malloc(sizeof(struct Node));
    rootNode3->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode3->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode3->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode3->constraint[length] = 0;
    rootNode3->numChild = 0;
    rootNode3->length = length;
    rootNode3->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode3->secStructure[i]=0;
        rootNode3->fixedPairs[i]=0;
        //rootNode3->constraint[i]='.';
    }
    strcpy(rootNode3->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode3,rootNode3, MaxHotspots, CC2006a);  // ADDED: for CC2006 energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile3, "%s%s%s%s", outpath, outpath3, fileName, "_CCa.ct");  // ADDED: for CCa energy model (filename)
    FILE *cfile3 = fopen(ctFile3, "w");
    if (cfile3 == NULL) {
        printf("can't open %s \n", ctFile3);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i,"_CCa.ps");  // ADDED: for CCa energy model (filename)
            sprintf(bpseqFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i, "_CCa.bpseq");  // ADDED: for CCa energy model (filename)
        printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
        printRnaStruct(listOfNodes[i]->secStructure,length);
        printf("\n");
        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile3, listOfNodes[i]->score);
            FILE *bfile3 = fopen(bpseqFile3, "w");
            if (bfile3 == NULL) {
                printf("can't open %s \n", bpseqFile3);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile3, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
            }
        }
        fprintf(cfile3, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile3, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }

    /////////////////////////// CAO and CHEN (b) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode4=(struct Node *)malloc(sizeof(struct Node));
    rootNode4->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode4->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode4->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode4->constraint[length] = 0;
    rootNode4->numChild = 0;
    rootNode4->length = length;
    rootNode4->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode4->secStructure[i]=0;
        rootNode4->fixedPairs[i]=0;
        //rootNode4->constraint[i]='.';
    }
    strcpy(rootNode4->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode4,rootNode4, MaxHotspots, CC2006b);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile4, "%s%s%s%s", outpath, outpath4, fileName, "_CCb.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile4 = fopen(ctFile4, "w");
    if (cfile4 == NULL) {
        printf("can't open %s \n", ctFile4);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i,"_CCb.ps");  // ADDED: for CC2006b energy model (filename)
            sprintf(bpseqFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i, "_CCb.bpseq");  // ADDED: for CCb energy model (filename)
        printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
        printRnaStruct(listOfNodes[i]->secStructure,length);
        printf("\n");
        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile4, listOfNodes[i]->score);
            FILE *bfile4 = fopen(bpseqFile4, "w");
            if (bfile4 == NULL) {
                printf("can't open %s \n", bpseqFile4);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile4, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
            }
        }
        fprintf(cfile4, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile4, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }

    /////////////////////////// CAO and CHEN (c) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode5=(struct Node *)malloc(sizeof(struct Node));
    rootNode5->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode5->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode5->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode5->constraint[length] = 0;
    rootNode5->numChild = 0;
    rootNode5->length = length;
    rootNode5->score = 0;
    for(i=0;i<length;i++)
    {
        rootNode5->secStructure[i]=0;
        rootNode5->fixedPairs[i]=0;
        //rootNode4->constraint[i]='.';
    }
    strcpy(rootNode5->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode5,rootNode5, MaxHotspots, CC2006c);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile5, "%s%s%s%s", outpath, outpath5, fileName, "_CCc.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile5 = fopen(ctFile5, "w");
    if (cfile5 == NULL) {
        printf("can't open %s \n", ctFile5);
        exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
        sprintf(outPSFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i,"_CCc.ps");  // ADDED: for CC2006b energy model (filename)
            sprintf(bpseqFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i, "_CCc.bpseq");  // ADDED: for CCb energy model (filename)
        printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
        printRnaStruct(listOfNodes[i]->secStructure,length);
        printf("\n");
        if ((first == 1 && i==0) || (first == 0)) {
            if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile5, listOfNodes[i]->score);
            FILE *bfile5 = fopen(bpseqFile5, "w");
            if (bfile5 == NULL) {
                printf("can't open %s \n", bpseqFile5);
                exit(1);
            }
            for (int b = 0; b < length; b++) {
                fprintf(bfile5, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
            }
        }
        fprintf(cfile5, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
            fprintf(cfile5, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/
    /////////////// *** Add new energy model code here (similar to above section) /////////////


    //(void) fflush(stdout);
    //free(string);
    free(cstructure);

    // TODO
    // fill out structure with the predicted structure
    
    /*    
    //  Mirela: commented out this part
    char dpFile_temp[1000];
    sprintf(dpFile_temp, "%s%s", seqName, "0_DP.dp");

    // make sure the bpseq file is closed
 
    char bpseqFileTemp[1000];
    sprintf(bpseqFileTemp, "%s%s", seqName, "0_DP.bpseq");
    
//  FILE* ftemp = fopen (bpseqFileTemp, "r");
//  if (ftemp == NULL)
//  {
//      printf("Before system call, cannot open file %s\n", bpseqFileTemp);
//      exit(1);
//  }
//  else
//      fclose(ftemp);
  

    char command[200];
    // assumes bpseq2dp.pl and Analyser (needed for pl script) are in same folder as parameter tuning executables (ie. tools)
    sprintf(command, "/cs/local/bin/perl ./bpseq2dp.pl %s%s %s%s", seqName, "0_DP.bpseq", seqName, "0_DP.dp");

        if (DEBUG)
                printf("System call: %s\n", command);

    // Assumes bpseq2dp.pl is in the same folder as the executable Hotknots
    system(command);

    if (DEBUG)
        printf("Done system call to bpseq2dp.pl\n");

    FILE *dpfile_temp;
        char buffer[MaxN];
        char  v1[MaxN];
//  char dpFile_temp[] = "RNASeq_temp.dp";

    //printf ("FILENAME: %s\n", filename);
    if ((dpfile_temp = fopen (dpFile_temp, "r")) == NULL)
    {
        printf ("Cannot open file %s\n", dpFile_temp);
        exit(1);
    }

    fgets (buffer, sizeof(buffer), dpfile_temp);
    sscanf (buffer, "%s", v1);
    strcpy(structure, v1);
    fclose(dpfile_temp);

    if (DEBUGH)
        printf("Structure returned: %s\n", structure);

    // Mirela: added the following 3 lines to clean up the mess
    remove(dpFile_temp);
    remove(bpseqFileTemp);
    //remove (ctFile2);
*/

	num_suboptimals = numRnaStruct;
	for (int i=0; i < numRnaStruct; i++)
	{
		energies[i] = -listOfNodes[i]->score /1000.0;
	    bpseq2dp (length, listOfNodes[i]->secStructure, structures[i]);

		if (DEBUGH)
			printf("Structure returned: %s\n", structures[i]);
	}

    return -listOfNodes[0]->score /1000.0;   //return best suboptimal energy
}


double hotknots_CC2006b_name(char *seq, char *structure, char* fileName)
// fileName: input - does not include path, just the name of the file without extension (sequence name)
// sequence: input
// structure: output
// By Mirela: THE NAME iS NOT USED ANY MORE, I wrote a function called bpseq2dp in utils.c
//            TO CLEAN UP THE MESS IN THIS FILE.
{
	if (DEBUG)
		printf("Input parameter seq = %s, strlen = %d\n", seq, (int) strlen(seq));

	int length = strlen(seq);
	char* sequence = new char[length+1];
	for (int i = 0; i < length; i++)
		sequence[i] = toupper(seq[i]);
	sequence[length] = '\0';

	if (DEBUG)
		printf("\nInput to hotknots_energy: len=%d, seq = %s, fileName = %s\n", length, sequence, fileName);

  	char *string;
  	char line[5000];
  	char *cstructure=NULL, *cstruc=NULL;  // constraint structure
  	char  fname[13], ffname[20], gfname[20];

  	char  *ns_bases=NULL, *c;
  	int   i, l, sym, r;
  	double energy, min_en;
  	double kT, sfact=1.07;
  	int   pf=0, istty;
  	int noconv=0;
  	int InputFile = 1;
  	int StrCon = 0;    
  	int endFlag;
  	char seqName[1000];  // includes path, used for PlotRNA
	char inFile[1000];
//	char ctFile[100], outPSFile[100], bpseqFile[100];
//	char ctFile2[100], outPSFile2[100], bpseqFile2[100];
//	char ctFile3[100], outPSFile3[100], bpseqFile3[100];
	char ctFile4[1000], outPSFile4[1000], bpseqFile4[1000];
//	char ctFile5[100], outPSFile5[100], bpseqFile5[100];
//	// *** Add new energy model code here

// 	struct Node *rootNode;
//  	struct Node *rootNode2;
//  	struct Node *rootNode3;
  	struct Node *rootNode4;
//  	struct Node *rootNode5;
//	// *** Add new energy model code here

	int len = 0;
    int index = -1;  
    char separator;
//    char fileName[200];  // filename of input file sequence, without the path
    
    // *** Add new energy model code here (also change output path as desired)

    //Mirela: These paths are not used any more, I'm not writing any file any more, just using the function bpseq2dp, defined in utils.c   
    char outpath[30] = "../HotKnots/bin/TestSeq/";  // starting part of path for output files assuming directory of the parameter tuning files are in ../tools
//	char outpath[30] = "../bin/TestSeq/";  // for use with drivertemp.c
//    char outpath1[20] = "RivasEddy/";  // directory to put files generated using Rivas&Eddy (preceeded by outpath)
//    char outpath2[20] = "DirksPierce/";  // see above
    //char outpath3[20] = "CaoChen_a/";  // see above
    char outpath4[20] = "CaoChen_b/";  // see above
    //char outpath5[20] = "CaoChen_c/";  // see above

//	sprintf(seqName, "%s%s%s", outpath, outpath2, fileName);
	sprintf(seqName, "%s%s%s", outpath, outpath4, fileName);

	char tempFile[10]="hello";
  	FILE* input_file;
  	int MaxHotspots = 200;

//	string=NULL;
//   	length = (int) strlen(sequence);

	cstructure = (char *) malloc((length+1)*sizeof(char));

   	for (i = 0; i < length; i++) 
    	cstructure[i] = '.';
	cstructure[length] = 0;

	int first = 1;  // don't print suboptimals
	int noPS = 1;   // don't do the PS file
	int noGU = 0;   // include GU base pairs

//	char bpseqFile_temp[] = "RNASeq_temp.bpseq";

//        cout << "after init, COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;


	/////////////////////////// RIVAS and EDDY ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode=(struct Node *)malloc(sizeof(struct Node));
    rootNode->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode->constraint[length] = 0;
    rootNode->numChild = 0;
    rootNode->length = length;
    rootNode->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode->secStructure[i]=0;
		rootNode->fixedPairs[i]=0;
		//rootNode->constraint[i]='.';
    }
    strcpy(rootNode->constraint, cstructure);
    //===end of initialization    
 
	if (DEBUGH)
	{
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", string);    
	}
    
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);
    GenerateStemList(length, string, cstructure);
    
    endFlag=secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, RE);  // ADDED: for RE energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;

	if (DEBUGH)
	{
    printf("In total, %d nodes created. \n", count);
    printf(" total number of Rna structures: %d \n", numRnaStruct);
	}

    sprintf(ctFile, "%s%s%s%s", outpath, outpath1, fileName, "_RE.ct");  // ADDED: for RE energy model
    FILE *cfile = fopen(ctFile, "w");
    if (cfile == NULL) {
		printf("can't open %s \n", ctFile);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.ps");  // ADDED: for RE energy model
        	sprintf(bpseqFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.bpseq");  // ADDED: for RE energy model
		
		if (DEBUGH)
		{
			printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			printRnaStruct(listOfNodes[i]->secStructure,length);
			printf("\n");
		}

		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile, listOfNodes[i]->score);
	        FILE *bfile = fopen(bpseqFile, "w");
			if (bfile == NULL) {
		    	printf("can't open %s \n", bpseqFile);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);            
		  	}          
		}
	    fprintf(cfile, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/

/*
	/////////////////////////// DIRKS and PIERCE ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode2=(struct Node *)malloc(sizeof(struct Node));
    rootNode2->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode2->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode2->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode2->constraint[length] = 0;
    rootNode2->numChild = 0;
    rootNode2->length = length;
    rootNode2->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode2->secStructure[i]=0;
		rootNode2->fixedPairs[i]=0;
		//rootNode2->constraint[i]='.';
    }
    strcpy(rootNode2->constraint, cstructure);
    //===end of initialization    

	if (DEBUGH)
	{
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", sequence); 
	}
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section

//	printf("Input to Generate Stem List: len = %d, seq = %s, scruct = %s", length, sequence, cstructure);

    GenerateStemList(length, sequence, cstructure);  // depends only on hairpin, stacked, internal loop parameters

//	printf("After generate stem list\n");
//              cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//		cout << &dangle_top[2][1][0] << endl;

    endFlag=secondaryStruct(sequence,length,rootNode2,rootNode2,MaxHotspots, DP);  // ADDED: for DP energy model

//	printf("After secondaryStruct\n");

    ClearHotspots(MaxHotspots);

//                cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//		cout << &dangle_top[2][1][0] << endl;

//	printf("After ClearHotspots\n");

    min_en=endFlag;

	if (DEBUGH)
	{
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
	}

    sprintf(ctFile2, "%s%s%s%s", outpath, outpath2, fileName, "_DP.ct");
    FILE *cfile2 = fopen(ctFile2, "w");
    if (cfile2 == NULL) {
		printf("can't open %s \n", ctFile2);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i,"_DP.ps");  // ADDED: for DP energy model
        	sprintf(bpseqFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i, "_DP.bpseq");  // ADDED: for DP energy model

		if (DEBUGH)
		{
			printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			printRnaStruct(listOfNodes[i]->secStructure,length);
			printf("\n");
		}

		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, sequence, listOfNodes[i]->secStructure, outPSFile2, listOfNodes[i]->score);
	        FILE *bfile2 = fopen(bpseqFile2, "w");

		// Used later:
//		FILE *bfile_temp = fopen(bpseqFile_temp, "w");  // erase previous contents

			if (bfile2 == NULL) {
			    	printf("can't open %s \n", bpseqFile2);
		    		exit(1);
		  	}

//			if (bfile_temp == NULL) {
//				printf("can't open %s \n", bpseqFile_temp);
//			    	exit(1);
//		  	}

		  	for (int b = 0; b < length; b++) {
			    	fprintf(bfile2, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
//				fprintf(bfile_temp, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
		  	}          
//			fclose(bfile_temp);
			fclose(bfile2);
		}
	    fprintf(cfile2, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
	        for (int b = 0; b < length; b++) {
        		fprintf(cfile2, "%5d %c    %4d %4d %4d %4d\n", b+1, sequence[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
	        }
    	}
	fclose(cfile2);
*/
	/////////////////////////// CAO and CHEN (a) ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode3=(struct Node *)malloc(sizeof(struct Node));
    rootNode3->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode3->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode3->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode3->constraint[length] = 0;
    rootNode3->numChild = 0;
    rootNode3->length = length;
    rootNode3->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode3->secStructure[i]=0;
		rootNode3->fixedPairs[i]=0;
		//rootNode3->constraint[i]='.';
    }
    strcpy(rootNode3->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode3,rootNode3, MaxHotspots, CC2006a);  // ADDED: for CC2006 energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile3, "%s%s%s%s", outpath, outpath3, fileName, "_CCa.ct");  // ADDED: for CCa energy model (filename)
    FILE *cfile3 = fopen(ctFile3, "w");
    if (cfile3 == NULL) {
		printf("can't open %s \n", ctFile3);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i,"_CCa.ps");  // ADDED: for CCa energy model (filename)
        	sprintf(bpseqFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i, "_CCa.bpseq");  // ADDED: for CCa energy model (filename)
		printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
		printRnaStruct(listOfNodes[i]->secStructure,length);
		printf("\n");
		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile3, listOfNodes[i]->score);
	        FILE *bfile3 = fopen(bpseqFile3, "w");
			if (bfile3 == NULL) {
		    	printf("can't open %s \n", bpseqFile3);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile3, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
		  	}
		}
	    fprintf(cfile3, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile3, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/

	/////////////////////////// CAO and CHEN (b) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode4=(struct Node *)malloc(sizeof(struct Node));
    rootNode4->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode4->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode4->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode4->constraint[length] = 0;
    rootNode4->numChild = 0;
    rootNode4->length = length;
    rootNode4->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode4->secStructure[i]=0;
		rootNode4->fixedPairs[i]=0;
		//rootNode4->constraint[i]='.';
    }
    strcpy(rootNode4->constraint, cstructure);
    //===end of initialization    
 
	if (DEBUGH)
	{
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", sequence); 
	}

    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, sequence, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    // Mirela: added comment
    if (DEBUGH)  printf ("Stem list generated.\n");
    
    endFlag=secondaryStruct(sequence,length,rootNode4,rootNode4, MaxHotspots, CC2006b);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;

	if (DEBUGH)
	{
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    for (int k=0; k < numRnaStruct; k++)
    {
        char structure2[2000];
        bpseq2dp (length, listOfNodes[k]->secStructure, structure2);    
	printf("Structure %d: %s, energy=%.2lf\n", k, structure2, -listOfNodes[k]->score /1000.0);
    }
	}

    // Mirela: commented this out, I'm not writing any files any more, just calling the function bpse2dp defined in utils.c

    /*
    sprintf(ctFile4, "%s%s%s%s", outpath, outpath4, fileName, "_CCb.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile4 = fopen(ctFile4, "w");
    if (cfile4 == NULL) {
		printf("can't open %s \n", ctFile4);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i,"_CCb.ps");  // ADDED: for CC2006b energy model (filename)
        	sprintf(bpseqFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i, "_CCb.bpseq");  // ADDED: for CCb energy model (filename)

		if (DEBUGH)
		{
			printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			printRnaStruct(listOfNodes[i]->secStructure,length);
			printf("\n");
		}

		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, sequence, listOfNodes[i]->secStructure, outPSFile4, listOfNodes[i]->score);
	        FILE *bfile4 = fopen(bpseqFile4, "w");
			if (bfile4 == NULL) {
		    	printf("can't open %s \n", bpseqFile4);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile4, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
		  	}

			fclose(bfile4);
		}
	    fprintf(cfile4, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile4, "%5d %c    %4d %4d %4d %4d\n", b+1, sequence[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
	fclose(cfile4);
    */   


/*
	/////////////////////////// CAO and CHEN (c) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode5=(struct Node *)malloc(sizeof(struct Node));
    rootNode5->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode5->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode5->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode5->constraint[length] = 0;
    rootNode5->numChild = 0;
    rootNode5->length = length;
    rootNode5->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode5->secStructure[i]=0;
		rootNode5->fixedPairs[i]=0;
		//rootNode4->constraint[i]='.';
    }
    strcpy(rootNode5->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode5,rootNode5, MaxHotspots, CC2006c);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile5, "%s%s%s%s", outpath, outpath5, fileName, "_CCc.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile5 = fopen(ctFile5, "w");
    if (cfile5 == NULL) {
		printf("can't open %s \n", ctFile5);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i,"_CCc.ps");  // ADDED: for CC2006b energy model (filename)
        	sprintf(bpseqFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i, "_CCc.bpseq");  // ADDED: for CCb energy model (filename)
		printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
		printRnaStruct(listOfNodes[i]->secStructure,length);
		printf("\n");
		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile5, listOfNodes[i]->score);
	        FILE *bfile5 = fopen(bpseqFile5, "w");
			if (bfile5 == NULL) {
		    	printf("can't open %s \n", bpseqFile5);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile5, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
		  	}
		}
	    fprintf(cfile5, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile5, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/
	/////////////// *** Add new energy model code here (similar to above section) /////////////


    //(void) fflush(stdout);
    //free(string);
	free(cstructure);

    // Mirela: I'm not writing in files, I'm just calling bpseq2dp defined in utils.c
	// fill out structure with the predicted structure

    /*
	char dpFile_temp[150];
	sprintf(dpFile_temp, "%s%s", seqName, "0_CCb.dp");

//	printf("%s\n", dpFile_temp);

	// make sure the bpseq file is closed
	char bpseqFileTemp[1000];
	sprintf(bpseqFileTemp, "%s%s", seqName, "0_CCb.bpseq");
	FILE* ftemp = fopen (bpseqFileTemp, "r");
	if (ftemp == NULL)
	{
		printf("Before system call, cannot open file %s\n", bpseqFileTemp);
		exit(1);
	}
	else
		fclose(ftemp);

	char command[200];
	// assumes bpseq2dp.pl and Analyser (needed for pl script) are in same folder as parameter tuning executables (ie. tools)
	sprintf(command, "/cs/local/bin/perl ./bpseq2dp.pl %s%s %s%s", seqName, "0_CCb.bpseq", seqName, "0_CCb.dp");

        if (DEBUG)
                printf("System call: %s\n", command);

	// Assumes bpseq2dp.pl is in the same folder as the executable Hotknots
	system(command);

	if (DEBUG)
		printf("Done system call to bpseq2dp.pl\n");

	FILE *dpfile_temp;
    	char buffer[MaxN];
    	char  v1[MaxN];
//	char dpFile_temp[] = "RNASeq_temp.dp";

	//printf ("FILENAME: %s\n", filename);
	if ((dpfile_temp = fopen (dpFile_temp, "r")) == NULL)
	{
	    printf ("Cannot open file %s\n", dpFile_temp);
		exit(1);
	}

	fgets (buffer, sizeof(buffer), dpfile_temp);
	sscanf (buffer, "%s", v1);
	strcpy(structure, v1);
	fclose(dpfile_temp);
    
    // Mirela: added the following 3 lines to clean up the mess
    remove(dpFile_temp);
    remove(bpseqFileTemp);
    remove (ctFile4);    
    */

    bpseq2dp (length, listOfNodes[0]->secStructure, structure);
    if (DEBUGH)
	printf("Structure returned: %s, energy=%.2lf\n", structure, -listOfNodes[0]->score /1000.0);

    return -listOfNodes[0]->score /1000.0;   //return best suboptimal energy
}

// Prototype: double hotknots_CC2006b (char *sequence, char *structure)
//      - sequence is the given RNA sequence to be folded (input parameter)
//      - structure is the predicted secondary structure (output parameter)
//      - the function may return the energy
// Uses the Cao&Chen (b) model
// Defaults: no PS is output, GU is included (if want to change this, set it manually in ? TODO )
double hotknots_CC2006b(char *sequence, char *structure)
{       
        return hotknots_CC2006b_name (sequence, structure, (char *)  "Default");
}

// Prototype: double hotknots_CC2006b_suboptimals (char *sequence, char **structures, double* energies)
double hotknots_CC2006b_suboptimals(char *seq, char **structures, double* energies, int& num_suboptimals)
// fileName: input - does not include path, just the name of the file without extension (sequence name)
// sequence: input
// structure: output
{
	if (DEBUG)
		printf("Input parameter seq = %s, strlen = %d\n", seq, (int) strlen(seq));

	int length = strlen(seq);
	char* sequence = new char[length+1];
	for (int i = 0; i < length; i++)
		sequence[i] = toupper(seq[i]);
	sequence[length] = '\0';

	if (DEBUG)
		printf("\nInput to hotknots_energy: len=%d, seq = %s\n", length, sequence);

  	char *string;
  	char line[5000];
  	char *cstructure=NULL, *cstruc=NULL;  // constraint structure
  	char  fname[13], ffname[20], gfname[20];

  	char  *ns_bases=NULL, *c;
  	int   i, l, sym, r;
  	double energy, min_en;
  	double kT, sfact=1.07;
  	int   pf=0, istty;
  	int noconv=0;
  	int InputFile = 1;
  	int StrCon = 0;    
  	int endFlag;
  	char seqName[1000];  // includes path, used for PlotRNA
	char inFile[1000];
//	char ctFile[100], outPSFile[100], bpseqFile[100];
//	char ctFile2[100], outPSFile2[100], bpseqFile2[100];
//	char ctFile3[100], outPSFile3[100], bpseqFile3[100];
	char ctFile4[1000], outPSFile4[1000], bpseqFile4[1000];
//	char ctFile5[100], outPSFile5[100], bpseqFile5[100];
//	// *** Add new energy model code here

// 	struct Node *rootNode;
//  	struct Node *rootNode2;
//  	struct Node *rootNode3;
  	struct Node *rootNode4;
//  	struct Node *rootNode5;
//	// *** Add new energy model code here

	int len = 0;
    int index = -1;  
    char separator;
//    char fileName[200];  // filename of input file sequence, without the path
    
    // *** Add new energy model code here (also change output path as desired)

    //Mirela: These paths are not used any more, I'm not writing any file any more, just using the function bpseq2dp, defined in utils.c   
    char outpath[30] = "../HotKnots/bin/TestSeq/";  // starting part of path for output files assuming directory of the parameter tuning files are in ../tools
//	char outpath[30] = "../bin/TestSeq/";  // for use with drivertemp.c
//    char outpath1[20] = "RivasEddy/";  // directory to put files generated using Rivas&Eddy (preceeded by outpath)
//    char outpath2[20] = "DirksPierce/";  // see above
    //char outpath3[20] = "CaoChen_a/";  // see above
    char outpath4[20] = "CaoChen_b/";  // see above
    //char outpath5[20] = "CaoChen_c/";  // see above

//	sprintf(seqName, "%s%s%s", outpath, outpath2, fileName);
	sprintf(seqName, "%s", "default");

	char tempFile[10]="hello";
  	FILE* input_file;
  	int MaxHotspots = 200;

//	string=NULL;
//   	length = (int) strlen(sequence);

	cstructure = (char *) malloc((length+1)*sizeof(char));

   	for (i = 0; i < length; i++) 
    	cstructure[i] = '.';
	cstructure[length] = 0;

	int first = 1;  // don't print suboptimals
	int noPS = 1;   // don't do the PS file
	int noGU = 0;   // include GU base pairs

//	char bpseqFile_temp[] = "RNASeq_temp.bpseq";

//        cout << "after init, COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;


	/////////////////////////// RIVAS and EDDY ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode=(struct Node *)malloc(sizeof(struct Node));
    rootNode->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode->constraint[length] = 0;
    rootNode->numChild = 0;
    rootNode->length = length;
    rootNode->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode->secStructure[i]=0;
		rootNode->fixedPairs[i]=0;
		//rootNode->constraint[i]='.';
    }
    strcpy(rootNode->constraint, cstructure);
    //===end of initialization    
 
	if (DEBUGH)
	{
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", string);    
	}
    
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);
    GenerateStemList(length, string, cstructure);
    
    endFlag=secondaryStruct(string,length,rootNode,rootNode, MaxHotspots, RE);  // ADDED: for RE energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;

	if (DEBUGH)
	{
    printf("In total, %d nodes created. \n", count);
    printf(" total number of Rna structures: %d \n", numRnaStruct);
	}

    sprintf(ctFile, "%s%s%s%s", outpath, outpath1, fileName, "_RE.ct");  // ADDED: for RE energy model
    FILE *cfile = fopen(ctFile, "w");
    if (cfile == NULL) {
		printf("can't open %s \n", ctFile);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.ps");  // ADDED: for RE energy model
        	sprintf(bpseqFile, "%s%s%s%d%s", outpath, outpath1, fileName, i, "_RE.bpseq");  // ADDED: for RE energy model
		
		if (DEBUGH)
		{
			printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			printRnaStruct(listOfNodes[i]->secStructure,length);
			printf("\n");
		}

		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile, listOfNodes[i]->score);
	        FILE *bfile = fopen(bpseqFile, "w");
			if (bfile == NULL) {
		    	printf("can't open %s \n", bpseqFile);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);            
		  	}          
		}
	    fprintf(cfile, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/

/*
	/////////////////////////// DIRKS and PIERCE ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode2=(struct Node *)malloc(sizeof(struct Node));
    rootNode2->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode2->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode2->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode2->constraint[length] = 0;
    rootNode2->numChild = 0;
    rootNode2->length = length;
    rootNode2->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode2->secStructure[i]=0;
		rootNode2->fixedPairs[i]=0;
		//rootNode2->constraint[i]='.';
    }
    strcpy(rootNode2->constraint, cstructure);
    //===end of initialization    

	if (DEBUGH)
	{
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", sequence); 
	}
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section

//	printf("Input to Generate Stem List: len = %d, seq = %s, scruct = %s", length, sequence, cstructure);

    GenerateStemList(length, sequence, cstructure);  // depends only on hairpin, stacked, internal loop parameters

//	printf("After generate stem list\n");
//              cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//		cout << &dangle_top[2][1][0] << endl;

    endFlag=secondaryStruct(sequence,length,rootNode2,rootNode2,MaxHotspots, DP);  // ADDED: for DP energy model

//	printf("After secondaryStruct\n");

    ClearHotspots(MaxHotspots);

//                cout << "COUT: dangle_top[2][1][0] " << dangle_top[2][1][0] << endl;
//		cout << &dangle_top[2][1][0] << endl;

//	printf("After ClearHotspots\n");

    min_en=endFlag;

	if (DEBUGH)
	{
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
	}

    sprintf(ctFile2, "%s%s%s%s", outpath, outpath2, fileName, "_DP.ct");
    FILE *cfile2 = fopen(ctFile2, "w");
    if (cfile2 == NULL) {
		printf("can't open %s \n", ctFile2);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i,"_DP.ps");  // ADDED: for DP energy model
        	sprintf(bpseqFile2, "%s%s%s%d%s", outpath, outpath2, fileName, i, "_DP.bpseq");  // ADDED: for DP energy model

		if (DEBUGH)
		{
			printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			printRnaStruct(listOfNodes[i]->secStructure,length);
			printf("\n");
		}

		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, sequence, listOfNodes[i]->secStructure, outPSFile2, listOfNodes[i]->score);
	        FILE *bfile2 = fopen(bpseqFile2, "w");

		// Used later:
//		FILE *bfile_temp = fopen(bpseqFile_temp, "w");  // erase previous contents

			if (bfile2 == NULL) {
			    	printf("can't open %s \n", bpseqFile2);
		    		exit(1);
		  	}

//			if (bfile_temp == NULL) {
//				printf("can't open %s \n", bpseqFile_temp);
//			    	exit(1);
//		  	}

		  	for (int b = 0; b < length; b++) {
			    	fprintf(bfile2, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
//				fprintf(bfile_temp, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
		  	}          
//			fclose(bfile_temp);
			fclose(bfile2);
		}
	    fprintf(cfile2, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
	        for (int b = 0; b < length; b++) {
        		fprintf(cfile2, "%5d %c    %4d %4d %4d %4d\n", b+1, sequence[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
	        }
    	}
	fclose(cfile2);
*/
	/////////////////////////// CAO and CHEN (a) ENERGY MODEL /////////////////////////
/*
    //-----initialization of rootNode
    rootNode3=(struct Node *)malloc(sizeof(struct Node));
    rootNode3->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode3->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode3->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode3->constraint[length] = 0;
    rootNode3->numChild = 0;
    rootNode3->length = length;
    rootNode3->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode3->secStructure[i]=0;
		rootNode3->fixedPairs[i]=0;
		//rootNode3->constraint[i]='.';
    }
    strcpy(rootNode3->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode3,rootNode3, MaxHotspots, CC2006a);  // ADDED: for CC2006 energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile3, "%s%s%s%s", outpath, outpath3, fileName, "_CCa.ct");  // ADDED: for CCa energy model (filename)
    FILE *cfile3 = fopen(ctFile3, "w");
    if (cfile3 == NULL) {
		printf("can't open %s \n", ctFile3);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i,"_CCa.ps");  // ADDED: for CCa energy model (filename)
        	sprintf(bpseqFile3, "%s%s%s%d%s", outpath, outpath3, fileName, i, "_CCa.bpseq");  // ADDED: for CCa energy model (filename)
		printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
		printRnaStruct(listOfNodes[i]->secStructure,length);
		printf("\n");
		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile3, listOfNodes[i]->score);
	        FILE *bfile3 = fopen(bpseqFile3, "w");
			if (bfile3 == NULL) {
		    	printf("can't open %s \n", bpseqFile3);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile3, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
		  	}
		}
	    fprintf(cfile3, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile3, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/

	/////////////////////////// CAO and CHEN (b) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode4=(struct Node *)malloc(sizeof(struct Node));
    rootNode4->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode4->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode4->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode4->constraint[length] = 0;
    rootNode4->numChild = 0;
    rootNode4->length = length;
    rootNode4->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode4->secStructure[i]=0;
		rootNode4->fixedPairs[i]=0;
		//rootNode4->constraint[i]='.';
    }
    strcpy(rootNode4->constraint, cstructure);
    //===end of initialization    
 
	if (DEBUGH)
	{
    printf("LENGTH OF THE RNA IS %d.\n",length);
    printf("%s\n", sequence); 
	}

    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, sequence, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(sequence,length,rootNode4,rootNode4, MaxHotspots, CC2006b);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;

	if (DEBUGH)
	{
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
	}

    // Mirela: I'm not writing in files, I'm just calling bpseq2dp defined in utils.c
	// fill out structure with the predicted structure

    /*
    sprintf(ctFile4, "%s%s%s%s", outpath, outpath4, fileName, "_CCb.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile4 = fopen(ctFile4, "w");
    if (cfile4 == NULL) {
		printf("can't open %s \n", ctFile4);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i,"_CCb.ps");  // ADDED: for CC2006b energy model (filename)
        	sprintf(bpseqFile4, "%s%s%s%d%s", outpath, outpath4, fileName, i, "_CCb.bpseq");  // ADDED: for CCb energy model (filename)

		if (DEBUGH)
		{
			printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
			printRnaStruct(listOfNodes[i]->secStructure,length);
			printf("\n");
		}

		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, sequence, listOfNodes[i]->secStructure, outPSFile4, listOfNodes[i]->score);
	        FILE *bfile4 = fopen(bpseqFile4, "w");
			if (bfile4 == NULL) {
		    	printf("can't open %s \n", bpseqFile4);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile4, "%d %c %d\n", b+1, sequence[b], listOfNodes[i]->secStructure[b]);
		  	}

			fclose(bfile4);
		}
	    fprintf(cfile4, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile4, "%5d %c    %4d %4d %4d %4d\n", b+1, sequence[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
	fclose(cfile4);
    */   


/*
	/////////////////////////// CAO and CHEN (c) ENERGY MODEL /////////////////////////

    //-----initialization of rootNode
    rootNode5=(struct Node *)malloc(sizeof(struct Node));
    rootNode5->secStructure=(short *)malloc((length)*sizeof(short));
    rootNode5->fixedPairs=(short *)malloc((length)*sizeof(short));
    rootNode5->constraint=(char *)malloc((length+1)*sizeof(char));
    rootNode5->constraint[length] = 0;
    rootNode5->numChild = 0;
    rootNode5->length = length;
    rootNode5->score = 0;
    for(i=0;i<length;i++)
    {
		rootNode5->secStructure[i]=0;
		rootNode5->fixedPairs[i]=0;
		//rootNode4->constraint[i]='.';
    }
    strcpy(rootNode5->constraint, cstructure);
    //===end of initialization    
 
    numRnaStruct = 0;
    InitHotspots(MaxHotspots,length);  // it was cleared in the previous section
    GenerateStemList(length, string, cstructure);  // depends only on hairpin, stacked, internal loop parameters
    
    endFlag=secondaryStruct(string,length,rootNode5,rootNode5, MaxHotspots, CC2006c);  // ADDED: for CC2006withDP energy model
    ClearHotspots(MaxHotspots);
    min_en=endFlag;
    printf("In total, %d nodes created. \n", count);  // 'count' is set in function: secondaryStruct()
    printf(" total number of Rna structures: %d \n", numRnaStruct);
    sprintf(ctFile5, "%s%s%s%s", outpath, outpath5, fileName, "_CCc.ct");  // ADDED: for CCb energy model (filename)
    FILE *cfile5 = fopen(ctFile5, "w");
    if (cfile5 == NULL) {
		printf("can't open %s \n", ctFile5);
	  	exit(1);
    }

    for (i=0; i < numRnaStruct; i++) {
		sprintf(outPSFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i,"_CCc.ps");  // ADDED: for CC2006b energy model (filename)
        	sprintf(bpseqFile5, "%s%s%s%d%s", outpath, outpath5, fileName, i, "_CCc.bpseq");  // ADDED: for CCb energy model (filename)
		printf(" sub optimal energy  %f kal/mol\n", -listOfNodes[i]->score /1000.0);
		printRnaStruct(listOfNodes[i]->secStructure,length);
		printf("\n");
		if ((first == 1 && i==0) || (first == 0)) {
			if (noPS == 0) PlotRna(seqName, string, listOfNodes[i]->secStructure, outPSFile5, listOfNodes[i]->score);
	        FILE *bfile5 = fopen(bpseqFile5, "w");
			if (bfile5 == NULL) {
		    	printf("can't open %s \n", bpseqFile5);
		    	exit(1);
		  	}
		  	for (int b = 0; b < length; b++) {
		    	fprintf(bfile5, "%d %c %d\n", b+1, string[b], listOfNodes[i]->secStructure[b]);
		  	}
		}
	    fprintf(cfile5, "%5d   ENERGY = %.2f %s suboptimal structure %d\n", length, -listOfNodes[i]->score/1000.0, seqName, i);
        for (int b = 0; b < length; b++) {
        	fprintf(cfile5, "%5d %c    %4d %4d %4d %4d\n", b+1, string[b], b, b+2, listOfNodes[i]->secStructure[b], b+1);
        }
    }
*/
	/////////////// *** Add new energy model code here (similar to above section) /////////////


    //(void) fflush(stdout);
    //free(string);
	free(cstructure);

    // Mirela: I'm not writing in files, I'm just calling bpseq2dp defined in utils.c
	// fill out structure with the predicted structure

    /*
	char dpFile_temp[150];
	sprintf(dpFile_temp, "%s%s", seqName, "0_CCb.dp");

//	printf("%s\n", dpFile_temp);

	// make sure the bpseq file is closed
	char bpseqFileTemp[1000];
	sprintf(bpseqFileTemp, "%s%s", seqName, "0_CCb.bpseq");
	FILE* ftemp = fopen (bpseqFileTemp, "r");
	if (ftemp == NULL)
	{
		printf("Before system call, cannot open file %s\n", bpseqFileTemp);
		exit(1);
	}
	else
		fclose(ftemp);

	char command[200];
	// assumes bpseq2dp.pl and Analyser (needed for pl script) are in same folder as parameter tuning executables (ie. tools)
	sprintf(command, "/cs/local/bin/perl ./bpseq2dp.pl %s%s %s%s", seqName, "0_CCb.bpseq", seqName, "0_CCb.dp");

        if (DEBUG)
                printf("System call: %s\n", command);

	// Assumes bpseq2dp.pl is in the same folder as the executable Hotknots
	system(command);

	if (DEBUG)
		printf("Done system call to bpseq2dp.pl\n");

	FILE *dpfile_temp;
    	char buffer[MaxN];
    	char  v1[MaxN];
//	char dpFile_temp[] = "RNASeq_temp.dp";

	//printf ("FILENAME: %s\n", filename);
	if ((dpfile_temp = fopen (dpFile_temp, "r")) == NULL)
	{
	    printf ("Cannot open file %s\n", dpFile_temp);
		exit(1);
	}

	fgets (buffer, sizeof(buffer), dpfile_temp);
	sscanf (buffer, "%s", v1);
	strcpy(structure, v1);
	fclose(dpfile_temp);
    
    // Mirela: added the following 3 lines to clean up the mess
    remove(dpFile_temp);
    remove(bpseqFileTemp);
    remove (ctFile4);    
    */

	num_suboptimals = numRnaStruct;
	for (int i=0; i < numRnaStruct; i++)
	{
		energies[i] = -listOfNodes[i]->score /1000.0;
	    bpseq2dp (length, listOfNodes[i]->secStructure, structures[i]);

		if (DEBUGH)
			printf("Structure returned: %s\n", structures[i]);
	}

    return -listOfNodes[0]->score /1000.0;   //return best suboptimal energy
}
