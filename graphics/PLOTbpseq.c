#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "graphics.h"

extern void PlotRna(char* seqname, char *sequence, short *structure, char *filename,float score);

int main(int argc, char ** argv) {
  char outPSFile[100], ilmOutFile[100], bpseqFile[100], seqFile[100], mwmFile[100];
  char prefix[100];

  if (argc < 3) {
    printf("Usage: Plot seqFile bpseqFile\n");
    exit(1);
  }
  
  strcpy(seqFile, argv[1]);
  strcpy(bpseqFile, argv[2]);
  
  strcpy(prefix, argv[2]);
  int len;
  len = strlen(prefix);
  prefix[len-6] = '\0';
  strcpy(outPSFile, prefix);
  strcat(outPSFile,".ps");


  int MAX = 1000;
  int MaxLen = 5000;
  short structure[MaxLen];  
  int j;



  FILE *bpfile = fopen(bpseqFile, "r+");
  int i = 0;
  while(!feof(bpfile)) {
    fscanf(bpfile, "%d %d", &i, &j); 
    structure[i-1] = j;
  }

  fclose(bpfile);

  FILE* seqfile = fopen(seqFile, "r+");
  char* seq = (char*)malloc(MaxLen*sizeof(char));
  if (seqfile == NULL) {
    printf("Can't open file %s \n", seqFile);
    exit(1);
  }

  fscanf(seqfile,"%s", seq);
  fclose(seqfile);
  
    
  
  float energy = 0;
  
  PlotRna(prefix, seq, structure, outPSFile, energy*1000);
  
}
 
