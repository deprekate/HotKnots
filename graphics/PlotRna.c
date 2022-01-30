#include <stdio.h>
#include "graphics.h"
#include <stdlib.h>
#include <string.h>

#define CELLSIZE 100


static void drawInit(int sizex, int sizey, char* seqName, char* filename, float score);
static void button_press (float x, float y);
static void drawscreen (void);
void PlotRna(char* seqName, char* sequence, short* structure, char* filename, float score);
static void plotBP(char* sequence, short* structure);
char* fname;

/*
int main() {
  char sequence[] = "UUAGAGUGUUCAAAGCCAGGCGUUUGCCGUGAAUACAUUAGCAUGGAA";
  short structure[48] = {0, 0, 0, 0, 0, 36, 35, 0, 33, 32, 31, 30, 0, 0, 0, 46, 45, 44, 28, 27, 26, 0, 0, 0, 0, 21, 20, 19, 0, 12, 11, 10, 9, 0, 7, 6, 0, 0, 0, 0, 0, 0, 0, 18, 17, 16, 0, 0};
  char filename[] = "test.ps";

  PlotRna(sequence, structure, filename, int score);
}*/

static void drawInit( int sizex, int sizey, char* seqName, char* filename, float score) {

 int i,k=0;
 int x,y;
 char* c;
 c =(char*)malloc(sizeof(char)*1000);
 


/* initialize display */
 init_graphics("RNA");


/* still picture drawing allows user to zoom, etc. */
 init_world (0.,0.,CELLSIZE*sizex,CELLSIZE*sizey);
 //event_loop(button_press, drawscreen);   

/* animation section */
 clearscreen();
 init_postscript(filename);
 clearscreen();
 
 sprintf(c, "%s%s%s%d%s%f%s", "RNA: ", seqName, ", length = ", sizex, ", energy = ",-score/1000.0, " kal/mol");
 //printf(" output string is %s\n",c);
 update_message("RNA secondary structure");
 flushinput();
 
 drawtext (sizex*CELLSIZE/2,sizey*CELLSIZE/10,c,1.0e6);
 flushinput();
 setcolor (BLACK);
 setlinewidth(1);
 setlinestyle (SOLID);
 for (i=0;i<=sizex;i++) {
    drawline (i*CELLSIZE,(sizey/2-1)*CELLSIZE,i*CELLSIZE,CELLSIZE*sizey/2);
    flushinput();
 }
 drawline(0,(sizey/2-1)*CELLSIZE, sizex*CELLSIZE, CELLSIZE*(sizey/2-1));
 drawline(0,sizey*CELLSIZE/2, sizex*CELLSIZE, CELLSIZE*sizey/2);
 flushinput();
 free(c);
}

void PlotRna(char* seqName, char* sequence, short* structure, char* filename, float score)
{
  /*A:red, U: Green, G:Blue, C:Yellow*/
  int i, j;
  int length;
  int sizex, sizey;
  
  
  fname = filename;
  length = (int)strlen(sequence) ;
  sizex = length;
  sizey = length;
  /*printf("here......length %d\n",length);*/
  drawInit(length, length, seqName, filename, score);
  /* print out stems */
  int st, stp;
  char* c = (char*)malloc(sizeof(char)*1000);


  i = 0; j = 0;
  while(i < length) {
    while(i < length && (structure[i] == 0 || structure[i] < i)) i++;
    st = i;    
    while(i < length && (structure[i]-1 == structure[i+1])) i++; 
    stp = i;
    if (i >= length) break;
    
    sprintf(c, "%d-%d; %d-%d ",st+1, stp+1, structure[stp],structure[st]);
    i++; j++;
    if ((j%2) == 1) drawtext (sizex*CELLSIZE/4,sizey*CELLSIZE/10.0+((int)(j/2)+1)*CELLSIZE*sizey/30,c,1.0e6);
    else drawtext (sizex*CELLSIZE/2,sizey*CELLSIZE/10.0+((int)(j/2))*CELLSIZE*sizey/30,c,1.0e6);
    flushinput();

  }

 
 
  
  for(i=0; i<length; i++) {
    if(sequence[i] == 'A' || sequence[i] == 'a') 
      setcolor(RED);
    else if(sequence[i] == 'U'|| sequence[i] == 'u' || sequence[i] == 'T' || sequence[i] == 't')
      setcolor(GREEN);
    else if(sequence[i] == 'G' || sequence[i] == 'g')
      setcolor(BLUE);
    else if(sequence[i] == 'C' || sequence[i] == 'c')
      setcolor(YELLOW);
    else
      setcolor(BLACK);
    fillrect((i+0.2)*CELLSIZE, (sizey/2 - 0.8)*CELLSIZE, (i+0.8)*CELLSIZE, (sizey/2 - 0.2)*CELLSIZE);
    flushinput();
  }
  plotBP(sequence, structure);
  close_postscript();
  close_graphics();
}

static void plotBP(char* sequence, short* structure)
{
  int i;
  int length;
  float lx, ly;
  float arcx, arcy, rad;
  length = (int)strlen(sequence) ;

  setlinestyle(SOLID);
  setcolor(BLACK);

  for(i = 0; i<length; i++) {
    if(structure[i] > i) {
      lx = structure[i] - i + 1;
      
      //drawline((i+0.5)*CELLSIZE, length*CELLSIZE/2 - 1, (i + 0.5)*CELLSIZE, (length+lx)*CELLSIZE/2 -1);
      flushinput();
    
      //drawline((structure[i]-0.5)*CELLSIZE, length*CELLSIZE/2-1,(structure[i]-0.5)*CELLSIZE, (length+lx)*CELLSIZE/2 -1);
      flushinput();
      //drawline((i+0.5)*CELLSIZE, (length+lx)*CELLSIZE/2-1, (structure[i]-0.5)*CELLSIZE, (length+lx)*CELLSIZE/2-1);
      arcx = (i + structure[i] )/2.0 * CELLSIZE;
      arcy = length*CELLSIZE/2 -1;
      rad = (structure[i] - 1 - i )/2.0 *CELLSIZE;
      drawarc(arcx, arcy,  rad, 0, -180.0);
      flushinput();
    }
  }
}
      
  
  
  


static void drawscreen(void) {

  return;
}

static void button_press (float x, float y) {

/* Called whenever event_loop gets a button press in the graphics *
 * area.  Allows the user to do whatever he/she wants with button *
 * clicks.                                                        */

 printf("User clicked a button at coordinates (%f, %f)\n", x, y);
}


 































