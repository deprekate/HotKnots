#ifndef SCORE_H
#define SCORE_H

class Score{
  int sc;
  int F1( int len,  char* csequence, int* sequence, int i, int j);
  int intloop( int len,  int* s,
     int i, int k, int l, int j);
 public:
  Score();
  ~Score();
  long score( int len, char* s, short* p,int TRACE, int energyModel);
  int F2( int len,  int* sequence,  int i, int k, int l, int j);
  int sPair( int len,  int* sequence,  int i,  int j);
};

#endif
