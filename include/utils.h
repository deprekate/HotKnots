short *make_pair_table(char *structure);
void nrerror(char *message);
void printRnaStruct(short* structure, int length);   

void bpseq2dp (int length, short int *pair, char *structure);
// Function that transforms from "bpseq" format to dot-parentheses format
// Input parameters: 
//      length = the length of the sequence
//      pair = an int array with the pairings. unpaired is 0 (not -1!), so all indeces are shifted by 1
// Output: structure
// It support 30 types of symbols (see symbol_left and symbol_right below.
// It uses the first available symbol, and when needed it goes to the next symbol. 
// Note: there's not relationship between the symbols used and the 
//      minimum number of base pairs to break in order to remove the pseudoknots!!
// Written by Mirela on August 2, 2008.
