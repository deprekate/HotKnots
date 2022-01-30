//#include "Sequence.h"
#ifndef LOOPLIST_H
#define LOOPLIST_H

#include "Defines.h"
#include "Input.h"
#include "params.h"  // Cristina: May 7
#include "paramsPK.h"  // Cristina: May 16

class LoopList{

    public:

	//FUNCTIONS	
	LoopList(ReadInput* R, int b1, int b2);
	~LoopList();
	float 	interiorPseudoEnergyRE();
	float	interiorPseudoEnergyDP();
	float	interiorPseudoEnergyCC2006a();
	float	multiPseudoEnergyRE();
	float	multiPseudoEnergyDP();
	float	multiPseudoEnergyDP(double *c, double &f, int reset_c, int ignore_dangles);  // for parameter tuning//	float	multiPseudoEnergyCC2006();  // Unecessary
	float   stackPseudoEnergyRE();
	float   stackPseudoEnergyDP();
	float	stackPseudoEnergyCC2006a();
	//	float	LoopList::stackEnergy();  // line gives error in Eclipse 3.0.0 (this version is 3.1.0)
	void FindChildren();

	//ATTRIBUTES

	int base1, base2;
	int Size;
	ReadInput* Input;

};

#endif
