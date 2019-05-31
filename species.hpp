#ifndef __species_hpp__
#define __species_hpp__
#include<string>
#include<vector>
#include<fstream>

struct species 
{
	//	Species name
	std::string name;
	//	Species mass scaled by the smallest mass starting from 1
	double mass;
	//	Species diameter scaled by the smallest diameter starting from 1
	double diameter;
	//	Species Collision Energy
	double energy;
};

#endif /* __species_hpp__*/
