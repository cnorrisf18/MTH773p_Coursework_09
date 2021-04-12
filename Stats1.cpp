#include "Stats1.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <random>
#include <cmath>
using namespace std;

//Stats1 class definitions

Stats1::Stats1()
{
	sumofnumbers = 0;
	sumofsquare = 0;
	length = 0;
	currentstDev = 0;
	currentmean = 0;
	currentdevsum = 0;
}
;
void Stats1::add(double num) { //add numbers in one at a time and sum them
	sumofnumbers += num;
	sumofsquare += pow(num, 2);
	length += 1.0;


}
double Stats1::mean() {
	return sumofnumbers / length;

}
double Stats1::stDev() {
	return sqrt((sumofsquare / length) - pow((sumofnumbers / length), 2));

}


//NormalRandomGenerator class definitions
NormalRandomGenerator::NormalRandomGenerator()
{
	usedZ1 = false;
	innerZ2 = 0;
};
double NormalRandomGenerator::generate()
{
	if (usedZ1 == true)
	{
		usedZ1 = false;
		return innerZ2;
	}
	else
	{
		double U1 = rand() / double(RAND_MAX);
		bool done = false;
		while (not done)
		{
			if (U1 == 0)
			{
				U1 = rand() / double(RAND_MAX);
			}
			else {
				done = true;
			}
		}
		double U2 = rand() / double(RAND_MAX);
		double Z1 = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);
		innerZ2 = sqrt(-2.0 * log(U1)) * sin(2.0 * M_PI * U2);
		usedZ1 = true;
		return Z1;
	}
}

