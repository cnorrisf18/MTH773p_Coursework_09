#pragma once
class Stats1
{
public:
	Stats1();
	void add(double num);
	double mean();
	double stDev();
	double length;
	double sumofnumbers;
private:


	double sumofsquare;
	double currentstDev;
	double currentmean;
	double currentdevsum;
};

class NormalRandomGenerator
{
public:
	NormalRandomGenerator();
	double generate();

private:
	bool usedZ1;
	double innerZ2;
};