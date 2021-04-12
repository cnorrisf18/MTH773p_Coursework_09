#pragma once
#include <vector>
using namespace std;

struct Results
{
	double price;
	double error;
	int N;
};

class EurCall
{
public:
	EurCall(double T, double K) : m_T(T), m_K(K) {};

	Results HestonMonteCarloPricewithN(double S0, double r, double a, double vsquare, double vol_of_vol, double v0, double p, int time_steps, int sample_paths) const;
	Results HestonMonteCarloPricewithError(double S0, double r, double a, double vsquare, double vol_of_vol, double v0, double p, int time_steps, double req_error) const;

	double PriceByBSFormula(double S0, double sigma, double r) const;
	double N(double x) const;
	double d_plus(double S0, double sigma, double r) const;
	double d_minus(double S0, double sigma, double r) const;

	double AnalyticalPriceByJumpDiffusion(double S0, double r, double sigma, double lambda, double m, double s) const;
	Results MCPriceByJumpDiffusion(double S0, double r, double sigma, double lambda, double m, double s, double time_steps, double req_error) const;

private:
	const double m_T;
	const double m_K;
	double payoff(double S) const;
};

class Function
{
public:
	virtual double eval(double x) const = 0;
};


double BisectionSolver(const Function* f, double c, double left, double right, double acc);

class Adapter : public Function
{
public:
	Adapter(const EurCall* call, double S0, double r)
		: m_call(call), m_S0(S0), m_r(r) {};
	double eval(double x) const
	{
		return m_call->PriceByBSFormula(m_S0,x,m_r);
	}
private:
	const EurCall* const m_call;
	const double m_S0;
	const double m_r;
};
