#include "Options.h"
#include "Stats1.h"
#include <vector>
#include <iostream>
using namespace std;

Results EurCall::HestonMonteCarloPricewithN(double S0, double r, double a, double vsquare, double vol_of_vol, double v0, double p, int time_steps, int sample_paths) const
{
	const double dt = m_T / time_steps;

	Stats1 paths_stats;
	NormalRandomGenerator nrg;

	for (unsigned int i = 0; i < sample_paths; ++i)
	{
		double S = S0;
		double v = v0;
		for (unsigned int j = 1; j <= time_steps; ++j)
		{
			double e1 = nrg.generate();
			double e2 = nrg.generate();

			double Z1 = e1;
			double Z2 = (p * e1) + (sqrt(1 - pow(p, 2)) * e2);

			S *= exp(((r - (v / 2)) * dt) + (sqrt(v * dt) * Z1));
			v = pow((sqrt(v) + ((vol_of_vol / 2.0) * sqrt(dt) * Z2)), 2.0) - (a * (v - vsquare) * dt) - ((pow(vol_of_vol, 2.0) / 4.0) * dt);
		}
		const double po = payoff(S);
		paths_stats.add(po);
	}
	Results results;
	results.price = exp(-r * m_T) * paths_stats.mean();
	results.error = exp(-r * m_T) * paths_stats.stDev() / sqrt(sample_paths - 1.0);
	results.N = sample_paths;
	return results;
}
Results EurCall::HestonMonteCarloPricewithError(double S0, double r, double a, double vsquare, double vol_of_vol, double v0, double p, int time_steps, double req_error) const
{
	const double dt = m_T / time_steps;

	Stats1 paths_stats;
	NormalRandomGenerator nrg;

	int sample_paths = 1;

	bool done = false;
	while (not done)
	{
		double S = S0;
		double v = v0;
		for (unsigned int j = 1; j <= time_steps; ++j)
		{
			double e1 = nrg.generate();
			double e2 = nrg.generate();

			double Z1 = e1;
			double Z2 = (p * e1) + (sqrt(1 - pow(p, 2)) * e2);

			S *=  exp(((r - (v/ 2)) * dt) + (sqrt(v * dt) * Z1));
			v = pow((sqrt(v) + ((vol_of_vol / 2.0) * sqrt(dt) * Z2)), 2.0) - (a * (v - vsquare) * dt) - ((pow(vol_of_vol, 2.0) / 4.0) * dt);
		}
		const double po = payoff(S);
		paths_stats.add(po);
		sample_paths += 1;
		double error = exp(-r * m_T) * paths_stats.stDev() / sqrt(sample_paths - 1.0);
		if (error <= req_error && sample_paths > 10)
		{
			done = true;
		} 
	}
	Results results;
	results.price = exp(-r * m_T) * paths_stats.mean();
	results.error = exp(-r * m_T) * paths_stats.stDev() / sqrt(sample_paths - 1);
	results.N = sample_paths;
	return results;
}

double EurCall::payoff(double S) const
{
	return S - m_K > 0 ? S - m_K : 0.0;
}

double EurCall::PriceByBSFormula(double S0, double sigma, double r) const
{
	return S0 * N(d_plus(S0, sigma, r))
		- m_K * exp(-r * m_T) * N(d_minus(S0, sigma, r));
}

double EurCall::N(double x) const
{
	double gamma = 0.2316419;     double a1 = 0.319381530;
	double a2 = -0.356563782;   double a3 = 1.781477937;
	double a4 = -1.821255978;   double a5 = 1.330274429;
	double pi = 4.0 * atan(1.0); double k = 1.0 / (1.0 + gamma * x);
	if (x >= 0.0)
	{
		return 1.0 - ((((a5 * k + a4) * k + a3) * k + a2) * k + a1)
			* k * exp(-x * x / 2.0) / sqrt(2.0 * pi);
	}
	else return 1.0 - N(-x);
}
double EurCall::d_plus(double S0, double sigma, double r) const
{
	return (log(S0 / m_K) +
		(r + 0.5 * pow(sigma, 2.0)) * m_T)
		/ (sigma * sqrt(m_T));
}
double EurCall::d_minus(double S0, double sigma, double r) const
{
	return d_plus(S0, sigma, r) - sigma * sqrt(m_T);
}

double BisectionSolver(const Function* f, double c, double left, double right, double acc) 
{
	double fl = f->eval(left);
	double fr = f->eval(right);
		while (true)
	{
		const double mid = (left + right) / 2.0;
		if (right - left <= acc)
			return mid;
		const double fm = f->eval(mid);
		if (fl <= c && c <= fm || fm <= c && c <= fl)
		{
			right = mid;
			fr = fm;
		}
		else
		{
			left = mid;
			fl = fm;
		}
	}

	return 0.0;
}

long long fact(int n)
{
	if ((n == 0) || (n == 1))
		return 1;
	else
		return n * fact(n - 1);
}

double EurCall::AnalyticalPriceByJumpDiffusion(double S0, double r, double sigma, double lambda, double m, double s) const
{
	const double k = exp(m + (pow(s, 2) / 2)) - 1;
	const double lambda_prime = lambda * (1 + k);
	
	double sum = 0;
	double n = 0;

	while (true)
	{
		double r_n = r - (lambda * k) + ((n * log(1 + k)) / m_T);
		double sigma_n = sqrt(pow(sigma, 2) + ((n * pow(s, 2)) / m_T));

		double mult_factor = (exp(-lambda_prime * m_T) * pow(lambda_prime * m_T, n)) / fact(n);

		if (mult_factor < .00000000000001)
		{
			return sum;
		}

		sum += mult_factor * PriceByBSFormula(S0, sigma_n, r_n);
		n += 1;
	}
} 


double Pdf(double x, double stdev, double mean, double strike)
{
	//pdf of the function
	if (x == 0)
	{
		return 0;
	}
	double pi = 4.0 * atan(1.0);
	double first_part = 1.0 / (x * stdev * sqrt(2.0 * pi));
	double numer = pow((log(x) - mean), 2.0);
	double denom = 2 * stdev * stdev;
	double second_part = exp(-numer / denom);
	double pdf = first_part * second_part;
	//payoff of the function
	double payoff = 0;
	if (x - strike > 0)
	{
		payoff = x - strike;
	}
	return payoff * pdf;
}

double Cdf(double stdev, double mean, double strike, double a, double b, double N)
{
	double h = (b - a) / N;
	double sum = 0;
	for (double i = 0; i < N + 1; ++i)
	{
		double multiplier = 2.0;
		if (i == 0 || i == N)
		{
			multiplier = 1.0;
		}
		sum += multiplier * Pdf((a + (i * h)), stdev, mean, strike);
	}
	return (h / 2.0) * sum;
}

double InverseCdf(double stdev, double mean, double strike, double a, double b, double N)
{
	return 1 / Cdf(stdev, mean, strike, a, b, N);
}

double Nprime(double x)
{
	double pi = 4.0 * atan(1.0);
	return exp(-pow(x, 2.0) / 2.0) / sqrt(2.0 * pi);
}

double AbsoluteValue(double x)
{
	return x > 0 ? x : -x;
}

int random(int a, int b)
{
	srand(time(NULL));
	int r = rand() % 2;
	if (r == 0)
		return a;
	else
		return b;
}

Results EurCall::MCPriceByJumpDiffusion(double S0, double r, double sigma, double lambda, double m, double s, double time_steps, double req_error) const
{
	double dt = m_T / time_steps;
	double drift = r - pow(sigma, 2) / 2 - (lambda * (exp(m + pow(s, 2) / 2) - 1));

	NormalRandomGenerator nrg;



	int Nsteps = 1;

	double mean = log(S0) + ((r - (.5 * pow(sigma, 2.0))) * m_T);
	double stdev = sigma * sqrt(m_T);



	Stats1 avgpo;

	Results results;
	while (true)
	{
		double S = S0;
		double J = 1;
		for (int i = 1; i <= time_steps; ++i)
		{
			double u = nrg.generate();
			double z1 = nrg.generate();
			if (u <= lambda * dt) //we have a jump
			{
				//J = AbsoluteValue(u);
				//J = Nprime(u);
				//S = J>.20? S * (1 + J) : S * (1-J);
				int a = .5;
				int b = 1.5;
				int direction = random(a, b);
				S = S * direction;
			}
			///cout << J << endl;
			//cout << S << "\t" << J << endl;
			//cout << exp(drift * dt + sigma * sqrt(dt) * z1) << endl;
			S *= exp(drift * dt + sigma * sqrt(dt) * z1);
		}
		cout << S << "\t" << payoff(S) << endl;
		avgpo.add(payoff(S));
		++Nsteps;
		double error = exp(-r * m_T) * avgpo.stDev() / sqrt(Nsteps - 1);
		cout << error << "\t" << Nsteps << endl;
		if (error <= req_error && Nsteps > 10)
		{
			results.price = avgpo.mean();
			results.error = error;
			results.N = Nsteps;
			return results;
		}
	}
}
