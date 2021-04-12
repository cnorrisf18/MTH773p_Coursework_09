#include "Options.h"
#include <iostream>
using namespace std;

int main()
{
	const double T = 1.0;
	const double K = 105.0;

	EurCall example_call(T, K);

	const double S0 = 100.0;
	const double r = 0.05;
	const double v0 = pow(0.2, 2);
	const double vsquare = pow(0.3, 2);
	const double a = 1.25;
	const double vol_of_vol = 0.3;
	const double p = 0.1;
	const int time_steps = 365;

	const double req_error = 0.1;
	int N = 200'000;

	//Results example_results1 = example_call.HestonMonteCarloPricewithN(S0, r, a, vsquare, vol_of_vol, v0, p, 365, N);
	cout << "Done with part 1" << endl;
	//Results example_results2 = example_call.HestonMonteCarloPricewithError(S0, r, a, vsquare, vol_of_vol, v0, p, 365, req_error);

	//cout << example_results1.price << "\t" << example_results1.error << "\t" << example_results1.N << endl;
	//cout << example_results2.price << "\t" << example_results2.error << "\t" << example_results2.N << endl;

	//const Adapter adapter(&example_call, S0, r);
	//const double impliedVol = BisectionSolver(&adapter, example_results1.price, 0.0, 1.0, .01);

	//cout << "Implied vol: " << impliedVol << endl;

	EurCall call2(0.5, 105);
	Results example_results3 = call2.MCPriceByJumpDiffusion(100, 0.05, 0.15, 1, 0.05, 0.3, time_steps, 0.1);

	cout << call2.AnalyticalPriceByJumpDiffusion(100, 0.05, 0.15, 1, 0.05, 0.3) << endl;
	cout << example_results3.price << "\t" << example_results3.error << "\t" << example_results3.N << endl;
	
	return 0;
}