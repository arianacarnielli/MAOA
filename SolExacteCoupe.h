#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_EXACTE_COUPE_
#define SOL_EXACTE_COUPE_

using namespace std;

class SolExacteCoupe {
public:
	PRP* instance;

	Solution solution; 

	vector<vector<double>> w_sol;

	SolExacteCoupe(PRP*);

	void solve(Solution* sol_init = nullptr, double tolerance = 1e-4, double time_limit = -1, string coupe = "G1", bool verbose = false);
};

#endif
