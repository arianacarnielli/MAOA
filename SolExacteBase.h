#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_EXACTE_BASE_
#define SOL_EXACTE_BASE_

using namespace std;

class SolExacteBase {
public:
	PRP* instance;

	Solution solution; 

	vector<vector<double>> w_sol;

	SolExacteBase(PRP*);

	void solve(Solution* sol_init = nullptr, double tolerance = 0.01, double time_limit = -1, bool verbose = false);
};

#endif
