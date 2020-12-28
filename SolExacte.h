#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_EXACTE_
#define SOL_EXACTE_

using namespace std;

class SolExacte {
public:
	PRP* instance;

	Solution solution; 

	vector<vector<double>> w_sol;

	SolExacte(PRP*);

	void solve(bool verbose = false);
};

#endif
