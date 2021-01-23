#pragma once

#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_EXACTE_V2_
#define SOL_EXACTE_V2_

using namespace std;

class SolExacteV2 {
public:
	PRP* instance;

	Solution solution;

	SolExacteV2(PRP*);

	void solve(Solution* sol_init = nullptr, double tolerance = 0.01, double time_limit = -1, bool verbose = false);
};

#endif
