#pragma once

#include<vector>

#include "PRP.h"
#include "SolutionV2.h"

#ifndef SOL_EXACTE_V2_
#define SOL_EXACTE_V2_

using namespace std;

class SolExacteV2 {
public:
	PRP* instance;

	SolutionV2 solution;

	vector<vector<double>> w_sol;

	SolExacteV2(PRP*);

	void solve(SolutionV2* sol_init = nullptr, double tolerance = 0.01, bool verbose = false);
};

#endif
