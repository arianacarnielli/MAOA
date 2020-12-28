#include<vector>

#include "PRP.h"

#ifndef SOL_EXACTE_
#define SOL_EXACTE_

using namespace std;

class SolExacte {
public:
	PRP* instance;

	// Valeurs des variables de décision du problème LSP à la fin
	vector<double> p_sol;
	vector<bool> y_sol;
	vector<vector<double>> I_sol;
	vector<vector<double>> q_sol;
	vector<vector<bool>> z_sol;
	vector<vector<double>> w_sol;
	vector<vector<vector<bool>>> x_sol;

	SolExacte(PRP*);

	void solve(bool verbose = false);
};

#endif#pragma once
