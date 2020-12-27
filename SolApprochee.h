#include<vector>

#include "PRP.h"

#ifndef SOL_APPROCHEE_
#define SOL_APPROCHEE_

using namespace std;

class SolApprochee {
public:
	PRP* instance;
	vector<vector<double>> SC;

	// Valeurs des variables de décision du problème LSP à la fin
	vector<double> p_sol;
	vector<bool> y_sol;
	vector<vector<double>> I_sol;
	vector<vector<double>> q_sol;
	vector<vector<bool>> z_sol;

	SolApprochee(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, bool verbose = false);
};

#endif