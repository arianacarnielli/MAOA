#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_APPROCHEE_COUPE_
#define SOL_APPROCHEE_COUPE_

using namespace std;

class SolApprocheeCoupe {
public:
	PRP* instance; // instance du probl�me

	Solution meilleure, courante; // solutions : meilleure trouv�e jusqu'� une �tape et courante

	vector<vector<double>> SC; // cout heuristique de la visite du client i � l'instant t

	SolApprocheeCoupe(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, double time_limit = -1, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void solve(int max_iter, double VSP_time_limit = -1, bool verbose = false);
};

#endif