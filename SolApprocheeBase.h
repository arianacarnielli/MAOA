#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_APPROCHEE_BASE_
#define SOL_APPROCHEE_BASE_

using namespace std;

class SolApprocheeBase {
public:
	PRP* instance; // instance du problème

	Solution meilleure, courante; // solutions : meilleure trouvée jusqu'à une étape et courante

	vector<vector<double>> SC; // cout heuristique de la visite du client i à l'instant t

	SolApprocheeBase(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void solve(int max_iter, bool verbose = false);
};

#endif