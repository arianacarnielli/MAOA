#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_APPROCHEE_
#define SOL_APPROCHEE_

using namespace std;

class SolApprochee {
public:
	PRP* instance; // instance du problème

	Solution meilleure, courante; // solutions : meilleure trouvée jusqu'à une étape et courante

	vector<vector<double>> SC; // cout heuristique de la visite du client i à l'instant t

	SolApprochee(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void main_loop(int max_iter, bool verbose = false);
};

#endif