#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_APPROCHEE_HEURISTIQUE_
#define SOL_APPROCHEE_HEURISTIQUE_

using namespace std;

class SolApprocheeHeuristique {
public:
	PRP* instance; // instance du problème

	Solution meilleure, courante; // solutions : meilleure trouvée jusqu'à une étape et courante

	vector<vector<double>> SC; // cout heuristique de la visite du client i à l'instant t

	SolApprocheeHeuristique(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_LSP_bis(bool verbose = false);
	void solve_VRP_heuristique(int t, int nb_steps_optim, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void solve(int max_iter, int nb_steps_optim = 1000, bool verbose = false);
};

#endif