#include<vector>

#include "PRP.h"
#include "Solution.h"

#ifndef SOL_APPROCHEE_
#define SOL_APPROCHEE_

using namespace std;

class SolApprochee {
public:
	PRP* instance; // instance du probl�me

	Solution meilleure, courante; // solutions : meilleure trouv�e jusqu'� une �tape et courante

	vector<vector<double>> SC; // cout heuristique de la visite du client i � l'instant t

	SolApprochee(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void main_loop(int max_iter, bool verbose = false);
};

#endif