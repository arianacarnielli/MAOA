#include<vector>

#include "PRP.h"

#ifndef SOL_APPROCHEE_
#define SOL_APPROCHEE_

typedef vector<vector<int>> Graph; // Graphe repr�sent� par liste de successeurs

using namespace std;

class SolApprochee {
public:
	PRP* instance;
	vector<vector<double>> SC; // cout heuristique de la visite du client i � l'instant t

	// Valeurs des variables de d�cision du probl�me LSP � la fin
	vector<double> p_sol; // quantit� produite � chaque instant (vecteur de temps)
	vector<bool> y_sol; // bool�an indiquant si la production a �t� lanc� � chaque instant (vecteur de temps)
	vector<vector<double>> I_sol; // �tat du stock : d�pend de l'indice du client / fournisseur et du temps
	vector<vector<double>> q_sol; // quantit� envoy� � chaque client : d�pend de l'indice du client et du temps
	vector<vector<bool>> z_sol; // bool�an indiquant si un envoi est fait au client i � l'instant t

	// Valeurs des variables de d�cision du probl�me VRP
	vector<Graph> x_sol; // graphes avec les tourn�es de v�hicules (vecteur de temps)

	SolApprochee(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void main_loop(int max_iter, bool verbose = false);
};

#endif