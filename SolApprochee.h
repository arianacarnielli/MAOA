#include<vector>

#include "PRP.h"

#ifndef SOL_APPROCHEE_
#define SOL_APPROCHEE_

typedef vector<vector<int>> Graph; // Graphe représenté par liste de successeurs

using namespace std;

class SolApprochee {
public:
	PRP* instance;
	vector<vector<double>> SC; // cout heuristique de la visite du client i à l'instant t

	// Valeurs des variables de décision du problème LSP à la fin
	vector<double> p_sol; // quantité produite à chaque instant (vecteur de temps)
	vector<bool> y_sol; // booléan indiquant si la production a été lancé à chaque instant (vecteur de temps)
	vector<vector<double>> I_sol; // état du stock : dépend de l'indice du client / fournisseur et du temps
	vector<vector<double>> q_sol; // quantité envoyé à chaque client : dépend de l'indice du client et du temps
	vector<vector<bool>> z_sol; // booléan indiquant si un envoi est fait au client i à l'instant t

	// Valeurs des variables de décision du problème VRP
	vector<Graph> x_sol; // graphes avec les tournées de véhicules (vecteur de temps)

	SolApprochee(PRP*);

	void init_SC();
	void solve_LSP(bool verbose = false);
	void solve_VRP_MTZ(int t, bool verbose = false);
	void calcul_SC(int t, bool verbose = false);
	void main_loop(int max_iter, bool verbose = false);
};

#endif