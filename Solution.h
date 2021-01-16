#include <vector>
#include "PRP.h"

#ifndef SOLUTION_
#define SOLUTION_

using namespace std;

typedef vector<vector<int>> Graph; // Graphe représenté par liste de successeurs

class Solution {
public:
	int n, l; // quantité de clients et nombre de pas de temps

	double valeur; // valeur de la solution

	vector<double> p; // quantité produite à chaque instant (vecteur de temps)
	vector<bool> y; // booléan indiquant si la production a été lancé à chaque instant (vecteur de temps)
	vector<vector<double>> I; // état du stock : dépend de l'indice du client / fournisseur et du temps
	vector<vector<double>> q; // quantité envoyé à chaque client : dépend de l'indice du client et du temps
	vector<vector<bool>> z; // booléan indiquant si un envoi est fait au client i à l'instant t
	vector<Graph> x; // graphes avec les tournées de véhicules (vecteur de temps)

	Solution(int nn, int ll);

	void calcul_valeur(PRP& inst);

	friend ostream& operator<<(ostream& os, const Solution& solution);
};

#endif