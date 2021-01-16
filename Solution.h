#include <vector>
#include "PRP.h"

#ifndef SOLUTION_
#define SOLUTION_

using namespace std;

typedef vector<vector<int>> Graph; // Graphe repr�sent� par liste de successeurs

class Solution {
public:
	int n, l; // quantit� de clients et nombre de pas de temps

	double valeur; // valeur de la solution

	vector<double> p; // quantit� produite � chaque instant (vecteur de temps)
	vector<bool> y; // bool�an indiquant si la production a �t� lanc� � chaque instant (vecteur de temps)
	vector<vector<double>> I; // �tat du stock : d�pend de l'indice du client / fournisseur et du temps
	vector<vector<double>> q; // quantit� envoy� � chaque client : d�pend de l'indice du client et du temps
	vector<vector<bool>> z; // bool�an indiquant si un envoi est fait au client i � l'instant t
	vector<Graph> x; // graphes avec les tourn�es de v�hicules (vecteur de temps)

	Solution(int nn, int ll);

	void calcul_valeur(PRP& inst);

	friend ostream& operator<<(ostream& os, const Solution& solution);
};

#endif