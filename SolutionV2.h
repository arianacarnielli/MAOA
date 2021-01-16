#include <vector>
#include "PRP.h"

#ifndef SOLUTIONV2_
#define SOLUTIONV2_

using namespace std;

typedef vector<vector<int>> Graph; // Graphe repr�sent� par liste de successeurs

class SolutionV2 {
public:
	int n, l, m; // quantit� de clients et nombre de pas de temps

	double valeur; // valeur de la solution

	vector<double> p; // quantit� produite � chaque instant (vecteur de temps)
	vector<bool> y; // bool�an indiquant si la production a �t� lanc� � chaque instant (vecteur de temps)
	vector<vector<double>> I; // �tat du stock : d�pend de l'indice du client / fournisseur et du temps
	vector<vector<vector<double>>> q; // quantit� envoy� � chaque client : d�pend de l'indice du client et du temps
	vector<vector<vector<bool>>> z; // bool�an indiquant si un envoi est fait au client i � l'instant t
	vector < vector < vector < vector<bool>>>> x; // graphes avec les tourn�es de v�hicules (vecteur de temps)

	Solution(int nn, int ll, int mm);

	void calcul_valeur(PRP& inst);
};

#endif
#pragma once
