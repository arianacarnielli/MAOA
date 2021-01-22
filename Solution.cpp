#include <limits>
#include <iostream>
#include <cstdio>

#include "PRP.h"
#include "Solution.h"

Solution::Solution() {}

Solution::Solution(int nn, int ll) {
	n = nn;
	l = ll;

	valeur = numeric_limits<double>::infinity();

	p.resize(l);
	y.resize(l);
	I.resize(n + 1);
	q.resize(n + 1);
	z.resize(n + 1);
	x.resize(l);

	for (int t = 0; t < l; t++) {
		x[t].resize(n + 1);
	}

	I[0].resize(l);
	z[0].resize(l);
	for (int i = 1; i <= n; i++) {
		I[i].resize(l);
		q[i].resize(l);
		z[i].resize(l);
	}
}

void Solution::calcul_valeur(PRP& inst) {
	valeur = 0;
	for (int t = 0; t < l; t++) {
		valeur += inst.u * p[t];
		valeur += inst.f * y[t];
		for (int i = 0; i <= n; i++) {
			valeur += inst.h[i] * I[i][t];
			vector<int> succs = x[t][i];
			for (int j = 0; j < succs.size(); j++) {
				valeur += inst.cost(i, succs[j]);
			}
		}
	}
}

ostream& operator<<(ostream& os, const Solution& solution) {
	char buffer[200];
	os << "Instance avec " << solution.n << " clients et " << solution.l << " pas de temps" << endl;

	for (int t = 0; t < solution.l; t++) {
		os << "Pas de temps " << t << " :" << endl;
		os << "  Production : p_" << t << " = " << solution.p[t] << endl;
		os << "  Production lancee ? y_" << t << " = " << solution.y[t] << endl;
		os << "  Stock a l'usine : I_0_" << t << " = " << solution.I[0][t] << endl;
		for (int i = 1; i <= solution.n; i++) {
			sprintf_s(buffer, 200, "  Client %03d : I_%03d_%02d = %10.3f ; q_%03d_%02d = %10.3f ; z_%03d_%02d = %1d",
				i, i, t, solution.I[i][t], i, t, solution.q[i][t], i, t, solution.z[i][t]);
			os << buffer << endl;
		}
		os << "  Tournees : ";
		if (!solution.x[t][0].size()) {
			os << "aucune" << endl;
		}
		else {
			os << endl;
			int v;
			for (int j = 0; j < solution.x[t][0].size(); j++) {
				os << "  0 --> ";
				v = solution.x[t][0][j];
				while (v) {
					os << v << " --> ";
					v = solution.x[t][v][0];
				}
				os << "0" << endl;
			}
		}
	}
	os << "Valeur de la solution : " << solution.valeur << endl;
	return os;
}