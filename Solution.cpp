#include <limits>
#include <iostream>

#include "PRP.h"
#include "Solution.h"

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