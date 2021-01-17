#include <limits>
#include <iostream>

#include "PRP.h"
#include "SolutionV2.h"

SolutionV2::Solution(int nn, int ll, int mm) {
	n = nn;
	l = ll;
	m = mm;

	valeur = numeric_limits<double>::infinity();

	p.resize(l);
	y.resize(l);
	I.resize(n + 1);
	q.resize(n + 1);
	z.resize(n + 1);

	x.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		x[i].resize(n + 1);
		for (int j = 0; j <= n; j++) {
			x[i][j].resize(m);
			for (int k = 0; k <= m; k++) {
				x[i][j][k].resize(l);
			}
		}
	}


	I[0].resize(l);
	z[0].resize(l);
	for (int i = 1; i <= n; i++) {
		I[i].resize(l);
		q[i].resize(m);
		z[i].resize(m);
		for (int j = 1; j <= m; j++) {
			q[i][j].resize(l);
			z[i][j].resize(l);
		}
	}
}

void SolutionV2::calcul_valeur(PRP& inst) {
	valeur = 0;
	for (int t = 0; t < l; t++) {
		valeur += inst.u * p[t];
		valeur += inst.f * y[t];
		for (int i = 0; i <= n; i++) {
			valeur += inst.h[i] * I[i][t];
			for (int j = 0; j <= n; j++) {
				if (j != i) {
					for (int k = 0; k < m; k++) {
						valeur += x[i][j][k][t] * inst.cost(i, j);
					}
				}
			}
		}
	}
}