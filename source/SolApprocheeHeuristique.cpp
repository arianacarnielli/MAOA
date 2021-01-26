#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <random>
#include <limits>

#include "SolApprocheeHeuristique.h"
#include "PRP.h"

typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

extern default_random_engine alea;

SolApprocheeHeuristique::SolApprocheeHeuristique(PRP* inst) :
	meilleure(inst->n, inst->l), 
	courante(inst->n, inst->l) {
	instance = inst;
	int n = instance->n;
	int l = instance->l;

	SC.resize(n + 1);
	for (int i = 1; i <= n; i++) {
		SC[i].resize(l);
	}
}

void SolApprocheeHeuristique::init_SC() {
	for (int i = 1; i <= instance->n; i++) {
		SC[i][0] = 2 * (instance->cost(0, i));
		for (int t = 1; t < instance->l; t++) {
			SC[i][t] = SC[i][0];
		}
	}
}

void SolApprocheeHeuristique::solve_LSP(bool verbose) {
	int n = instance->n;
	int l = instance->l;

	vector<double> M;
	vector<vector<double>> tildeM;
	M.resize(l);
	tildeM.resize(n + 1);
	for (int i = 1; i <= n ; i++) tildeM[i].resize(l);

	vector<double> part_sum;
	part_sum.resize(n + 1);

	// on initialise les valeurs de M et tildeM
	for (int i = 1; i <= n; i++) part_sum[i] = 0;
	double sum;
	for (int t = l - 1; t >= 0; t--) {
		sum = 0;
		for (int i = 1; i <= n; i++) {
			part_sum[i] += instance->d[i][t];
			tildeM[i][t] = min({ instance->L[i], instance->Q, part_sum[i] });
			sum += part_sum[i];
		}
		M[t] = min(instance->C, sum);
	}


	//////////////
	//////  CPLEX INITIALIZATION
	//////////////

	IloEnv   env;
	IloModel model(env);

	////////////////////////
	//////  VAR
	////////////////////////
	
	vector<IloNumVar> p;
	vector<IloNumVar> y;
	vector<vector<IloNumVar>> I;
	vector<vector<IloNumVar>> q;
	vector<vector<IloNumVar>> z;

	p.resize(l);
	y.resize(l);
	I.resize(n + 1);
	q.resize(n + 1);
	z.resize(n + 1);

	I[0].resize(l);
	for (int i = 1; i <= n; i++) {
		I[i].resize(l);
		q[i].resize(l);
		z[i].resize(l);
	}

	ostringstream varname;
	for (int t = 0; t < l; t++) {

		p[t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
		varname.str("");
		varname << "p_" << t;
		p[t].setName(varname.str().c_str());

		y[t] = IloNumVar(env, 0.0, 1.0, ILOINT);
		varname.str("");
		varname << "y_" << t;
		y[t].setName(varname.str().c_str());

		I[0][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
		varname.str("");
		varname << "I_0_" << t;
		I[0][t].setName(varname.str().c_str());

		for (int i = 1; i <= n; i++) {
			I[i][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
			varname.str("");
			varname << "I_" << i << "_" << t;
			I[i][t].setName(varname.str().c_str());

			q[i][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
			varname.str("");
			varname << "q_" << i << "_" << t;
			q[i][t].setName(varname.str().c_str());

			z[i][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
			varname.str("");
			varname << "z_" << i << "_" << t;
			z[i][t].setName(varname.str().c_str());
		}
	}
	//////////////
	//////  CONTRAINTES
	//////////////

	IloRangeArray CC(env);
	int nbcst = 0;
	ostringstream cstname;

	//Contrainte 1 : I0,t-1 + pt = sum{i}(qi,t) + I0,t
	// il faut traiter le cas t = 0 separement
	IloExpr cst(env);
	for (int i = 1; i <= n; i++) {
		cst += q[i][0];
	}
	cst += I[0][0] - p[0];
	CC.add(cst == instance->L0[0]);

	cstname.str("");
	cstname << "Cst_onecol_0";
	CC[nbcst].setName(cstname.str().c_str());
	nbcst++;

	for (int t = 1; t < l; t++) {
		IloExpr cst(env);
		for (int i = 1; i <= n; i++) {
			cst += q[i][t];
		}
		cst += I[0][t] - I[0][t - 1] - p[t];
		CC.add(cst == 0);

		cstname.str("");
		cstname << "Cst_onecol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 2 : 
	// il faut traiter le cas t = 0 separement
	for (int i = 1; i <= n; i++) {
		IloExpr cst(env);
		cst += q[i][0] - I[i][0];
		CC.add(cst == instance->d[i][0] - instance->L0[i]);
		cstname.str("");
		cstname << "Cst_twocol_" << i << "_line_0";
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	for (int t = 1; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			cst += I[i][t - 1] + q[i][t] - I[i][t];
			CC.add(cst == instance->d[i][t]);
			cstname.str("");
			cstname << "Cst_twocol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//contrainte 3 : 
	for (int t = 0; t < l; t++) {
		IloExpr cst(env);
		cst += M[t] * y[t] - p[t];
		CC.add(cst >= 0);
		cstname.str("");
		cstname << "Cst_threecol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 4 : 
	for (int t = 0; t < l; t++) {
		IloExpr cst(env);
		cst += I[0][t];
		CC.add(cst <= instance->L[0]);
		cstname.str("");
		cstname << "Cst_fourcol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 5 : 
	// il faut traiter le cas t = 0 separement
	for (int i = 1; i <= n; i++) {
		IloExpr cst(env);
		cst += q[i][0];
		CC.add(cst <= instance->L[i] - instance->L0[i]);
		cstname.str("");
		cstname << "Cst_fivecol_" << i << "_line_0";
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	for (int t = 1; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			cst += I[i][t - 1] + q[i][t];
			CC.add(cst <= instance->L[i]);
			cstname.str("");
			cstname << "Cst_fivecol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//contrainte 6 (les contraintes couplantes entre q et z): 
	for (int t = 0; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			cst += tildeM[i][t] * z[i][t] - q[i][t];
			CC.add(cst >= 0);
			cstname.str("");
			cstname << "Cst_sixcol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	// Contrainte de charge maximale
	for (int t = 0; t < l; t++) {
		IloExpr cst(env);
		for (int i = 1; i <= n; i++) {
			cst += q[i][t];
		}
		CC.add(cst <= instance->Q * instance->m);
		cstname.str("");
		cstname << "Cst_charge_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	model.add(CC);

	//////////////
	////// OBJ
	//////////////

	IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));

	for (int t = 0; t < l; t++) {

		obj.setLinearCoef(p[t], instance->u);
		obj.setLinearCoef(y[t], instance->f);
		obj.setLinearCoef(I[0][t], instance->h[0]);

		for (int i = 1; i <= n; i++) {
			obj.setLinearCoef(I[i][t], instance->h[i]);

			//avec le cout de visite pour l'heuristique de la partie 1 :
			obj.setLinearCoef(z[i][t], SC[i][t]);
		}
	}

	///////////
	//// RESOLUTION
	//////////

	IloCplex cplex(model);
	if (!verbose) {
		cplex.setOut(env.getNullStream());
	}
	cplex.solve();

	for (int t = 0; t < l; t++) {
		courante.p[t] = cplex.getValue(p[t]);
		courante.y[t] = cplex.getValue(y[t]);
		courante.I[0][t] = cplex.getValue(I[0][t]);
		for (int i = 1; i <= n; i++) {
			courante.I[i][t] = cplex.getValue(I[i][t]);
			courante.q[i][t] = cplex.getValue(q[i][t]);
			courante.z[i][t] = cplex.getValue(z[i][t]);
		}
	}

	env.end();
}

void SolApprocheeHeuristique::solve_LSP_bis(bool verbose) {
	int n = instance->n;
	int l = instance->l;
	int K = min(instance->n, instance->m);


	vector<double> M;
	vector<vector<double>> tildeM;
	M.resize(l);
	tildeM.resize(n + 1);
	for (int i = 1; i <= n; i++) tildeM[i].resize(l);

	vector<double> part_sum;
	part_sum.resize(n + 1);

	// on initialise les valeurs de M et tildeM
	for (int i = 1; i <= n; i++) part_sum[i] = 0;
	double sum;
	for (int t = l - 1; t >= 0; t--) {
		sum = 0;
		for (int i = 1; i <= n; i++) {
			part_sum[i] += instance->d[i][t];
			tildeM[i][t] = min({ instance->L[i], instance->Q, part_sum[i] });
			sum += part_sum[i];
		}
		M[t] = min(instance->C, sum);
	}


	//////////////
	//////  CPLEX INITIALIZATION
	//////////////

	IloEnv   env;
	IloModel model(env);

	////////////////////////
	//////  VAR
	////////////////////////

	vector<IloNumVar> p;
	vector<IloNumVar> y;
	vector<vector<IloNumVar>> I;
	vector<vector<vector<IloNumVar>>> q;
	vector<vector<vector<IloNumVar>>> z;

	p.resize(l);
	y.resize(l);
	I.resize(n + 1);
	q.resize(n + 1);
	z.resize(n + 1);

	I[0].resize(l);
	for (int i = 1; i <= n; i++) {
		I[i].resize(l);
		q[i].resize(K);
		z[i].resize(K);
		for (int k = 0; k < K; k++) {
			q[i][k].resize(l);
			z[i][k].resize(l);
		}
	}

	ostringstream varname;
	for (int t = 0; t < l; t++) {

		p[t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
		varname.str("");
		varname << "p_" << t;
		p[t].setName(varname.str().c_str());

		y[t] = IloNumVar(env, 0.0, 1.0, ILOINT);
		varname.str("");
		varname << "y_" << t;
		y[t].setName(varname.str().c_str());

		I[0][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
		varname.str("");
		varname << "I_0_" << t;
		I[0][t].setName(varname.str().c_str());

		for (int i = 1; i <= n; i++) {
			I[i][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
			varname.str("");
			varname << "I_" << i << "_" << t;
			I[i][t].setName(varname.str().c_str());

			for (int k = 0; k < K; k++) {
				q[i][k][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
				varname.str("");
				varname << "q_" << i << "_" << k << "_" << t;
				q[i][k][t].setName(varname.str().c_str());

				z[i][k][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
				varname.str("");
				varname << "z_" << i << "_" << k << "_" << t;
				z[i][k][t].setName(varname.str().c_str());
			}
		}
	}
	//////////////
	//////  CONTRAINTES
	//////////////

	IloRangeArray CC(env);
	int nbcst = 0;
	ostringstream cstname;

	//Contrainte 1 : I0,t-1 + pt = sum{i}(qi,t) + I0,t
	// il faut traiter le cas t = 0 separement
	IloExpr cst(env);
	for (int i = 1; i <= n; i++) {
		for (int k = 0; k < K; k++) {
			cst += q[i][k][0];
		}
	}
	cst += I[0][0] - p[0];
	CC.add(cst == instance->L0[0]);

	cstname.str("");
	cstname << "Cst_onecol_0";
	CC[nbcst].setName(cstname.str().c_str());
	nbcst++;

	for (int t = 1; t < l; t++) {
		IloExpr cst(env);
		for (int i = 1; i <= n; i++) {
			for (int k = 0; k < K; k++) {
				cst += q[i][k][t];
			}
		}
		cst += I[0][t] - I[0][t - 1] - p[t];
		CC.add(cst == 0);

		cstname.str("");
		cstname << "Cst_onecol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 2 : 
	// il faut traiter le cas t = 0 separement
	for (int i = 1; i <= n; i++) {
		IloExpr cst(env);
		for (int k = 0; k < K; k++) {
			cst += q[i][k][0];
		}
		cst -= I[i][0];
		CC.add(cst == instance->d[i][0] - instance->L0[i]);
		cstname.str("");
		cstname << "Cst_twocol_" << i << "_line_0";
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	for (int t = 1; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			for (int k = 0; k < K; k++) {
				cst += q[i][k][t];
			}
			cst += I[i][t - 1] - I[i][t];
			CC.add(cst == instance->d[i][t]);
			cstname.str("");
			cstname << "Cst_twocol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//contrainte 3 : 
	for (int t = 0; t < l; t++) {
		IloExpr cst(env);
		cst += M[t] * y[t] - p[t];
		CC.add(cst >= 0);
		cstname.str("");
		cstname << "Cst_threecol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 4 : 
	for (int t = 0; t < l; t++) {
		IloExpr cst(env);
		cst += I[0][t];
		CC.add(cst <= instance->L[0]);
		cstname.str("");
		cstname << "Cst_fourcol_" << t;
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 5 : 
	// il faut traiter le cas t = 0 separement
	for (int i = 1; i <= n; i++) {
		IloExpr cst(env);
		for (int k = 0; k < K; k++) {
			cst += q[i][k][0];
		}
		CC.add(cst <= instance->L[i] - instance->L0[i]);
		cstname.str("");
		cstname << "Cst_fivecol_" << i << "_line_0";
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	for (int t = 1; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			for (int k = 0; k < K; k++) {
				cst += q[i][k][t];
			}
			cst += I[i][t - 1];
			CC.add(cst <= instance->L[i]);
			cstname.str("");
			cstname << "Cst_fivecol_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	//contrainte 6 (les contraintes couplantes entre q et z): 
	for (int t = 0; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			for (int k = 0; k < K; k++) {
				IloExpr cst(env);
				cst += tildeM[i][t] * z[i][k][t] - q[i][k][t];
				CC.add(cst >= 0);
				cstname.str("");
				cstname << "Cst_sixcol_" << i << "_line_" << t;
				CC[nbcst].setName(cstname.str().c_str());
				nbcst++;
			}
		}
	}

	// Contraintes d'un vehicule par client
	for (int t = 0; t < l; t++) {
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			for (int k = 0; k < K; k++) {
				cst += z[i][k][t];
			}
			CC.add(cst <= 1);
			cstname.str("");
			cstname << "Cst_vehicule_" << i << "_line_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	// Contraintes de charge maximale d'un vehicule
	for (int t = 0; t < l; t++) {
		for (int k = 0; k < K; k++) {
			IloExpr cst(env);
			for (int i = 1; i <= n; i++) {
				cst += q[i][k][t];
			}
			CC.add(cst <= instance->Q);
			cstname.str("");
			cstname << "Cst_charge_" << t << "_" << k;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}
	}

	model.add(CC);

	//////////////
	////// OBJ
	//////////////

	IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));

	for (int t = 0; t < l; t++) {

		obj.setLinearCoef(p[t], instance->u);
		obj.setLinearCoef(y[t], instance->f);
		obj.setLinearCoef(I[0][t], instance->h[0]);

		for (int i = 1; i <= n; i++) {
			obj.setLinearCoef(I[i][t], instance->h[i]);

			//avec le cout de visite pour l'heuristique de la partie 1 :
			for (int k = 0; k < K; k++) {
				obj.setLinearCoef(z[i][k][t], SC[i][t]);
			}
		}
	}

	///////////
	//// RESOLUTION
	//////////

	IloCplex cplex(model);
	if (!verbose) {
		cplex.setOut(env.getNullStream());
	}
	cplex.solve();

	for (int t = 0; t < l; t++) {
		courante.p[t] = cplex.getValue(p[t]);
		courante.y[t] = cplex.getValue(y[t]);
		courante.I[0][t] = cplex.getValue(I[0][t]);
		for (int i = 1; i <= n; i++) {
			courante.I[i][t] = cplex.getValue(I[i][t]);
			courante.q[i][t] = 0;
			courante.z[i][t] = false;
			for (int k = 0; k < K; k++) {
				courante.q[i][t] += cplex.getValue(q[i][k][t]);
				courante.z[i][t] = courante.z[i][t] || (cplex.getValue(z[i][k][t]) > 0.5);
			}
		}
	}

	env.end();
}

vector<int> tournee_vois_2TSP(vector<int> tournee, int u, int v) {
	// Hypothèses : en notant n = tournee.size() :
	// n >= 3
	// -1 <= u <= n-3
	// u+2 <= v <= n-1

	vector<int> new_tournee;

	int n = tournee.size();
	int i;
	for (i = 0; i <= u; i++) {
		new_tournee.push_back(tournee[i]);
	}
	for (i = v; i > u; i--) {
		new_tournee.push_back(tournee[i]);
	}
	for (i = v + 1; i < n; i++) {
		new_tournee.push_back(tournee[i]);
	}

	return new_tournee;
}

void SolApprocheeHeuristique::solve_VRP_heuristique(int t, int nb_steps_optim, bool verbose) {
	if (verbose) {
		cout << "t = " << t << endl;
	}

	vector<vector<int>> tournees;
	unordered_set<int> a_visiter;
	int n = instance->n;

	// Nettoyage des valeurs courantes de x
	for (int i = 0; i <= instance->n; i++) {
		courante.x[t][i].clear();
	}

	// Construction des noeuds à visiter
	for (int i = 1; i <= n; i++) {
		if (courante.z[i][t]) {
			a_visiter.insert(i);
		}
	}

	if (a_visiter.size() == 0) {
		return;
	}

	// Boucle pour construire toutes les tournees
	while (a_visiter.size() > 0) {
		// Construction d'une tournée
		vector<int> new_tournee;
		double charge = 0;
		int u = 0;
		while (true) {
			double dist_min = numeric_limits<double>::infinity();
			int v_min = -1;
			for (int v : a_visiter) {
				double dist = instance->cost(u, v);
				if ((dist < dist_min) && (charge + courante.q[v][t] <= instance->Q)) {
					dist_min = dist;
					v_min = v;
				}
			}
			if (v_min == -1) {
				break;
			}
			a_visiter.erase(v_min);
			new_tournee.push_back(v_min);
			charge += courante.q[v_min][t];
			u = v_min;
		}
		tournees.push_back(new_tournee);
	}

	// Si l'approche gloutonne ci-dessus ne donne pas une solution réalisable,
	// (ce qui est rare) on essaie avec la résolution d'un PLNE
	if (tournees.size() > instance->m) {
		a_visiter.clear();
		for (int i = 1; i <= n; i++) {
			if (courante.z[i][t]) {
				a_visiter.insert(i);
			}
		}
		int l = instance->l;
		int K = instance->m;
		IloEnv env;
		IloModel model(env);

		vector<IloNumVar> z(K); // z[j] == 1 si camion j utilisé
		vector<vector<IloNumVar>> x(n + 1); // x[i][j] == 1 si client i dans le camion j
		for (int i = 1; i <= n; i++) {
			x[i].resize(K);
		}

		for (int j = 0; j < K; j++) {
			z[j] = IloNumVar(env, 0.0, 1.0, ILOINT);
			for (int i : a_visiter) {
				x[i][j] = IloNumVar(env, 0.0, 1.0, ILOINT);
			}
		}

		IloRangeArray CC(env);

		// Contrainte 1 : chaque client dans un seul camion
		for (int i : a_visiter) {
			IloExpr cst(env);
			for (int j = 0; j < K; j++) {
				cst += x[i][j];
			}
			CC.add(cst == 1);
		}

		// Contrainte 2 : capacité des camions pris
		for (int j = 0; j < K; j++) {
			IloExpr cst(env);
			for (int i : a_visiter) {
				cst += courante.q[i][t] * x[i][j];
			}
			cst -= (instance->Q) * z[j];
			CC.add(cst <= 0);
		}

		model.add(CC);

		IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));
		for (int j = 0; j < K; j++) {
			obj.setLinearCoef(z[j], 1.0);
		}

		IloCplex cplex(model);
		if (!verbose) {
			cplex.setOut(env.getNullStream());
		}
		//cplex.setParam(IloCplex::Param::MIP::Limits::Solutions, 1); // Arrêt à la première solution
		if (!cplex.solve()) {
			if (verbose) {
				cout << "Probleme infaisable ! On doit envoyer :" << endl;
				for (int i : a_visiter) {
					cout << "q_" << i << " = " << courante.q[i][t] << endl;
				}
				cout << "On a " << instance->m << " vehicules de capacite " << instance->Q << endl;
			}
			env.end();
			throw -1; // Si pas de solution retournée, alors on leve une exception (traitée dans notre solve)
		}

		tournees.clear();
		for (int j = 0; j < K; j++) {
			if (cplex.getValue(z[j]) > 0.5) {
				vector<int> new_tournee;
				for (int i : a_visiter) {
					if (cplex.getValue(x[i][j]) > 0.5) {
						new_tournee.push_back(i);
					}
				}
				tournees.push_back(new_tournee);
			}
		}

		env.end();
	}

	if (verbose) {
		cout << "  Tournees au debut :" << endl;
		for (auto tournee : tournees) {
			cout << "    0 ";
			for (auto u : tournee) {
				cout << u << " ";
			}
			cout << endl;
		}
	}

	uniform_int_distribution<int> what_to_do(0, 1);

	// Critère d'arrêt : nombre de pas maximal
	for (int counter = 0; counter < nb_steps_optim; counter++) {
		bool modified = false;

		if (what_to_do(alea)) { // améliorer l'une des tournées
			// choix aléatoire de la tournée à être améliorée
			uniform_int_distribution<int> quelle_tournee(0, tournees.size() - 1);
			int ind_tournee = quelle_tournee(alea);

			// amélioration possible uniquement si taille au moins 3
			if (tournees[ind_tournee].size() >= 3) {
				int best_u = -1;
				int best_v = -1;
				double best_cost = numeric_limits<double>::infinity();

				// On parcourt toutes les paires d'arcs non-adjacents
				for (int u = -1; u <= (int) tournees[ind_tournee].size() - 3; u++) {
					for (int v = u + 2; v <= (int) tournees[ind_tournee].size() - 1; v++) {
						if (u == -1 && v == (int) tournees[ind_tournee].size() - 1) { // Seul cas impossible
							continue;
						}
						int val_u = (u == -1) ? 0 : tournees[ind_tournee][u];
						int val_v = tournees[ind_tournee][v];
						int val_succ_u = tournees[ind_tournee][u+1];
						int val_succ_v = (v == (int)tournees[ind_tournee].size() - 1) ? 0 : tournees[ind_tournee][v+1];
						// Différence de coût
						double delta_cost = instance->cost(val_u, val_v)
							+ instance->cost(val_succ_u, val_succ_v)
							- instance->cost(val_u, val_succ_u)
							- instance->cost(val_v, val_succ_v);

						if (delta_cost < best_cost) {
							best_cost = delta_cost;
							best_u = u;
							best_v = v;
						}
					}
				}

				if (best_cost < 0) {
					tournees[ind_tournee] = tournee_vois_2TSP(tournees[ind_tournee], best_u, best_v);
					modified = true;
				}
			}
		} else { // essayer de changer quelqu'un de tournée
			// Choix aléatoire de la tournée d'où on sortira quelqu'un
			uniform_int_distribution<int> quelle_tournee(0, tournees.size() - 1);
			int ind_tournee = quelle_tournee(alea);

			// Choix de l'élément qui sortira de la tournée
			uniform_int_distribution<int> quel_element(0, tournees[ind_tournee].size() - 1);
			int ind_element = quel_element(alea);
			int elem = tournees[ind_tournee][ind_element];

			// Cout de le supprimer de sa tournée actuelle
			int pred = (ind_element == 0) ? 0 : tournees[ind_tournee][ind_element - 1];
			int succ = (ind_element == tournees[ind_tournee].size() - 1) ? 0 : tournees[ind_tournee][ind_element + 1];
			double cout_suppression = instance->cost(pred, succ) - instance->cost(pred, elem) - instance->cost(elem, succ);

			// Parcours des autres tournées pour savoir où l'insérer
			int best_i = -1;
			int best_j = -1;
			double best_cout_insertion = numeric_limits<double>::infinity();
			for (int i = 0; i < tournees.size(); i++) {
				if (i == ind_tournee) {
					continue; // On ne considère pas le cas où il rentre dans sa propre tournee
				}
				double charge = 0;
				for (int e : tournees[i]) {
					charge += courante.q[e][t];
				}
				if (charge + courante.q[elem][t] > instance->Q) {
					continue; // Il ne peut pas aller vers une tournée où il ne rentre pas
				}

				// Si on est ici, il peut rentrer, on cherche s'il y a une position interessante
				for (int j = 0; j <= tournees[i].size(); j++) {
					pred = (j == 0) ? 0 : tournees[i][j - 1];
					succ = (j == tournees[i].size()) ? 0 : tournees[i][j];
					double cout_insertion = instance->cost(pred, elem) + instance->cost(elem, succ) - instance->cost(pred, succ);
					if (cout_insertion < best_cout_insertion) {
						best_cout_insertion = cout_insertion;
						best_i = i;
						best_j = j;
					}
				}
			}

			// On voit s'il y a l'intérêt de le mettre dans une nouvelle tournée, si possible
			if ((tournees.size() < instance->m) && 
				(2 * instance->cost(elem, 0) < best_cout_insertion) && 
				(2 * instance->cost(elem, 0) + cout_suppression < 0)) {
				// Suppression
				tournees[ind_tournee].erase(tournees[ind_tournee].begin() + ind_element);
				if (tournees[ind_tournee].size() == 0) {
					tournees.erase(tournees.begin() + ind_tournee);
				}
				// Insertion
				vector<int> new_tournee = { elem };
				tournees.push_back(new_tournee);
				modified = true;
			} else {
				if (best_cout_insertion + cout_suppression < 0) {
					// Insertion
					tournees[best_i].emplace(tournees[best_i].begin() + best_j, elem);
					// Suppression
					tournees[ind_tournee].erase(tournees[ind_tournee].begin() + ind_element);
					if (tournees[ind_tournee].size() == 0) {
						tournees.erase(tournees.begin() + ind_tournee);
					}
					modified = true;
				}
			}
		}

		if (verbose && modified) {
			cout << "  Tournees apres " << counter + 1 << " pas :" << endl;
			for (auto tournee : tournees) {
				cout << "    0 ";
				for (auto u : tournee) {
					cout << u << " ";
				}
				cout << endl;
			}
		}
	}

	// Sauvegarde des tournees trouvees
	for (auto tournee : tournees) {
		int u = 0;
		int v;
		for (int i = 0; i < tournee.size(); i++) {
			v = tournee[i];
			courante.x[t][u].push_back(v);
			u = v;
		}
		courante.x[t][u].push_back(0);
	}
}

void SolApprocheeHeuristique::calcul_SC(int t, bool verbose) {
	Graph g = courante.x[t]; // récupération du graphe
	int n = g.size() - 1;
	
	// Si 0 n'a pas de successeurs, il n'y a aucune tournée à cet instant t
    // Dans ce cas, le coût de l'insertion est le coût d'aller-retour depuis 0.
	if (g[0].empty()) {
		for (int i = 1; i <= n; i++) {
			SC[i][t] = 2 * instance->cost(0, i);
		}

	}
	else {

		vector<int> non_visites; // contient les noeuds non visités par aucune tournée
		for (int i = 0; i <= n; i++) {
			if (g[i].empty()) {
				non_visites.push_back(i);
				SC[i][t] = 2 * instance->cost(0, i); // initialisation avec le cout de l'aller-retour direct
			}
		}

		int nb_tournees = g[0].size();
		for (int i = 0; i < nb_tournees; i++) { // parcours de chaque tournée
			int pred = 0;
			int cour = g[0][i];
			int succ;
			double cost_pred_cour;
			while (cour != 0) {
				// dans chaque sommet de chaque tournée, il faut faire deux choses :
				// 1. recalculer son cout
				// 2. passer par les sommets non-visités et voir si le cout de l'inserer dans cette tournée est plus petit
				//    que le meilleur cout d'insertion connu.
				succ = g[cour][0];
				cost_pred_cour = instance->cost(pred, cour);
				SC[cour][t] = max(0.0, cost_pred_cour + instance->cost(cour, succ) - instance->cost(pred, succ));

				for (int j = 0; j < non_visites.size(); j++) {
					double cost_insertion = max(0.0, instance->cost(pred, non_visites[j]) + instance->cost(non_visites[j], cour) - cost_pred_cour);
					SC[non_visites[j]][t] = min(SC[non_visites[j]][t], cost_insertion);
				}

				pred = cour;
				cour = succ;
			}
			// Il faut aussi tester le coût de l'insertion de chaque sommet non-visité entre le dernier sommet de la tournée (pred) et 0 (cour)
			cost_pred_cour = instance->cost(pred, cour);
			for (int j = 0; j < non_visites.size(); j++) {
				double cost_insertion = max(0.0, instance->cost(pred, non_visites[j]) + instance->cost(non_visites[j], cour) - cost_pred_cour);
				SC[non_visites[j]][t] = min(SC[non_visites[j]][t], cost_insertion);
			}
		}
	}
}

void SolApprocheeHeuristique::solve(int max_iter, int nb_steps_optim, bool verbose) {
	init_SC();
	for (int i = 0; i < max_iter; i++) {
		solve_LSP(verbose);
		try { // Peut lever une exception si infaisable (impossible de diviser en camions)
			for (int t = 0; t < instance->l; t++) {
				solve_VRP_heuristique(t, nb_steps_optim, verbose);
			}
		}
		catch (int e) { // Résoudre une version renforcée du LSP pour permettre la division en camions
			solve_LSP_bis(verbose);
			for (int t = 0; t < instance->l; t++) {
				solve_VRP_heuristique(t, nb_steps_optim, verbose);
			}
		}
		for (int t = 0; t < instance->l; t++) {
			calcul_SC(t, verbose);
		}
		courante.calcul_valeur(*instance);
		if (courante.valeur < meilleure.valeur) {
			meilleure = courante;
		}
		if (verbose) {
			cout << "Valeur de la solution courante = " << courante.valeur << endl;
			cout << "Valeur de la meilleure solution trouvee = " << meilleure.valeur << endl << endl;
		}
	}
}