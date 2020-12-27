#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>

#include "SolApprochee.h"
#include "PRP.h"

typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

SolApprochee::SolApprochee(PRP* inst) {
	instance = inst;
	int n = instance->n;
	int l = instance->l;

	SC.resize(n + 1);

	p_sol.resize(l);
	y_sol.resize(l);
	I_sol.resize(n + 1);
	q_sol.resize(n + 1);
	z_sol.resize(n + 1);

	I_sol[0].resize(l);
	for (int i = 1; i <= n; i++) {
		SC[i].resize(l);

		I_sol[i].resize(l);
		q_sol[i].resize(l);
		z_sol[i].resize(l);
	}
}

void SolApprochee::init_SC() {
	for (int i = 1; i <= instance->n; i++) {
		SC[i][0] = 2 * (instance->cost(0, i));
		for (int t = 1; t < instance->l; t++) {
			SC[i][t] = SC[i][0];
		}
	}
}

void SolApprochee::solve_LSP(bool verbose) {
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

	model.add(CC);

	//////////////
	////// OBJ
	//////////////

	IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));

	for (int t = 0; t < l; t++) {

		obj.setLinearCoef(p[t], instance->u);
		obj.setLinearCoef(y[t], instance->f);

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
		p_sol[t] = cplex.getValue(p[t]);
		y_sol[t] = cplex.getValue(y[t]);
		I_sol[0][t] = cplex.getValue(I[0][t]);
		for (int i = 1; i <= n; i++) {
			I_sol[i][t] = cplex.getValue(I[i][t]);
			q_sol[i][t] = cplex.getValue(q[i][t]);
			z_sol[i][t] = cplex.getValue(z[i][t]);
		}
	}

	if (verbose) {
		for (int t = 0; t < l; t++) {
			cout << "p_" << t << ": " << p_sol[t] << endl;
			cout << "y_" << t << ": " << y_sol[t] << endl;
			cout << "I_0_" << t << ": " << I_sol[0][t] << endl;
			for (int i = 1; i <= n; i++) {
				cout << "I_" << i << "_" << t << ": " << I_sol[i][t] << endl;
				cout << "q_" << i << "_" << t << ": " << q_sol[i][t] << endl;
				cout << "z_" << i << "_" << t << ": " << z_sol[i][t] << endl;
			}
			cout << endl;
		}
	}
}

void SolApprochee::solve_VRP_MTZ(int t, bool verbose) {
	// N contient les noeuds par lesquels on doit passer
	vector<int> N;
	N.push_back(0);
	for (int i = 1; i <= instance->n; i++) {
		if (z_sol[i][t]) {
			N.push_back(i);
		}
	}
	int n = N.size() - 1;

	//////////////
	//////  CPLEX INITIALIZATION
	//////////////

	IloEnv env;
	IloModel model(env);

	////////////////////////
	//////  VAR
	////////////////////////
	ostringstream varname;

	//tous les arcs entre tous les clients et le productuer
	vector<vector<IloNumVar>> x;
	x.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		x[i].resize(n + 1);
		for (int j = 0; j <= n; j++) {
			if (i != j) { // pas necessaire d'avoir des variables pour la diagonale
				x[i][j] = IloNumVar(env, 0.0, 1.0, ILOINT);
				varname.str("");
				varname << "x_" << N[i] << "_" << N[j];
				x[i][j].setName(varname.str().c_str());
			}
		}
	}

	vector<IloNumVar> w;
	w.resize(n + 1);
	for (int i = 1; i <= n; i++) {
		w[i] = IloNumVar(env, 0.0, instance->Q, ILOFLOAT);
		varname.str("");
		varname << "w_" << N[i];
		w[i].setName(varname.str().c_str());
	}

	//////////////
	//////  CST
	//////////////
	ostringstream cstname;

	IloRangeArray CC(env);
	int nbcst = 0;

	//contrainte 6 (numérotés par rapport au sujet)
	IloExpr c6(env);
	for (int j = 1; j <= n; j++)
		c6 += x[0][j];
	CC.add(c6 <= instance->m);
	cstname.str("");
	cstname << "Cst_6";
	CC[nbcst].setName(cstname.str().c_str());
	nbcst++;

	//contrainte 7
	IloExpr c7(env);
	for (int i = 1; i <= n; i++)
		c7 += x[i][0];
	CC.add(c7 <= instance->m);
	cstname.str("");
	cstname << "Cst_7";
	CC[nbcst].setName(cstname.str().c_str());
	nbcst++;

	//contrainte 8
	for (int i = 1; i <= n; i++) {
		IloExpr c8(env);
		for (int j = 0; j <= n; j++) {
			if (i != j) {
				c8 += x[i][j];
			}
		}
		CC.add(c8 == 1);
		cstname.str("");
		cstname << "Cst_8col_" << N[i];
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 9, presque pareil
	for (int j = 1; j <= n; j++) {
		IloExpr c9(env);
		for (int i = 0; i <= n; i++) {
			if (i != j) {
				c9 += x[i][j];
			}
		}
		CC.add(c9 == 1);
		cstname.str("");
		cstname << "Cst_9col_" << N[j];
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;
	}

	//contrainte 10 
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			if (i != j) {
				IloExpr c10(env);
				c10 += w[i] - w[j] - (instance->Q + q_sol[N[i]][t]) * x[i][j];
				CC.add(c10 >= -instance->Q);
				cstname.str("");
				cstname << "Cst_11col_" << N[i] << "_line_" << N[j];
				CC[nbcst].setName(cstname.str().c_str());
				nbcst++;
			}
		}
	}

	model.add(CC);

	//////////////
	////// OBJ
	//////////////

	IloObjective obj = IloAdd(model, IloMinimize(env, 0.0));

	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			if (i != j) {
				obj.setLinearCoef(x[i][j], instance->cost(N[i], N[j]));
			}
		}
	}

	///////////
	//// RESOLUTION
	//////////

	IloCplex cplex(model);

	/*
	IloNumVarArray vars(env);
	IloNumArray vals(env);
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			if (i != j) {
				x[i][j].setLB((i == 0 && j == 1) || (i == 1 && j == 4) || (i == 4 && j == 8) || (i == 8 && j == 0));
				x[i][j].setUB((i == 0 && j == 1) || (i == 1 && j == 4) || (i == 4 && j == 8) || (i == 8 && j == 0));
				//vars.add(x[i][j]);
				//vals.add((i == 0 && j == 1) || (i == 1 && j == 4) || (i == 4 && j == 8) || (i == 8 && j == 0));
			}
		}
	}
	for (int i = 1; i <= n; i++) {
		//vars.add(w[i]);
		switch (i) {
		case 1:
			w[i].setLB(60);
			w[i].setUB(60);
			//vals.add(60);
			break;
		case 4:
			w[i].setLB(40);
			w[i].setUB(40);
			//vals.add(40);
			break;
		case 8:
			w[i].setLB(26);
			w[i].setUB(26);
			//vals.add(26);
			break;
		default:
			w[i].setLB(0);
			w[i].setUB(0);
			//vals.add(0);
			break;
		}
	}
	//cplex.addMIPStart(vars, vals);

	//cout << cplex.isMIP() << endl;

	//cout << CC[0] << endl;

	cout << CC << endl;
	*/
	if (!verbose) {
		cplex.setOut(env.getNullStream());
	}
	cplex.solve();

	if (verbose) {
		for (int i = 1; i <= n; i++) {
			cout << "w_" << N[i] << ": " << cplex.getValue(w[i]) << endl;
		}

		for (int i = 0; i <= n; i++) {
			for (int j = 0; j <= n; j++) {
				if (i != j) {
					cout << "x_" << N[i] << "_" << N[j] << ": " << cplex.getValue(x[i][j]) << endl;
				}
			}
		}
	}
}