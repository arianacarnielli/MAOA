#include <vector>
#include <random>
#include <unordered_set>
#include <ilcplex/ilocplex.h>
#include <algorithm>

#include "SolExacteV2.h"
#include "PRP.h"
#include "SolutionV2.h"

//typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

default_random_engine alea;


SolExacteV2::SolExacteV2(PRP* inst) : solution(inst->n, inst->l) {
	instance = inst;
	int n = instance->n;
	int l = instance->l;
	int m = instance->m;

	w_sol.resize(n + 1);

	for (int i = 1; i <= n; i++) {
		w_sol[i].resize(l);
	}
}

void SolExacteV2::solve(SolutionV2* sol_init, double tolerance, bool verbose) {
	int n = instance->n;
	int l = instance->l;
	int m = instance->m;

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
	vector<vector<IloNumVar>> w;

	vector < vector<vector<IloNumVar>>> q;
	vector < vector<vector<IloNumVar>>> z;
	vector < vector<vector<vector<IloNumVar>>>> x;

	p.resize(l);
	y.resize(l);
	I.resize(n + 1);
	q.resize(n + 1);
	z.resize(n + 1);
	w.resize(n + 1);

	I[0].resize(l);
	z[0].resize(l);
	for (int i = 1; i <= n; i++) {
		I[i].resize(l);
		w[i].resize(l);
		q[i].resize(m);
		z[i].resize(m);
		for (int j = 1; j <= m; j++) {
			q[i][j].resize(l);
			z[i][j].resize(l);
		}
	}

	x.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		x[i].resize(n+1);
		for (int j = 0; j <= n; j++) {
			x[i][j].resize(m);
			for (int k = 0; k <= m; k++) {
				x[i][j][k].resize(l);
			}
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

		for (int k = 0; k < m; k++) {
			z[0][k][t] = IloNumVar(env, 0.0, 1.0, ILOINT); //inférieur au nombre de véhicules m (contrainte 10)
			varname.str("");
			varname << "z_0_" << k << "_" << t;
			z[0][k][t].setName(varname.str().c_str());
		}


		for (int i = 1; i <= n; i++) {
			I[i][t] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT);
			varname.str("");
			varname << "I_" << i << "_" << t;
			I[i][t].setName(varname.str().c_str());

			w[i][t] = IloNumVar(env, 0.0, instance->Q, ILOFLOAT);
			varname.str("");
			varname << "w_" << i << "_" << t;
			w[i][t].setName(varname.str().c_str());

			for (int k = 0; k < m; k++) {
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

	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			for (int k = 0; k < m; k++) {
				for (int t = 0; t < l; t++) {
					if (i != j) { // pas necessaire d'avoir des variables pour la diagonale
						x[i][j][k][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
						varname.str("");
						varname << "x_" << i << "_" << j << "_" << k << "_" << t;
						x[i][j][k][t].setName(varname.str().c_str());
					}
				}
			}
		}
	}

	//////////////
	//////  CONTRAINTES
	//////////////


	IloRangeArray CC(env);
	int nbcst = 0;

	ostringstream cstname;

	if (true) {
		//contrainte 21
		IloExpr cst(env);
		for (int i = 1; i <= n; i++) {
			for (int k = 0; k < m; k++) {
				cst += q[i][k][0];
			}
		}
		cst += I[0][0] - p[0];
		CC.add(cst == instance->L0[0]);
		

		cstname.str("");
		cstname << "Cst_onecol_0";
		CC[nbcst].setName(cstname.str().c_str());
		nbcst++;

		for (int t = 0; t < l; t++) {
			IloExpr cst(env);
			for (int i = 1; i <= n; i++) {
				for (int k = 0; k < m; k++) {
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

		//contrainte 22
		// il faut traiter le cas t = 0 separement
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			for (int k = 0; k < m; k++) {
				cst += q[i][k][0];
			}
			cst += -I[i][0];
			CC.add(cst == instance->d[i][0] - instance->L0[i]);
			cstname.str("");
			cstname << "Cst_twocol_" << i << "_line_0";
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}

		for (int t = 1; t < l; t++) {
			for (int i = 1; i <= n; i++) {
				IloExpr cst(env);
				for (int k = 0; k < m; k++) {
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

		//contrainte 23
		for (int t = 0; t < l; t++) {
			IloExpr cst(env);
			cst += M[t] * y[t] - p[t];
			CC.add(cst >= 0);
			cstname.str("");
			cstname << "Cst_threecol_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}

		//contrainte 24
		for (int t = 0; t < l; t++) {
			IloExpr cst(env);
			cst += I[0][t];
			CC.add(cst <= instance->L[0]);
			cstname.str("");
			cstname << "Cst_fourcol_" << t;
			CC[nbcst].setName(cstname.str().c_str());
			nbcst++;
		}

		//contrainte 25
		// il faut traiter le cas t = 0 separement
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			for (int k = 0; k < m; k++) {
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
				for (int k = 0; k < m; k++) {
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

		//contrainte 26 (les contraintes couplantes entre q et z): 
		for (int t = 0; t < l; t++) {
			for (int i = 1; i <= n; i++) {
				IloExpr cst(env);
				for (int k = 0; k < m; k++) {
					cst += z[i][k][t] * tildeM[i][t];
					cst -= q[i][k][t];
					CC.add(cst >= 0);
					cstname.str("");
					cstname << "Cst_sixcol_" << i << "_line_" << t << "_vehi_" << k;
					CC[nbcst].setName(cstname.str().c_str());
					nbcst++;
				}
			}
		}
	}
	
	if (true) {
		//contrainte 28
		for (int t = 0; t < l; t++) {
			for (int i = 0; i <= n; i++) {
				for (int k = 0; k < m; k++) {
					IloExpr cst(env);
					for (int j = 0; j <= n; j++) {
						cst += x[j][i][k][t] + x[i][j][k][t];
					}
					cst -= 2*z[i][k][t];
					CC.add(cst >= 0);
					cstname.str("");
					cstname << "Cst_8col_" << i << "_line_" << t << "_vehi_" << k;
					CC[nbcst].setName(cstname.str().c_str());
					nbcst++;
				}
			}
		}

		//contrainte 29 exponentielle
		//CPXaddindconstraints()

		//contrainte 30

		for (int t = 0; t < l; t++) {
			for (int k = 0; k < m; k++) {
				IloExpr cst(env);

				for (int i = 1; i <= n; i++) {
					cst += q[i][k][t];
				}
				cst -= z[0][k][t];
				CC.add(cst <= 0);
				cstname.str("");
				cstname << "Cst_10col_" << t  << "_vehi_" << k;
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

	for (int t = 0; t < l; t++) {

		obj.setLinearCoef(p[t], instance->u);
		obj.setLinearCoef(y[t], instance->f);

		for (int i = 0; i <= n; i++) {
			obj.setLinearCoef(I[i][t], instance->h[i]);
			for (int j = 0; j <= n; j++) {
				if (j != i) {
					for (int k = 0; k < m; k++) {
						obj.setLinearCoef(x[i][j][k][t], instance->cost(i, j));
					}
					
				}
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

	// Si une solution iniale a été passée en argument, on la charge
	if (sol_init) {

		IloNumVarArray vars(env);
		IloNumArray vals(env);

		for (int t = 0; t < l; t++) {
			vars.add(p[t]);
			vals.add(sol_init->p[t]);
			//p[t].setBounds(sol_init->p[t], sol_init->p[t]);

			vars.add(y[t]);
			vals.add(sol_init->y[t]);
			//y[t].setBounds(sol_init->y[t], sol_init->y[t]);

			vars.add(I[0][t]);
			vals.add(sol_init->I[0][t]);
			//I[0][t].setBounds(sol_init->I[0][t], sol_init->I[0][t]);

			for (int i = 1; i <= n; i++) {
				vars.add(I[i][t]);
				vals.add(sol_init->I[i][t]);
				//I[i][t].setBounds(sol_init->I[i][t], sol_init->I[i][t]);

				int qtemp = 0;
				int ztemp = 0;
				for (int k = 0; k < m; k++) {
					vars.add(q[i][k][t]);
					vals.add(sol_init->q[i][k][t]);
					//q[i][t].setBounds(sol_init->q[i][t], sol_init->q[i][t]);

					vars.add(z[i][k][t]);
					vals.add(sol_init->z[i][k][t]);
					//z[i][t].setBounds(sol_init->z[i][t], sol_init->z[i][t]);
				}
			}

			for (int i = 0; i <= n; i++) {
				for (int j = 0; j <= n; j++) {
					if (i != j) {
						for (int k = 0; k < m; k++) {
							vars.add(x[i][j][k][t]);
							vals.add(sol_init->x[i][j][k][t]);
							//x[i][j][t].setBounds(x_tab[i][j], x_tab[i][j]);
						}
					}
				}
			}
		}

		// Chargement de la condition initiale passee en argument
		cplex.addMIPStart(vars, vals);
	}

	cplex.solve();

	// Récupération des valeurs retrouvées dans solution
	for (int t = 0; t < l; t++) {
		solution.p[t] = cplex.getValue(p[t]);
		solution.y[t] = cplex.getValue(y[t]);
		solution.I[0][t] = cplex.getValue(I[0][t]);
		for (int k = 0; k < m; k++) solution.z[0][k][t] = cplex.getValue(z[0][k][t]);
		for (int i = 1; i <= n; i++) {
			solution.I[i][t] = cplex.getValue(I[i][t]);
			for (int k = 0; k < m; k++) {
				solution.q[i][k][t] = cplex.getValue(q[i][k][t]);
				solution.z[i][k][t] = cplex.getValue(z[i][k][t]);
			}

			w_sol[i][t] = cplex.getValue(w[i][t]);
		}
		for (int i = 0; i <= n; i++) {
			for (int j = 0; j <= n; j++) {
				for (int k = 0; k < m; k++) {
					if (i != j && cplex.getValue(x[i][j][k][t]) > 0.5) {
						solution.x[i][j][k][t] = 1;
					}
				}
			}
		}
	}

	// Calcul de la valeur de la solution par l'algorithme calcul_valeur
	// Doit donner la même chose que cplex.getObjValue()
	solution.calcul_valeur(*instance);

	// Affichages
	if (verbose) {
		for (int t = 0; t < l; t++) {
			cout << "p_" << t << ": " << solution.p[t] << endl;
			cout << "y_" << t << ": " << solution.y[t] << endl;
			cout << "I_0_" << t << ": " << solution.I[0][t] << endl;
			for (int k = 0; k < m; k++) cout << "z_0_" << t << "_" << k << ": " << solution.z[0][k][t] << endl;
			for (int i = 1; i <= n; i++) {
				cout << "I_" << i << "_" << t << ": " << solution.I[i][t] << endl;
				for (int k = 0; k < m; k++) {
					cout << "q_" << i << "_" << t << "_" << k << ": " << solution.q[i][k][t] << endl;
					cout << "z_" << i << "_" << t << "_" << k << ": " << solution.z[i][k][t] << endl;
				}

				cout << "w_" << i << "_" << t << ": " << w_sol[i][t] << endl;
			}

			/*
			cout << "tournee a l'instant t : " << t << endl;
			for (int i = 0; i <= n; i++) {
				cout << "Successeurs de " << i << " : ";
				for (int k2 = 0; k2 < solution.x[t][i].size(); k2++) {
					cout << solution.x[t][i][k2] << " ";
				}
				cout << endl;
			}
			cout << endl;
			*/

			for (int i = 0; i <= n; i++) {
				cout << "i = " << i << ": ";
				for (int j = 0; j <= n; j++) {
					cout << "j = " << i << ": ";
					if (i != j) {
						for (int k = 0; k < m; k++) {
							cout << "k = " << solution.x[i][j][k][t] << " ";
						}
					}
				}
				cout << endl;
			}
		}
		cout << "valeur de la solution cplex : " << cplex.getObjValue() << endl;
		cout << "valeur de la solution : " << solution.valeur << endl;
	}
}