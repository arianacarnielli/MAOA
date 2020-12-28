#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>

#include "SolExacte.h"
#include "PRP.h"

typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

SolExacte::SolExacte(PRP* inst) {
	instance = inst;
	int n = instance->n;
	int l = instance->l;

	p_sol.resize(l);
	y_sol.resize(l);
	I_sol.resize(n + 1);
	q_sol.resize(n + 1);
	z_sol.resize(n + 1);
	w_sol.resize(n + 1);

	I_sol[0].resize(l);
	z_sol[0].resize(l);
	for (int i = 1; i <= n; i++) {
		I_sol[i].resize(l);
		q_sol[i].resize(l);
		z_sol[i].resize(l);
		w_sol[i].resize(l);
	}

	x_sol.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		x_sol[i].resize(n + 1);
		for (int j = 0; j <= n; j++) {
			x_sol[i][j].resize(l);
		}
	}
}


void SolExacte::solve(bool verbose) {
	
	int n = instance->n;
	int l = instance->l;
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
	vector<vector<IloNumVar>> q;
	vector<vector<IloNumVar>> z;
	vector<vector<IloNumVar>> w;
	vector<vector<vector<IloNumVar>>> x;

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
		q[i].resize(l);
		z[i].resize(l);
		w[i].resize(l);
	}

	x.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		x[i].resize(n + 1);
		for (int j = 0; j <= n; j++) {
			x[i][j].resize(l);
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

		z[0][t] = IloNumVar(env, 0.0, instance ->m, ILOINT); //inf�rieur au nombre de v�hicules m (contrainte 10)
		varname.str("");
		varname << "z_0_" << t;
		z[0][t].setName(varname.str().c_str());


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

			w[i][t] = IloNumVar(env, 0.0, instance -> Q, ILOINT);
			varname.str("");
			varname << "w_" << i << "_" << t;
			w[i][t].setName(varname.str().c_str());
		}
	}

	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			for (int t = 0; t < l; t++) {
				x[i][j][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
				varname.str("");
				varname << "x_" << i << "_" << j << "_" << t;
				x[i][j][t].setName(varname.str().c_str());
			}
		}
	}


	//////////////
	//////  CONTRAINTES
	//////////////

	IloRangeArray CC(env);
	int nbcst = 0;
	ostringstream cstname;


	//m�me code que les contraintes du LSP pour les 6 premi�res contraintes, c'est dans un if pour pouvoir le r�tr�cir et ne pas se perdre dans le code
	if (true) {
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
	}
	//fin de la copie

	//contraintes 8 � 16
	if (true) {
		//contrainte 8 (sur le papier) :
		for (int i = 1; i <= n; i++) {
			for (int t = 0; t < l; t++) {
				IloExpr cst(env);
				for (int j = 0; j <= n; j++) {
					cst += x[i][j][t];
				}

				cst -= z[i][t];
				CC.add(cst == 0);
				cstname.str("");
				cstname << "Cst_8col_" << i << "_line_" << t;
				CC[nbcst].setName(cstname.str().c_str());
				nbcst++;
			}
		}

		//contrainte 9
		for (int i = 0; i <= n; i++) {
			for (int t = 0; t < l; t++) {
				IloExpr cst(env);
				for (int j = 0; j <= n; j++) {
					cst += x[i][j][t] + x[j][i][t];
				}

				cst -= 2 * z[i][t];
				CC.add(cst == 0);
				cstname.str("");
				cstname << "Cst_9col_" << i << "_line_" << t;
				CC[nbcst].setName(cstname.str().c_str());
				nbcst++;
			}
		}
		
		//l'usine n'est pas pris en compte dedans (remarque je suis ce qu'il y a dans le papier)
		//contrainte 11
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				for (int t = 0; t < l; t++) {
					if (i != j) {
						IloExpr cst(env);
						cst += w[i][t] - w[j][t] - q[i][t] - tildeM[i][t] * x[i][j][t];
						CC.add(cst >= -tildeM[i][t]);
						cstname.str("");
						cstname << "Cst_11col_" << i << "_line_" << j << "_case_" << t;
						CC[nbcst].setName(cstname.str().c_str());
						nbcst++;
					}

				}
			}
		}
		
		//contrainte 12 (all�g�e)
		for (int i = 1; i <= n; i++) {
			for (int t = 0; t < l; t++) {
				IloExpr cst(env);
				cst += w[i][t] - instance->Q * z[i][t];
				CC.add(cst <= 0);
				cstname.str("");
				cstname << "Cst_12col_" << i << "_line_" << t;
				CC[nbcst].setName(cstname.str().c_str());
				nbcst++;
			}
		}
	}

	//TODO : ajouter les autres contraintes (17 et +)

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
					int cos = instance->cost(i, j);
					obj.setLinearCoef(x[i][j][t], cos);
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
	cplex.solve(); //�a merde l� et jsp pourquoi

	for (int t = 0; t < l; t++) {
		p_sol[t] = cplex.getValue(p[t]);
		y_sol[t] = cplex.getValue(y[t]);
		I_sol[0][t] = cplex.getValue(I[0][t]);
		z_sol[0][t] = cplex.getValue(z[0][t]);
		for (int i = 1; i <= n; i++) {
			I_sol[i][t] = cplex.getValue(I[i][t]);
			q_sol[i][t] = cplex.getValue(q[i][t]);
			z_sol[i][t] = cplex.getValue(z[i][t]);
			w_sol[i][t] = cplex.getValue(w[i][t]);
		}

		//j'aurai pu factoriser mais je pr�f�re privil�gier la clart�
		for(int i = 0; i <= n; i++){
			for (int j = 0; j <= n; j++) {
				if (i = j) {
					x_sol[i][j][t] = 0;
				}
				else {
					x_sol[i][j][t] = cplex.getValue(x[i][j][t]);
				}
			}
		}
	}

	if (verbose) {
		for (int t = 0; t < l; t++) {
			cout << "p_" << t << ": " << p_sol[t] << endl;
			cout << "y_" << t << ": " << y_sol[t] << endl;
			cout << "I_0_" << t << ": " << I_sol[0][t] << endl;
			cout << "z_0_" << t << ": " << z_sol[0][t] << endl;
			for (int i = 1; i <= n; i++) {
				cout << "I_" << i << "_" << t << ": " << I_sol[i][t] << endl;
				cout << "q_" << i << "_" << t << ": " << q_sol[i][t] << endl;
				cout << "z_" << i << "_" << t << ": " << z_sol[i][t] << endl;
				cout << "w_" << i << "_" << t << ": " << w_sol[i][t] << endl;
			}

			for (int i = 0; i <= n; i++) {
				for (int j = 0; j <= n; j++) {
					cout << "x_" << i << "_" << j << "_" << t << ": " << x_sol[i][j][t] << endl;
				}
			}
			cout << endl;
		}
	}
}