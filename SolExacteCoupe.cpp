#include <vector>
#include <random>
#include <unordered_set>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <ctime>

#include "SolExacteCoupe.h"
#include "PRP.h"
#include "Solution.h"

//typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

default_random_engine alea(time(0));

// Coupe pour les contraintes de capacité fractionnaire
ILOUSERCUTCALLBACK4(CoupeFractionalCapacity, PRP*, instance, vector<vector<vector<IloNumVar>>>, x, vector<vector<IloNumVar>>, q, IloNum, eps) {
	int n = instance->n;
	int l = instance->l;
	
	uniform_int_distribution<int> size(1, (n+1)/2);
	uniform_int_distribution<int> element(1, n);

	for (int t = 0; t < l; t++) {
		// S est un ensemble aléatoire de taille aléatoire aussi
		// S ne contient jamais l'usine
		// Sbar = {0, 1, ..., n} \ S
		unordered_set<int> S, Sbar;
		unordered_set<int> best_S, best_Sbar;
		int r_size = size(alea);
		for (int i = 0; i < r_size; i++) {
			S.insert(element(alea));
		}
		for (int i = 0; i <= n; i++) {
			if (!S.count(i)) {
				Sbar.insert(i);
			}
		}

		// Calcul de la valeur de la contrainte sur l'ensemble S à l'instant t
		double value = 0;
		double best_value;
		for (int i : S) {
			for (int j : Sbar) {
				value += getValue(x[i][j][t]);
			}
			value -= getValue(q[i][t]) / (instance->Q);
		}
		best_value = value;
		best_S = S;
		best_Sbar = Sbar;

		// Arrêt lorsque Sbar = {0}, S = {1, ..., n}
		while (Sbar.size() > 1) {
			// Choix du sommet v à ajouter à S
			// Celui qui maximise Somme_{i in S} x[i][v][t]
			// Heuristique de Augerat et al., EJOR, 1997
			double xv_max = -1;
			int v_max = 0;
			for (int v : Sbar) {
				if (v != 0) {
					double xv_value = 0;
					for (int i : S) {
						xv_value += getValue(x[i][v][t]);
					}
					if (xv_value > xv_max) {
						xv_max = xv_value;
						v_max = v;
					}
				}
			}

			S.insert(v_max);
			Sbar.erase(v_max);

			// Calcul de la valeur de la contrainte sur l'ensemble S à l'instant t
			value = 0;
			for (int i : S) {
				for (int j : Sbar) {
					value += getValue(x[i][j][t]);
				}
				value -= getValue(q[i][t]) / (instance->Q);
			}
			// On garde la meilleure déjà trouvée
			if (value < best_value) {
				best_value = value;
				best_S = S;
				best_Sbar = Sbar;
			}
		}

		// Si négatif : contrainte violée
		// On ajoute la contrainte au PLNE
		if (best_value < -eps) {
			IloExpr cst(getEnv());
			for (int i : best_S) {
				for (int j : best_Sbar) {
					cst += x[i][j][t];
				}
				cst -= q[i][t] / (instance -> Q);
			}
			add(cst >= 0).end();
		}
	}
}


SolExacteCoupe::SolExacteCoupe(PRP* inst): solution(inst->n, inst->l) {
	instance = inst;
	int n = instance->n;
	int l = instance->l;

	w_sol.resize(n + 1);

	for (int i = 1; i <= n; i++) {
		w_sol[i].resize(l);
	}
}

void SolExacteCoupe::solve(Solution* sol_init, double tolerance, bool verbose) {
	
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

		z[0][t] = IloNumVar(env, 0.0, instance ->m, ILOINT); //inférieur au nombre de véhicules m (contrainte 10)
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

			w[i][t] = IloNumVar(env, 0.0, instance -> Q, ILOFLOAT);
			varname.str("");
			varname << "w_" << i << "_" << t;
			w[i][t].setName(varname.str().c_str());
		}
	}

	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			for (int t = 0; t < l; t++) {
				if (i != j) { // pas necessaire d'avoir des variables pour la diagonale
					x[i][j][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
					varname.str("");
					varname << "x_" << i << "_" << j << "_" << t;
					x[i][j][t].setName(varname.str().c_str());
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

	//même code que les contraintes du LSP pour les 6 premières contraintes, c'est dans un if pour pouvoir le rétrécir et ne pas se perdre dans le code
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

	//contraintes 8 à 16
	if (true) {
		//contrainte 8 (sur le papier) :
		for (int i = 1; i <= n; i++) {
			for (int t = 0; t < l; t++) {
				IloExpr cst(env);
				for (int j = 0; j <= n; j++) {
					if (i != j) {
						cst += x[i][j][t];
					}
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
					if (i != j) {
						cst += x[i][j][t] + x[j][i][t];
					}
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
		// il y a une erreur dans la contrainte du papier, on fait une version corrigée ici.
		//contrainte 11
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				for (int t = 0; t < l; t++) {
					if (i != j) {
						IloExpr cst(env);
						cst += w[i][t] - w[j][t] - q[i][t] - (instance->Q + tildeM[i][t]) * x[i][j][t];
						CC.add(cst >= -instance->Q - tildeM[i][t]);
						cstname.str("");
						cstname << "Cst_11col_" << i << "_line_" << j << "_case_" << t;
						CC[nbcst].setName(cstname.str().c_str());
						nbcst++;
					}

				}
			}
		}
		
		//contrainte 12 (allégée)
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
					obj.setLinearCoef(x[i][j][t], instance->cost(i, j));
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
	// Définition de la tolérance
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, tolerance);
	// Coupe
	cplex.use(CoupeFractionalCapacity(env, instance, x, q, cplex.getParam(IloCplex::EpRHS)));

	// Si une solution iniale a été passée en argument, on la charge
	if (sol_init){

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

				vars.add(q[i][t]);
				vals.add(sol_init->q[i][t]);
				//q[i][t].setBounds(sol_init->q[i][t], sol_init->q[i][t]);

				vars.add(z[i][t]);
				vals.add(sol_init->z[i][t]);
				//z[i][t].setBounds(sol_init->z[i][t], sol_init->z[i][t]);

			}

			vector<vector<bool>> x_tab;
			x_tab.resize(n + 1);
			for (int i = 0; i <= n; i++) {
				x_tab[i].resize(n + 1);
				for (int j = 0; j <= n; j++) {
					x_tab[i][j] = false;
				}
			}
			for (int i = 0; i <= n; i++) {
				for (int j = 0; j < sol_init->x[t][i].size(); j++) {
					x_tab[i][sol_init->x[t][i][j]] = true;
				}
			}
			for (int i = 0; i <= n; i++) {
				for (int j = 0; j <= n; j++) {
					if (i != j) {
						vars.add(x[i][j][t]);
						vals.add(x_tab[i][j]);
						//x[i][j][t].setBounds(x_tab[i][j], x_tab[i][j]);
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
		solution.z[0][t] = cplex.getValue(z[0][t]);
		for (int i = 1; i <= n; i++) {
			solution.I[i][t] = cplex.getValue(I[i][t]);
			solution.q[i][t] = cplex.getValue(q[i][t]);
			solution.z[i][t] = cplex.getValue(z[i][t]);
			w_sol[i][t] = cplex.getValue(w[i][t]);
		}
		for (int i = 0; i <= n; i++) {
			for (int j = 0; j <= n; j++) {
				if (i != j && cplex.getValue(x[i][j][t]) > 0.5) {
					solution.x[t][i].push_back(j);
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
			cout << "z_0_" << t << ": " << solution.z[0][t] << endl;
			for (int i = 1; i <= n; i++) {
				cout << "I_" << i << "_" << t << ": " << solution.I[i][t] << endl;
				cout << "q_" << i << "_" << t << ": " << solution.q[i][t] << endl;
				cout << "z_" << i << "_" << t << ": " << solution.z[i][t] << endl;
				cout << "w_" << i << "_" << t << ": " << w_sol[i][t] << endl;
			}

			cout << "tournee a l'instant t : " << t << endl;
			for (int i = 0; i <= n; i++) {
				cout << "Successeurs de " << i << " : ";
				for (int k = 0; k < solution.x[t][i].size(); k++) {
					cout << solution.x[t][i][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
			
			/*for (int i = 0; i <= n; i++) {
				cout << "i = " << i << ": ";
				for (int j = 0; j <= n; j++) {
					if (i != j) {
						cout << cplex.getValue(x[i][j][t]) << " ";
					}
				}
				cout << endl;
			}*/
		}
		cout << "valeur de la solution cplex : " << cplex.getObjValue() << endl;
		cout << "valeur de la solution : " << solution.valeur << endl;
	}
}