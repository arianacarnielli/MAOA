#include <vector>
#include <random>
#include <unordered_set>
#include <ilcplex/ilocplex.h>
#include <algorithm>

#include "SolExacteV2.h"
#include "PRP.h"
#include "Solution.h"

ILOSTLBEGIN

extern default_random_engine alea;

// Contraintes de coupes pour renforcement de la formulation
class Coupes2I : public IloCplex::UserCutCallbackI {
	PRP* instance;
	vector<vector<vector<vector<IloNumVar>>>> x;
	vector<vector<vector<IloNumVar>>> q;
	IloNum eps;
	int repeat;
public:
	Coupes2I(IloEnv env, PRP* instance1, vector<vector<vector<vector<IloNumVar>>>> x1,
		vector<vector<vector<IloNumVar>>> q1, IloNum eps1, int repeat1) :
		IloCplex::UserCutCallbackI(env),
		instance(instance1),
		x(x1),
		q(q1),
		eps(eps1),
		repeat(repeat1) {}

	IloCplex::CallbackI* duplicateCallback() const ILO_OVERRIDE {
		return (new (getEnv()) Coupes2I(*this));
	}

	double valeurCapaciteFractionnaire(int t, int k, unordered_set<int>& S, unordered_set<int>& Sbar);
	int coupeGloutonCapaciteFractionnaire(int n, int l, int K);
	void main() ILO_OVERRIDE;
};

IloCplex::Callback Coupes2(IloEnv env, PRP* instance1, vector<vector<vector<vector<IloNumVar>>>> x1,
	vector<vector<vector<IloNumVar>>> q1, IloNum eps1, int repeat1) {
	return (IloCplex::Callback(new (env) Coupes2I(env, instance1, x1, q1, eps1, repeat1)));
}

double Coupes2I::valeurCapaciteFractionnaire(int t, int k, unordered_set<int>& S, unordered_set<int>& Sbar) {
	double value = 0;
	for (int i : S) {
		for (int j : Sbar) {
			value += getValue(x[i][j][k][t]);
		}
		value -= getValue(q[i][k][t]) / (instance->Q);
	}
	return value;
}

// Coupe pour les contraintes de capacité fractionnaire
int Coupes2I::coupeGloutonCapaciteFractionnaire(int n, int l, int K) {
	int cpt_constraints = 0;

	uniform_int_distribution<int> size(1, (n + 1) / 2);
	uniform_int_distribution<int> element(1, n);

	for (int t = 0; t < l; t++) {
		unordered_set<int> best_S, best_Sbar;
		double best_value = 0;

		// Pour éviter trop de calculs, k est choisi de façon aléatoire
		// Proba d'un k proportionnelle à sum_i q[i][k][t]
		vector<double> weight_k(K, 0);
		double sum_weight = 0;
		for (int k = 0; k < K; k++) {
			for (int i = 1; i <= n; i++) {
				weight_k[k] += max(0, (int) getValue(q[i][k][t]));
			}
			sum_weight += weight_k[k];
		}
		if (sum_weight < eps) {
			// Pas de contrainte violée à cet instant si il n'y a rien à transporter
			// On passe au prochain t
			continue;
		}

		discrete_distribution<int> random_k(weight_k.begin(), weight_k.end());
		int k = random_k(alea);

		for (int counter = 0; counter < repeat; counter++) {
			// S est un ensemble aléatoire de taille aléatoire aussi
			// S ne contient jamais l'usine
			// Sbar = {0, 1, ..., n} \ S
			unordered_set<int> S, Sbar;
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
			double value = valeurCapaciteFractionnaire(t, k, S, Sbar);
			if (value < best_value) {
				best_value = value;
				best_S = S;
				best_Sbar = Sbar;
			}

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
							xv_value += getValue(x[i][v][k][t]);
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
				value = valeurCapaciteFractionnaire(t, k, S, Sbar);

				// On garde la meilleure déjà trouvée
				if (value < best_value) {
					best_value = value;
					best_S = S;
					best_Sbar = Sbar;
				}
			}
		}
		
		// Si négatif : contrainte violée
		// On ajoute la contrainte au PLNE
		if (best_value < -eps) {
			IloExpr cst(getEnv());
			for (int i : best_S) {
				for (int j : best_Sbar) {
					cst += x[i][j][k][t];
				}
				cst -= q[i][k][t] / (instance->Q);
			}
			add(cst >= 0).end();
			cpt_constraints++;
		}
	}

	return cpt_constraints;
}

// Appel des coupes
void Coupes2I::main() {
	//cout << "UserCut" << endl;
	int n = instance->n;
	int l = instance->l;
	int K = min(instance->m, instance->n);

	coupeGloutonCapaciteFractionnaire(n, l, K);
	//cout << coupeGloutonCapaciteFractionnaire(n, l, K) << " contraintes ajoutees" << endl;
}

// Contrainte 29
ILOLAZYCONSTRAINTCALLBACK5(BriserCycle, int, n, int, l, int, K, vector<vector<vector<vector<IloNumVar>>>>, x, vector<vector<vector<IloNumVar>>>, z) {
	for (int t = 0; t < l; t++) {
		for (int k = 0; k < K; k++) {
			unordered_set<int> S;
			// On initialise S avec tous les sommets qui doivent être visités
			for (int i = 1; i <= n; i++) {
				if (getValue(z[i][k][t]) > 0.5) {
					S.insert(i);
				}
			}

			if (S.size() == 0) {
				continue;
			}

			// On parcourt l'unique cycle contenant l'usine et on les supprime de S
			int u;
			for (u = 1; (u <= n) && (getValue(x[0][u][k][t]) < 0.5); u++);
			if (u == n + 1) {
				// Dans ce cas, rien ne sort de l'usine, on veut rajouter la contrainte avec S entier
				// pour sauteur la prochaine boucle, on met u = 0
				u = 0;
			}
			// Parcours du cycle de l'usine
			while (u) {
				S.erase(u);
				int v;
				for (v = 0; (v==u) || (getValue(x[u][v][k][t]) < 0.5); v++);
				u = v;
			}

			// S'il reste des éléments dans S après ce parcours, alors forcément S viole la contrainte
			// On la rajoute pour tous les k
			if (S.size() > 0) {
				for (int k2 = 0; k2 < K; k2++) {
					IloExpr cst(getEnv());
					for (int i : S) {
						for (int j : S) {
							if (i != j) {
								cst += x[i][j][k2][t];
							}
						}
					}
					add(cst <= (int)S.size() - 1).end();
				}
			}
		}
	}
}

SolExacteV2::SolExacteV2(PRP* inst) : solution(inst->n, inst->l) {
	instance = inst;
}

void SolExacteV2::solve(Solution* sol_init, double tolerance, double time_limit, bool verbose) {
	int n = instance->n;
	int l = instance->l;
	int K = min(instance->m, instance->n);

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
	vector<vector<vector<vector<IloNumVar>>>> x;

	p.resize(l);
	y.resize(l);
	I.resize(n + 1);
	q.resize(n + 1);
	z.resize(n + 1);

	I[0].resize(l);
	z[0].resize(K);
	for (int k = 0; k < K; k++) {
		z[0][k].resize(l);
	}

	for (int i = 1; i <= n; i++) {
		I[i].resize(l);
		q[i].resize(K);
		z[i].resize(K);
		for (int k = 0; k < K; k++) {
			q[i][k].resize(l);
			z[i][k].resize(l);
		}
	}

	x.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		x[i].resize(n+1);
		for (int j = 0; j <= n; j++) {
			x[i][j].resize(K);
			for (int k = 0; k < K; k++) {
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

		for (int k = 0; k < K; k++) {
			z[0][k][t] = IloNumVar(env, 0.0, 1.0, ILOINT);
			varname.str("");
			varname << "z_0_" << k << "_" << t;
			z[0][k][t].setName(varname.str().c_str());
		}


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

	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			if (i != j) { // pas necessaire d'avoir des variables pour la diagonale
				for (int k = 0; k < K; k++) {
					for (int t = 0; t < l; t++) {
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
		//Contrainte 21 : I0,t-1 + pt = sum{i, k}(qi,k,t) + I0,t
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

		//contrainte 22
		// il faut traiter le cas t = 0 separement
		for (int i = 1; i <= n; i++) {
			IloExpr cst(env);
			for (int k = 0; k < K; k++) {
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

		//contrainte 26 (les contraintes couplantes entre q et z): 
		for (int t = 0; t < l; t++) {
			for (int i = 1; i <= n; i++) {
				for (int k = 0; k < K; k++) {
					IloExpr cst(env);
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

		//contrainte 27
		for (int t = 0; t < l; t++) {
			for (int i = 1; i <= n; i++) {
				IloExpr cst(env);
				for (int k = 0; k < K; k++) {
					cst += z[i][k][t];
				}
				CC.add(cst <= 1);
				cstname.str("");
				cstname << "Cst_sevencol_" << i << "_line_" << t;
				CC[nbcst].setName(cstname.str().c_str());
				nbcst++;
			}
		}

		//contrainte 28 (découpée en deux car formulation incorrecte sinon)
		for (int t = 0; t < l; t++) {
			for (int i = 0; i <= n; i++) {
				for (int k = 0; k < K; k++) {
					IloExpr cst_a(env);
					IloExpr cst_b(env);
					for (int j = 0; j <= n; j++) {
						if (i != j) {
							cst_a += x[j][i][k][t];
							cst_b += x[i][j][k][t];
						}
					}
					cst_a -= z[i][k][t];
					cst_b -= z[i][k][t];
					
					CC.add(cst_a == 0);
					cstname.str("");
					cstname << "Cst_8a_col_" << i << "_line_" << t << "_vehi_" << k;
					CC[nbcst].setName(cstname.str().c_str());
					nbcst++;

					CC.add(cst_b == 0);
					cstname.str("");
					cstname << "Cst_8b_col_" << i << "_line_" << t << "_vehi_" << k;
					CC[nbcst].setName(cstname.str().c_str());
					nbcst++;
				}
			}
		}

		//contrainte 29 exponentielle
		//traitée plus tard

		//contrainte 30
		for (int t = 0; t < l; t++) {
			for (int k = 0; k < K; k++) {
				IloExpr cst(env);
				for (int i = 1; i <= n; i++) {
					cst += q[i][k][t];
				}
				cst -= (instance->Q) * z[0][k][t];
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
					double cost = instance->cost(i, j);
					for (int k = 0; k < K; k++) {
						obj.setLinearCoef(x[i][j][k][t], cost);
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
		cplex.setWarning(env.getNullStream());
	}
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, tolerance);

	//contrainte 29 comme LazyConstraint
	cplex.use(BriserCycle(env, n, l, K, x, z));

	//UserCuts pour les contraintes de renforcement
	cplex.use(Coupes2(env, instance, x, q, cplex.getParam(IloCplex::EpRHS), 1));

	// Temps maximal de résolution
	if (time_limit > 0) {
		cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
	}

	// Si une solution initiale a été passée en argument, on la charge
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

			for (int i = 0; i <= n; i++) {
				vars.add(I[i][t]);
				vals.add(sol_init->I[i][t]);
				//I[i][t].setBounds(sol_init->I[i][t], sol_init->I[i][t]);
			}

			// Parcours des tournées
			int qtt_camions = sol_init->x[t][0].size();
			for (int k = 0; k < qtt_camions; k++) {
				// Ces variables vont contenir les valeurs de x, q, z pour la tournée k
				vector<vector<bool>> x_tab;
				vector<double> q_tab;
				vector<bool> z_tab;
				x_tab.resize(n + 1);
				q_tab.resize(n + 1);
				z_tab.resize(n + 1);
				for (int i = 0; i <= n; i++) {
					x_tab[i].resize(n + 1);
					for (int j = 0; j <= n; j++) {
						x_tab[i][j] = false;
					}
				}
				z_tab[0] = true;
				for (int i = 1; i <= n; i++) {
					q_tab[i] = 0;
					z_tab[i] = false;
				}

				// Parcours de la tornée et stockage dans les bonnes variables
				int u = 0;
				int v = sol_init->x[t][0][k];
				x_tab[u][v] = true;
				while (v) {
					q_tab[v] = sol_init->q[v][t];
					z_tab[v] = sol_init->z[v][t];

					u = v;
					v = sol_init->x[t][v][0];
					x_tab[u][v] = true;
				}

				// Ajout dans les listes de vars et vals
				for (int i = 0; i <= n; i++) {
					for (int j = 0; j <= n; j++) {
						if (i != j) {
							vars.add(x[i][j][k][t]);
							vals.add(x_tab[i][j]);
							//x[i][j][t].setBounds(x_tab[i][j], x_tab[i][j]);
						}
					}
				}

				vars.add(z[0][k][t]);
				vals.add(z_tab[0]);

				for (int i = 1; i <= n; i++) {
					vars.add(q[i][k][t]);
					vals.add(q_tab[i]);

					vars.add(z[i][k][t]);
					vals.add(z_tab[i]);
				}
			}

			// Les autres camions ne sont pas utilisés : initialisation à 0
			for (int k = qtt_camions; k < K; k++) {
				// Initialisation à 0
				vars.add(z[0][k][t]);
				vals.add(0);

				for (int i = 1; i <= n; i++) {
					vars.add(q[i][k][t]);
					vals.add(0);

					vars.add(z[i][k][t]);
					vals.add(0);
				}

				for (int i = 0; i <= n; i++) {
					for (int j = 0; j <= n; j++) {
						if (i != j) {
							vars.add(x[i][j][k][t]);
							vals.add(0);
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

		solution.z[0][t] = false;
		for (int k = 0; k < K; k++) {
			if (cplex.getValue(z[0][k][t]) > 0.5) {
				solution.z[0][t] = true;
			}
		}

		solution.I[0][t] = cplex.getValue(I[0][t]);

		for (int i = 1; i <= n; i++) {
			solution.z[i][t] = false;
			solution.I[i][t] = cplex.getValue(I[i][t]);
			solution.q[i][t] = 0;

			for (int k = 0; k < K; k++) {
				if (cplex.getValue(z[i][k][t]) > 0.5) {
					solution.z[i][t] = true;
				}
				solution.q[i][t] += cplex.getValue(q[i][k][t]);
			}
		}

		for (int i = 0; i <= n; i++) {
			for (int j = 0; j <= n; j++) {
				if (i != j) {
					for (int k = 0; k < K; k++) {
						if (cplex.getValue(x[i][j][k][t]) > 0.5) {
							solution.x[t][i].push_back(j);
						}
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
			cout << "Pas de temps " << t << endl;

			for (int k = 0; k < K; k++) {
				cout << "  Camion " << k << endl;

				for (int i = 0; i <= n; i++) {
					for (int j = 0; j <= n; j++) {
						if (i != j) {
							cout << cplex.getValue(x[i][j][k][t]) << " ";
						} else {
							cout << "0 ";
						}
					}
					cout << endl;
				}
			}
		}

		cout << "valeur de la solution cplex : " << cplex.getObjValue() << endl;
		cout << "valeur de la solution : " << solution.valeur << endl;
	}
}