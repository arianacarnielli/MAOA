#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <random>
#include <limits>

#include "SolApprocheeCoupe.h"
#include "PRP.h"

typedef IloArray<IloNumVarArray> NumVarMatrix;

ILOSTLBEGIN

extern default_random_engine alea;

class CoupesVRPI : public IloCplex::UserCutCallbackI {
	int n;
	PRP* instance;
	vector<vector<IloNumVar>> x;
	vector<double> d;
	IloNum eps;
	int repeat;
public:
	CoupesVRPI(IloEnv env, int n1, PRP* instance1, vector<vector<IloNumVar>> x1,
		vector<double> d1, IloNum eps1, int repeat1) :
		IloCplex::UserCutCallbackI(env),
		n(n1),
		instance(instance1),
		x(x1),
		d(d1),
		eps(eps1),
		repeat(repeat1) {}

	IloCplex::CallbackI* duplicateCallback() const ILO_OVERRIDE {
		return (new (getEnv()) CoupesVRPI(*this));
	}

	double valeurContrainteCoupes(unordered_set<int>& S, unordered_set<int>& Sbar);
	int coupeGloutonContrainteCoupes();
	void main() ILO_OVERRIDE;
};

IloCplex::Callback CoupesVRP(IloEnv env, int n1, PRP* instance1, vector<vector<IloNumVar>> x1,
	vector<double> d1, IloNum eps1, int repeat1) {
	return (IloCplex::Callback(new (env) CoupesVRPI(env, n1, instance1, x1, d1, eps1, repeat1)));
}

double CoupesVRPI::valeurContrainteCoupes(unordered_set<int>& S, unordered_set<int>& Sbar) {
	double lhs = 0;
	double rhs = 0;

	for (int i : S) {
		for (int j : Sbar) {
			lhs += getValue(x[i][j]);
		}
		rhs += d[i];
	}
	rhs /= instance->Q;
	return lhs - ceil(rhs);
}

int CoupesVRPI::coupeGloutonContrainteCoupes() {
	int cpt_constraints = 0;

	uniform_int_distribution<int> size(1, (n + 1) / 2);
	uniform_int_distribution<int> element(1, n);
	
	unordered_set<int> best_S, best_Sbar;
	double best_value = numeric_limits<double>::infinity();

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
		double value = valeurContrainteCoupes(S, Sbar);
		if (value < best_value) {
			best_value = value;
			best_S = S;
			best_Sbar = Sbar;
		}

		// Arrêt lorsque Sbar = {0}, S = {1, ..., n}
		while (Sbar.size() > 1) {
			// Choix du sommet v à ajouter à S
			// Celui qui maximise Somme_{i in S} x[i][v]
			// Heuristique de Augerat et al., EJOR, 1997
			double xv_max = -1;
			int v_max = 0;
			for (int v : Sbar) {
				if (v != 0) {
					double xv_value = 0;
					for (int i : S) {
						xv_value += getValue(x[i][v]);
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
			value = valeurContrainteCoupes(S, Sbar);

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
		double rhs = 0;
		for (int i : best_S) {
			for (int j : best_Sbar) {
				cst += x[i][j];
			}
			rhs += d[i];
		}
		rhs /= (instance->Q);
		add(cst >= ceil(rhs)).end();
		cpt_constraints++;
	}
	
	return cpt_constraints;
}

void CoupesVRPI::main() {
	coupeGloutonContrainteCoupes();
}


SolApprocheeCoupe::SolApprocheeCoupe(PRP* inst) :
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

void SolApprocheeCoupe::init_SC() {
	for (int i = 1; i <= instance->n; i++) {
		SC[i][0] = 2 * (instance->cost(0, i));
		for (int t = 1; t < instance->l; t++) {
			SC[i][t] = SC[i][0];
		}
	}
}

void SolApprocheeCoupe::solve_LSP(bool verbose) {
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
		courante.p[t] = cplex.getValue(p[t]);
		courante.y[t] = cplex.getValue(y[t]);
		courante.I[0][t] = cplex.getValue(I[0][t]);
		for (int i = 1; i <= n; i++) {
			courante.I[i][t] = cplex.getValue(I[i][t]);
			courante.q[i][t] = cplex.getValue(q[i][t]);
			courante.z[i][t] = cplex.getValue(z[i][t]);
		}
	}

	//if (verbose) {
	//	for (int t = 0; t < l; t++) {
	//		cout << "p_" << t << ": " << courante.p[t] << endl;
	//		cout << "y_" << t << ": " << courante.y[t] << endl;
	//		cout << "I_0_" << t << ": " << courante.I[0][t] << endl;
	//		for (int i = 1; i <= n; i++) {
	//			cout << "I_" << i << "_" << t << ": " << courante.I[i][t] << endl;
	//			cout << "q_" << i << "_" << t << ": " << courante.q[i][t] << endl;
	//			cout << "z_" << i << "_" << t << ": " << courante.z[i][t] << endl;
	//		}
	//		cout << endl;
	//	}
	//}
}

void SolApprocheeCoupe::solve_VRP_MTZ(int t, double time_limit, bool verbose) {
	// N contient les noeuds par lesquels on doit passer
	vector<int> N;
	N.push_back(0);
	for (int i = 1; i <= instance->n; i++) {
		if (courante.z[i][t] > 0.5) {
			N.push_back(i);
		}
	}
	int n = N.size() - 1;

	// Tableau de demandes
	vector<double> d;
	d.resize(n + 1);
	for (int i = 1; i <= n; i++) {
		d[i] = courante.q[N[i]][t];
	}

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
				c10 += w[i] - w[j] - (instance->Q + d[i]) * x[i][j];
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

	if (!verbose) {
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
	}
	// Coupe
	cplex.use(CoupesVRP(env, n, instance, x, d, cplex.getParam(IloCplex::EpRHS), 1));

	// Temps maximal de résolution
	if (time_limit > 0) {
		cplex.setParam(IloCplex::Param::TimeLimit, time_limit);
	}

	cplex.solve();

	for (int i = 0; i <= instance->n; i++) {
		courante.x[t][i].clear();
	}
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= n; j++) {
			if (i != j && cplex.getValue(x[i][j]) > 0.5) {
				courante.x[t][N[i]].push_back(N[j]);
			}
		}
	}

	//if (verbose) {
	//	if (n > 1) { // Quand n == 1, la variable w[1] n'apparaît sur aucune contrainte ni sur la fonction objectif, dont elle n'est pas calculée
	//		for (int i = 1; i <= n; i++) {
	//			cout << "w_" << N[i] << ": " << cplex.getValue(w[i]) << endl;
	//		}
	//	}

	//	for (int i = 0; i <= n; i++) {
	//		cout << "Successeurs de " << N[i] << ": ";
	//		for (int k = 0; k < courante.x[t][N[i]].size(); k++) {
	//			cout << courante.x[t][N[i]][k] << " ";
	//		}
	//		cout << endl;
	//	}
	//}
}

void SolApprocheeCoupe::calcul_SC(int t, bool verbose) {
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

	// Affichage des résultats si verbose
	/*if (verbose) {
		for (int i = 1; i <= instance->n; i++) {
			cout << "SC[" << i << "][" << t << "] = " << SC[i][t] << endl;
		}

		cout << "tableau des distances" << endl;

		cout << "   | ";
		for (int i = 0; i <= instance->n; i++) printf("%4d | ", i);
		cout << endl;
		for (int i = 0; i <= instance->n; i++) {
			printf("%2d | ", i);
			for (int j = 0; j <= instance->n; j++) {
				printf("%4.0f | ", instance->cost(i, j));
			}
			cout << endl;
		}
	}*/
}

void SolApprocheeCoupe::solve(int max_iter, double VSP_time_limit, bool verbose) {
	init_SC();
	for (int i = 0; i < max_iter; i++) {
		solve_LSP(verbose);
		for (int t = 0; t < instance->l; t++) {
			solve_VRP_MTZ(t, VSP_time_limit, verbose);
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