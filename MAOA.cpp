// MAOA.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <ctime>
#include <chrono>

#include "PRP.h"
#include "SolApprocheeBase.h"
#include "SolApprocheeCoupe.h"
#include "SolApprocheeHeuristique.h"
#include "SolExacteBase.h"
#include "SolExacteCoupe.h"
#include "SolExacteV2.h"

using namespace std;

default_random_engine alea(time(0));

void tester_A_14(int indice = 1, int repeat = 10) {
    ofstream sortie("resultats_A_14.csv");
    sortie << "Instances,Approx Base Valeur,Approx Base Temps,Approx Coupe Valeur,Approx Coupe Temps,Approx Heur Valeur,Approx Heur Temps,Exacte Coupe Valeur,Exacte Coupe Temps" << endl;
    sortie.flush();

    chrono::time_point<chrono::high_resolution_clock> start;
    chrono::time_point<chrono::high_resolution_clock> end;

    for (int i = indice; i <= 96; i++) {
        for (int j = 1; j <= 5; j++) {
            ostringstream nom("");
            nom << "A_014_ABS" << i << "_15_" << j;

            sortie << nom.str() << ",";

            ifstream fic("PRP_instances/" + nom.str() + ".prp");

            PRP I(fic);

            fic.close();

            cout << "Instance " << nom.str() << " :" << endl;
            cout << "calcul approche base... ";

            SolApprocheeBase sol_app_base = SolApprocheeBase(&I);
            start = chrono::high_resolution_clock::now();
            sol_app_base.solve(repeat, false);
            end = chrono::high_resolution_clock::now();
            auto duration_app_base = chrono::duration_cast<chrono::microseconds>(end - start);
            
            sortie << sol_app_base.meilleure.valeur << "," << duration_app_base.count() / 1000.0 << ",";
            sortie.flush();
            cout << "ok " << duration_app_base.count() / 1000.0 << endl;

            cout << "calcul approche coupe... ";

            SolApprocheeCoupe sol_app_coupe = SolApprocheeCoupe(&I);
            start = chrono::high_resolution_clock::now();
            sol_app_coupe.solve(repeat, -1, false);
            end = chrono::high_resolution_clock::now();
            auto duration_app_coupe = chrono::duration_cast<chrono::microseconds>(end - start);

            sortie << sol_app_coupe.meilleure.valeur << "," << duration_app_coupe.count() / 1000.0 << ",";
            sortie.flush();
            cout << "ok " << duration_app_coupe.count() / 1000.0 << endl;

            cout << "calcul approche heuristique... ";

            SolApprocheeHeuristique sol_app_heur = SolApprocheeHeuristique(&I);
            start = chrono::high_resolution_clock::now();
            sol_app_heur.solve(repeat, 1000, false);
            end = chrono::high_resolution_clock::now();
            auto duration_app_heur = chrono::duration_cast<chrono::microseconds>(end - start);

            sortie << sol_app_heur.meilleure.valeur << "," << duration_app_heur.count() / 1000.0 << ",";
            sortie.flush();
            cout << "ok " << duration_app_heur.count() / 1000.0 << endl;

            cout << "calcul exacte coupe... ";

            Solution* best;
            if (sol_app_base.meilleure.valeur < sol_app_coupe.meilleure.valeur){
                best = &sol_app_base.meilleure;
            } else {
                best = &sol_app_coupe.meilleure;
            }
            if (best->valeur > sol_app_heur.meilleure.valeur) {
                best = &sol_app_heur.meilleure;
            }

            SolExacteCoupe sol_exa_coupe = SolExacteCoupe(&I);
            start = chrono::high_resolution_clock::now();
            sol_exa_coupe.solve(best, 1e-4, false);
            end = chrono::high_resolution_clock::now();
            auto duration_exa_coupe = chrono::duration_cast<chrono::microseconds>(end - start);

            sortie << sol_exa_coupe.solution.valeur << "," << duration_exa_coupe.count() / 1000.0 << endl;
            sortie.flush();
            cout << "ok " << duration_exa_coupe.count() / 1000.0 << endl;
        }
    }
    sortie.close();
}


int main(int argc, char** argv)
{
    int s = 1;
    if (argc == 2) {
        s = stoi(argv[1]);
    }

    tester_A_14(s);

    return 0;

    //instance de test (type A)
    //ifstream fic("PRP_instances/Instance_0.prp");

    //type A
    ifstream fic("PRP_instances/A_014_ABS1_15_1.prp");

    //type B
    //ifstream fic("../PRP_instances/B_200_instance1.prp");

    if (!fic.is_open()) {
        cerr << "impossible d'ouvrir le fichier";
        return 1;
    }

    PRP I(fic);

    //I.write_screen_txt();
    
    chrono::time_point<chrono::high_resolution_clock> start;
    chrono::time_point<chrono::high_resolution_clock> end;
    
    //SolApprocheeBase sol_app_base = SolApprocheeBase(&I);
    //start = chrono::high_resolution_clock::now();
    //sol_app_base.solve(10, false);
    //end = chrono::high_resolution_clock::now();
    //auto duration_app_base = chrono::duration_cast<chrono::microseconds>(end - start);

    //SolApprocheeCoupe sol_app_coupe = SolApprocheeCoupe(&I);
    //start = chrono::high_resolution_clock::now();
    //sol_app_coupe.solve(10, -1, false);
    //end = chrono::high_resolution_clock::now();
    //auto duration_app_coupe = chrono::duration_cast<chrono::microseconds>(end - start);

    SolApprocheeHeuristique sol_app_heur = SolApprocheeHeuristique(&I);
    start = chrono::high_resolution_clock::now();
    sol_app_heur.solve(1, 1000, false);
    end = chrono::high_resolution_clock::now();
    auto duration_app_heur = chrono::duration_cast<chrono::microseconds>(end - start);

    /*SolExacteBase sol_ex_base = SolExacteBase(&I);
    sol_ex_base.solve(nullptr, 1e-6, false);*/
    
    SolExacteCoupe sol_ex_coupe = SolExacteCoupe(&I);
	sol_ex_coupe.solve(&(sol_app_heur.meilleure), 1e-6, true);

    /*SolExacteV2 sol_ex_v2 = SolExacteV2(&I);
    sol_ex_v2.solve(&(sol_app_coupe.meilleure), 1e-6, true);*/

    //cout << sol_app_base.meilleure << endl;
    //cout << sol_app_coupe.meilleure << endl;
    //cout << sol_app_heur.meilleure << endl;
    //cout << sol_ex_base.solution << endl;
    //cout << sol_ex_coupe.solution << endl;
    //cout << sol_ex_v2.solution << endl;

    cout << "Valeurs des solutions :" << endl;
    //cout << "  Approchee base : " << sol_app_base.meilleure.valeur << endl;
    //cout << "  |-- Temps de calcul : " << duration_app_base.count() / 1000 << " ms" << endl;
    //cout << "  Approchee coupe : " << sol_app_coupe.meilleure.valeur << endl;
    //cout << "  |-- Temps de calcul : " << duration_app_coupe.count() / 1000 << " ms" << endl;
    cout << "  Approchee heuristique : " << sol_app_heur.meilleure.valeur << endl;
    cout << "  |-- Temps de calcul : " << duration_app_heur.count() / 1000 << " ms" << endl;
    //cout << "  Exacte base : " << sol_ex_base.solution.valeur << endl;
    cout << "  Exacte coupe : " << sol_ex_coupe.solution.valeur << endl;
    //cout << "  Exacte v2 : " << sol_ex_v2.solution.valeur << endl;

    return 0;

}
