// MAOA.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <chrono>

#include "PRP.h"
#include "SolApprocheeBase.h"
#include "SolApprocheeCoupe.h"
#include "SolExacteBase.h"
#include "SolExacteCoupe.h"

using namespace std;

default_random_engine alea(time(0));

int main()
{
    //instance de test (type A)
    //ifstream fic("PRP_instances/Instance_0.prp");

    //type A
    //ifstream fic("PRP_instances/A_014_ABS1_15_1.prp");

    //type B
    ifstream fic("../PRP_instances/B_050_instance1.prp");

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
    //sol_app_base.solve(1, true);
    //end = chrono::high_resolution_clock::now();
    //auto duration_app_base = chrono::duration_cast<chrono::microseconds>(end - start);

    SolApprocheeCoupe sol_app_coupe = SolApprocheeCoupe(&I);
    start = chrono::high_resolution_clock::now();
    sol_app_coupe.solve(1, 30, true);
    end = chrono::high_resolution_clock::now();
    auto duration_app_coupe = chrono::duration_cast<chrono::microseconds>(end - start);

    //SolExacteBase sol_ex_base = SolExacteBase(&I);
    //sol_ex_base.solve(&(sol_app.meilleure), 0.02, true);
    
    //SolExacteCoupe sol_ex_coupe = SolExacteCoupe(&I);
	//sol_ex_coupe.solve(&(sol_app_coupe.meilleure), 1e-4, true);

    //cout << sol_app_base.meilleure << endl;
    cout << sol_app_coupe.meilleure << endl;
    //cout << sol_ex_base.solution << endl;
    //cout << sol_ex_coupe.solution << endl;

    cout << "Valeurs des solutions :" << endl;
    //cout << "  Approchee base : " << sol_app_base.meilleure.valeur << endl;
    //cout << "  Temps de calcul : " << duration_app_base.count() / 1000 << " ms" << endl;
    cout << "  Approchee coupe : " << sol_app_coupe.meilleure.valeur << endl;
    cout << "  Temps de calcul : " << duration_app_coupe.count() / 1000 << " ms" << endl;
    //cout << "  Exacte base : " << sol_ex_base.solution.valeur << endl;
    //cout << "  Exacte coupe : " << sol_ex_coupe.solution.valeur << endl;

    

    return 0;

}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
