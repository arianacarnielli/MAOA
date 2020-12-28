// MAOA.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <fstream>


#include <cstdio> // à effacer après, c'est juste pour printf

#include "PRP.h"
#include "SolApprochee.h"
#include "SolExacte.h"

using namespace std;

int main()
{
    //instance de test (type A)
    ifstream fic("PRP_instances/Instance_0.prp");

    //type A
    //ifstream fic("PRP_instances/A_014_ABS1_15_1.prp");

    //type B
    //ifstream fic("../PRP_instances/B_050_instance1.prp");

    if (!fic.is_open()) {
        cerr << "impossible d'ouvrir le fichier";
        return 1;
    }

    PRP I(fic);

    //I.write_screen_txt();
    
    //cout << I.mc << endl;
    //cout << I.xy[0].first << " " << I.xy[0].second << endl;
    //cout << I.xy[1].first << " " << I.xy[1].second << endl;
    //cout << I.cost(0, 1) <<  endl;
    
    SolApprochee sol_app = SolApprochee(&I);
    /*
    sol_app.init_SC();
    sol_app.solve_LSP(false);
    sol_app.solve_VRP_MTZ(0, true);
    sol_app.calcul_SC(0, true);
	*/
    sol_app.solve(10, true);

	SolExacte sol_ex = SolExacte(&I);
	sol_ex.solve(true);

    cout << "approchee : " << sol_app.meilleure.valeur << endl;
    cout << "exacte : " << sol_ex.solution.valeur << endl;

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
