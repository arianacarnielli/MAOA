// MAOA.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//

#include <iostream>
#include <fstream>

#include "PRP.h"
#include "SolApprochee.h"
#include "SolExacte.h"

using namespace std;

int main()
{
    //type A
    ifstream fic("PRP_instances/A_014_ABS1_15_1.prp");

    //type B
    //ifstream fic("../PRP_instances/B_050_instance1.prp");

    PRP I(fic);

    //I.write_screen_txt();
    
    //cout << I.mc << endl;
    //cout << I.xy[0].first << " " << I.xy[0].second << endl;
    //cout << I.xy[1].first << " " << I.xy[1].second << endl;
    //cout << I.cost(0, 1) <<  endl;
    
	/*
    SolApprochee sol = SolApprochee(&I); 
    sol.init_SC();
    sol.solve_LSP(true);

    sol.solve_VRP_MTZ(2, true);
	*/
	

	
	SolExacte sol = SolExacte(&I);
	sol.solve(true);
	


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
