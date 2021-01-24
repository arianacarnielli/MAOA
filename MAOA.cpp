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

void tester_A_14(int indice = 1, int indice2 = 1, int repeat = 10) {
    ofstream sortie("resultats_A_14.csv");
    sortie << "Instances,Approx Base Valeur,Approx Base Temps,Approx Coupe Valeur,Approx Coupe Temps,Approx Heur Valeur,Approx Heur Temps,Exacte Coupe Valeur,Exacte Coupe Temps" << endl;
    sortie.flush();

    chrono::time_point<chrono::high_resolution_clock> start;
    chrono::time_point<chrono::high_resolution_clock> end;

    for (int i = indice; i <= 96; i++) {
        for (int j = indice2; j <= 5; j++) {
            indice2 = 1;
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
            sol_exa_coupe.solve(best, 1e-4, -1, "G1", false);
            end = chrono::high_resolution_clock::now();
            auto duration_exa_coupe = chrono::duration_cast<chrono::microseconds>(end - start);

            sortie << sol_exa_coupe.solution.valeur << "," << duration_exa_coupe.count() / 1000.0 << endl;
            sortie.flush();
            cout << "ok " << duration_exa_coupe.count() / 1000.0 << endl;
        }
    }
    sortie.close();
}

void tester_coupe_vs_v2(int repeat) {
    ifstream fic("PRP_instances/A_014_ABS1_15_1.prp");
    if (!fic.is_open()) {
        cerr << "impossible d'ouvrir le fichier";
        return;
    }

    PRP I(fic);

    fic.close();

    ofstream sortie("resultsts_coupe_vs_v2.csv");

    sortie << "Iteration,Coupe Valeur,Coupe Temps,V2 Valeur,V2 Temps" << endl;

    chrono::time_point<chrono::high_resolution_clock> start;
    chrono::time_point<chrono::high_resolution_clock> end;

    SolApprocheeHeuristique sol_app = SolApprocheeHeuristique(&I);
    start = chrono::high_resolution_clock::now();
    sol_app.solve(10, 1000, false);
    end = chrono::high_resolution_clock::now();
    double duration_app = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
    double duration;

    for (int i = 0; i < repeat; i++) {
        sortie << i << ",";
        cout << "i = " << i << endl;

        cout << "  Coupe : ";

        SolExacteCoupe sol_coupe = SolExacteCoupe(&I);
        start = chrono::high_resolution_clock::now();
        sol_coupe.solve(&(sol_app.meilleure), 1e-6, -1, "G1", false);
        end = chrono::high_resolution_clock::now();
        duration = duration_app + chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        sortie << sol_coupe.solution.valeur << "," << duration << ",";

        cout << duration << endl;

        cout << "  V2 : ";

        SolExacteV2 sol_v2 = SolExacteV2(&I);
        start = chrono::high_resolution_clock::now();
        sol_v2.solve(&(sol_app.meilleure), 1e-6, -1, false);
        end = chrono::high_resolution_clock::now();
        duration = duration_app + chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        sortie << sol_v2.solution.valeur << "," << duration << endl;

        cout << duration << endl;
    }

    sortie.close();
}

void tester_instance_0(int repeat) {
    ifstream fic("PRP_instances/Instance_1.prp");
    if (!fic.is_open()) {
        cerr << "impossible d'ouvrir le fichier";
        return;
    }

    PRP I(fic);

    double mean_time = 0;
    double std_time = 0;
    double mean_val = 0;
    double std_val = 0;

    chrono::time_point<chrono::high_resolution_clock> start;
    chrono::time_point<chrono::high_resolution_clock> end;

    cout << "Resolution approchee base";
    for (int i = 0; i < repeat; i++) {
        SolApprocheeBase sol = SolApprocheeBase(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(10, false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.meilleure.valeur;
        std_val += sol.meilleure.valeur * sol.meilleure.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time/(repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val/(repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution approchee coupe";
    for (int i = 0; i < repeat; i++) {
        SolApprocheeCoupe sol = SolApprocheeCoupe(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(10, -1, false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.meilleure.valeur;
        std_val += sol.meilleure.valeur * sol.meilleure.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution approchee heuristique";
    for (int i = 0; i < repeat; i++) {
        SolApprocheeHeuristique sol = SolApprocheeHeuristique(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(10, 1000, false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.meilleure.valeur;
        std_val += sol.meilleure.valeur * sol.meilleure.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution exacte base";
    for (int i = 0; i < repeat; i++) {
        SolExacteBase sol = SolExacteBase(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(nullptr, 1e-6, -1, false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.solution.valeur;
        std_val += sol.solution.valeur * sol.solution.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution exacte coupe gloutonne 1";
    for (int i = 0; i < repeat; i++) {
        SolExacteCoupe sol = SolExacteCoupe(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(nullptr, 1e-6, -1, "G1", false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.solution.valeur;
        std_val += sol.solution.valeur * sol.solution.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution exacte coupe gloutonne 2";
    for (int i = 0; i < repeat; i++) {
        SolExacteCoupe sol = SolExacteCoupe(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(nullptr, 1e-6, -1, "G2", false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.solution.valeur;
        std_val += sol.solution.valeur * sol.solution.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution exacte coupe tabou 1";
    for (int i = 0; i < repeat; i++) {
        SolExacteCoupe sol = SolExacteCoupe(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(nullptr, 1e-6, -1, "T1", false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.solution.valeur;
        std_val += sol.solution.valeur * sol.solution.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution exacte coupe tabou 2";
    for (int i = 0; i < repeat; i++) {
        SolExacteCoupe sol = SolExacteCoupe(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(nullptr, 1e-6, -1, "T2", false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.solution.valeur;
        std_val += sol.solution.valeur * sol.solution.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;

    cout << "Resolution exacte v2";
    for (int i = 0; i < repeat; i++) {
        SolExacteV2 sol = SolExacteV2(&I);
        start = chrono::high_resolution_clock::now();
        sol.solve(nullptr, 1e-6, -1, false);
        end = chrono::high_resolution_clock::now();
        double duration = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        mean_time += duration;
        std_time += duration * duration;

        mean_val += sol.solution.valeur;
        std_val += sol.solution.valeur * sol.solution.valeur;

        cout << ".";
    }
    cout << " ok!" << endl;
    mean_time /= repeat;
    std_time = sqrt(std_time / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_time * mean_time);
    mean_val /= repeat;
    std_val = sqrt(std_val / (repeat - 1.0) - repeat / (repeat - 1.0) * mean_val * mean_val);
    cout << "Temps moyen : " << mean_time << " +- " << std_time << " ms" << endl;
    cout << "Valeur moyenne : " << mean_val << " +- " << std_val << endl << endl;
    mean_time = std_time = mean_val = std_val = 0;
}

void tester_app_heuristique_all(int repeat) {
    ofstream sortie("resultats_app_heuristique.csv");
    sortie << "Instances,Valeur,Temps" << endl;
    sortie.flush();

    vector<string> filenames;
    vector<pair<string, string>> A_sizes = { {"014", "15"}, {"050", "50"}, {"100", "100"} };
    vector<string> B_sizes = { "050", "100", "200" };

    chrono::time_point<chrono::high_resolution_clock> start;
    chrono::time_point<chrono::high_resolution_clock> end;

    for (auto size : A_sizes) {
        for (int i = 1; i <= 96; i++) {
            for (int j = 1; j <= 5; j++) {
                ostringstream nom("");
                nom << "A_" << size.first << "_ABS" << i << "_" << size.second << "_" << j;
                filenames.push_back(nom.str());
            }
        }
    }
    for (string size : B_sizes) {
        for (int i = 1; i <= 30; i++) {
            ostringstream nom("");
            nom << "B_" << size << "_instance" << i;
            filenames.push_back(nom.str());
        }
    }

    for (string filename : filenames) {
        sortie << filename << ",";

        ifstream fic("PRP_instances/" + filename + ".prp");
        PRP I(fic);
        fic.close();

        cout << "Instance " << filename << " : ";

        SolApprocheeHeuristique sol_app_heur = SolApprocheeHeuristique(&I);
        start = chrono::high_resolution_clock::now();
        sol_app_heur.solve(repeat, 1000, false);
        end = chrono::high_resolution_clock::now();
        double duration_app_heur = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;

        sortie << sol_app_heur.meilleure.valeur << "," << duration_app_heur << endl;
        sortie.flush();
        cout << duration_app_heur << endl;
    }

    sortie.close();
}

void cli(int argc, char** argv) {
    string mode = "EC";
    string start = "AH";
    string coupe = "G1";
    int verbose = 1; 
    int repeat = 10;
    double timeout = -1;
    int heur = 1000;
    double tol = 1e-4;

    PRP I;
    ofstream out_file;
    string filename;

    ostringstream msg("");
    msg << "un probleme est survenu dans le programme MAOA : " << endl;
    try {
        if (argc < 2) {
            msg << "le programme a besoin d'un fichier .prp avec l'instance du probleme a resoudre.";
            throw msg.str();
        }
        filename = argv[1];
        ifstream fic(filename);
        if (!fic.is_open()) {
            msg << "impossible d'ouvrir le fichier, il doit etre un fichier du type .prp.";
            throw msg.str();
        }
        try {
            PRP inst(fic);
            I = inst;
        }
        catch (...){
            msg << "impossible d'ouvrir le fichier, il doit etre un fichier du type .prp.";
            throw msg.str();
        }

        for (int i = 2; i < argc; i+=2) {
            string opt = argv[i];
            
            if (opt == "--mode") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                mode = argv[i + 1];
                if (mode != "AB" && mode != "AC" && mode != "AH" && mode != "EB" && mode != "EC" && mode != "E2") {
                    msg << "Mode de calcul non reconnu." << endl;
                    msg << "Le mode de calcul doit etre l'un des suivants :" << endl;
                    msg << "AB : pour le mode approche de base" << endl;
                    msg << "AC : pour le mode approche avec coupe" << endl;
                    msg << "AH : pour le mode approche heuristique" << endl;
                    msg << "EB : pour le mode exact (premiere formulation) de base" << endl;
                    msg << "EC : pour le mode exact (premiere formulation) avec coupe (mode par defaut)" << endl;
                    msg << "E2 : pour le mode exact (deuxieme formulation)" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--start") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                start = argv[i + 1];
                if (start != "AB" && start != "AC" && start != "AH" && start != "none") {
                    msg << "Mode de calcul de la solution de debut non reconnu. " << endl;
                    msg << "Le mode de calcul doit etre l'un des suivants :" << endl;
                    msg << "AB : pour le mode approche de base" << endl;
                    msg << "AC : pour le mode approche avec coupe" << endl;
                    msg << "AH : pour le mode approche heuristique (start par defaut)" << endl;
                    msg << "none : sans calcul de solution initiale" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--coupe") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                coupe = argv[i + 1];
                if (coupe != "G1" && coupe != "G2" && coupe != "T1" && coupe != "T2") {
                    msg << "Type de coupe pour la resolution exacte avec coupe non reconnu. " << endl;
                    msg << "Le type de coupe doit etre l'un des suivants :" << endl;
                    msg << "G1 : heuristique gloutonne 1 (coupe par defaut)" << endl;
                    msg << "G2 : heuristique gloutonne 2" << endl;
                    msg << "T1 : recherche tabou 1" << endl;
                    msg << "T2 : recherche tabou 2" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--verbose") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                try {
                    verbose = stoi(argv[i + 1]);
                }
                catch (...) {
                    verbose = -1;
                }
                if (verbose != 0 && verbose != 1 && verbose != 2) {
                    msg << "Verbose doit etre 0, 1 et 2." << endl;
                    msg << "0 : pas de messages affiches" << endl;
                    msg << "1 : messages synthetiques synthetiques " << endl;
                    msg << "2 : tous les messages affiches" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--repeat_approx") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                try {
                    repeat = stoi(argv[i + 1]);
                }
                catch (...) {
                    repeat = -1;
                }
                if (repeat <= 0) {
                    msg << "repeat_approx doit etre un nombre entier plus grand que zero" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--timeout") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                try {
                    timeout = stod(argv[i + 1]);
                }
                catch (...) {
                    timeout = -1;
                }
                if (timeout <= 0) {
                    msg << "timeout doit etre un nombre plus grand que zero" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--heuristic_steps") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                try {
                    heur = stoi(argv[i + 1]);
                }
                catch (...) {
                    heur = -1;
                }
                if (heur <= 0) {
                    msg << "heuristic_steps doit etre un nombre entier plus grand que zero" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--tol_exact") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                try {
                    tol = stod(argv[i + 1]);
                }
                catch (...) {
                    tol = -1;
                }
                if (tol <= 0 || tol >= 1) {
                    msg << "tol_exact doit etre un nombre strictement compris entre 0 et 1" << endl;
                    throw msg.str();
                }
            }
            else if (opt == "--output") {
                if (i + 1 >= argc) {
                    msg << "Syntaxe invalide. L'option " << opt << " doit avoir une valeur." << endl;
                    throw msg.str();
                }
                try {
                    out_file.open(argv[i + 1]);
                }
                catch (...) {
                    msg << "Impossible d'ouvrir le fichier de sortie " << argv[i + 1] << endl;
                    throw msg.str();
                }
            }
            else {
                msg << "Option " << argv[i] << " non reconnue. Les options sont :" << endl;
                msg << "--mode : mode de resolution" << endl;
                msg << "--start : algorithme approche pour la solution initiale" << endl;
                msg << "--coupe : algorithme de coupe utilise pour le mode exact avec coupe" << endl;
                msg << "--verbose : niveau de messages affiches" << endl;
                msg << "--repeat_approx : nombre de tours de boucle de l'algorithme approche" << endl;
                msg << "--timeout : temps maximal de resolution des algorithmes exacts et de l'algorithe approche avec coupe" << endl;
                msg << "--heuristic_steps : nombre de pas heuristiques de l'algorithme approche heuristique" << endl;
                msg << "--tol_exact : tolerance pour l'algorithme exact" << endl;
                msg << "--output : nom du fichier de sortie (sans espaces)" << endl;
                throw msg.str();
            }
        }

    } catch (string s) {
        cout << s << endl;
        return;
    }

    if (!out_file.is_open()) {
        size_t pos_point = filename.rfind(".");
        if (pos_point == string::npos) {
            out_file.open(filename + ".txt");
        }
        else {
            out_file.open(filename.substr(0, pos_point) + ".txt");
        }
    }

    if (mode == "AB") {
        if (verbose >= 1) {
            cout << "Resolution de l'instance " << filename << " en cours..." << endl;
        }
        
        SolApprocheeBase sol = SolApprocheeBase(&I);
        auto start_time = chrono::high_resolution_clock::now();
        sol.solve(repeat, verbose >= 2);
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

        if (verbose >= 1) {
            cout << "Solution :" << endl;
            cout << sol.meilleure << endl;

            cout << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
        }
        out_file << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
        out_file << sol.meilleure << endl;
    }
    else if (mode == "AC") {
        if (verbose >= 1) {
            cout << "Resolution de l'instance " << filename << " en cours..." << endl;
        }

        SolApprocheeCoupe sol = SolApprocheeCoupe(&I);
        auto start_time = chrono::high_resolution_clock::now();
        sol.solve(repeat, timeout, verbose >= 2);
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

        if (verbose >= 1) {
            cout << "Solution :" << endl;
            cout << sol.meilleure << endl;

            cout << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
        }
        out_file << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
        out_file << sol.meilleure << endl;
    }
    else if (mode == "AH") {
        if (verbose >= 1) {
            cout << "Resolution de l'instance " << filename << " en cours..." << endl;
        }

        SolApprocheeHeuristique sol = SolApprocheeHeuristique(&I);
        auto start_time = chrono::high_resolution_clock::now();
        sol.solve(repeat, heur, verbose >= 2);
        auto end_time = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

        if (verbose >= 1) {
            cout << "Solution :" << endl;
            cout << sol.meilleure << endl;

            cout << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
        }
        out_file << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
        out_file << sol.meilleure << endl;
    }
    else if (mode == "EB" || mode == "EC" || mode == "E2") {
        Solution s;
        Solution* ps = nullptr;

        if (verbose >= 1) {
            cout << "Resolution de l'instance " << filename << " en cours..." << endl;
        }

        chrono::time_point<chrono::high_resolution_clock> start_app;
        chrono::time_point<chrono::high_resolution_clock> end_app;
        chrono::time_point<chrono::high_resolution_clock> start_ex;
        chrono::time_point<chrono::high_resolution_clock> end_ex;
        chrono::microseconds duration;

        if (start == "AB") {
            if (verbose >= 1) {
                cout << "Pre-resolution par l'algorithme approche de base... ";
            }
            SolApprocheeBase sol = SolApprocheeBase(&I);
            start_app = chrono::high_resolution_clock::now();
            sol.solve(repeat, verbose >= 2);
            end_app = chrono::high_resolution_clock::now();
            s = sol.meilleure;
            ps = &s;
            if (verbose >= 1) {
                cout << "ok!" << endl;
            }
        }
        else if (start == "AC") {
            if (verbose >= 1) {
                cout << "Pre-resolution par l'algorithme approche avec coupe... ";
            }
            SolApprocheeCoupe sol = SolApprocheeCoupe(&I);
            start_app = chrono::high_resolution_clock::now();
            sol.solve(repeat, timeout, verbose >= 2);
            end_app = chrono::high_resolution_clock::now();
            s = sol.meilleure;
            ps = &s;
            if (verbose >= 1) {
                cout << "ok!" << endl;
            }
        }
        else if (start == "AH") {
            if (verbose >= 1) {
                cout << "Pre-resolution par l'algorithme approche heuristique... ";
            }
            SolApprocheeHeuristique sol = SolApprocheeHeuristique(&I);
            start_app = chrono::high_resolution_clock::now();
            sol.solve(repeat, heur, verbose >= 2);
            end_app = chrono::high_resolution_clock::now();
            s = sol.meilleure;
            ps = &s;
            if (verbose >= 1) {
                cout << "ok!" << endl;
            }
        }

        if (mode == "EB") {
            SolExacteBase sol = SolExacteBase(&I);
            start_ex = chrono::high_resolution_clock::now();
            sol.solve(ps, tol, timeout, verbose >= 2);
            end_ex = chrono::high_resolution_clock::now();
            
            duration = (start != "none") ? chrono::duration_cast<chrono::microseconds>(end_ex - start_ex + end_app - start_app) : 
                chrono::duration_cast<chrono::microseconds>(end_ex - start_ex);

            if (verbose >= 1) {
                cout << "Solution :" << endl;
                cout << sol.solution << endl;

                cout << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
            }
            out_file << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
            out_file << sol.solution << endl;
        }
        else if (mode == "EC") {
            SolExacteCoupe sol = SolExacteCoupe(&I);
            start_ex = chrono::high_resolution_clock::now();
            sol.solve(ps, tol, timeout, coupe, verbose >= 2);
            end_ex = chrono::high_resolution_clock::now();

            duration = (start != "none") ? chrono::duration_cast<chrono::microseconds>(end_ex - start_ex + end_app - start_app) :
                chrono::duration_cast<chrono::microseconds>(end_ex - start_ex);

            if (verbose >= 1) {
                cout << "Solution :" << endl;
                cout << sol.solution << endl;

                cout << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
            }
            out_file << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
            out_file << sol.solution << endl;
        }
        else if (mode == "E2") {
            SolExacteV2 sol = SolExacteV2(&I);
            start_ex = chrono::high_resolution_clock::now();
            sol.solve(ps, tol, timeout, verbose >= 2);
            end_ex = chrono::high_resolution_clock::now();

            duration = (start != "none") ? chrono::duration_cast<chrono::microseconds>(end_ex - start_ex + end_app - start_app) :
                chrono::duration_cast<chrono::microseconds>(end_ex - start_ex);

            if (verbose >= 1) {
                cout << "Solution :" << endl;
                cout << sol.solution << endl;

                cout << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
            }
            out_file << "Temps de calcul : " << duration.count() / 1000 << " ms" << endl;
            out_file << sol.solution << endl;
        }

    }
}

int main(int argc, char** argv)
{
    /*cli(argc, argv);
    return 0;*/

    /*int indice = 1;
    int indice2 = 1;
    if (argc >= 2) {
        indice = stoi(argv[1]);
    }
    if (argc >= 3) {
        indice2 = stoi(argv[2]);
    }
    tester_A_14(indice, indice2);
    return 0;*/

    tester_coupe_vs_v2(10);
    //tester_instance_0(10);
    return 0;

    //instance de test (type A)
    //ifstream fic("PRP_instances/Instance_1.prp");

    //type A
    ifstream fic("PRP_instances/A_014_ABS1_15_1.prp");

    //type B
    //ifstream fic("../PRP_instances/B_050_instance1.prp");

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
    sol_app_heur.solve(10, 1000, false);
    end = chrono::high_resolution_clock::now();
    auto duration_app_heur = chrono::duration_cast<chrono::microseconds>(end - start);

    /*SolExacteBase sol_ex_base = SolExacteBase(&I);
    sol_ex_base.solve(nullptr, 1e-6, false);*/
    
    /*SolExacteCoupe sol_ex_coupe = SolExacteCoupe(&I);
	sol_ex_coupe.solve(&(sol_app_heur.meilleure), 1e-6, true);*/

    SolExacteV2 sol_ex_v2 = SolExacteV2(&I);
    sol_ex_v2.solve(&(sol_app_heur.meilleure), 1e-6, -1, true);

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
    //cout << "  Exacte coupe : " << sol_ex_coupe.solution.valeur << endl;
    cout << "  Exacte v2 : " << sol_ex_v2.solution.valeur << endl;

    return 0;
}
