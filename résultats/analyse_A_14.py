# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 01:25:52 2021

@author: arian
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("resultats_A_14.csv")

plt.close("all")

# Histogrammes de temps de calcul
L = [(df["Approx Base Temps"] / 1000, "SolApprocheeBase"),\
     (df["Approx Coupe Temps"] / 1000, "SolApprocheeCoupe"),\
     (df["Approx Heur Temps"] / 1000, "SolApprocheeHeuristique"),\
     (df["Exacte Coupe Temps"] / 1000, "SolExacteCoupe")]
    
minTime = 10**np.floor(np.log10(min(col.min() for col, _ in L)))
maxTime = 10**np.ceil(np.log10(max(col.max() for col, _ in L)))

totalTime = sum(col.sum() for col, _ in L)
totalTimeD = int(totalTime/(3600*24))
totalTimeH = int((totalTime - totalTimeD*3600*24)/3600)
totalTimeM = int((totalTime - totalTimeD*3600*24 - 3600*totalTimeH)/60)
totalTimeS = totalTime - totalTimeD*3600*24 - 3600*totalTimeH - 60*totalTimeM
print("Temps total de calcul : ", totalTimeD, "j", totalTimeH, "h", totalTimeM, "min", totalTimeS, "s", sep="")

totalTimeExa = L[-1][0].sum()
totalTimeExaD = int(totalTimeExa/(3600*24))
totalTimeExaH = int((totalTimeExa - totalTimeExaD*3600*24)/3600)
totalTimeExaM = int((totalTimeExa - totalTimeExaD*3600*24 - 3600*totalTimeExaH)/60)
totalTimeExaS = totalTimeExa - totalTimeExaD*3600*24 - 3600*totalTimeExaH - 60*totalTimeExaM
print("Temps de calcul de la solution exacte : ", totalTimeExaD, "j", totalTimeExaH, "h", totalTimeExaM, "min", totalTimeExaS, "s", sep="")
print("Proportion du temps passé sur la solution exacte : {:.1f}%".format(100*totalTimeExa/totalTime))

print()

for vals, name in L:
    print(name + " : ", end="")
    print(vals.mean(), vals.std(), vals.quantile(0.25), vals.median(), vals.quantile(0.75))

for vals, name in L:
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.set_xscale("log")
    ax.hist(vals, bins=np.logspace(np.log10(vals.min()), np.log10(vals.max()), 10))
    ax.set_xlim([minTime, maxTime])
    ax.set_title(name)
    ax.set_xlabel("Temps (en secondes)")
    fig.savefig(name + "_temps.pdf")

print()
    
# Histogrammes d'écarts
L = [(df["Approx Base Valeur"], "SolApprocheeBase"),\
     (df["Approx Coupe Valeur"], "SolApprocheeCoupe"),\
     (df["Approx Heur Valeur"], "SolApprocheeHeuristique")]
 
L = [(a / df["Exacte Coupe Valeur"] - 1, b) for a, b in L]
    
minEcart = min(col.min() for col, _ in L)
maxEcart = max(col.max() for col, _ in L)

for vals, name in L:
    print(name + " : ", end="")
    print(100*vals.mean(), 100*vals.std(), 100*vals.quantile(0.25), 100*vals.median(), 100*vals.quantile(0.75))

for vals, name in L:
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.grid(True)
    ax.set_axisbelow(True)
    ax.hist(vals, bins=15)
    ax.set_xlim([minEcart-0.01, maxEcart+0.01])
    ax.set_title(name)
    ax.set_xlabel("Écart relatif par rapport à la solution exacte")
    ax.set_xticklabels(["{:.0f}%".format(100*val) for val in ax.get_xticks()])
    fig.savefig(name + "_ecart.pdf")
    
# Graphes de meilleure résolution
L = [(df["Approx Base Valeur"], "SolApprocheeBase"),\
     (df["Approx Coupe Valeur"], "SolApprocheeCoupe"),\
     (df["Approx Heur Valeur"], "SolApprocheeHeuristique")]
minVals = np.minimum(np.minimum(L[0][0].values, L[1][0].values), L[2][0].values)

print()
for vals, name in L:
    print(name, (vals.values == minVals).sum())
    print(name, (vals.values == df["Exacte Coupe Valeur"].values).sum())
