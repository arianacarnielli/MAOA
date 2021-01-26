# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 23:39:17 2021

@author: arian
"""

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("resultats_app_heuristique.csv", index_col=0);

classes = ["A_014", "A_050", "A_100", "B_050", "B_100"]

data = [df.filter(axis = 0, regex = c)["Temps"] / 1000 for c in classes]

plt.close("all")

fig, ax = plt.subplots(figsize = (5, 5))
ax.grid(True)
ax.set_axisbelow(True)
ax.set_yscale("log")
ax.boxplot(data, whis = "range", widths = 1/3, patch_artist = True,\
           labels = classes, showmeans = True,\
           meanprops = dict(markerfacecolor="k", marker="o", markeredgewidth=0),\
           medianprops = dict(linewidth=2, color="C3"),\
           boxprops = dict(color = "C9", facecolor="C9"))
ax.set_xlabel("Catégorie d'instance")
ax.set_ylabel("Temps de résolution (en s)")
ax.set_title("Temps de résolution de SolApprocheeHeuristique")
fig.savefig("resultats_app_heuristique.pdf")