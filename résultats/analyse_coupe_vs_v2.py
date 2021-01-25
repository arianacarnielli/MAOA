# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 00:44:22 2021

@author: arian
"""

import pandas as pd
import matplotlib.pyplot as plt
#import numpy as np

df = pd.read_csv("resultats_coupe_vs_v2.csv")

print("SolExacteCoupe : ", end="")
print(df["Coupe Temps"].mean() / 1000, df["Coupe Temps"].std() / 1000)
print("SolExacteV2 : ", end="")
print(df["V2 Temps"].mean() / 1000, df["V2 Temps"].std() / 1000)

plt.close("all")

fig, ax = plt.subplots(figsize=(5, 2))
ax.grid(True)
ax.set_axisbelow(True)
ax.plot(df["Coupe Temps"].values / 1000, [0]*10, 'o')
ax.plot(df["V2 Temps"].values / 1000, [1]*10, 'o')
ax.set_ylim([-0.5, 1.5])
# ax.set_xlim([10**int(np.log10(min(df["Coupe Temps"].values.min(), df["V2 Temps"].values.min()) / 1000)),\
#              10**int(1+np.log10(max(df["Coupe Temps"].values.max(), df["V2 Temps"].values.max()) / 1000))])
#ax.set_xscale("log")