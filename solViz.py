# -*- coding: utf-8 -*-

import graphviz as gv
import re
import sys
from matplotlib import cm

def extractNodesFromPRP(inst):
    """
    :param inst: chaine de caractères représentant le nom de l'instance
    :return: la liste des positions des noeuds (liste de paires)
    """
    l = []
    with open(inst + ".prp", "r") as file:
        for line in file:
            newLine = line.split(" ")
            if newLine[0].isdigit():
                l.append((int(newLine[1]), int(newLine[2])))
            elif newLine[0] == "d\n":
                break
    return l

def show(inst, scale = 10) :
    """
    :param inst: chaine de caractères représentant le nom de l'instance
    :return: rien
    affiche le graphe de la solution de l'instance (fichier txt) et le sauve
    en format pdf
    """
    g = gv.Digraph(name = inst +"graph", engine = 'neato', format = "pdf")

    l = extractNodesFromPRP(inst)

    lx = [x for x, _ in l]
    ly = [y for _, y in l]
    a_x = max(lx) - min(lx)
    a_y = max(ly) - min(ly)

    for n, (x, y) in enumerate(l):
        g.node(str(n), pos = str(x / a_x * scale) + ',' + str(y / a_y * scale) + '!')
        
    t = -1
    tMax = 1
    infoClient = [None] * len(l)

    with open(inst + ".txt", "r") as file:
        for line in file:
            if line.startswith("Instance"):
                match = re.match("Instance avec (\d+) clients et (\d+) pas de temps", line)
                nbClients, tMax = match.groups()
                nbClients = int(nbClients)
                tMax = int(tMax)
            
            if line.startswith("Pas de temps "):
                t += 1
    
            if line.startswith("  0"):
                color_tuple = cm.plasma(int(t/(tMax-1) * cm.plasma.N))[:-1]
                color_rgba = "#" + "".join(["{:2x}".format(int(c*255)) for c in color_tuple])
                color_text = "#" + "".join(["{:2x}".format(int(c*255*0.7)) for c in color_tuple])
                
                prevS = 0
                nline = line.split(" ")
                for s in nline :
                    if s.isdigit() and s != "0":
                        g.edge(str(prevS), str(s),\
                            color = color_rgba,\
                            headlabel = "t=" + str(t) + ", I=" + str(infoClient[int(s)][0]) + ", q=" +str(infoClient[int(s)][1]),\
                            fontcolor = color_text,\
                            fontsize = "10")
                        prevS = s
                g.edge(str(prevS), "0", color = color_rgba)
    
    
            if line.startswith("  Client") :
                nline = line.split(' ')
                nline = [i for i in nline if i != '']
                infoClient[int(nline[1])] = (float(nline[5]), float(nline[9]))

    g.view()
    g.save(inst+"graph")
    return

if __name__ == "__main__":
    msg = "Un probleme est survenu dans le programme " + sys.argv[0] + " :\n"
    try:
        # Vérification du nombre minimal d'arguments
        if len(sys.argv)==1:
            msg += "le programme a besoin du nom du fichier prp (sans extension ni espaces)\n"
            raise Exception(msg)
            
        filename = sys.argv[1]
        scale = 10
        
        # Vérification des arguments optionnels
        i = 2
        while i < len(sys.argv):
            if sys.argv[i] == "--scale":
                i += 1
                try:
                    scale = int(sys.argv[i])
                    if scale <= 0:
                        raise Exception()
                except:
                    msg += "Argument scale invalide\n"
                    msg += "scale doit être un entier strictement positif\n"
                    raise Exception(msg)
                    
            else:
                msg += "Option " + sys.argv[i] + " non reconnue\n"
                raise Exception(msg)
                
            i += 1
        
        # Exécution
        try:
            show(filename, scale)
        except Exception as e:
            msg += str(e)
            raise Exception(msg)
            
    except Exception as e:
        print(e)