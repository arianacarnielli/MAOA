import numpy as np
import graphviz as gv

def extractNodesFromPRP(inst):
    """
    :param inst: chaine de caractères représentant le nom de l'instance
    :return: la liste des positions des noeuds (liste de paires)
    """
    l = []
    file = open(inst + ".prp", "r")
    line = "rien"

    while line != "d\n" :
        line = file.readline()

        newLine = line.split(" ")
        if( newLine[0].isdigit() ) :
            l.append((int(newLine[1]), int(newLine[2])))
    return l

def show(inst) :
    """

    :param inst: chaine de caractères représentant le nom de l'instance
    :return: rien
    affiche le graphe de la solution de l'instance (fichier txt) et le sauve en format pdf
    """
    g = gv.Digraph(name= inst +"graph", engine ='neato', format = "pdf")

    l= extractNodesFromPRP(inst)

    for n in range(len(l)) :
        g.node(str(n), pos = str(l[n][0]/50) + ',' + str(l[n][1]/50)+ '!')

    file = open(inst + ".txt", "r")
    line = "rien"
    t = -1

    infoClient = [None] * len(l)

    while line != '' :
        line = file.readline()

        if line.startswith( "Pas de temps ") :
            t += 1

        if line.startswith("  0")  :
            prevS = 0
            nline = line.split(" ")
            for s in nline :
                if s.isdigit() and s != "0":
                    g.edge(str(prevS), str(s), color = str(t/len(l))+" 1 1", headlabel = "t=" + str(t) + ", I=" + str(infoClient[int(s)][0]) + ", q=" +str(infoClient[int(s)][1]), fontcolor=str(t/len(l)*0.7)+" 0.7 0.7", fontsize = "10")
                    prevS = s
            g.edge(str(prevS), "0", color=str(t / len(l)) + " 1 1")


        if line.startswith("  Client") :
            nline = line.split(' ')
            nline = [i for i in nline if i != '']
            infoClient[int(nline[1])] = (int(float(nline[5])), int(float(nline[9])))

    g.view()
    g.save(inst+"graph")
    return

show("A_014_ABS1_15_1")