import pandas as pd
import numpy as np
import dendropy
import sys
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import rf_distance

#def num_taxa_in_trees(pathToTrees):

pre1 = "../../data/evaluating-fast-ml-based/singlegene_trees/single-gene_trees/"
pathList = [pre1+"BoroA6/Best_observed/BoroA6-Best_observed.trees",
            pre1+"ChenA4/Best_observed/ChenA4-Best_observed.trees",
            pre1+"JarvD5a/Best_observed/Ja00:00rvD5a-Best_observed.trees",
            pre1+"JarvD5b/Best_observed/JarvD5b-Best_observed.trees",
            pre1+"MisoA2/Best_observed/MisoA2-Best_observed.trees",
            pre1+"MisoD2a/Best_observed/MisoD2a-Best_observed.trees",
            pre1+"MisoD2b/Best_observed/MisoD2b-Best_observed.trees",
            pre1+"NagyA1/Best_observed/NagyA1-Best_observed.trees",
            pre1+"PrumD6/Best_observed/PrumD6-Best_observed.trees",
            pre1+"ShenA9/Best_observed/ShenA9-Best_observed.trees",
            pre1+"SongD1/Best_observed/SongD1-Best_observed.trees",
            pre1+"StruA5/Best_observed/StruA5-Best_observed.trees",
            pre1+"TarvD7/Best_observed/TarvD7-Best_observed.trees",
            pre1+"WhelA7/Best_observed/WhelA7-Best_observed.trees",
            pre1+"WickA3/Best_observed/WickA3-Best_observed.trees",
            pre1+"WickD3a/Best_observed/WickD3a-Best_observed.trees",
            pre1+"WickD3b/Best_observed/WickD3b-Best_observed.trees",
            pre1+"XiD4/Best_observed/XiD4-Best_observed.trees",
            pre1+"YangA8/Best_observed/YangA8-Best_observed.trees"]
"""
#pathToTrees = sys.argv[1]
for pathToTrees in pathList:
    print(pathToTrees[pathToTrees.rfind("/")+1:])
    tree_list = dendropy.TreeList.get(path=pathToTrees, schema="newick")

    count_taxa = [0] * (len(tree_list.taxon_namespace) + 1)

    for tree in tree_list:
        tree.prune_leaves_without_taxa()
        count_taxa[len(tree.leaf_nodes())] = count_taxa[len(tree.leaf_nodes())] + 1


    df = pd.DataFrame({"Number_of_Trees": count_taxa})
    df.index.name = 'Number_of_Taxa'

    #print(df)
    df=df[(df["Number_of_Trees"]).cumsum() > 0]
    print(df)
    #ax = df.plot.bar()
    ax=df.plot(kind='bar',alpha=0.75, rot=45, figsize=(8.3, 11.7))
    #ax=df.plot(kind='bar',alpha=0.75)
    plt.savefig("out/"+pathToTrees[pathToTrees.rfind("/")+1:]+'.png')
"""

import os
folder = "../../data/evaluating-fast-ml-based/singlegene_trees/single-gene_trees/WhelA7/Best_observed/"
files = []

for file in os.listdir(folder):
    #print(file)
    if file.endswith(".tre"):
        files.append(os.path.join(folder, file))
        #print(os.path.join(folder, file))
print(len(files))

rfs = []
error = 0
for tree1 in files:

    for tree2 in files:
        try:
            rf = rf_distance.normalized_RF_distance(tree1, tree2)
            print(rf)
            rfs.append(rf)
        except:
            print("Error computing RF")
            error = error + 1

rfs = np.array(rfs, dtype=np.float64)
#print(rfs)
print("number of errors (probably multfurcating tree): " + str(error))
print("mean norm. RF distance: " + str(np.mean(rfs)))
print("var: " + str(np.var(rfs)))
