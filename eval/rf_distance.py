#! /usr/env/python
import subprocess
import sys
import os
import dendropy
"""
def taxa_to_prune(tree1, tree2):
    prune_labels = []
    for taxon1 in tree1.taxon_namespace:
        found = False
        for taxon2 in tree2.taxon_namespace:
            if str(taxon1) == str(taxon2):
                found = True
                break
        if not found:
            prune_labels.append(str(taxon1))
    for taxon2 in tree2.taxon_namespace:
        found = False
        for taxon1 in tree1.taxon_namespace:
            if str(taxon2) == str(taxon1):
                found = True
                break
        if not found:
            prune_labels.append(str(taxon2))
    return prune_labels
"""
def taxa_to_prune(tree1, tree2):
    prune = []
    for taxon1 in tree1.taxon_namespace:
        found = False
        for taxon2 in tree2.taxon_namespace:
            if str(taxon1) == str(taxon2):
                found = True
                break
        if not found:
            prune.append(taxon1)
    for taxon2 in tree2.taxon_namespace:
        found = False
        for taxon1 in tree1.taxon_namespace:
            if str(taxon2) == str(taxon1):
                found = True
                break
        if not found:
            prune.append(taxon2)
    return prune


def normalized_RF_distance(tree_1_path, tree_2_path):
    print("tree_1_path: " + tree_1_path)
    print("tree_2_path: " + tree_2_path)

    # remove previous intermediary files
    if os.path.exists('RAxML_info.TEST'):
        os.remove('RAxML_info.TEST')
    if os.path.exists('RAxML_RF-Distances.TEST'):
        os.remove('RAxML_RF-Distances.TEST')
    if os.path.exists('trees'):
        os.remove('trees')

    # unroot the trees
    tree_1_unrooted = dendropy.Tree.get(
        path=tree_1_path,
        schema='newick',
        rooting='default-rooted')
    tree_1_unrooted.deroot()
    tree_2_unrooted = dendropy.Tree.get(
        path=tree_2_path,
        schema='newick',
        rooting='default-rooted')
    tree_2_unrooted.deroot()

    # Retry unrooting if the previous didn't work
    if tree_1_unrooted.is_rooted:
        tree_1_unrooted.is_rooted = False
        tree_1_unrooted.update_bipartitions()
    if tree_2_unrooted.is_rooted:
        tree_2_unrooted.is_rooted = False
        tree_2_unrooted.update_bipartitions()
    if tree_1_unrooted.is_rooted or tree_2_unrooted.is_rooted:
        print("WHAT THE FUCK")
        sys.exit()

    # prune taxa that are only present in one of the trees, but not in the other one
    pruneme = taxa_to_prune(tree_1_unrooted, tree_2_unrooted)
    tree_1_unrooted.prune_taxa(pruneme)
    tree_2_unrooted.prune_taxa(pruneme)

    tree_2_unrooted.migrate_taxon_namespace(tree_1_unrooted.taxon_namespace)

    internodes=0
    for edge in tree_1_unrooted.edges():
        if edge.is_internal():
            internodes = internodes + 1
    for edge in tree_2_unrooted.edges():
        if edge.is_internal():
            internodes = internodes + 1
    rf = dendropy.calculate.treecompare.symmetric_difference(tree_1_unrooted, tree_2_unrooted)
    rf_normalized=rf/internodes
    return rf_normalized
