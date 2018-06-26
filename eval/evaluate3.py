#! /usr/env/python
import subprocess
import sys
import os
import timeit
import pandas as pd
import json
from compare_rf import compare_rf
from eval_utils import parse_score, parse_times

def run(algo, pathToEvaluationTrees, pathToReferenceTree):
    a = ["../build/treesearch"] + algo + ["-e", pathToEvaluationTrees]
    print("----------------------------------------")
    print(" ".join(a))
    process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    lqic = parse_score(out, "LQIC", False)
    qpic = parse_score(out, "QPIC", False)
    eqpic= parse_score(out, "EQPIC", False)
    (rf_plain, rf_normalized) = compare_rf('../../out/out.tre', pathToReferenceTree)
    return (rf_normalized, lqic, qpic, eqpic)

def referenceTree(pathToEvaluationTrees, pathToReferenceTree):
    return run(["-a", "no", "--starttree", pathToReferenceTree, "-t", "2"], pathToEvaluationTrees, pathToReferenceTree)

def clustered_lqic(pathToEvaluationTrees, pathToReferenceTree):
    return run(["-a", "nni", "-s", "random", "--clustering", "--treesearchAlgorithmClustered", "simann", "-c", "-x", "-t", "2"], pathToEvaluationTrees, pathToReferenceTree)


def clustered_qpic(pathToEvaluationTrees, pathToReferenceTree):
    return run(["-a", "nni", "-s", "random", "--clustering", "--treesearchAlgorithmClustered", "simann", "-c", "-x", "--objectiveFunction", "qpic", "-t", "2"], pathToEvaluationTrees, pathToReferenceTree)

def clustered_eqpic(pathToEvaluationTrees, pathToReferenceTree):
    return run(["-a", "nni", "-s", "random", "--clustering", "--treesearchAlgorithmClustered", "simann", "-c", "-x", "--objectiveFunction", "eqpic", "-t", "2"], pathToEvaluationTrees, pathToReferenceTree)

data = [["../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre",
         "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre"],
        ["../../data/simulated_datasets_for_sarah/model102M6E.n01.estimated_gene.trees",
         "../../data/simulated_datasets_for_sarah/model102M6E.n01.true_species.tre"],
        ["../../data/ICTC-master/data/Empirical/Avian/avian_all.tre",
         "../../data/ICTC-master/data/Empirical/Avian/avian_reference.tre"],
        ["../../data/evaluating-fast-ml-based/singlegene_trees/single-gene_trees/BoroA6/Best_observed/BoroA6-Best_observed.trees",
         "../../data/evaluating-fast-ml-based/coalescent_species_trees/BoroA6/BoroA6.Best_observed.ASTRAL.tre"],
        ["../../data/evaluating-fast-ml-based/singlegene_trees/single-gene_trees/TarvD7/Best_observed/TarvD7-Best_observed.trees",
         "../../data/evaluating-fast-ml-based/coalescent_species_trees/TarvD7/TarvD7.Best_observed.ASTRAL.tre"],
        ["../../data/evaluating-fast-ml-based/singlegene_trees/single-gene_trees/ChenA4/Best_observed/ChenA4-Best_observed.trees",
         "../../data/evaluating-fast-ml-based/coalescent_species_trees/ChenA4/ChenA4.Best_observed.ASTRAL.tre"],
        ["../../data/evaluating-fast-ml-based/singlegene_trees/single-gene_trees/NagyA1/Best_observed/NagyA1-Best_observed.trees",
         "../../data/evaluating-fast-ml-based/coalescent_species_trees/NagyA1/NagyA1.Best_observed.ASTRAL.tre"]]

col_data = []
col_ref_rf = []
col_ref_lqic = []
col_ref_qpic = []
col_ref_eqpic = []
col_lqic_rf = []
col_lqic_lqic = []
col_lqic_qpic = []
col_lqic_eqpic = []
col_qpic_rf = []
col_qpic_lqic = []
col_qpic_qpic = []
col_qpic_eqpic = []
col_eqpic_rf = []
col_eqpic_lqic = []
col_eqpic_qpic = []
col_eqpic_eqpic = []

for d in data:
    col_data.append(d[0][d[0].rfind('/')+1:d[0].rfind('.')])
    (rf_normalized, lqic, qpic, eqpic) = referenceTree(d[0], d[1])
    col_ref_rf.append(rf_normalized)
    col_ref_lqic.append(lqic)
    col_ref_qpic.append(qpic)
    col_ref_eqpic.append(eqpic)
    os.rename('../../out/out.tre', 'out/'+d[0][d[0].rfind('/')+1:d[0].rfind('.')]+'.ref.tre')
    (rf_normalized, lqic, qpic, eqpic) = clustered_lqic(d[0], d[1])
    col_lqic_rf.append(rf_normalized)
    col_lqic_lqic.append(lqic)
    col_lqic_qpic.append(qpic)
    col_lqic_eqpic.append(eqpic)
    os.rename('../../out/out.tre', 'out/'+d[0][d[0].rfind('/')+1:d[0].rfind('.')]+'.lqic.tre')
    (rf_normalized, lqic, qpic, eqpic) = clustered_qpic(d[0], d[1])
    col_qpic_rf.append(rf_normalized)
    col_qpic_lqic.append(lqic)
    col_qpic_qpic.append(qpic)
    col_qpic_eqpic.append(eqpic)
    os.rename('../../out/out.tre', 'out/'+d[0][d[0].rfind('/')+1:d[0].rfind('.')]+'.qpic.tre')
    (rf_normalized, lqic, qpic, eqpic) = clustered_eqpic(d[0], d[1])
    col_eqpic_rf.append(rf_normalized)
    col_eqpic_lqic.append(lqic)
    col_eqpic_qpic.append(qpic)
    col_eqpic_eqpic.append(eqpic)
    os.rename('../../out/out.tre', 'out/'+d[0][d[0].rfind('/')+1:d[0].rfind('.')]+'.eqpic.tre')

df=pd.DataFrame({'Dataset':col_data, 'Reference Tree RF':col_ref_rf, 'Reference Tree LQIC':col_ref_lqic, 'Reference Tree QPIC':col_ref_qpic, 'Reference Tree EQPIC':col_ref_eqpic, 'LQIC Tree RF':col_lqic_rf, 'LQIC Tree LQIC':col_lqic_lqic, 'LQIC Tree QPIC':col_lqic_qpic, 'LQIC Tree EQPIC':col_lqic_eqpic, 'QPIC Tree RF':col_qpic_rf, 'QPIC Tree LQIC':col_qpic_lqic, 'QPIC Tree QPIC':col_qpic_qpic, 'QPIC Tree EQPIC':col_qpic_eqpic, 'EQPIC Tree RF':col_eqpic_rf, 'EQPIC Tree LQIC':col_eqpic_lqic, 'EQPIC Tree QPIC':col_eqpic_qpic, 'EQPIC Tree EQPIC':col_eqpic_eqpic})
df=df[['Dataset', 'Reference Tree RF', 'LQIC Tree RF', 'QPIC Tree RF', 'EQPIC Tree RF', 'Reference Tree LQIC', 'LQIC Tree LQIC', 'QPIC Tree LQIC', 'EQPIC Tree LQIC', 'Reference Tree QPIC', 'LQIC Tree QPIC', 'QPIC Tree QPIC', 'EQPIC Tree QPIC', 'Reference Tree EQPIC', 'LQIC Tree EQPIC', 'QPIC Tree EQPIC', 'EQPIC Tree EQPIC']]


print(df)
df.to_csv('df.csv')
