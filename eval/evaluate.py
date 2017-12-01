#! /usr/env/python
import subprocess
import sys
import os
import timeit
from compare_rf import compare_rf


TREESEARCH_EXCTBL = '../build/treesearch'


data = [('../../data/ICTC-master/data/Empirical/Avian/avian_all.tre',
         '../../data/ICTC-master/data/Empirical/Avian/avian_reference.tre'),
        ('../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre',
         '../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre'),
        ('../../data/simulated_datasets_for_sarah/model102M6E.n01.estimated_gene.trees',
         '../../data/simulated_datasets_for_sarah/model102M6E.n01.true_species.tre'),
        ('../../data/simulated_datasets_for_sarah/model502M6E.n01.estimated_gene.trees',
         '../../data/simulated_datasets_for_sarah/model502M6E.n01.true_species.tre'),
        ('../../data/simulated_datasets_for_sarah/model1002M6E.n01.estimated_gene.trees',
         '../../data/simulated_datasets_for_sarah/model1002M6E.n01.true_species.tre')]

data = data[1:3]


args1 = [TREESEARCH_EXCTBL, '-s', 'stepwiseaddition', '-a', 'nni', '-l', 'Info']
args2 = [TREESEARCH_EXCTBL, '-s', 'stepwiseaddition', '-a', 'spr', '-l', 'Info']
args3 = [TREESEARCH_EXCTBL, '-s', 'stepwiseaddition', '-a', 'combo', '-l', 'Info']
args4 = [TREESEARCH_EXCTBL, '-s', 'random', '-a', 'nni', '-l', 'Info']
args5 = [TREESEARCH_EXCTBL, '-s', 'random', '-a', 'spr', '-l', 'Info']
args6 = [TREESEARCH_EXCTBL, '-s', 'random', '-a', 'combo', '-l', 'Info']

configs = [args1, args2, args3, args4, args5, args6]

for d in data:
    for arguments in configs:
        a = arguments + ['-e', d[0], '-r', d[1]]
        print("----------------------------------------")
        print(" ".join(a))
        start = timeit.default_timer()
        process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = process.communicate()
        end = timeit.default_timer()
        print( "took " + str(end - start) + "s")
        outString = out.decode("UTF-8")
        str_pos = outString.find("INFO Sum lqic final Tree:",0)
        if str_pos != -1:
            start = str_pos + 25
            end = outString.find("\n",str_pos)
            lqic = float(outString[start:end])
            print("sum lqic: " + str(lqic))
        (rf_plain, rf_normalized) = compare_rf('../../out/out.tre', d[1])
        print("plain RF distance: " + rf_plain)
        print("normalized RF distance: " + rf_normalized)

