#! /usr/env/python
import subprocess
import sys
import os
import pandas as pd
import json
from compare_rf import compare_rf
import matplotlib.pyplot as plt

from eval_utils import parse_lqic, lqic

EVAL_TREES = "../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre"
REF_TREE = "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre"

ALGORITHMS = [["../build/treesearch", "-s", "random", "-a", "nni", "-x", "-c", "-t 4"], ["../build/treesearch", "-s", "random", "-a", "simann", "--simannfactor", "0.0001", "-c", "-t 4"], ["../build/treesearch", "-a", "no", "-x", "-c", "-t 4"], ["../build/treesearch", "-s", "random", "-a", "no", "--clustering", "-t 4"], ["java", "-jar", "../../ASTRAL/Astral/astral.5.6.1.jar", "-t", "0"], ["java", "-jar", "../../ASTRAL/Astral/astral.5.6.1.jar", "-t", "0"]]

REPEAT = 10

df_col_rfnorm = []
df_col_lqic = []

for alg in ALGORITHMS:
    for seed in range(REPEAT):
        try:
            os.remove("out.tre")
        except:
            pass
        if "astral" in alg[2]:
            a = alg + ['-i', EVAL_TREES, "--seed", str(seed)] + ["-o", "out.tre"]
        else:
            a = alg + ['-e', EVAL_TREES, "--seed", str(seed)] + ["-o", "out.tre"]
        print(" ".join(a))
        process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
        out, err = process.communicate()
        (rf_plain, rf_normalized) = compare_rf('out.tre', REF_TREE)
        if "astral" in alg[2]:
            astral_lqic = lqic("out.tre", EVAL_TREES)
            print(astral_lqic)
            df_col_lqic.append(astral_lqic)
        else:
            df_col_lqic.append(parse_lqic(out))
        df_col_rfnorm.append(rf_normalized)
        df = pd.DataFrame({'RF_normalized':df_col_rfnorm, 'LQIC':df_col_lqic})
        print(df)

df = pd.DataFrame({'RF':df_col_rfnorm, 'LQIC':df_col_lqic})
print(df.to_string())
print(df.corr())
df.plot.scatter('LQIC','RF')
plt.savefig('lqic-rf-scatter.png')
