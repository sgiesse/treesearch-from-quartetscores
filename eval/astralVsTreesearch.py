#! /usr/env/python
import subprocess
import sys
import os
import timeit
import pandas as pd
import json
from compare_rf import compare_rf

from eval_utils import parse_lqic, lqic

# Provide a json-config file as first argument. Sample:
'''
{
    "data": [["../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre",
              "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre"],
             ["../../data/simulated_datasets_for_sarah/model102M6E.n01.estimated_gene.trees",
              "../../data/simulated_datasets_for_sarah/model102M6E.n01.true_species.tre"]],
    "tool": [["../build/treesearch", "-s", "random", "-a", "nni", "-x", "-c", "--clustering", "--treesearchAlgorithmClustered", "simann", "-t 4"], ["java", "-jar", "../../ASTRAL/Astral/astral.5.6.1.jar"]],
    "repeat": 1,
}
'''
with open(sys.argv[1]) as json_config_file:
    config = json.load(json_config_file)

df_col_runtime = []
df_col_dataset = []
df_col_tool = []
df_col_rf = []
df_col_rfnorm = []
df_col_lqic = []

for d in config["data"]:
    for tool in config["tool"]:
        for seed in range(config["repeat"]):
            try:
                os.remove("out.tre")
            except:
                pass
            a = []
            if "treesearch" in tool[0]:
                a = tool + ['-e', d[0], "--seed", str(seed)] + ["-o", "out.tre"]
            elif "astral" in tool[2]:
                a = tool + ['-i', d[0], "--seed", str(seed)] + ["-o", "out.tre"]
            else:
                print("Unknown Tool")
                sys.exit()
            print("----------------------------------------")
            print(" ".join(a))
            start = timeit.default_timer()
            process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
            out, err = process.communicate()
            end = timeit.default_timer()
            runtime = end - start
            (rf_plain, rf_normalized) = compare_rf('out.tre', d[1])
            #df_col_tool.append(" ".join(a))
            if "treesearch" in tool[0]:
                df_col_tool.append("treesearch")
                df_col_lqic.append(parse_lqic(out))
            elif "astral" in tool[2]:
                df_col_tool.append("astral")
                df_col_lqic.append(lqic("out.tre", d[0]))
            df_col_dataset.append(d[0][d[0].rfind('/')+1:d[0].rfind('.')])
            df_col_runtime.append(runtime)
            df_col_rf.append(rf_plain)
            df_col_rfnorm.append(rf_normalized)
            df = pd.DataFrame({'Dataset':df_col_dataset, 'Tool':df_col_tool, 'runtime':df_col_runtime,'RF':df_col_rf, 'RF_normalized':df_col_rfnorm, 'LQIC':df_col_lqic})
            print(df)

df = pd.DataFrame({'Dataset':df_col_dataset, 'Tool':df_col_tool,  'runtime':df_col_runtime,'RF':df_col_rf, 'RF_normalized':df_col_rfnorm, 'LQIC':df_col_lqic})

print(df.to_string())

#print(df.groupby(['Dataset', 'StartTree']).mean())
#print(df.groupby(['Dataset', 'StartTree']).var())

#df_stats = df.groupby(['Dataset', 'StartTree', 'Algorithm']).agg(['min', 'max', 'mean', 'var', 'std'])
df_stats = df.groupby(['Dataset', 'Tool']).agg(['mean', 'var'])

print(df_stats)
df_stats.to_csv('df.csv')

#os.remove(progfile)
