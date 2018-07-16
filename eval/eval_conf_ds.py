#! /usr/env/python
import subprocess
import sys
import os
import timeit
import pandas as pd
import json
from compare_rf import compare_rf
from eval_utils import parse_score, parse_times

# Provide a json-config file as first argument. Sample:
'''
{
    "data": [["../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre",
              "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre"],
             ["../../data/simulated_datasets_for_sarah/model102M6E.n01.estimated_gene.trees",
              "../../data/simulated_datasets_for_sarah/model102M6E.n01.true_species.tre"]],
    "exe": "../build/uquest",
    "algorithms": [["custom", "-s", "random", "-a", "combo"],["custom", "-s", "random", "-a", "combo", "-x"], ["custom", "-s", "random", "-a", "combo", "-x", "-c"], ["custom", "-s", "random", "-a", "combo", "-x", "--clustering"], ["custom", "-s", "random", "-a", "combo", "-x", "--clustering", "-c"]]]
    "repeat": 1,
    "rename": {"runtime_total":"time","RF_normalized":"RF_norm"},
    "replace": {"model102M6E.n01.estimated_gene":"model102M6E.n01"},
    "columns":["time","LQIC","RF_norm"]
}
'''
with open(sys.argv[1]) as json_config_file:
    config = json.load(json_config_file)

def make_treesearch(d, algo, seed):
    a = [config["exe"]] + algo + ["-e", d[0], "-o", "../../out/out.tre", "--seed", str(seed)]
    print("----------------------------------------")
    print(" ".join(a))
    process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    lqic = parse_score(out, "LQIC", False)
    qpic = parse_score(out, "QPIC", False)
    eqpic= parse_score(out, "EQPIC", False)
    (t_count, t_start, t_search1, t_search2, t_total) = parse_times(out)
    (rf_plain, rf_normalized) = compare_rf('../../out/out.tre', d[1])
    return (t_count, t_start, t_search1, t_search2, t_total, lqic, qpic, eqpic, rf_plain, rf_normalized)

df_col_runtime_start = []
df_col_runtime_count = []
df_col_runtime_search = []
df_col_lqic = []
df_col_qpic = []
df_col_eqpic = []
df_col_dataset = []
df_col_algo = []
df_col_rf = []
df_col_rfnorm = []


progfile = sys.argv[1]+".progress.csv"

try:
    if os.path.getmtime(config["exe"]) < os.path.getmtime(progfile) and os.path.getmtime(sys.argv[1]) < os.path.getmtime(progfile):
        dfp = pd.read_csv(progfile)
        df_col_runtime_start = dfp["runtime_start"].tolist()
        df_col_runtime_count = dfp["runtime_count"].tolist()
        df_col_runtime_search = dfp["runtime_search"].tolist()
        df_col_lqic = dfp["LQIC"].tolist()
        df_col_qpic = dfp["QPIC"].tolist()
        df_col_eqpic = dfp["EQPIC"].tolist()
        df_col_dataset = dfp["Dataset"].tolist()
        df_col_algo = dfp["Algorithm"].tolist()
        df_col_rf = dfp["RF"].tolist()
        df_col_rfnorm = dfp["RF_normalized"].tolist()
except FileNotFoundError:
    pass

i = 0
for d in config["data"]:
    for algo in config["algorithms"]:
        for seed in range(config["repeat"]):
            if i == len(df_col_runtime_search):
                (t_count, t_start, t_search1, t_search2, t_total, lqic, qpic, eqpic, rf_plain, rf_normalized) = make_treesearch(d, algo, seed)
                print("took " + str(t_total) + "s")
                print("sum lqic: " + str(lqic))
                print("plain RF distance: " + str(rf_plain))
                print("normalized RF distance: " + str(rf_normalized))
                df_col_dataset.append(d[0][d[0].rfind('/')+1:d[0].rfind('.')])
                df_col_algo.append(" ".join(algo[1:]))
                df_col_runtime_start.append(t_start)
                df_col_runtime_count.append(t_count)
                df_col_runtime_search.append(t_search1+t_search2)
                df_col_lqic.append(lqic)
                df_col_qpic.append(qpic)
                df_col_eqpic.append(eqpic)
                df_col_rf.append(rf_plain)
                df_col_rfnorm.append(rf_normalized)
                df = pd.DataFrame({'Dataset':df_col_dataset, 'Algorithm':df_col_algo, 'runtime_start':df_col_runtime_start, 'runtime_count':df_col_runtime_count, 'runtime_search':df_col_runtime_search, 'LQIC': df_col_lqic, 'QPIC': df_col_qpic, 'EQPIC': df_col_eqpic, 'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
                df.to_csv(progfile)
            else:
                print("skip, because already computed")
            i = i+1

            df = pd.DataFrame({'Dataset':df_col_dataset, 'Algorithm':df_col_algo, 'runtime_start':df_col_runtime_start, 'runtime_count':df_col_runtime_count, 'runtime_search':df_col_runtime_search, 'LQIC': df_col_lqic, 'QPIC': df_col_qpic, 'EQPIC': df_col_eqpic, 'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
            df['runtime_total'] = df['runtime_start']+df['runtime_count']+df['runtime_search']
print(df.to_string())

#print(df.groupby(['Dataset', 'StartTree']).mean())
#print(df.groupby(['Dataset', 'StartTree']).var())

#df_stats = df.groupby(['Dataset', 'StartTree', 'Algorithm']).agg(['min', 'max', 'mean', 'var', 'std'])
df_stats = df.groupby(['Dataset', 'Algorithm']).agg(['mean', 'var'])

print(df_stats)
df_stats.to_csv('df.csv')

#dfp=df[["Algorithm", "Dataset", "runtime_total", "LQIC", "RF_normalized"]].groupby(['Algorithm','Dataset']).mean().unstack().swaplevel(axis=1)
dfp=df
dfp=dfp.rename(config["rename"], axis="columns")
dfp=dfp.replace(config["replace"])
dfp=dfp[["Algorithm", "Dataset"]+config["columns"]]
dfp=dfp.round(2)
dfp=dfp.groupby(['Algorithm','Dataset']).mean().unstack().swaplevel(axis=1)
dfp.sort_index(level=0, axis=1, inplace=True)
print(dfp)
dfp.to_csv('dfp.csv')

print(dfp.to_latex())
