#! /usr/env/python
import subprocess
import sys
import os
import timeit
import pandas as pd
import json
from compare_rf import compare_rf

# Provide a json-config file as first argument. Sample:
'''
{
    "data": [["../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre",
              "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre"],
             ["../../data/simulated_datasets_for_sarah/model102M6E.n01.estimated_gene.trees",
              "../../data/simulated_datasets_for_sarah/model102M6E.n01.true_species.tre"]],
    "exe": "../build/treesearch",
    "args": [ "-l", "Info"],
    "starttrees": [["-s", "random"],["-s","stepwiseaddition"]],
    "file_starttree": "../../out/starttree.tre",
    "algorithms": [["-a", "nni"], ["-a", "combo", "-x"]],
    "repeat": 1,
    "repeatStartTreeMethod": 1
}
'''
with open(sys.argv[1]) as json_config_file:
    config = json.load(json_config_file)


def parse_lqic(out):
    outString = out.decode("UTF-8")
    str_pos = outString.find("INFO Sum lqic final Tree:",0)
    if str_pos != -1:
        start = str_pos + 25
        end = outString.find("\n",str_pos)
        lqic = float(outString[start:end])
        return lqic
    else:
        print(out.decode("UTF-8"))
        raise RuntimeError("no LQIC found")

def make_starttree(d, file_st, args_st):
    process = subprocess.Popen([config["exe"],'-e', d[0], '-r', d[1],  '-a', 'no', '-o', file_st]+args_st,
                                   stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    if (process.returncode != 0): 
        print(out.decode("UTF-8"))
        raise RuntimeError("")
    lqic = parse_lqic(out)
    print("----------------------------------------")
    print("sum lqic of startTree: " + str(lqic))

def make_treesearch(d, file_st, algo):
    a = [config["exe"]] + config["args"] + ['-e', d[0], '-r', d[1], '-t', file_st] + algo
    print("----------------------------------------")
    print(" ".join(a))
    start = timeit.default_timer()
    process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    end = timeit.default_timer()
    runtime = end - start
    lqic = parse_lqic(out)
    (rf_plain, rf_normalized) = compare_rf('../../out/out.tre', d[1])
    return (runtime, lqic, rf_plain, rf_normalized)



df_col_runtime = []
df_col_lqic = []
df_col_dataset = []
df_col_starttree = []
df_col_algo = []
df_col_rf = []
df_col_rfnorm = []

progfile = sys.argv[1]+".progress.csv"

try:
    dfp = pd.read_csv(progfile)
    df_col_runtime = dfp["runtime"].tolist()
    df_col_lqic = dfp["LQIC"].tolist()
    df_col_dataset = dfp["Dataset"].tolist()
    df_col_starttree = dfp["StartTree"].tolist()
    df_col_algo = dfp["Algorithm"].tolist()
    df_col_rf = dfp["RF"].tolist()
    df_col_rfnorm = dfp["RF_normalized"].tolist()
except FileNotFoundError:
    pass

i = 0
for d in config["data"]:
    for args_st in config["starttrees"]:
        for _i in range(config["repeatStartTreeMethod"]):
            make_starttree(d, config["file_starttree"], args_st)
            for algo in config["algorithms"]:
                for _j in range(config["repeat"]):
                    if i == len(df_col_runtime):
                        (runtime, lqic, rf_plain, rf_normalized) = make_treesearch(d, config["file_starttree"], algo)
                        print("took " + str(runtime) + "s")
                        print("sum lqic: " + str(lqic))
                        print("plain RF distance: " + str(rf_plain))
                        print("normalized RF distance: " + str(rf_normalized))
                        df_col_dataset.append(d[0][d[0].rfind('/')+1:d[0].rfind('.')])
                        df_col_starttree.append(args_st[1])
                        df_col_algo.append(" ".join(algo[1:]))
                        df_col_runtime.append(runtime)
                        df_col_lqic.append(lqic)
                        df_col_rf.append(rf_plain)
                        df_col_rfnorm.append(rf_normalized)
                        df = pd.DataFrame({'Dataset':df_col_dataset, 'StartTree':df_col_starttree, 'Algorithm':df_col_algo, 'runtime':df_col_runtime, 'LQIC': df_col_lqic,'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
                        df.to_csv(progfile)
                    else:
                        print("skip, because already computed")
                    i = i+1

df = pd.DataFrame({'Dataset':df_col_dataset, 'StartTree':df_col_starttree, 'Algorithm':df_col_algo, 'runtime':df_col_runtime, 'LQIC': df_col_lqic,'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
print(df.to_string())

#print(df.groupby(['Dataset', 'StartTree']).mean())
#print(df.groupby(['Dataset', 'StartTree']).var())

#df_stats = df.groupby(['Dataset', 'StartTree', 'Algorithm']).agg(['min', 'max', 'mean', 'var', 'std'])
df_stats = df.groupby(['Dataset', 'StartTree', 'Algorithm']).agg(['mean', 'var'])

print(df_stats)
df_stats.to_csv('df.csv')

#os.remove(progfile)
