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

def parse_times(out):
    STR_START = "Finished computing start tree. It took: "
    STR_COUNT = "Finished counting quartets.\nIt took: "
    STR_SEARCH = "Finished computing final tree. It took: "
    outString = out.decode("UTF-8")
    starttreetime = countquartetstime = treesearchtime = 0
    str_pos = outString.find(STR_START, 0)
    if str_pos != -1:
        start = str_pos + len(STR_START)
        end = outString.find(" seconds",str_pos)
        starttreetime = float(outString[start:end])
    str_pos = outString.find(STR_COUNT, str_pos)
    if str_pos != -1:
        start = str_pos + len(STR_COUNT)
        end = outString.find(" seconds",str_pos)
        countquartetstime = float(outString[start:end])
    str_pos = outString.find(STR_SEARCH, str_pos)
    if str_pos != -1:
        start = str_pos + len(STR_SEARCH)
        end = outString.find(" seconds",str_pos)
        treesearchtime = float(outString[start:end])

    return (starttreetime, countquartetstime, treesearchtime)

def make_treesearch(d, starttreemethod, algo, seed):
    a = [config["exe"]] + config["args"] + ['-e', d[0], "--seed", str(seed)] + algo + starttreemethod
    print("----------------------------------------")
    print(" ".join(a))
    process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    lqic = parse_lqic(out)
    (starttreetime, countquartetstime, treesearchtime) = parse_times(out)
    (rf_plain, rf_normalized) = compare_rf('../../out/out.tre', d[1])
    return (starttreetime, countquartetstime, treesearchtime, lqic, rf_plain, rf_normalized)

df_col_runtime_start = []
df_col_runtime_count = []
df_col_runtime_search = []
df_col_lqic = []
df_col_dataset = []
df_col_starttree = []
df_col_algo = []
df_col_rf = []
df_col_rfnorm = []

"""
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
"""
i = 0
for d in config["data"]:
    for args_st in config["starttrees"]:
        for _i in range(config["repeatStartTreeMethod"]):
            for algo in config["algorithms"]:
                for _j in range(config["repeat"]):
                    seed = _i
                    if i == len(df_col_runtime_search):
                        (starttreetime, countquartetstime, treesearchtime, lqic, rf_plain, rf_normalized) = make_treesearch(d, args_st, algo, seed)
                        print("took " + str(starttreetime+countquartetstime+treesearchtime) + "s")
                        print("sum lqic: " + str(lqic))
                        print("plain RF distance: " + str(rf_plain))
                        print("normalized RF distance: " + str(rf_normalized))
                        df_col_dataset.append(d[0][d[0].rfind('/')+1:d[0].rfind('.')])
                        df_col_starttree.append(args_st[1])
                        df_col_algo.append(" ".join(algo[1:]))
                        df_col_runtime_start.append(starttreetime)
                        df_col_runtime_count.append(countquartetstime)
                        df_col_runtime_search.append(treesearchtime)
                        df_col_lqic.append(lqic)
                        df_col_rf.append(rf_plain)
                        df_col_rfnorm.append(rf_normalized)
                        df = pd.DataFrame({'Dataset':df_col_dataset, 'StartTree':df_col_starttree, 'Algorithm':df_col_algo, 'runtime_start':df_col_runtime_start, 'runtime_count':df_col_runtime_count, 'runtime_search':df_col_runtime_search, 'LQIC': df_col_lqic,'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
                        #df.to_csv(progfile)
                    else:
                        print("skip, because already computed")
                    i = i+1

df = pd.DataFrame({'Dataset':df_col_dataset, 'StartTree':df_col_starttree, 'Algorithm':df_col_algo, 'runtime_start':df_col_runtime_start, 'runtime_count':df_col_runtime_count, 'runtime_search':df_col_runtime_search, 'LQIC': df_col_lqic,'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
df['runtime_total'] = df['runtime_start']+df['runtime_count']+df['runtime_search']
print(df.to_string())

#print(df.groupby(['Dataset', 'StartTree']).mean())
#print(df.groupby(['Dataset', 'StartTree']).var())

#df_stats = df.groupby(['Dataset', 'StartTree', 'Algorithm']).agg(['min', 'max', 'mean', 'var', 'std'])
df_stats = df.groupby(['Dataset', 'StartTree', 'Algorithm']).agg(['mean', 'var'])

print(df_stats)
df_stats.to_csv('df.csv')

#os.remove(progfile)
