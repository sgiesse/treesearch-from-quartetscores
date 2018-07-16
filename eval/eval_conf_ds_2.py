#! /usr/env/python
import sys
import os
import pandas as pd
import json
from compare_rf import compare_rf
from eval_utils import parse_score, parse_times
from eval_conf_ds_1 import ds_name, outdir, out_info, out_tree

# Provide a json-config file as first argument.
# Provide a directory with out and tree files as second argument.

with open(sys.argv[1]) as json_config_file:
    config = json.load(json_config_file)

def analyze_output(d, algo, seed):
    with open(out_info(d, algo, seed)) as f:
        out = f.read()
    lqic = parse_score(out, "LQIC", False)
    qpic = parse_score(out, "QPIC", False)
    eqpic= parse_score(out, "EQPIC", False)
    (t_count, t_start, t_search1, t_search2, t_total) = parse_times(out)
    (rf_plain, rf_normalized) = compare_rf(out_tree(d, algo, seed), d[1])
    return (t_count, t_start, t_search1, t_search2, t_total, lqic, qpic, eqpic, rf_plain, rf_normalized)


if __name__ == "__main__":
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

    for d in config["data"]:
        for algo in config["algorithms"]:
            for seed in range(config["repeat"]):
                (t_count, t_start, t_search1, t_search2, t_total, lqic, qpic, eqpic, rf_plain, rf_normalized) = analyze_output(d, algo, seed)
                print("took " + str(t_total) + "s")
                print("sum lqic: " + str(lqic))
                print("plain RF distance: " + str(rf_plain))
                print("normalized RF distance: " + str(rf_normalized))
                df_col_dataset.append(ds_name(d))
                #df_col_algo.append(" ".join(algo[1:]))
                df_col_algo.append(algo[0])
                df_col_runtime_start.append(t_start)
                df_col_runtime_count.append(t_count)
                df_col_runtime_search.append(t_search1+t_search2)
                df_col_lqic.append(lqic)
                df_col_qpic.append(qpic)
                df_col_eqpic.append(eqpic)
                df_col_rf.append(rf_plain)
                df_col_rfnorm.append(rf_normalized)

    df = pd.DataFrame({'Dataset':df_col_dataset, 'Algorithm':df_col_algo, 'runtime_start':df_col_runtime_start, 'runtime_count':df_col_runtime_count, 'runtime_search':df_col_runtime_search, 'LQIC': df_col_lqic, 'QPIC': df_col_qpic, 'EQPIC': df_col_eqpic, 'RF':df_col_rf, 'RF_normalized':df_col_rfnorm})
    df['runtime_total'] = df['runtime_start']+df['runtime_count']+df['runtime_search']
    print(df.to_string())
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
