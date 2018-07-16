#! /usr/env/python
import sys
import os
import json

# Provide a json-config file as first argument. Sample:
'''
{
    "data": [["../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre",
              "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre"],
             ["../../data/simulated_datasets_for_sarah/model102M6E.n01.estimated_gene.trees",
              "../../data/simulated_datasets_for_sarah/model102M6E.n01.true_species.tre"]],
    "exe": "../build/uquest",
    "algorithms": [["combo", ["custom", "-s", "random", "-a", "combo"]],["combox", ["custom", "-s", "random", "-a", "combo", "-x"]], ["comboxc", ["custom", "-s", "random", "-a", "combo", "-x", "-c"]], ["comboxcl", ["custom", "-s", "random", "-a", "combo", "-x", "--clustering"]], ["comboxccl", ["custom", "-s", "random", "-a", "combo", "-x", "--clustering", "-c"]]],
    "repeat": 1,
    "rename": {"runtime_total":"time","RF_normalized":"RF_norm"},
    "replace": {"model102M6E.n01.estimated_gene":"model102M6E.n01"},
    "columns":["time","LQIC","RF_norm"],
    "script":["#!/bin/bash"]
}
'''
with open(sys.argv[1]) as json_config_file:
    config = json.load(json_config_file)

def ds_name(d):
    #return d[0][d[0].rfind('/')+1:d[0].rfind('.')]
    return d[2]

def outdir():
    return "out/"+sys.argv[1][sys.argv[1].rfind("/")+1:sys.argv[1].find(".config.json")]+"/"

def out_info(d, algo, seed):
    return outdir()+ds_name(d)+"_"+algo[0]+"_"+str(seed)+".out"

def out_tree(d, algo, seed):
    return outdir()+ds_name(d)+"_"+algo[0]+"_"+str(seed)+".tree"

def make_treesearch(d, algo, seed):
    a = [config["exe"]] + algo[1] + ["-e", d[0], "-o", out_tree(d, algo, seed), "--seed", str(seed), "|", "tee", out_info(d, algo, seed)]
    return " ".join(a)

if __name__ == "__main__":
    with open(sys.argv[1][sys.argv[1].rfind("/")+1:sys.argv[1].find(".config.json")]+".sh", "w") as f:
        for s in config["script"]:
            f.write(s+"\n")
        #f.write("#!/bin/sh\n")
        f.write("mkdir out\n")
        f.write("mkdir " + outdir() + "\n")
        for d in config["data"]:
            for algo in config["algorithms"]:
                for seed in range(config["repeat"]):
                    a = make_treesearch(d, algo, seed)
                    print(a)
                    f.write(a+"\n")
