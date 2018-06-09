import subprocess
import sys
import os
import timeit
import pandas as pd
import json
from compare_rf import compare_rf

def parse_score(out, scoreName, raiseError):
    outString = out.decode("UTF-8")
    searchStr = "Sum " + scoreName + " final Tree:"
    str_pos = outString.find(searchStr, 0)
    if str_pos != -1:
        start = str_pos + len(searchStr)
        end = outString.find("\n",str_pos)
        score = float(outString[start:end])
        return score
    else:
        if raiseError:
            print(out.decode("UTF-8"))
            raise RuntimeError("no " + scoreName + " found")
        else:
            return float('nan')

def lqic(treePath, evalTreesPath, raiseError):
    a = ["../build/treesearch", "-a", "no", "-e", evalTreesPath, "--starttree", treePath]
    process = subprocess.Popen(a, stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    lqic = parse_lqic(out, "LQIC", raiseError)
    return lqic

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
