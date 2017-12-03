#! /usr/env/python
import subprocess
import sys
import os
import dendropy

def compare_rf(tree_1_path, tree_2_path):
    # remove previous intermediary files
    if os.path.exists('RAxML_info.TEST'):
        os.remove('RAxML_info.TEST')
    if os.path.exists('RAxML_RF-Distances.TEST'):
        os.remove('RAxML_RF-Distances.TEST')
    if os.path.exists('trees'):
        os.remove('trees')

    # unroot the trees
    tree_1_unrooted = dendropy.Tree.get(
        path=tree_1_path,
        schema='newick',
        rooting='default-rooted')
    tree_1_unrooted.deroot()
    tree_2_unrooted = dendropy.Tree.get(
        path=tree_2_path,
        schema='newick',
        rooting='default-rooted')
    tree_2_unrooted.deroot()

    # merge the trees into one file
    outfile = open('trees', 'w')
    outfile.write(tree_1_unrooted.as_string(schema='newick'))
    outfile.write(tree_2_unrooted.as_string(schema='newick'))
    outfile.close()

    # compute RF distance using RAxML
    pathToRaxml = '../../raxml/standard-RAxML/'
    process = subprocess.Popen([pathToRaxml+'raxmlHPC-PTHREADS-AVX', '-m', 'GTRCAT', '-z', 'trees', '-f', 'r', '-n', 'TEST'], stdout = subprocess.PIPE, stderr=subprocess.STDOUT)
    out, err = process.communicate()
    if not os.path.exists('RAxML_RF-Distances.TEST'):
	    print("Something went wrong")
	    sys.exit()

    rf_line_splitted = open('RAxML_RF-Distances.TEST', 'r').read().split(' ')
    rf_plain = float(rf_line_splitted[2])
    rf_normalized = float(rf_line_splitted[3])
    return (rf_plain, rf_normalized)
