#ifndef UQUEST_WRAPPER_HPP
#define UQUEST_WRAPPER_HPP

#include "genesis/genesis.hpp"
#include "QuartetScoreComputer.hpp"
#include "QuartetLookupTable.hpp"

#include "utils.hpp"
#include "random.hpp"
#include "objective_function.hpp"
#include "reduce_tree.hpp"
#include "greedy.hpp"
#include "simulated_annealing.hpp"
#include "starttree.hpp"

using namespace genesis;
using namespace genesis::tree;

Tree inferTree(const std::vector<std::string>& taxonLabels, QuartetLookupTable& table) {
    QuartetScoreComputer<CINT> qsc =
        QuartetScoreComputer<CINT>(rand_tree, table, true);
    Tree start_tree = random_tree_from_leaves(leaves);
    Tree final_tree = simulated_annealing<CINT>(start_tree, qsc, false, LQIC);
    return final_tree;
}



#endif
