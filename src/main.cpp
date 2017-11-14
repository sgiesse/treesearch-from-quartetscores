#include "genesis/genesis.hpp"
#include <genesis/tree/printer/compact.hpp>
#include <genesis/tree/printer/table.hpp>
#include <genesis/tree/function/operators.hpp>
#include <genesis/utils/core/logging.hpp>

#include "QuartetScoreComputer.hpp"

#include <string>
#include <limits>
#include <random>

#include "treesearch.hpp"
#include "tree_operations.hpp"

using namespace genesis;
using namespace genesis::tree;

std::mt19937 mt;
std::uniform_int_distribution<int> distribution_ab;
std::uniform_int_distribution<int> distribution_edges;

int main() {
    Logging::log_to_stdout ();
    Tree ref_tree = DefaultTreeNewickReader().from_file(
        "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre");

    std::random_device rd;
    mt = std::mt19937(rd());
    distribution_ab = std::uniform_int_distribution<int>(0,1);
    distribution_edges = std::uniform_int_distribution<int>(0,ref_tree.edge_count()-1);

    // Print the tree topology including node names and branch lengths.
    std::cout << PrinterCompact().print( ref_tree );
    //std::cout << PrinterTable().print( tree );

    Tree sa_tree = stepwise_addition_tree(
        "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", mt);
    std::cout << PrinterCompact().print(sa_tree);

    Tree rand_tree = random_tree(
        "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", mt);
    std::cout << PrinterCompact().print( rand_tree );


    QuartetScoreComputer<uint64_t> qsc(ref_tree, "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", 1218, true, true); //TODO int type
    std::cout << "Sum lqic Reference Tree: " << sum_lqic_scores(qsc) << std::endl;

    qsc = QuartetScoreComputer<uint64_t>(rand_tree, "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", 1218, true, true); //TODO int type
    std::cout << "Sum lqic pure random Tree: " << sum_lqic_scores(qsc) << std::endl;

    qsc = QuartetScoreComputer<uint64_t>(sa_tree, "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", 1218, true, true); //TODO int type
    std::cout << "Sum lqic stepwise addition Tree: " << sum_lqic_scores(qsc) << std::endl;

    Tree tree = tree_search(sa_tree, qsc, mt);
/*
    Tree tree = rand_tree;
    for (int i = 0; i < 5; ++i) {
        tree = make_random_nni_moves(tree, 49, distribution_edges, distribution_ab, mt);
        //Tree tree = rand_tree;
        qsc.recomputeScores(tree, false);
        std::cout << "Sum lqic random start tree: " << sum_lqic_scores(qsc) << std::endl;
        tree = tree_search(tree, qsc, mt);
        std::cout << std::endl;
    }
*/
}

