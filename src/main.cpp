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

using namespace genesis;
using namespace genesis::tree;

std::mt19937 mt;
std::uniform_int_distribution<int> distribution_ab;
std::uniform_int_distribution<int> distribution_edges;

int main() {
    Logging::log_to_stdout ();
    Tree tree = DefaultTreeNewickReader().from_file(
        "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre");

    std::random_device rd;
    mt = std::mt19937(rd());
    distribution_ab = std::uniform_int_distribution<int>(0,1);
    distribution_edges = std::uniform_int_distribution<int>(0,tree.edge_count()-1);

    // Print the tree topology including node names and branch lengths.
    std::cout << PrinterCompact().print( tree );
    std::cout << PrinterTable().print( tree );

    QuartetScoreComputer<uint64_t> qsc(tree, "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", 1218, true, true); //TODO int type

    std::cout << "Sum lqic Reference Tree: " << sum_lqic_scores(qsc) << std::endl;

    for (int i = 0; i < 5; ++i) {
        tree = make_random_nni_moves(tree, 49, distribution_edges, distribution_ab, mt);
        qsc.recomputeScores(tree, false);
        std::cout << "Sum lqic random start tree: " << sum_lqic_scores(qsc) << std::endl;
        tree = tree_search(tree, qsc);
        std::cout << std::endl;
    }
}

