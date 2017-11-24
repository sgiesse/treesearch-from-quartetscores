#include "genesis/genesis.hpp"
#include <genesis/tree/printer/compact.hpp>
#include <genesis/tree/printer/table.hpp>
#include <genesis/tree/function/operators.hpp>
#include <genesis/utils/core/logging.hpp>

#include "QuartetScoreComputer.hpp"

#include <string>
#include <limits>

#include "treesearch.hpp"
#include "tree_operations.hpp"
#include "utils.hpp"

using namespace genesis;
using namespace genesis::tree;

template<typename CINT>
void doStuff(std::string pathToEvaluationTrees, int m) {
    //Tree sa_tree = stepwise_addition_tree<CINT>(pathToEvaluationTrees, mt, m);
    Tree sa_tree = random_tree(pathToEvaluationTrees);
    std::cout << PrinterCompact().print(sa_tree, print_help);
    QuartetScoreComputer<CINT> qsc =
        QuartetScoreComputer<CINT>(sa_tree, pathToEvaluationTrees, m, true, true);
    std::cout << "Sum lqic stepwise addition Tree: " << sum_lqic_scores(qsc) << std::endl;
    //Tree tree = tree_search_with_spr<CINT>(sa_tree, qsc);
    Tree tree = tree_search<CINT>(sa_tree, qsc);
}

int main2() {
    Logging::log_to_stdout ();
    std::string newick = "(((A,B),C),D,((E1,E2),(F,G)));";
    //std::string newick = "(((A,B),C),((D,E),(F,G)),(H,(I,((J,K),L))));";
    Tree tree = DefaultTreeNewickReader().from_string(newick);

    std::cout << PrinterCompact().print(tree, print_help);

    size_t edge1 = find_node( tree, "D" )->link().outer().next().next().edge().index();
    size_t edge2 = find_node( tree, "G" )->link().edge().index();
    //edge1 = 2;//20;
    //edge2 = 8;//4;
    std::cout << edge1 << " " << edge2 << std::endl;
    spr(tree, edge1, edge2);
    std::cout << "valid:  " << validate_topology(tree) << std::endl;
    std::cout << PrinterCompact().print(tree, print_help);
    return 0;
}

int main() {
    Logging::log_to_stdout ();

    std::string pathToEvaluationTrees =
        "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre";
    std::string pathToReferenceTree =
        "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre";

    Tree ref_tree = DefaultTreeNewickReader().from_file(pathToReferenceTree);

    std::cout << PrinterCompact().print( ref_tree );

    size_t m = countEvalTrees(pathToEvaluationTrees);
    if (m < (size_t(1) << 8)) doStuff<uint8_t>(pathToEvaluationTrees, m);
    else if (m < (size_t(1) << 16)) doStuff<uint16_t>(pathToEvaluationTrees, m);
    else if (m < (size_t(1) << 32)) doStuff<uint32_t>(pathToEvaluationTrees, m);
    else doStuff<uint64_t>(pathToEvaluationTrees, m);

    return 0;
}

