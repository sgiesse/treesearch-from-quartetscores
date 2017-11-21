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


size_t countEvalTrees(const std::string &evalTreesPath) {
    size_t count = 0;
    utils::InputStream instream(utils::make_unique<utils::FileInputSource>(evalTreesPath));
    auto it = NewickInputIterator(instream);
    while (it) {
        count++;
        ++it;
    }
    return count;
}

std::string print_help(TreeNode const& node,TreeEdge const& edge) {
    std::string out = node.data<DefaultNodeData>().name;
    out += " ";
    out += std::to_string(node.index());
    out += " --- (";
    out += std::to_string(edge.index());
    out += ": ";
    out += std::to_string(edge.primary_link().index());
    out += ",";
    out += std::to_string(edge.secondary_link().index());
    out += ")";
    return out;
}

template<typename CINT>
void doStuff(std::string pathToEvaluationTrees, int m, std::mt19937 mt) {
    Tree sa_tree = stepwise_addition_tree<CINT>(pathToEvaluationTrees, mt, m);
    std::cout << PrinterCompact().print(sa_tree, print_help);
    QuartetScoreComputer<CINT> qsc =
        QuartetScoreComputer<CINT>(sa_tree, pathToEvaluationTrees, m, true, true);
    std::cout << "Sum lqic stepwise addition Tree: " << sum_lqic_scores(qsc) << std::endl;
    Tree tree = tree_search_with_spr<CINT>(sa_tree, qsc);
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

    std::random_device rd;
    mt = std::mt19937(rd());
    distribution_ab = std::uniform_int_distribution<int>(0,1);
    distribution_edges = std::uniform_int_distribution<int>(0,ref_tree.edge_count()-1);

    std::cout << PrinterCompact().print( ref_tree );

    size_t m = countEvalTrees(pathToEvaluationTrees);
    if (m < (size_t(1) << 8)) doStuff<uint8_t>(pathToEvaluationTrees, m, mt);
    else if (m < (size_t(1) << 16)) doStuff<uint16_t>(pathToEvaluationTrees, m, mt);
    else if (m < (size_t(1) << 32)) doStuff<uint32_t>(pathToEvaluationTrees, m, mt);
    else doStuff<uint64_t>(pathToEvaluationTrees, m, mt);

    return 0;
}

