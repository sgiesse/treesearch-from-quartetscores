#ifndef UTILS_HPP
#define UTILS_HPP

#define DO_PRAGMA(x) _Pragma (#x)
#define TODO(x) DO_PRAGMA(message ("TODO - " #x))

#include "QuartetScoreComputer.hpp"
#include "nni.hpp"
#include "spr.hpp"

size_t countEvalTrees(const std::string &evalTreesPath);
std::string print_help(TreeNode const& node,TreeEdge const& edge);
std::string print_data(TreeNode const& node,TreeEdge const& edge);
void print_tree_with_lqic(Tree& tree, const std::vector<double>& lqic);

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

std::string print_data(TreeNode const& node,TreeEdge const& edge) {
    std::string out = node.data<DefaultNodeData>().name;
    out += " --- LQIC: ";
    out += std::to_string(edge.data<DefaultEdgeData>().branch_length);
    return out;
}

void print_tree_with_lqic(Tree& tree, const std::vector<double>& lqic) {
    for (size_t i = 0; i < tree.edge_count(); ++i) {
        tree.edge_at(i).data<DefaultEdgeData>().branch_length = lqic[i];
    }
    LOG_INFO << PrinterCompact().print(tree, print_data);
}


#endif
