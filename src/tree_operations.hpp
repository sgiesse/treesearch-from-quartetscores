#ifndef TREEOPERATIONS_HPP
#define TREEOPERATIONS_HPP

#include <algorithm>
#include <cmath>
#include <memory>
#include "genesis/tree/function/manipulation.hpp"
#include "utils.hpp"

// --------- Forward Declarations
std::vector<std::string> leafNames(const std::string &evalTreesPath);
bool verify_leaf_ids_match(Tree tree1, Tree tree2, bool verbose);

// -----------------------------

void reconnect_node_primary(Tree& tree, size_t edge, size_t new_primary_link) {
    size_t secondary_link = tree.edge_at(edge).secondary_link().index();
    tree.edge_at(edge).reset_primary_link(&tree.link_at(new_primary_link));
    tree.link_at(new_primary_link).reset_edge(&tree.edge_at(edge));
    tree.link_at(new_primary_link).reset_outer(&tree.link_at(secondary_link));
    tree.link_at(secondary_link).reset_outer(&tree.link_at(new_primary_link));
}

void reconnect_node_secondary(Tree& tree, size_t edge, size_t new_secondary_link) {
    size_t primary_link = tree.edge_at(edge).primary_link().index();
    tree.edge_at(edge).reset_secondary_link(&tree.link_at(new_secondary_link));
    tree.link_at(new_secondary_link).reset_edge(&tree.edge_at(edge));
    tree.link_at(new_secondary_link).reset_outer(&tree.link_at(primary_link));
    tree.link_at(primary_link).reset_outer(&tree.link_at(new_secondary_link));
}

void swap_edges(Tree& tree, size_t edge1, size_t edge2) {
    size_t e1pl = tree.edge_at(edge1).primary_link().index();
    size_t e1sl = tree.edge_at(edge1).secondary_link().index();
    size_t e2pl = tree.edge_at(edge2).primary_link().index();
    size_t e2sl = tree.edge_at(edge2).secondary_link().index();

    tree.edge_at(edge1).reset_primary_link(&tree.link_at(e2pl));
    tree.edge_at(edge1).reset_secondary_link(&tree.link_at(e2sl));
    tree.edge_at(edge2).reset_primary_link(&tree.link_at(e1pl));
    tree.edge_at(edge2).reset_secondary_link(&tree.link_at(e1sl));
    tree.link_at(e1pl).reset_edge(&tree.edge_at(edge2));
    tree.link_at(e1sl).reset_edge(&tree.edge_at(edge2));
    tree.link_at(e2pl).reset_edge(&tree.edge_at(edge1));
    tree.link_at(e2sl).reset_edge(&tree.edge_at(edge1));
}

std::vector<std::string> leafNames(const std::string &evalTreesPath) {
    std::set<std::string> leaf_names;
    utils::InputStream instream(utils::make_unique<utils::FileInputSource>(evalTreesPath));
    auto itTree = NewickInputIterator(instream, DefaultTreeNewickReader());
    while (itTree) {
        Tree const& tree = *itTree;
        for(auto const& node : tree.nodes()) {
            if (node->is_leaf()) {
                auto const& name = node->data<DefaultNodeData>().name;
                leaf_names.insert(name);
            }
        }
        ++itTree;
    }

    return std::vector<std::string>(leaf_names.begin(), leaf_names.end());
}

bool verify_leaf_ids_match(Tree tree1, Tree tree2, bool verbose = false) {
    std::map<std::string, size_t> leafToID1;
    std::map<std::string, size_t> leafToID2;

    for (size_t i = 0; i < tree1.node_count(); ++i) {
        if (tree1.node_at(i).is_leaf())
            leafToID1[tree1.node_at(i).data<DefaultNodeData>().name] = i;
    }
    for (size_t i = 0; i < tree2.node_count(); ++i) {
        if (tree2.node_at(i).is_leaf())
            leafToID2[tree2.node_at(i).data<DefaultNodeData>().name] = i;
    }

    bool ok = true;
    if (leafToID1.size() < leafToID2.size()) {
        for (auto it = leafToID1.begin(); it != leafToID1.end(); ++it) {
            if (leafToID2.find(it->first) == leafToID2.end() or
                leafToID1[it->first] != leafToID2[it->first]) {
                if (verbose) {
                    std::cout << "leafToID1[" << it->first << "] = " << leafToID1[it->first] << " != leafToID2[" << it->first << "] = "  << leafToID2[it->first] << std::endl;
                } else {
                    return false;
                }
            }
        }
    } else {
        for (auto it = leafToID2.begin(); it != leafToID2.end(); ++it) {
            if (leafToID1.find(it->first) == leafToID1.end() or
                leafToID1[it->first] != leafToID2[it->first]) {
                if (verbose) {
                    std::cout << "leafToID1[" << it->first << "] = " << leafToID1[it->first] << " != leafToID2[" << it->first << "] = "  << leafToID2[it->first] << std::endl;
                } else {
                    return false;
                }
            }
        }
    }

    return ok;
}

#endif
