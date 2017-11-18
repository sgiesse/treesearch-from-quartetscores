#ifndef TREEOPERATIONS_HPP
#define TREEOPERATIONS_HPP

#include <algorithm>
#include <cmath>
#include <memory>
#include "genesis/tree/function/manipulation.hpp"

// --------- Forward Declarations
std::vector<Tree> nni(Tree& tree);
Tree nni_a(Tree& tree, int i);
Tree nni_b(Tree& tree, int i);
Tree make_random_nni_moves(Tree& tree, int n, std::uniform_int_distribution<int> distribution_edges, std::uniform_int_distribution<int> distribution_ab, std::mt19937 mt);
std::vector<std::string> leafNames(const std::string &evalTreesPath);
bool verify_leaf_ids_match(Tree tree1, Tree tree2, bool verbose);
void add_leaf(Tree& tree, size_t i, std::string lname);

// -----------------------------

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

void add_leaf(Tree& tree, size_t i, std::string lname) {
    add_new_node(tree, tree.edge_at(i).secondary_link().node());
    add_new_node(tree, tree.edge_at(i).secondary_link().node());
    tree.edge_at(i).secondary_link().next().outer().node().data_cast<DefaultNodeData>()->name = lname;

    size_t l2 = tree.edge_at(i).primary_link().index();
    size_t l3 = tree.edge_at(i).secondary_link().index();
    size_t l6 = tree.link_at(l3).next().index();
    size_t l8 = tree.link_at(l6).next().index();
    size_t l9 = tree.link_at(l8).outer().index();

    tree.link_at(l9).reset_next(&tree.link_at(l6));
    tree.link_at(l8).reset_next(&tree.link_at(l9));
    tree.link_at(l3).reset_next(&tree.link_at(l3));

    tree.link_at(l2).reset_outer(&tree.link_at(l9));
    tree.link_at(l9).reset_outer(&tree.link_at(l2));
    tree.link_at(l3).reset_outer(&tree.link_at(l8));
    tree.link_at(l8).reset_outer(&tree.link_at(l3));

    tree.link_at(l6).reset_node(&tree.link_at(l9).node());
    tree.link_at(l8).reset_node(&tree.link_at(l9).node());

    tree.link_at(l2).edge().reset_secondary_link(&tree.link_at(l9));
    tree.link_at(l8).edge().reset_secondary_link(&tree.link_at(l3));

    tree.link_at(l3).reset_edge(&tree.link_at(l8).edge());
    tree.link_at(l9).reset_edge(&tree.link_at(l2).edge());
}




std::vector<Tree> nni(Tree& tree) {
    std::vector<Tree> trees;
    for(size_t i = 0; i < tree.edge_count(); ++i ) {
        auto const& edge = tree.edge_at( i );
        if (!(edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
            continue; //edge is no internode

        Tree tnew = nni_a(tree, i);
        //std::cout << "valid tnew1:  " << validate_topology(tnew) << std::endl;
        trees.push_back(tnew);

        Tree tnew2 = nni_b(tree, i);
        //std::cout << "valid tnew2:  " << validate_topology(tnew2) << std::endl;
        trees.push_back(tnew2);
    }
    return trees;
}

Tree nni_a(Tree& tree, int i) {
    Tree tnew(tree);

    if (tnew.edge_at(i).primary_link().next().edge().secondary_link().index() ==
        tnew.edge_at(i).primary_link().next().index()) {

        tnew.edge_at( i ).primary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next());
        tnew.edge_at(i).secondary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next().next());
        auto l1 = &tnew.edge_at( i ).primary_link().next().next().outer();
        tnew.edge_at( i ).primary_link().next().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().reset_outer(l1);

        tnew.edge_at(i).primary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().outer());

        tnew.edge_at(i).primary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().next().edge());
        tnew.edge_at(i).secondary_link().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().edge());
    } else {

        tnew.edge_at( i ).primary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next());
        tnew.edge_at(i).secondary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next());
        auto l1 = &tnew.edge_at( i ).primary_link().next().outer();
        tnew.edge_at( i ).primary_link().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().reset_outer(l1);

        tnew.edge_at(i).primary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().outer());

        tnew.edge_at(i).primary_link().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().edge());
        tnew.edge_at(i).secondary_link().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().edge());
    }

    return tnew;
}

Tree nni_b(Tree& tree, int i) {
    Tree tnew(tree);
    if (tnew.edge_at(i).primary_link().next().edge().secondary_link().index() ==
        tnew.edge_at(i).primary_link().next().index()) {

        tnew.edge_at( i ).primary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next().next());
        auto l2 = &tnew.edge_at( i ).primary_link().next().next().outer();
        tnew.edge_at( i ).primary_link().next().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().next().reset_outer(l2);

        tnew.edge_at(i).primary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().next().outer());

        tnew.edge_at(i).primary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().next().edge());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().next().edge());
    } else {
        tnew.edge_at( i ).primary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next());
        auto l2 = &tnew.edge_at( i ).primary_link().next().outer();
        tnew.edge_at( i ).primary_link().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().next().reset_outer(l2);

        tnew.edge_at(i).primary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().next().outer());

        tnew.edge_at(i).primary_link().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().edge());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().next().edge());
    }
    return tnew;
}

Tree make_random_nni_moves(Tree& tree, int n, std::uniform_int_distribution<int> distribution_edges, std::uniform_int_distribution<int> distribution_ab, std::mt19937 mt) {
    Tree tnew = tree;
    for (int j = 0; j < n; ++j) {
        int i = distribution_edges(mt);
        while (!(tnew.edge_at(i).primary_link().node().is_inner() && tnew.edge_at(i).secondary_link().node().is_inner())) {
            i = distribution_edges(mt);
        }

        int ab = distribution_ab(mt);
        if (ab == 0)
            tnew = nni_a(tnew, i);
        else
            tnew = nni_b(tnew, i);
    }
    return tnew;
}

#endif
