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

Tree random_tree_from_leaves(std::vector<std::string> leaves) {
    // Define simple Tree structure
    struct SimpleNode {
        std::vector<SimpleNode> children;
        std::string name;
    };

    // Root, 3xchilds, 2x childs per internal node recursively
    int L = leaves.size(); // number of leaves
    int N = 2*L-2;
    int D = log2((N-1)/3)+1;
    SimpleNode root;
    root.children.push_back(SimpleNode());
    root.children.push_back(SimpleNode());
    root.children.push_back(SimpleNode());
    N -= 4;

    std::function<void(SimpleNode&, int)> buildTree;
    buildTree = [&](SimpleNode& node, int d) {
        for (size_t i = 0; i < node.children.size(); ++i) {
            // Decide if this is leaf or inner node
            bool leaf = L==N or d==D;
            if (leaf) {
                std::string _name = leaves.back();
                leaves.pop_back();
                node.children[i].name = _name;
                L--;
            } else {
                if (N >= 2) {
                    N-=2;
                    node.children[i].children.push_back(SimpleNode());
                    node.children[i].children.push_back(SimpleNode());
                    buildTree(node.children[i],d+1);
                }
                else {
                    L--;
                    node.name = leaves.back();
                    leaves.pop_back();
                }
            }
        }
    };

    buildTree(root, 1);

    // Print tree to newick format string
    std::string newick;
    std::function<void(SimpleNode&, std::string&)> rec_make_newick;
    rec_make_newick = [&](SimpleNode& node, std::string& str) {
        str += '(';
        for (size_t i = 0; i < node.children.size(); ++i) {
            if (node.children[i].children.size() == 0)
                str += node.children[i].name;
            else
                rec_make_newick(node.children[i], str);
            if (i < node.children.size()-1) str += ',';
        }
        str += node.name;
        str += ')';
    };
    rec_make_newick(root, newick);
    newick += ';';
    std::cout << newick << std::endl;

    // Genesis Tree from newick string
    Tree tree = DefaultTreeNewickReader().from_string(newick);
    return tree;
/*
    // random NNI moves
    std::uniform_int_distribution<int> distribution_ab = std::uniform_int_distribution<int>(0,1);
    std::uniform_int_distribution<int> distribution_edges = std::uniform_int_distribution<int>(0,tree.edge_count()-1);
    return make_random_nni_moves(tree, 10, distribution_edges, distribution_ab, mt);*/
}

Tree random_tree(const std::string &evalTreesPath, std::mt19937 mt) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);

    std::shuffle(leaves.begin(), leaves.end(), mt);
    return random_tree_from_leaves(leaves);
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


Tree stepwise_addition_tree(const std::string &evalTreesPath, std::mt19937 mt) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);
    std::shuffle(leaves.begin(), leaves.end(), mt);

    std::string newick = "(" + leaves[leaves.size()-1] + "," + leaves[leaves.size()-2] + "," + leaves[leaves.size()-3] + ");";
    leaves.pop_back(); leaves.pop_back(); leaves.pop_back();
    Tree tree = DefaultTreeNewickReader().from_string(newick);

    Tree precalc_tree(tree);
    for (int i = leaves.size()-1; i >= 0; --i) {
        std::cout << leaves[i] << std::endl;
        for (size_t j = 0; j < precalc_tree.edge_count(); ++j) {
            auto const& edge = precalc_tree.edge_at(j);
            if ((edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
                continue;
            else {
                add_leaf(precalc_tree, j, leaves[i]);
                break;
            }
        }
    }
    QuartetScoreComputer<uint64_t> qsc(precalc_tree, evalTreesPath, 1218, true, true);

    while (!leaves.empty()) {
        std::string lname = leaves.back();
        leaves.pop_back();
        std::cout << "insert " << lname << std::endl;

        Tree best;
        double max = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < tree.edge_count(); ++i) {
            std::cout << "leaves left: " << leaves.size() << " edge " << i << "/" << tree.edge_count() << std::endl;

            auto const& edge = tree.edge_at(i);
            if ((edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
                continue; //edge is internode TODO how?

            Tree tnew = Tree(tree);
            add_leaf(tnew, i, lname);

            std::cout << "valid tnew:  " << validate_topology(tnew) << std::endl;
            //std::cout << PrinterTable().print(tnew);
            //std::cout << PrinterCompact().print(tnew);

            qsc.recomputeScores(tnew, false);
            double sum = sum_lqic_scores(qsc);
            std::cout << "sum lqic " << sum << std::endl;
            if (sum > max) {
                best = tnew;
                max = sum;
            }

            if (!verify_leaf_ids_match(tnew, precalc_tree)) {
                verify_leaf_ids_match(tnew, precalc_tree, true);
                throw std::runtime_error("leaf names don't match");
            }

        }
        std::cout << "choose tree to continue with sum lqic " << max << std::endl;
        tree = best;
    }

    return tree;
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
