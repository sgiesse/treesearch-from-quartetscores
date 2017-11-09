#ifndef TREEOPERATIONS_HPP
#define TREEOPERATIONS_HPP

#include <algorithm>
#include <cmath>

// --------- Forward Declarations
std::vector<Tree> nni(Tree& tree);
Tree nni_a(Tree& tree, int i);
Tree nni_b(Tree& tree, int i);
Tree make_random_nni_moves(Tree& tree, int n, std::uniform_int_distribution<int> distribution_edges, std::uniform_int_distribution<int> distribution_ab, std::mt19937 mt);

// -----------------------------

Tree random_tree(const std::string &evalTreesPath, std::mt19937 mt) {
    // Get set of node names

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

    std::vector<std::string> leaves(leaf_names.begin(), leaf_names.end());
    //std::random_shuffle(leaves.begin(), leaves.end());
    std::shuffle(leaves.begin(), leaves.end(), mt);

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

    // random NNI moves
    std::uniform_int_distribution<int> distribution_ab = std::uniform_int_distribution<int>(0,1);
    std::uniform_int_distribution<int> distribution_edges = std::uniform_int_distribution<int>(0,tree.edge_count()-1);
    return make_random_nni_moves(tree, 10, distribution_edges, distribution_ab, mt);
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
