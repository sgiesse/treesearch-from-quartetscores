#ifndef STARTTREE_HPP
#define STARTTREE_HPP

#include "objective_function.hpp"

Tree random_tree_from_leaves(std::vector<std::string>& leaves) {
    // Define simple Tree structure
    struct SimpleNode {
        std::vector<SimpleNode> children;
        std::string name;
    };

    // Root, 3xchilds, 2x childs per internal node recursively
    SimpleNode root;
    root.children.push_back(SimpleNode());
    root.children.push_back(SimpleNode());
    root.children.push_back(SimpleNode());

    std::function<void(SimpleNode&, int, int)> buildTree;
    buildTree = [&](SimpleNode& node, int a, int b) {
        if (b-a == 1) {
            // leaf
            node.name = leaves[a];
        } else {
            node.children.push_back(SimpleNode());
            node.children.push_back(SimpleNode());
            buildTree(node.children[0],a,(b+a)/2);
            buildTree(node.children[1],(b+a)/2,b);
        }
    };

    int a = leaves.size()/3;
    int b = 2*a;
    buildTree(root.children[0], 0, a);
    buildTree(root.children[1], a, b);
    buildTree(root.children[2], b, leaves.size());

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
    LOG_DBG << newick << std::endl;

    // Genesis Tree from newick string
    Tree tree = DefaultTreeNewickReader().from_string(newick);
    return tree;
}

Tree random_tree(const std::string &evalTreesPath) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);

    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());
    return random_tree_from_leaves(leaves);
}


template<typename CINT>
Tree stepwise_addition_tree_from_leaves(QuartetScoreComputer<CINT>& qsc, std::vector<std::string>& leaves, size_t m, ObjectiveFunction objective) {
    Functions<CINT> functions = Functions<CINT>(objective);

    std::string newick = "(" + leaves[leaves.size()-1] + "," + leaves[leaves.size()-2] + "," + leaves[leaves.size()-3] + ");";
    leaves.pop_back(); leaves.pop_back(); leaves.pop_back();
    Tree tree = DefaultTreeNewickReader().from_string(newick);

    while (!leaves.empty()) {
        std::string lname = leaves.back();
        leaves.pop_back();
        LOG_DBG << "insert " << lname << std::endl;

        Tree best;
        double max = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < tree.edge_count(); ++i) {
            LOG_DBG << "leaves left: " << leaves.size() << " edge " << i << "/" << tree.edge_count() << std::endl;
            Tree tnew = Tree(tree);
            add_new_node(tnew, tnew.edge_at(i)).
                secondary_link().node().data_cast<DefaultNodeData>()->name = lname;

            LOG_DBG << "valid tnew:  " << validate_topology(tnew) << std::endl;

            qsc.recomputeScores(tnew, false);
            //double sum = sum_lqic_scores(qsc);
            double sum = functions.obj_fun(qsc);
            LOG_DBG << "sum scores " << sum << std::endl;
            if (sum > max) {
                best = tnew;
                max = sum;
            }
        }
        tree = best;
    }

    return tree;
}

template<typename CINT>
Tree stepwise_addition_tree(const std::string &evalTreesPath, size_t m, ObjectiveFunction objective) {
    std::vector<std::string> leaves = leafNames(evalTreesPath);
    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());
    return stepwise_addition_tree_from_leaves<CINT>-(evalTreesPath, leaves, m, objective);
}


uint64_t treecount(uint64_t n) {
    if (n == 3) return 1;
    return treecount(n-1) * (2*n-5);
}

template<typename CINT>
void _rec_exhaustive_search(Tree& tree, Tree& best, std::vector<std::string>& leaves, int li, QuartetScoreComputer<CINT>& qsc, double& max, int& c, ObjectiveFunction objective) {
    Functions<CINT> functions = Functions<CINT>(objective);
    if (li < 0) {
        qsc.recomputeScores(tree, false);
        //double sum = sum_lqic_scores(qsc);
        double sum = functions.obj_fun(qsc);
        if (sum > max) {
            max = sum;
            best = tree;
        } 
        c++;
        auto tc = treecount(leaves.size()+3);
        auto oneP = tc/100;
        if (c % oneP == 0) std::cout << c << "/" << treecount(leaves.size()+3) << "\n";
        return;
    }

    for (size_t i = 0; i < tree.edge_count(); ++i) {
        Tree tnew = Tree(tree);
        add_new_node(tnew, tnew.edge_at(i)).
            secondary_link().node().data_cast<DefaultNodeData>()->name = leaves[li];
        _rec_exhaustive_search(tnew, best, leaves, li-1, qsc, max, c, objective);
    }
}

template<typename CINT>
Tree exhaustive_search_from_leaves(const std::string &evalTreesPath, std::vector<std::string>& leaves, size_t m, ObjectiveFunction objective) {

    std::string newick = "(" + leaves[leaves.size()-1] + "," + leaves[leaves.size()-2] + "," + leaves[leaves.size()-3] + ");";
    leaves.pop_back(); leaves.pop_back(); leaves.pop_back();
    Tree tree = DefaultTreeNewickReader().from_string(newick);

    Tree precalc_tree(tree);
    for (int i = leaves.size()-1; i >= 0; --i) {
        LOG_DBG << leaves[i] << std::endl;
        add_new_node(precalc_tree,
                     precalc_tree.edge_at(Random::get_rand_int(0, precalc_tree.edge_count()-1))).
            secondary_link().node().data_cast<DefaultNodeData>()->name = leaves[i];
    }
    QuartetScoreComputer<CINT> qsc(precalc_tree, evalTreesPath, m, true, true);

    double max = -2.0 * leaves.size();
    Tree best = tree;
    int c = 0;

    for (size_t i = 0; i < tree.edge_count(); ++i) {
        Tree tnew = Tree(tree);
        add_new_node(tnew, tnew.edge_at(i)).
            secondary_link().node().data_cast<DefaultNodeData>()->name = leaves[leaves.size()-1];
        _rec_exhaustive_search(tnew, best, leaves, leaves.size()-2, qsc, max, c, objective);
    }

    return best;
}


template<typename CINT>
Tree exhaustive_search(const std::string &evalTreesPath, size_t m, ObjectiveFunction objective) {
    std::vector<std::string> leaves = leafNames(evalTreesPath);
    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());
    return exhaustive_search_from_leaves<CINT>(evalTreesPath, leaves, m, objective);
}




#endif
