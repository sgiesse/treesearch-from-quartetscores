#ifndef TREESEARCH_HPP
#define TREESEARCH_HPP

#include "utils.hpp"
#include "tree_operations.hpp"

template<typename CINT>
Tree tree_search(Tree& tree, QuartetScoreComputer<CINT>& qsc) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);
    int restarts = 10;

    Tree global_best = tnew;
    double global_max = oldscore;

    while (true) {
        double max = std::numeric_limits<double>::lowest();
        int best = 0;
        std::vector<Tree> nb_trees = nni(tnew);
        for (size_t i = 0; i < nb_trees.size(); ++i) {

            qsc.recomputeScores(nb_trees[i], false);
            double sum = sum_lqic_scores(qsc);
            if (sum > max) {
                max = sum;
                best = i;
            }
        }
        if (max > oldscore) {
            tnew = nb_trees[best];
            oldscore = max;
            std::cout << "best: " << max << std::endl;
            if (max > global_max) {
                global_max = max;
                global_best = tnew;
            }
        } else if (restarts > 0) {
            tnew = make_random_nni_moves(nb_trees[best], 10);
            restarts--;
        } else {
            break;
        }
    }

    qsc.recomputeScores(global_best, false);
    std::cout << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

    return global_best;
}

template<typename CINT>
Tree tree_search_with_spr(Tree& tree, QuartetScoreComputer<CINT>& qsc) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);

    Tree global_best = tnew;
    double global_max = oldscore;

    while (true) {
        double max = std::numeric_limits<double>::lowest();
        int best = 0;
        std::vector<Tree> nb_trees;
        for (size_t i = 0; i < tnew.edge_count(); ++i) {
            std::vector<bool> spr_ok(tnew.edge_count(), true);
            for (auto it : eulertour(tnew.edge_at(i).primary_link())) {
                if (it.edge().index() == i and it.link().index() == it.edge().secondary_link().index()) break;
                spr_ok[it.edge().index()] = false;
            }
            size_t j = i;
            do {
                spr_ok[j] = false;
                j = tnew.edge_at(j).primary_link().node().link().edge().index();
            } while (!tnew.edge_at(j).primary_link().node().is_root());
            spr_ok[j] = false;

            spr_ok[tnew.edge_at(i).primary_link().next().edge().index()] = false;
            spr_ok[tnew.edge_at(i).primary_link().next().next().edge().index()] = false;

            for(size_t j = 0; j < tnew.edge_count(); ++j) {
                if (!spr_ok[j]) continue;
                Tree t(tnew);
                if (spr(t, i, j))
                    nb_trees.push_back(t);
            }
        }
        std::cout << nb_trees.size() << " trees\n";
        for (size_t i = 0; i < nb_trees.size(); ++i) {
            if (!validate_topology(nb_trees[i])) continue;
            qsc.recomputeScores(nb_trees[i], false);
            double sum = sum_lqic_scores(qsc);
            if (sum > max) {
                max = sum;
                best = i;
            }
        }
        if (max > oldscore) {
            tnew = nb_trees[best];
            oldscore = max;
            std::cout << "best: " << max << std::endl;
            if (max > global_max) {
                global_max = max;
                global_best = tnew;
            }
        } else {
            break;
        }
    }

    qsc.recomputeScores(global_best, false);
    std::cout << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

    return global_best;
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
}

Tree random_tree(const std::string &evalTreesPath) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);

    std::shuffle(leaves.begin(), leaves.end(), Random::mt);
    return random_tree_from_leaves(leaves);
}

template<typename CINT>
Tree stepwise_addition_tree(const std::string &evalTreesPath, size_t m) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);
    std::shuffle(leaves.begin(), leaves.end(), Random::mt);

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
    QuartetScoreComputer<CINT> qsc(precalc_tree, evalTreesPath, m, true, true);

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



#endif
