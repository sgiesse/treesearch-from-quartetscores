#ifndef TREESEARCH_HPP
#define TREESEARCH_HPP

#include "utils.hpp"
#include "tree_operations.hpp"
#include "spr_iterator.hpp"

template<typename CINT>
Tree tree_search(Tree& tree, QuartetScoreComputer<CINT>& qsc) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);

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
            LOG_INFO << "best: " << max << std::endl;
            if (max > global_max) {
                global_max = max;
                global_best = tnew;
            }
        } else {
            break;
        }
    }

    qsc.recomputeScores(global_best, false);
    //LOG_INFO << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

    return global_best;
}

template<typename CINT>
Tree tree_search_combo(Tree& tree, QuartetScoreComputer<CINT>& qsc) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);

    Tree best = tnew;
    double max = oldscore;

    while (true) {
        SPRtree sprs(best);
        bool found_tree = false;
        LOG_INFO << "SPR";
        for (auto it : sprs) {
            Tree t(it.get());
            qsc.recomputeScores(t, false);
            double sum = sum_lqic_scores(qsc);
            if (sum > max) {
                max = sum;
                best = t;
                found_tree = true;
                break;
            }
        }
        if (found_tree) {
            best = tree_search(best, qsc);
            qsc.recomputeScores(best, false);
            max = sum_lqic_scores(qsc);
        } else break;
    }

    return best;
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
        LOG_DBG << nb_trees.size() << " trees\n";
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
            LOG_INFO << "best: " << max << std::endl;
            if (max > global_max) {
                global_max = max;
                global_best = tnew;
            }
        } else {
            break;
        }
    }

    qsc.recomputeScores(global_best, false);
    LOG_INFO << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

    return global_best;
}



Tree random_tree_from_leaves(std::vector<std::string> leaves) {
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
Tree stepwise_addition_tree(const std::string &evalTreesPath, size_t m) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);
    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());

    std::string newick = "(" + leaves[leaves.size()-1] + "," + leaves[leaves.size()-2] + "," + leaves[leaves.size()-3] + ");";
    leaves.pop_back(); leaves.pop_back(); leaves.pop_back();
    Tree tree = DefaultTreeNewickReader().from_string(newick);

    Tree precalc_tree(tree);
    for (int i = leaves.size()-1; i >= 0; --i) {
        LOG_DBG << leaves[i] << std::endl;
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
        LOG_DBG << "insert " << lname << std::endl;

        Tree best;
        double max = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < tree.edge_count(); ++i) {
            LOG_DBG << "leaves left: " << leaves.size() << " edge " << i << "/" << tree.edge_count() << std::endl;

            auto const& edge = tree.edge_at(i);
            if ((edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
                continue; //edge is internode TODO how?

            Tree tnew = Tree(tree);
            add_leaf(tnew, i, lname);

            LOG_DBG << "valid tnew:  " << validate_topology(tnew) << std::endl;

            qsc.recomputeScores(tnew, false);
            double sum = sum_lqic_scores(qsc);
            LOG_DBG << "sum lqic " << sum << std::endl;
            if (sum > max) {
                best = tnew;
                max = sum;
            }

            if (!verify_leaf_ids_match(tnew, precalc_tree)) {
                verify_leaf_ids_match(tnew, precalc_tree, true);
                throw std::runtime_error("leaf names don't match");
            }

        }
        LOG_DBG << "choose tree to continue with sum lqic " << max << std::endl;
        tree = best;
    }

    return tree;
}



#endif
