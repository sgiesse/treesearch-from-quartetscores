#ifndef TREESEARCH_HPP
#define TREESEARCH_HPP

#include <cmath>

#include "genesis/tree/function/manipulation.hpp"

#include "utils.hpp"
#include "tree_operations.hpp"
#include "spr_iterator.hpp"

template<typename CINT>
Tree tree_search(Tree& tree, QuartetScoreComputer<CINT>& qsc, bool restrict_by_lqic = false) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);

    Tree global_best = tnew;
    double global_max = oldscore;

    while (true) {
        double max = std::numeric_limits<double>::lowest();
        int best = 0;

        std::vector<Tree> nb_trees;
        if (restrict_by_lqic) {
            qsc.recomputeScores(tnew, false);
            nb_trees = nni_only_negative_lqic(tnew, qsc.getLQICScores());
        }
        else
            nb_trees = nni(tnew);

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
Tree tree_search_combo(Tree& tree, QuartetScoreComputer<CINT>& qsc, bool restrict_by_lqic = false) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);

    Tree best = tnew;
    double max = oldscore;

    while (true) {
        SPRtree sprs(best);
        if (restrict_by_lqic) {
            qsc.recomputeScores(best, false);
            std::vector<double> lqic = qsc.getLQICScores();
            sprs.restrict_by_lqic(lqic);
        }
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
            best = tree_search(best, qsc, restrict_by_lqic);
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
        add_new_node(precalc_tree,
                     precalc_tree.edge_at(Random::get_rand_int(0, precalc_tree.edge_count()-1))).
            secondary_link().node().data_cast<DefaultNodeData>()->name = leaves[i];
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
            Tree tnew = Tree(tree);
            add_new_node(tnew, tnew.edge_at(i)).
                secondary_link().node().data_cast<DefaultNodeData>()->name = lname;

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


template<typename CINT>
void mcmc_helper(Tree& current, Tree& candidate, double T, double& score_min, double& score_max, double& score_curr, QuartetScoreComputer<CINT>& qsc) {
    const int m = Random::get_rand_int(0, 1);
    if (m == 0) { //NNI
        const int ab = Random::get_rand_int(0, 1);
        int e = Random::get_rand_int(0, current.edge_count()-1);
        while (current.edge_at(e).secondary_link().node().is_leaf())
            e = Random::get_rand_int(0, current.edge_count()-1);
        if (ab == 0)
            candidate = nni_a(current, e);
        else
            candidate = nni_b(current, e);
    } else { //SPR
        candidate = Tree(current);
        size_t p = Random::get_rand_int(0, current.edge_count()-1);
        size_t r = Random::get_rand_int(0, current.edge_count()-1);
        while (!spr(candidate, p, r)) {
            p = Random::get_rand_int(0, current.edge_count()-1);
            r = Random::get_rand_int(0, current.edge_count()-1);
        }
    }
    qsc.recomputeScores(candidate, false);
    double score = sum_lqic_scores(qsc);
    double score_interp;
    double score_curr_interp;

    if (score < score_min) score_min = score;
    if (score > score_max) score_max = score;

    if (score < 0) score_interp = (score_min - score) / score_min * 0.5;
    else score_interp = 0.5 + score / score_max * 0.5;

    if (score_curr < 0) score_curr_interp = (score_min - score_curr) / score_min * 0.5;
    else score_curr_interp = 0.5 + score_curr / score_max * 0.5;

    score_interp = pow(score_interp, 4);
    score_curr_interp = pow(score_curr_interp, 4);

    double R = score_interp / score_curr_interp;
    //std::cout << T << " " << R << " " << pow(R, 1/T) <<  "\n";
    //std::cout << score << " " << score_curr << " " << R << std::endl;
    if (R > 1) {
        current = candidate;
        score_curr = score;
    }
    else {
        if (Random::get_rand_float(0.0, 1.0) < pow(R, 1/T)) {
            current = candidate;
            score_curr = score;
        }
    }
}


template<typename CINT>
Tree mcmc(Tree& tree, QuartetScoreComputer<CINT>& qsc) {
    Tree current(tree);
    qsc.recomputeScores(current, false);
    double score_start = sum_lqic_scores(qsc);
    double score_min = score_start;
    double score_max = score_start;
    if (score_start < 0) score_max = 0;
    else score_min = 0;

    const int N = 2000;

    for (int i = 0; i < N; ++i) {
        qsc.recomputeScores(current, false);
        double score_curr = sum_lqic_scores(qsc);

        LOG_INFO << score_curr << std::endl;

        double T = pow(2.0, (N-i)/(N/10)-7);
        Tree candidate;
        mcmc_helper(current, candidate, T, score_min, score_max,  score_curr, qsc);
    }

    return current;
}



template<typename CINT>
Tree mcmcmc(Tree& tree, QuartetScoreComputer<CINT>& qsc) {
    Tree current(tree);
    qsc.recomputeScores(current, false);
    double score_start = sum_lqic_scores(qsc);
    double score_min = score_start;
    double score_max = score_start;
    if (score_start < 0) score_max = 0;
    else score_min = 0;

    const int N = 2000;
    const int M = 4;
    std::vector<Tree> currents;
    std::vector<double> Ts;
    std::vector<double> scores(M);

    for (int i = 0; i < M; ++i) {
        currents.push_back(Tree(current));
        //Ts.push_back(pow(2.0, i-(M/2)-2));
        Ts.push_back(pow(2.0, i-M+2));
        std::cout << pow(2.0, i-M+2) << "\n";

        qsc.recomputeScores(currents[i], false);
        scores[i] = sum_lqic_scores(qsc);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            qsc.recomputeScores(currents[j], false);
            scores[j] = sum_lqic_scores(qsc);

            if (j == 0) {LOG_INFO << scores[j] << "\n";}

            Tree candidate;
            mcmc_helper(currents[j], candidate, Ts[j], score_min, score_max,  scores[j], qsc);
        }

        size_t c1 = Random::get_rand_int(0, M-1);
        size_t c2 = Random::get_rand_int(0, M-1);
        while (c1 == c2) c2 = Random::get_rand_int(0, M-1);
        double p = (pow(scores[c1], 1/Ts[c1]) * pow(scores[c2], 1/Ts[c2]))/(pow(scores[c1], 1/Ts[c2]) * pow(scores[c2], 1/Ts[c1]));
        //if (Random::get_rand_float(0.0, 1.0) < p) {
        if (p>=1) { // TODO
            Tree tmp(currents[c1]);
            currents[c1] = currents[c2];
            currents[c2] = tmp;
            if (c1 == 0 or c2 == 0) std::cout << "swapped " << c1 << " <-> " << c2 << " " << p <<  "\n";
        }

        /*for (int j = 0; j < M; ++j) {
            std::cout << scores[j] << " ";
        } std::cout << std::endl;*/
    }

    return currents[0];
}



uint64_t treecount(uint64_t n) {
    if (n == 3) return 1;
    return treecount(n-1) * (2*n-5);
}


template<typename CINT>
void _rec_exhaustive_search(Tree& tree, Tree& best, std::vector<std::string>& leaves, int li, QuartetScoreComputer<CINT>& qsc, double& max, int& c) {
    if (li < 0) {
        qsc.recomputeScores(tree, false);
        double sum = sum_lqic_scores(qsc);
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
        _rec_exhaustive_search(tnew, best, leaves, li-1, qsc, max, c);
    }
}

template<typename CINT>
Tree exhaustive_search(const std::string &evalTreesPath, size_t m) {
    // Get set of node names
    std::vector<std::string> leaves = leafNames(evalTreesPath);
    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());

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
        _rec_exhaustive_search(tnew, best, leaves, leaves.size()-2, qsc, max, c);
    }

    return best;
}


#endif
