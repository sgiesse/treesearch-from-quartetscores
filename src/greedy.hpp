#ifndef GREEDY_HPP
#define GREEDY_HPP

#include "nni.hpp"
#include "spr.hpp"
#include "objective_function.hpp"

template<typename CINT>
Tree treesearch_nni(Tree& tree,
                    QuartetScoreComputer<CINT>& qsc,
                    ObjectiveFunction objective,
                    bool restricted) {
    Functions<CINT> functions = Functions<CINT>(objective);

    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = functions.obj_fun(qsc);

    Tree global_best = tnew;
    double global_max = oldscore;

    while (true) {
        double max = std::numeric_limits<double>::lowest();
        Tree best;
        qsc.recomputeScores(tnew, false);


        for (size_t i = 0; i < tnew.edge_count(); i++) {
            if (!(tnew.edge_at(i).primary_link().node().is_inner() && tnew.edge_at(i).secondary_link().node().is_inner()))
                continue; //edge is no internode
            if (functions.nni_restrict_edge(tnew, i, qsc, restricted)) continue;

            functions.nni_a(tnew, i, qsc);
            double sum = functions.obj_fun(qsc);
            if (sum > max) {
                max = sum;
                best = Tree(tnew);
            }
            functions.nni_a(tnew, i, qsc);

            functions.nni_b(tnew, i, qsc);
            sum = functions.obj_fun(qsc);
            if (sum > max) {
                max = sum;
                best = Tree(tnew);
            }
            functions.nni_b(tnew, i, qsc);
        }

        if (max > oldscore) {
            tnew = best;
            oldscore = max;
            LOG_INFO << "NNI best: " << max << std::endl;
            if (max > global_max) {
                global_max = max;
                global_best = tnew;
            }
        } else {
            break;
        }
    }
    qsc.recomputeScores(global_best, false);

    return global_best;
}


template<typename CINT>
Tree treesearch_combo(Tree& tree,
                      QuartetScoreComputer<CINT>& qsc,
                      ObjectiveFunction objective,
                      bool restricted) {
    Functions<CINT> functions = Functions<CINT>(objective);

    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);

    Tree best = tnew;
    double max = oldscore;

    while (true) {
        bool found_tree = false;
        tnew = best;
        qsc.recomputeScores(tnew, false);

        for (size_t i = 0; i < tnew.edge_count() and !found_tree; ++i) {
            for (size_t j = 0; j < tnew.edge_count() and !found_tree; ++j) {
                if (!validSprMove(tnew, i, j)) continue;
                if (functions.spr_restrict_edgepair(tnew, i, j, qsc, restricted)) continue;

                spr(tnew, i, j);
                functions.spr_score_update(tnew, i, j, qsc);

                double sum = functions.obj_fun(qsc);
                if (sum > max) {
                    max = sum;
                    best = tnew;
                    LOG_INFO << "best: " << max << std::endl;
                    found_tree = true;
                    break;
                }

                spr(tnew, i, j);
                functions.spr_score_update(tnew, i, j, qsc);
            }
        }
        if (found_tree) {
            tnew = treesearch_nni(best, qsc, objective, restricted);
            qsc.recomputeScores(tnew, false);
            double sum = functions.obj_fun(qsc);
            if (sum > max) {
                max = sum;
                best = tnew;
            }
        } else break;
    }

    return best;
}


#endif
