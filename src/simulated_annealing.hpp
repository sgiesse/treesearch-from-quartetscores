#ifndef SIMULATED_ANNEALING_HPP
#define SIMULATED_ANNEALING_HPP

#include "objective_function.hpp"
#include "nni.hpp"
#include "spr.hpp"


template<typename CINT>
void simulated_annealing_helper(Tree& tree, QuartetScoreComputer<CINT>& qsc, ObjectiveFunction objective) {
    Functions<CINT> functions = Functions<CINT>(objective);

    const float m = Random::get_rand_float(0,1);
    const int ab = Random::get_rand_int(0, 1);
    size_t e, p, r;
    if (m < 0.8) { //NNI
        e = Random::get_rand_int(0, tree.edge_count()-1);
        while (tree.edge_at(e).secondary_link().node().is_leaf())
            e = Random::get_rand_int(0, tree.edge_count()-1);
        if (ab == 0)
            functions.nni_a(tree, e, qsc);
        else
            functions.nni_b(tree, e, qsc);
    }
    else { //SPR
        p = Random::get_rand_int(0, tree.edge_count()-1);
        r = Random::get_rand_int(0, tree.edge_count()-1);
        while (!validSprMove(tree, p, r)) {
            p = Random::get_rand_int(0, tree.edge_count()-1);
            r = Random::get_rand_int(0, tree.edge_count()-1);
        }
        spr(tree, p, r);
        functions.spr_score_update(tree, p, r, qsc);
    }
}

template<typename CINT>
Tree simulated_annealing(Tree& tree, QuartetScoreComputer<CINT>& qsc, bool lowtemp, ObjectiveFunction objective, float factor = 0.005) {
    Functions<CINT> functions = Functions<CINT>(objective);
    Tree current(tree);
    const size_t M = tree.edge_count();
    const size_t Ntrial = 100;
    const size_t MAX_EPOCH_LENGTH = std::max((int)(factor*tree.edge_count()*tree.edge_count()), 10);

    double trial_sum_downhill = 0;
    size_t trial_count_downhill = 0;
    for (size_t i = 0; i < Ntrial or trial_count_downhill < 2; ++i) {
        double score_curr = functions.obj_fun(qsc);
        Tree candidate(current);
        simulated_annealing_helper(candidate, qsc, objective);
        double score = functions.obj_fun(qsc);
        //LOG_INFO << score << " " << score_curr << std::endl;
        if (score < score_curr) {
            trial_sum_downhill += (score - score_curr);
            trial_count_downhill++;
        }
        current = candidate;
    }
    current = Tree(tree);
    qsc.recomputeScores(current, false);

    const double P0 = lowtemp ? 0.002 : 0.2;
    const double T0 = (trial_sum_downhill/trial_count_downhill)/log(P0);
    const double TM = 0.001;
    const double alpha = pow(TM/T0, 1.0/(M-1));
    std::cout << T0 << " " << alpha << std::endl;
    double T = T0;

    size_t C = 0;
    const size_t MAX_NO_CHANGE = 2;
    const double P_ACCEPT = std::max(1.0/MAX_EPOCH_LENGTH, 0.02);

    Tree best = current;
    double max = functions.obj_fun(qsc);
    while (C < MAX_NO_CHANGE) {
        size_t accepted = 0;
        LOG_INFO << C << "/" << MAX_NO_CHANGE << " --  T:" << T << "  --  current: " <<  functions.obj_fun(qsc) << std::endl;
        for (size_t i = 0; i < MAX_EPOCH_LENGTH; ++i) {
            std::vector<double> scores = functions.getScores(qsc);
            double score_curr = functions.obj_fun(qsc);

            Tree candidate(current);
            simulated_annealing_helper(candidate, qsc, objective);
            //qsc.recomputeScores(candidate, false);
            double score = functions.obj_fun(qsc);
            double R = exp((score-score_curr)/T);
            //std::cout << R << " (" << score << ") ";

            if (R > 1) {
                current = candidate;
                accepted++;
                if (score > max) {
                    max = score;
                    best = Tree(candidate);
                }
            } else {
                if (Random::get_rand_float(0.0, 1.0) < R) {
                    current = candidate;
                    accepted++;
                }
                else {
                    for (size_t e = 0; e < current.edge_count(); ++e) {
                        functions.setScore(qsc, e, scores[e]);
                    }
                }
            }
        }
        if (accepted/(double)MAX_EPOCH_LENGTH < P_ACCEPT) C++;
        else C = 0;
        //std::cout << "  " << accepted << std::endl;

        T = alpha * T;
    }
    return best;
}





#endif
