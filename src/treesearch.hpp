#ifndef TREESEARCH_HPP
#define TREESEARCH_HPP

#include "utils.hpp"
#include "tree_operations.hpp"

Tree tree_search(Tree& tree, QuartetScoreComputer<uint64_t>& qsc) {
    Tree tnew = tree;
    qsc.recomputeScores(tnew, false);
    double oldscore = sum_lqic_scores(qsc);
    while (true) {
        double max = std::numeric_limits<double>::lowest();
        int best = 0;
        std::vector<Tree> nb_trees = nni(tnew);
        for (size_t i = 0; i < nb_trees.size(); ++i) {
            //std::cout << "NNI-tree #" << i << " /" << nb_trees.size() << std::endl;
            if (false) {
                //std::cout << PrinterCompact().print( nb_trees[i] );
                std::cout << PrinterTable().print( nb_trees[i] );
            }

            qsc.recomputeScores(nb_trees[i], false);
            double sum = sum_lqic_scores(qsc);
            //std::cout << sum << std::endl;
            if (sum > max) {
                max = sum;
                best = i;
            }
        }
        if (max > oldscore) {
            tnew = nb_trees[best];
            oldscore = max;
            std::cout << "best: " << max << std::endl;
        } else {
            break;
        }
    }

    qsc.recomputeScores(tnew, false);
    std::cout << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

    return tnew;
}

#endif
