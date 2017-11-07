#ifndef UTILS_HPP
#define UTILS_HPP

double sum_lqic_scores(QuartetScoreComputer<uint64_t>& qsc) {
    std::vector<double> lqic = qsc.getLQICScores();
    double sum = 0;
    for (uint64_t j = 0; j < lqic.size(); ++j)
        if (lqic[j] <= 1 && lqic[j] >= -1) sum += lqic[j];
    return sum;
}

#endif
