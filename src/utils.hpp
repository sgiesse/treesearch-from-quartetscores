#ifndef UTILS_HPP
#define UTILS_HPP

template<typename CINT>
double sum_lqic_scores(QuartetScoreComputer<CINT>& qsc) {
    std::vector<double> lqic = qsc.getLQICScores();
    double sum = 0;
    for (size_t j = 0; j < lqic.size(); ++j)
        if (lqic[j] <= 1 && lqic[j] >= -1) sum += lqic[j];
    return sum;
}

#endif
