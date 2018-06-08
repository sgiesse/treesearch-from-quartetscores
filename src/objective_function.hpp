#ifndef OBJECTIVE_FUNCTION_HPP
#define OBJECTIVE_FUNCTION_HPP


enum ObjectiveFunction { LQIC, QPIC, EQPIC };

template<typename CINT>
double sum_lqic_scores(QuartetScoreComputer<CINT>& qsc) {
    std::vector<double> lqic = qsc.getLQICScores();
    double sum = 0;
    for (size_t j = 0; j < lqic.size(); ++j)
        if (lqic[j] <= 1 && lqic[j] >= -1) sum += lqic[j];
    return sum;
}

template<typename CINT>
double sum_qpic_scores(QuartetScoreComputer<CINT>& qsc) {
    std::vector<double> qpic = qsc.getQPICScores();
    double sum = 0;
    for (size_t j = 0; j < qpic.size(); ++j)
        if (qpic[j] <= 1 && qpic[j] >= -1) sum += qpic[j];
    return sum;
}

template<typename CINT>
double sum_eqpic_scores(QuartetScoreComputer<CINT>& qsc) {
    std::vector<double> eqpic = qsc.getEQPICScores();
    double sum = 0;
    for (size_t j = 0; j < eqpic.size(); ++j)
        if (eqpic[j] <= 1 && eqpic[j] >= -1) sum += eqpic[j];
    return sum;
}

template<typename CINT>
struct Functions {
    double (*obj_fun)(QuartetScoreComputer<CINT>&);
    void (*nni_a)(Tree&, size_t, QuartetScoreComputer<CINT>&);
    void (*nni_b)(Tree&, size_t, QuartetScoreComputer<CINT>&);
    void (*spr_score_update)(Tree&, size_t, size_t, QuartetScoreComputer<CINT>&);
    bool (*nni_restrict_edge)(Tree&, size_t, QuartetScoreComputer<CINT>&, bool);
    bool (*spr_restrict_edgepair)(Tree&, size_t, size_t, QuartetScoreComputer<CINT>&, bool);
    std::vector<double> (*getScores)(QuartetScoreComputer<CINT>&);
    void (*setScore)(QuartetScoreComputer<CINT>&, size_t, double);

    Functions(ObjectiveFunction objective);
};

template<typename CINT>
Functions<CINT>::Functions(ObjectiveFunction objective) {
    switch (objective) {
    case LQIC :
        obj_fun = sum_lqic_scores;
        nni_a = nni_a_with_lqic_update;
        nni_b = nni_b_with_lqic_update;
        spr_score_update = spr_lqic_update;
        nni_restrict_edge = [](Tree& tree, size_t e, QuartetScoreComputer<CINT>& qsc, bool restricted) {
            (void)tree;
            return restricted and qsc.getLQICScores()[e] > 0; };
        spr_restrict_edgepair =
            [](Tree& tree, size_t p, size_t r, QuartetScoreComputer<CINT>& qsc, bool restricted) {
            return restricted and !has_negative_lqic_on_spr_path(tree, p, r, qsc.getLQICScores()); };
        getScores = [](QuartetScoreComputer<CINT>& qsc) { return qsc.getLQICScores(); };
        setScore = [](QuartetScoreComputer<CINT>& qsc, size_t e, double val) { qsc.setLQIC(e, val); };

        break;
    case QPIC:
        obj_fun = sum_qpic_scores;
        nni_a = nni_a_with_qpic_update;
        nni_b = nni_b_with_qpic_update;
        spr_score_update = spr_qpic_update;
        TODO(Restrict Edges for QPIC)
        nni_restrict_edge = [](Tree& tree, size_t e, QuartetScoreComputer<CINT>& qsc, bool restricted) {
            (void)tree; (void)e; (void)qsc; (void)restricted;
            return false; };
        spr_restrict_edgepair = 
            [](Tree& tree, size_t p, size_t r, QuartetScoreComputer<CINT>& qsc, bool restricted) {
            (void)tree; (void)p; (void)r; (void)qsc; (void)restricted;
            return false; };
        getScores = [](QuartetScoreComputer<CINT>& qsc) { return qsc.getQPICScores(); };
        setScore = [](QuartetScoreComputer<CINT>& qsc, size_t e, double val) { qsc.setQPIC(e, val); };
        //throw std::runtime_error("Not implemented");
        break;
    case EQPIC:
        TODO(EQPIC)
        obj_fun = sum_eqpic_scores;
        nni_a = nni_a_with_eqpic_update;
        nni_b = nni_b_with_eqpic_update;
        spr_score_update = spr_eqpic_update;
        TODO(Restrict Edges for EQPIC)
            nni_restrict_edge = [](Tree& tree, size_t e, QuartetScoreComputer<CINT>& qsc, bool restricted) {
                                    (void)tree; (void)e; (void)qsc; (void)restricted;
                                    return false; };
        spr_restrict_edgepair = 
            [](Tree& tree, size_t p, size_t r, QuartetScoreComputer<CINT>& qsc, bool restricted) {
                (void)tree; (void)p; (void)r; (void)qsc; (void)restricted;
                return false; };
        getScores = [](QuartetScoreComputer<CINT>& qsc) { return qsc.getEQPICScores(); };
        setScore = [](QuartetScoreComputer<CINT>& qsc, size_t e, double val) { qsc.setEQPIC(e, val); };
        //throw std::runtime_error("Not implemented");
        break;
    }
}




#endif
