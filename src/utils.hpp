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


size_t countEvalTrees(const std::string &evalTreesPath) {
    size_t count = 0;
    utils::InputStream instream(utils::make_unique<utils::FileInputSource>(evalTreesPath));
    auto it = NewickInputIterator(instream);
    while (it) {
        count++;
        ++it;
    }
    return count;
}

std::string print_help(TreeNode const& node,TreeEdge const& edge) {
    std::string out = node.data<DefaultNodeData>().name;
    out += " ";
    out += std::to_string(node.index());
    out += " --- (";
    out += std::to_string(edge.index());
    out += ": ";
    out += std::to_string(edge.primary_link().index());
    out += ",";
    out += std::to_string(edge.secondary_link().index());
    out += ")";
    return out;
}

#endif
