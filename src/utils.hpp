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


namespace {
    std::mt19937 mt;
    bool initialized = false;
}

namespace Random {

    void init() {
        std::random_device rd;
        mt = std::mt19937(rd());
        initialized = true;
    }

    void seed(int s) {
        mt.seed(s);
        initialized = true;
    }

    int get_rand_int(int a, int b) {
        if (!initialized) {
            init();
        }
        std::uniform_int_distribution<int> distr(a,b);
        return distr(mt);
    }

    std::mt19937 getMT() {
        if (!initialized) init();
        return mt;
    }
}

#endif
