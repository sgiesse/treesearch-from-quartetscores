#include "genesis/genesis.hpp"
#include <genesis/tree/printer/compact.hpp>
#include <genesis/tree/printer/table.hpp>
#include <genesis/tree/function/operators.hpp>
#include <genesis/utils/core/logging.hpp>
#include <genesis/utils/core/options.hpp>

#include "QuartetScoreComputer.hpp"

#include <string>
#include <limits>
#include <cstdio>

//#include "treesearch.hpp"
//#include "tree_operations.hpp"
#include "utils.hpp"
#include "random.hpp"
#include "objective_function.hpp"
#include "reduce_tree.hpp"
#include "greedy.hpp"
#include "simulated_annealing.hpp"
#include "starttree.hpp"

#include "../externals/cli11/CLI11.hpp"

using namespace genesis;
using namespace genesis::tree;

struct ResultsAndStats {
    double timeClustering;
    double timeCountingQuartets;
    double timeStartTree;
    double timeFirstTreesearch;
    double timeExpandCluster;
    double timeFinalTreesearch;

    Tree finalTree;

    ResultsAndStats() {
        timeClustering = timeCountingQuartets = timeStartTree = timeFirstTreesearch = timeExpandCluster = timeFinalTreesearch = 0;
    }

    double totalTime() {
        return timeClustering + timeCountingQuartets + timeStartTree + timeFirstTreesearch + timeExpandCluster + timeFinalTreesearch; }
};

template<typename CINT>
void doStuff(std::string pathToEvaluationTrees, int m, std::string startTreeMethod, std::string algorithm, std::string pathToOutput, std::string pathToStartTree, bool restrictByLqic, bool cached, float simannfactor, bool clustering, std::string treesearchAlgorithmClustered, ObjectiveFunction objectiveFunction) {

    ResultsAndStats res;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::vector<std::string> leaves;
    std::vector<std::vector<std::string> > leafSets;
    if (clustering) {
        begin = std::chrono::steady_clock::now();
        leafSets = leaf_sets(pathToEvaluationTrees);
        for (auto x : leafSets) leaves.push_back(x[0]);
        end = std::chrono::steady_clock::now();
        res.timeClustering =
            std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;
    } else {
        leaves = leafNames(pathToEvaluationTrees);
    }
    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());

    begin = std::chrono::steady_clock::now();
    Tree start_tree;
    if (pathToStartTree == "") {
        if (startTreeMethod == "stepwiseaddition")
            start_tree = stepwise_addition_tree_from_leaves<CINT>(pathToEvaluationTrees, leaves, m, objectiveFunction);
        else if (startTreeMethod == "random")
            start_tree = random_tree_from_leaves(leaves);
        else if (startTreeMethod == "exhaustive")
            start_tree = exhaustive_search_from_leaves<CINT>(pathToEvaluationTrees, leaves, m, objectiveFunction);
        else { LOG_ERR << startTreeMethod << " is unknown start tree method"; }
    } else {
        LOG_INFO << "Read start tree from file";
        start_tree = DefaultTreeNewickReader().from_file(pathToStartTree);
        if (start_tree.root_node().rank() == 1) {
            std::cout << start_tree.node_count() << " " << start_tree.edge_count() << std::endl;
            std::string newick = DefaultTreeNewickWriter().to_string(start_tree);
            std::cout << newick << std::endl;

            int c = 0; size_t a = 0; size_t b = newick.size(); size_t d = newick.size();
            for (size_t i = 0; i < newick.size() and b == newick.size(); ++i) {
                if (newick[i] == '(') c++;
                if (newick[i] == ')') c--;
                if (c == 2 and a == 0) a = i;
                if (c == 1 and a > 0) {
                    b = i;
                    d = b;
                    if (newick[i+1] == ':') {
                        d++;
                        while (newick[d] != ',' and newick[d] != ')') d++;
                    }
                }
            }
            std::string newick2 = newick.substr(0, a) + newick.substr(a+1, b-a-1) + newick.substr(d, newick.size() - d);
            start_tree = DefaultTreeNewickReader().from_string(newick2);
        }
    }

    end = std::chrono::steady_clock::now();
    end = std::chrono::steady_clock::now();
    res.timeStartTree =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;

    LOG_INFO << "Finished computing start tree. It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001 << " seconds." << std::endl;

    LOG_INFO << PrinterCompact().print(start_tree);

    if (!validate_topology(start_tree)) {
        LOG_WARN << "Topology of start tree is not valid!";
    } else { LOG_INFO << "Topology of start tree is ok!"; }

    begin = std::chrono::steady_clock::now();
    Tree rand_tree = random_tree(pathToEvaluationTrees);
    QuartetScoreComputer<CINT> qsc =
        QuartetScoreComputer<CINT>(rand_tree, pathToEvaluationTrees, m, true, true);

    qsc.recomputeScores(start_tree, false);
    switch (objectiveFunction) {
    case LQIC:
        LOG_INFO << "Sum LQIC start Tree: " << sum_lqic_scores(qsc) << std::endl;
        break;
    case QPIC:
        LOG_INFO << "Sum QPIC start Tree: " << sum_qpic_scores(qsc) << std::endl;
        break;
    case EQPIC:
        LOG_INFO << "Sum EQPIC start Tree: " << sum_eqpic_scores(qsc) << std::endl;
        break;
    }

    end = std::chrono::steady_clock::now();
    res.timeCountingQuartets =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;

    if (cached) qsc.enableCache();
    else qsc.disableCache();

    TODO(SPR algorithm in 1st and 2nd algorithm)
    if (clustering) {
        begin = std::chrono::steady_clock::now();

        if (treesearchAlgorithmClustered == "nni")
            start_tree = treesearch_nni<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
        else if (treesearchAlgorithmClustered == "spr")
            throw std::runtime_error("Not implemented");
        else if (treesearchAlgorithmClustered == "combo")
            start_tree = treesearch_combo<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
        else if (treesearchAlgorithmClustered == "simann")
            start_tree = simulated_annealing<CINT>(start_tree, qsc, false, objectiveFunction, simannfactor);
        else if (treesearchAlgorithmClustered == "no")
            start_tree = start_tree;
        else  { LOG_ERR << treesearchAlgorithmClustered << " is unknown algorithm"; }

        end = std::chrono::steady_clock::now();
        res.timeFirstTreesearch =
            std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;

        qsc.recomputeScores(start_tree, false);
        switch (objectiveFunction) {
        case LQIC:
            LOG_INFO << "Sum LQIC cluster Tree: " << sum_lqic_scores(qsc) << std::endl; break;
        case QPIC:
            LOG_INFO << "Sum QPIC cluster Tree: " << sum_qpic_scores(qsc) << std::endl; break;
        case EQPIC:
            LOG_INFO << "Sum EQPIC cluster Tree: " << sum_eqpic_scores(qsc) << std::endl; break;
        }

        start_tree = expanded_cluster_tree(start_tree, leafSets);

        qsc.recomputeScores(start_tree, false);
        switch (objectiveFunction) {
        case LQIC:
            LOG_INFO << "Sum LQIC expanded Tree: " << sum_lqic_scores(qsc) << std::endl; break;
        case QPIC:
            LOG_INFO << "Sum QPIC expanded Tree: " << sum_qpic_scores(qsc) << std::endl; break;
        case EQPIC:
            LOG_INFO << "Sum EQPIC expanded Tree: " << sum_eqpic_scores(qsc) << std::endl; break;
        }
    }

    begin = std::chrono::steady_clock::now();
    Tree final_tree;
    if (algorithm == "nni")
        final_tree = treesearch_nni<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
    else if (algorithm == "spr")
        throw std::runtime_error("Not implemented");
    else if (algorithm == "combo")
        final_tree = treesearch_combo<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
    else if (algorithm == "simann")
        final_tree = simulated_annealing<CINT>(start_tree, qsc, clustering, objectiveFunction, simannfactor);
    else if (algorithm == "no")
        final_tree = start_tree;
    else  { LOG_ERR << algorithm << " is unknown algorithm"; }

    end = std::chrono::steady_clock::now();
    res.timeFinalTreesearch =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;
    LOG_INFO << "Finished computing final tree. It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001 << " seconds." << std::endl;
    qsc.recomputeScores(final_tree, false);

    LOG_INFO << "--------------------------------------------------" << std::endl;
    //LOG_INFO << "Sum LQIC final Tree: " << sum_lqic_scores(qsc) << std::endl;
    //LOG_INFO << "Sum QPIC final Tree: " << sum_qpic_scores(qsc) << std::endl;
    //LOG_INFO << "Sum EQPIC final Tree: " << sum_eqpic_scores(qsc) << std::endl;
    LOG_INFO << "Sum LQIC final Tree: " << mean_lqic_scores(qsc) << std::endl;
    LOG_INFO << "Sum QPIC final Tree: " << mean_qpic_scores(qsc) << std::endl;
    LOG_INFO << "Sum EQPIC final Tree: " << mean_eqpic_scores(qsc) << std::endl;

    LOG_INFO << "Time Clustering: " << std::fixed << res.timeClustering << " seconds" << std::endl;
    LOG_INFO << "Time CountingQuartets: " << std::fixed << res.timeCountingQuartets << " seconds" << std::endl;
    LOG_INFO << "Time StartTree: " << std::fixed << res.timeStartTree << " seconds" << std::endl;
    LOG_INFO << "Time FirstTreesearch: " << std::fixed << res.timeFirstTreesearch << " seconds" << std::endl;
    LOG_INFO << "Time ExpandCluster: " << std::fixed << res.timeExpandCluster << " seconds" << std::endl;
    LOG_INFO << "Time FinalTreesearch: " << std::fixed << res.timeFinalTreesearch << " seconds" << std::endl;
    LOG_INFO << "Time Total: " << std::fixed << res.totalTime() << " seconds" << std::endl;

    DefaultTreeNewickWriter().to_file(final_tree, pathToOutput);
}

struct VectorValidator : public CLI::Validator {
    VectorValidator(std::vector<std::string> accepted) {
        std::stringstream out;
        out << "<";
        for (size_t i = 0; i < accepted.size()-1; ++i) out << accepted[i] << "|";
        out << accepted[accepted.size()-1] << ">";
        tname = out.str();

        func = [accepted](std::string input) {
                   bool ok = false;
                   for (std::string x : accepted) {
                       if (input.compare(x) == 0) ok = true;
                   }

                   if (!ok) return std::string("Unknown input");

                   return std::string();
               };
    }
};

int main(int argc, char* argv[]) {
    Logging::log_to_stdout ();
    Logging::details.level = false;
    Options::get().allow_file_overwriting(true);

    LOG_BOLD << "Compute quartet score based Tree" << std::endl;

    std::string pathToEvaluationTrees;
    std::string startTreeMethod;
    std::string algorithm;
    std::string pathToOutput;
    std::string pathToStartTree;
    bool restrictByLqic = false;
    size_t numThreads = 1;
    bool cached = false;
    size_t seed = 0;
    float simannfactor = 0.005;
    bool clustering = false;
    std::string treesearchAlgorithmClustered = "same";
    std::string objectiveFunctionStr = "lqic";
    ObjectiveFunction objectiveFunction;
    std::string loglevel = "Info";


    // --- Global Options
    CLI::App app{"uQuEST: uncertainty Quartet Estimated SuperTree"};
    app.require_subcommand(1);
    app.fallthrough(true);
    app.add_option("-e, --eval", pathToEvaluationTrees, "Path to the evaluation trees")->required()->check(CLI::ExistingFile);
    app.add_option("-o, --outfile", pathToOutput, "Path to output file")->required();
    app.add_option("--starttree", pathToStartTree, "Path to start tree file");
    app.add_option("-t, --numThreads", numThreads, "Number of Threads", true);
    app.add_option("--seed", seed, "Random seed", true);
    app.add_option("--objectiveFunction", objectiveFunctionStr, "The objective function to maximize.")->check(VectorValidator({ "lqic", "qpic", "eqpic" }));

    CLI::App* custom = app.add_subcommand("custom", "");
    custom->add_option("-s, --startTreeMethod", startTreeMethod, "Method to generate start tree")->required()->check(VectorValidator({"random", "stepwiseaddition", "exhaustive"}));
    custom->add_option("-a, --algorithm", algorithm, "Algorithm to search tree")->required()->check(VectorValidator({"nni", "simann", "spr", "combo", "no"}));
    custom->add_flag("-x, --restricted", restrictByLqic, "Restrict NNI and SPR moves to edges with negative LQIC score");
    custom->add_flag("-c, --cached", cached, "Cache Scores");
    custom->add_flag("--clustering", clustering, "Cluster Taxa before Treesearch.");
    custom->add_option("--factor", simannfactor, "Factor for simulated_annealing.", true)->check(CLI::Range(0.001, 0.01));
    custom->add_option("--treesearchAlgorithmClustered", treesearchAlgorithmClustered, "")->check(VectorValidator({"nni", "simann", "spr", "combo", "no", "same"}));
    custom->add_option("-l, --loglevel", loglevel, "Log Level")->check(VectorValidator({"None","Error","Warning","Info","Progress","Debug","Debug1","Debug2","Debug3","Debug4"}));

    CLI::App* ccsa = app.add_subcommand("ccsa", "Cached, clustered Simulated Annealing");
    ccsa->add_option("--factor", simannfactor, "Factor for simulated_annealing.", true)->check(CLI::Range(0.001, 0.01));

    CLI::App* ccnni = app.add_subcommand("ccnni", "Cached, Clustered Greedy with NNI moves");

    CLI::App* cccombo = app.add_subcommand("cccombo", "Cached, Clustered Combination of Hill-Climbing SPR moves and greedy NNI moves");

    CLI11_PARSE(app, argc, argv);
    if (app.got_subcommand(custom)) {

    } else if (app.got_subcommand(ccsa)) {
        startTreeMethod = "random";
        algorithm = "nni";
        treesearchAlgorithmClustered = "simann";
        clustering = true;
        cached = true;
    } else if (app.got_subcommand(ccnni)) {
        startTreeMethod = "stepwiseaddition";
        algorithm = "nni";
        treesearchAlgorithmClustered = "nni";
        clustering = true;
        cached = true;
        restrictByLqic = true;
    } else if (app.got_subcommand(cccombo)) {
        startTreeMethod = "stepwiseaddition";
        algorithm = "combo";
        treesearchAlgorithmClustered = "combo";
        clustering = true;
        cached = true;
        restrictByLqic = true;
    } else {
        throw std::runtime_error("Unknown subcommand");
    }

    if (loglevel == "None") Logging::max_level(utils::Logging::kNone);
    else if (loglevel == "Error") Logging::max_level(utils::Logging::kError);
    else if (loglevel == "Warning") Logging::max_level(utils::Logging::kWarning);
    else if (loglevel == "Info") Logging::max_level(utils::Logging::kInfo);
    else if (loglevel == "Progress") Logging::max_level(utils::Logging::kProgress);
    else if (loglevel == "Debug") Logging::max_level(utils::Logging::kDebug);
    else if (loglevel == "Debug1") Logging::max_level(utils::Logging::kDebug1);
    else if (loglevel == "Debug2") Logging::max_level(utils::Logging::kDebug2);
    else if (loglevel == "Debug3") Logging::max_level(utils::Logging::kDebug3);
    else if (loglevel == "Debug4") Logging::max_level(utils::Logging::kDebug4);

    if (treesearchAlgorithmClustered == "same") treesearchAlgorithmClustered = algorithm;

    if (objectiveFunctionStr == "lqic") objectiveFunction = LQIC;
    if (objectiveFunctionStr == "qpic") objectiveFunction = QPIC;
    if (objectiveFunctionStr == "eqpic") objectiveFunction = EQPIC;


    FILE *fp = fopen(pathToOutput.c_str(), "w");
    if (fp == NULL) {
        LOG_ERR << "Cannot write output file: " << pathToOutput << std::endl;
        return 0;
    }
    fclose(fp);

    omp_set_num_threads(numThreads);
    Random::seed(seed);

    size_t m = countEvalTrees(pathToEvaluationTrees);
    if (m < (size_t(1) << 8))
        doStuff<uint8_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached, simannfactor, clustering, treesearchAlgorithmClustered, objectiveFunction);
    else if (m < (size_t(1) << 16))
        doStuff<uint16_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached, simannfactor, clustering, treesearchAlgorithmClustered, objectiveFunction);
    else if (m < (size_t(1) << 32))
        doStuff<uint32_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached, simannfactor, clustering, treesearchAlgorithmClustered, objectiveFunction);
    else
        doStuff<uint64_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached, simannfactor, clustering, treesearchAlgorithmClustered, objectiveFunction);

    LOG_BOLD << "Done" << std::endl;

    return 0;
}
