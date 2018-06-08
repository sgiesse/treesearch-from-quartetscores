#include "genesis/genesis.hpp"
#include <genesis/tree/printer/compact.hpp>
#include <genesis/tree/printer/table.hpp>
#include <genesis/tree/function/operators.hpp>
#include <genesis/utils/core/logging.hpp>
#include <genesis/utils/core/options.hpp>

#include "QuartetScoreComputer.hpp"

#include "tclap/CmdLine.h"
#include <tclap/ValuesConstraint.h>

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
            //start_tree = tree_search<CINT>(start_tree, qsc, restrictByLqic);
            start_tree = treesearch_nni<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
        else if (treesearchAlgorithmClustered == "spr")
            //start_tree = tree_search_with_spr<CINT>(start_tree, qsc);
            throw std::runtime_error("Not implemented");
        else if (treesearchAlgorithmClustered == "combo")
            //start_tree = tree_search_combo<CINT>(start_tree, qsc, restrictByLqic);
            start_tree = treesearch_combo<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
        else if (treesearchAlgorithmClustered == "simann")
            //start_tree = simulated_annealing<CINT>(start_tree, qsc, false, simannfactor);
            start_tree = simulated_annealing<CINT>(start_tree, qsc, false, objectiveFunction, simannfactor);
        else if (treesearchAlgorithmClustered == "no")
            start_tree = start_tree;
        else  { LOG_ERR << treesearchAlgorithmClustered << " is unknown algorithm"; }

        end = std::chrono::steady_clock::now();
        res.timeFirstTreesearch =
            std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;

        qsc.recomputeScores(start_tree, false);
        LOG_INFO << "Sum lqic cluster Tree: " << sum_lqic_scores(qsc) << std::endl;

        start_tree = expanded_cluster_tree(start_tree, leafSets);

        qsc.recomputeScores(start_tree, false);
        LOG_INFO << "Sum lqic expanded final Tree: " << sum_lqic_scores(qsc) << std::endl;
    }

    begin = std::chrono::steady_clock::now();
    Tree final_tree;
    if (algorithm == "nni")
        //final_tree = tree_search<CINT>(start_tree, qsc, restrictByLqic);
        final_tree = treesearch_nni<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
    else if (algorithm == "spr")
        //final_tree = tree_search_with_spr<CINT>(start_tree, qsc);
        throw std::runtime_error("Not implemented");
    else if (algorithm == "combo")
        //final_tree = tree_search_combo<CINT>(start_tree, qsc, restrictByLqic);
        final_tree = treesearch_combo<CINT>(start_tree, qsc, objectiveFunction, restrictByLqic);
    else if (algorithm == "simann")
        //final_tree = simulated_annealing<CINT>(start_tree, qsc, clustering, simannfactor);
        final_tree = simulated_annealing<CINT>(start_tree, qsc, clustering, objectiveFunction, simannfactor);
    else if (algorithm == "no")
        final_tree = start_tree;
    else  { LOG_ERR << algorithm << " is unknown algorithm"; }

    end = std::chrono::steady_clock::now();
    res.timeFinalTreesearch =
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001;
    LOG_INFO << "Finished computing final tree. It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001 << " seconds." << std::endl;

    qsc.recomputeScores(final_tree, false);
    //LOG_INFO << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

    switch (objectiveFunction) {
    case LQIC:
        LOG_INFO << "Sum LQIC final Tree: " << sum_lqic_scores(qsc) << std::endl;
        break;
    case QPIC:
        LOG_INFO << "Sum QPIC final Tree: " << sum_qpic_scores(qsc) << std::endl;
        break;
    case EQPIC:
        LOG_INFO << "Sum EQPIC final Tree: " << sum_eqpic_scores(qsc) << std::endl;
        break;
    }

    LOG_INFO << "time Clustering: " << std::fixed << res.timeClustering << " seconds" << std::endl;
    LOG_INFO << "time CountingQuartets: " << std::fixed << res.timeCountingQuartets << " seconds" << std::endl;
    LOG_INFO << "time StartTree: " << std::fixed << res.timeStartTree << " seconds" << std::endl;
    LOG_INFO << "time FirstTreesearch: " << std::fixed << res.timeFirstTreesearch << " seconds" << std::endl;
    LOG_INFO << "time ExpandCluster: " << std::fixed << res.timeExpandCluster << " seconds" << std::endl;
    LOG_INFO << "time FinalTreesearch: " << std::fixed << res.timeFinalTreesearch << " seconds" << std::endl;
    LOG_INFO << "time Total: " << std::fixed << res.totalTime() << " seconds" << std::endl;

    DefaultTreeNewickWriter().to_file(final_tree, pathToOutput);
}

int main(int argc, char* argv[]) {
    Logging::log_to_stdout ();
    Options::get().allow_file_overwriting(true);

    LOG_BOLD << "Compute quartet score based Tree" << std::endl;

    std::string pathToEvaluationTrees;
    std::string startTreeMethod;
    std::string algorithm;
    std::string pathToOutput;
    std::string pathToStartTree;
    bool restrictByLqic;
    size_t numThreads;
    bool cached;
    size_t seed;
    float simannfactor;
    bool clustering;
    std::string treesearchAlgorithmClustered;
    ObjectiveFunction objectiveFunction;

    try {
        TCLAP::CmdLine cmd("Compute quartet score based Tree", ' ', "1.0");

        TCLAP::ValueArg<std::string> evalArg("e", "eval", "Path to the evaluation trees", false, "../../data/ICTC-master/data/Empirical/Yeast/yeast_all.tre", "string");
        cmd.add(evalArg);

        std::vector<std::string> allowedLogLevels = { "None","Error","Warning","Info","Progress","Debug","Debug1","Debug2","Debug3","Debug4" };
        TCLAP::ValuesConstraint<std::string> constraintLogLevels(allowedLogLevels);
        TCLAP::ValueArg<std::string> logLevelArg("l", "loglevel", "Log Level", false, "Info" , &constraintLogLevels);
        cmd.add(logLevelArg);

        std::vector<std::string> allowedStartTreeMethods = { "random", "stepwiseaddition", "exhaustive" };
        TCLAP::ValuesConstraint<std::string> constraintStart(allowedStartTreeMethods);
        TCLAP::ValueArg<std::string> startTreeMethodArg("s", "startTreeMethod", "Method to generate start tree", false, "stepwiseaddition", &constraintStart);
        cmd.add(startTreeMethodArg);

        std::vector<std::string> allowedAlgorithms = { "nni", "spr", "combo", "no", "simann" };
        TCLAP::ValuesConstraint<std::string> constraintAlgorithm(allowedAlgorithms);
        TCLAP::ValueArg<std::string> algorithmArg("a", "algorithm", "Algorithm to search tree", false, "nni", &constraintAlgorithm);
        cmd.add(algorithmArg);

        TCLAP::ValueArg<std::string> outArg("o", "outfile", "Path to output file", false, "../../out/out.tre", "string");
        cmd.add(outArg);

        TCLAP::ValueArg<std::string> startTreeArg("", "starttree", "Path to start tree file", false, "", "string");
        cmd.add(startTreeArg);

        TCLAP::SwitchArg restrictByLqicArg("x", "restricted", "Restrict NNI and SPR moves to edges with negative LQIC score");
        cmd.add(restrictByLqicArg);

        TCLAP::ValueArg<size_t> numThreadsArg("t", "numThreads", "Number of Threads", false, 1, "int");
        cmd.add(numThreadsArg);

        TCLAP::SwitchArg cachedArg("c", "cached", "Cache LQIC Scores.");
        cmd.add(cachedArg);

        TCLAP::ValueArg<size_t> seedArg("", "seed", "Seed", false, 1, "int");
        cmd.add(seedArg);

        TCLAP::ValueArg<float> simannfactorArg("", "simannfactor", "Factor for simulated_annealing. Choose value in (0.001, 0.01).", false, 0.005, "float");
        cmd.add(simannfactorArg);

        TCLAP::SwitchArg clusteringArg("", "clustering", "Cluster Taxa before Treesearch.");
        cmd.add(clusteringArg);

        std::vector<std::string> allowedTreesearchAlgorithmClustered = { "same", "nni", "spr", "combo", "no", "simann" }; 
        TCLAP::ValuesConstraint<std::string> constraintTreesearchAlgorithmClustered(allowedTreesearchAlgorithmClustered);
        TCLAP::ValueArg<std::string> treesearchAlgorithmClusteredArg("", "treesearchAlgorithmClustered", "", false, "same", &constraintTreesearchAlgorithmClustered);
        cmd.add(treesearchAlgorithmClusteredArg);

        std::vector<std::string> allowedObjectiveFunction = { "lqic", "qpic", "eqpic" };
        TCLAP::ValuesConstraint<std::string> constraintObjectiveFunction(allowedObjectiveFunction);
        TCLAP::ValueArg<std::string> objectiveFunctionArg("", "objectiveFunction", "The objective function to maximize.", false, "lqic", &constraintObjectiveFunction);
        cmd.add(objectiveFunctionArg);

        cmd.parse(argc, argv);

        pathToEvaluationTrees = evalArg.getValue();
        pathToOutput = outArg.getValue();
        pathToStartTree = startTreeArg.getValue();
        restrictByLqic = restrictByLqicArg.getValue();
        numThreads = numThreadsArg.getValue();
        cached = cachedArg.getValue();
        seed = seedArg.getValue();
        simannfactor = simannfactorArg.getValue();
        clustering = clusteringArg.getValue();
        treesearchAlgorithmClustered = treesearchAlgorithmClusteredArg.getValue();
        if (treesearchAlgorithmClustered == "same") treesearchAlgorithmClustered = algorithm;

        if (logLevelArg.getValue() == "None") Logging::max_level(utils::Logging::kNone);
        else if (logLevelArg.getValue() == "Error") Logging::max_level(utils::Logging::kError);
        else if (logLevelArg.getValue() == "Warning") Logging::max_level(utils::Logging::kWarning);
        else if (logLevelArg.getValue() == "Info") Logging::max_level(utils::Logging::kInfo);
        else if (logLevelArg.getValue() == "Progress") Logging::max_level(utils::Logging::kProgress);
        else if (logLevelArg.getValue() == "Debug") Logging::max_level(utils::Logging::kDebug);
        else if (logLevelArg.getValue() == "Debug1") Logging::max_level(utils::Logging::kDebug1);
        else if (logLevelArg.getValue() == "Debug2") Logging::max_level(utils::Logging::kDebug2);
        else if (logLevelArg.getValue() == "Debug3") Logging::max_level(utils::Logging::kDebug3);
        else if (logLevelArg.getValue() == "Debug4") Logging::max_level(utils::Logging::kDebug4);

        startTreeMethod = startTreeMethodArg.getValue();
        algorithm = algorithmArg.getValue();

        if (objectiveFunctionArg.getValue() == "lqic") objectiveFunction = LQIC;
        if (objectiveFunctionArg.getValue() == "qpic") objectiveFunction = QPIC;
        if (objectiveFunctionArg.getValue() == "eqpic") objectiveFunction = EQPIC;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }

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
