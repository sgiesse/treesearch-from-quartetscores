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

#include "treesearch.hpp"
#include "tree_operations.hpp"
#include "utils.hpp"

using namespace genesis;
using namespace genesis::tree;

template<typename CINT>
void doStuff(std::string pathToEvaluationTrees, int m, std::string startTreeMethod, std::string algorithm, std::string pathToOutput, std::string pathToStartTree, bool restrictByLqic, bool cached) {

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    Tree start_tree;
    if (pathToStartTree == "") {
        if (startTreeMethod == "stepwiseaddition")
            start_tree = stepwise_addition_tree<CINT>(pathToEvaluationTrees, m);
        else if (startTreeMethod == "random")
            start_tree = random_tree(pathToEvaluationTrees);
        else if (startTreeMethod == "exhaustive")
            start_tree = exhaustive_search<CINT>(pathToEvaluationTrees, m);
        else { LOG_ERR << startTreeMethod << " is unknown start tree method"; }
    } else {
        LOG_INFO << "Read start tree from file";
        start_tree = DefaultTreeNewickReader().from_file(pathToStartTree);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //LOG_INFO << "Finished computing start tree. It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << " microseconds." << std::endl;
    LOG_INFO << "Finished computing start tree. It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001 << " seconds." << std::endl;

    LOG_INFO << PrinterCompact().print(start_tree);

    if (!validate_topology(start_tree)) {
        LOG_WARN << "Topology of start tree is not valid!";
    } else { LOG_INFO << "Topology of start tree is ok!"; }

    std::vector<std::string> st_leafNames;
    for (size_t i = 0; i < start_tree.node_count(); ++i) {
        if (start_tree.node_at(i).is_leaf())
            st_leafNames.push_back(start_tree.node_at(i).data<DefaultNodeData>().name);
    }
    std::sort(st_leafNames.begin(), st_leafNames.end());
    std::vector<std::string> leaves = leafNames(pathToEvaluationTrees);
    std::sort(leaves.begin(), leaves.end());
    if (leaves.size() != st_leafNames.size()) {
        throw std::runtime_error("size of set of leaf names is wrong");
    }
    for (size_t i = 0; i < leaves.size(); ++i) {
        if (leaves[i] != st_leafNames[i]){
            LOG_ERR << leaves[i] << " != " << st_leafNames[i] << std::endl;
            throw std::runtime_error("leaf name incorrect");
        }
    }

    QuartetScoreComputer<CINT> qsc =
        QuartetScoreComputer<CINT>(start_tree, pathToEvaluationTrees, m, true, true);
    LOG_INFO << "Sum lqic start Tree: " << sum_lqic_scores(qsc) << std::endl;

    if (cached) qsc.enableCache();
    else qsc.disableCache();

    begin = std::chrono::steady_clock::now();
    Tree final_tree;
    if (algorithm == "nni")
        final_tree = tree_search<CINT>(start_tree, qsc, restrictByLqic);
    else if (algorithm == "spr")
        final_tree = tree_search_with_spr<CINT>(start_tree, qsc);
    else if (algorithm == "combo")
        final_tree = tree_search_combo<CINT>(start_tree, qsc, restrictByLqic);
    else if (algorithm == "simann")
        final_tree = simulated_annealing<CINT>(start_tree, qsc);
    else if (algorithm == "no")
        final_tree = start_tree;
    else  { LOG_ERR << algorithm << " is unknown algorithm"; }

    end = std::chrono::steady_clock::now();
    LOG_INFO << "Finished computing final tree. It took: " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()*0.000001 << " seconds." << std::endl;

    qsc.recomputeScores(final_tree, false);
    LOG_INFO << "Sum lqic final Tree: " << sum_lqic_scores(qsc) << std::endl;

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

        cmd.parse(argc, argv);

        pathToEvaluationTrees = evalArg.getValue();
        pathToOutput = outArg.getValue();
        pathToStartTree = startTreeArg.getValue();
        restrictByLqic = restrictByLqicArg.getValue();
        numThreads = numThreadsArg.getValue();
        cached = cachedArg.getValue();
        seed = seedArg.getValue();

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
    } catch (TCLAP::ArgException &e) {
        std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }

    omp_set_num_threads(numThreads);
    Random::seed(seed);

    size_t m = countEvalTrees(pathToEvaluationTrees);
    if (m < (size_t(1) << 8))
        doStuff<uint8_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached);
    else if (m < (size_t(1) << 16))
        doStuff<uint16_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached);
    else if (m < (size_t(1) << 32))
        doStuff<uint32_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached);
    else
        doStuff<uint64_t>(pathToEvaluationTrees, m, startTreeMethod, algorithm, pathToOutput, pathToStartTree, restrictByLqic, cached);

    LOG_BOLD << "Done" << std::endl;

    return 0;
}
