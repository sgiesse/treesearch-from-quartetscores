#include "genesis/genesis.hpp"
#include <genesis/tree/printer/compact.hpp>
#include <genesis/tree/printer/table.hpp>
#include <genesis/tree/function/operators.hpp>
#include <genesis/utils/core/logging.hpp>

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
void doStuff(std::string pathToEvaluationTrees, int m) {
    Tree start_tree = stepwise_addition_tree<CINT>(pathToEvaluationTrees, m);
    //Tree start_tree = random_tree(pathToEvaluationTrees);
    LOG_INFO << PrinterCompact().print(start_tree, print_help);
    QuartetScoreComputer<CINT> qsc =
        QuartetScoreComputer<CINT>(start_tree, pathToEvaluationTrees, m, true, true);
    LOG_INFO << "Sum lqic stepwise addition Tree: " << sum_lqic_scores(qsc) << std::endl;
    //Tree tree = tree_search_with_spr<CINT>(start_tree, qsc);
    Tree tree = tree_search<CINT>(start_tree, qsc);
}

int main(int argc, char* argv[]) {
    Logging::log_to_stdout ();

    LOG_BOLD << "Compute quartet score based Tree" << std::endl;

    std::string pathToEvaluationTrees;
    std::string pathToReferenceTree;

    try {
        TCLAP::CmdLine cmd("Compute quartet score based Tree", ' ', "1.0");
        TCLAP::ValueArg<std::string> refArg("r", "ref", "Path to the reference tree", false, "../../data/ICTC-master/data/Empirical/Yeast/yeast_reference.tre", "string");
        cmd.add(refArg);
        TCLAP::ValueArg<std::string> evalArg("e", "eval", "Path to the evaluation trees", false, "../../data/ICTC-master/data/Empirical/Yeast/yeast_partial_only.tre", "string");
        cmd.add(evalArg);

        std::vector<std::string> allowedLogLevels = { "None","Error","Warning","Info","Progress","Debug","Debug1","Debug2","Debug3","Debug4" };
        TCLAP::ValuesConstraint<std::string> constraintLogLevels(allowedLogLevels);
        TCLAP::ValueArg<std::string> logLevelArg("l", "loglevel", "Log Level", false, "Info" , &constraintLogLevels);
        cmd.add(logLevelArg);

        cmd.parse(argc, argv);

        pathToReferenceTree = refArg.getValue();
        pathToEvaluationTrees = evalArg.getValue();

        if (logLevelArg.getValue() == "None") Logging::max_level(utils::Logging::kNone);
        if (logLevelArg.getValue() == "Error") Logging::max_level(utils::Logging::kError);
        if (logLevelArg.getValue() == "Warning") Logging::max_level(utils::Logging::kWarning);
        if (logLevelArg.getValue() == "Info") Logging::max_level(utils::Logging::kInfo);
        if (logLevelArg.getValue() == "Progress") Logging::max_level(utils::Logging::kProgress);
        if (logLevelArg.getValue() == "Debug") Logging::max_level(utils::Logging::kDebug);
        if (logLevelArg.getValue() == "Debug1") Logging::max_level(utils::Logging::kDebug1);
        if (logLevelArg.getValue() == "Debug2") Logging::max_level(utils::Logging::kDebug2);
        if (logLevelArg.getValue() == "Debug3") Logging::max_level(utils::Logging::kDebug3);
        if (logLevelArg.getValue() == "Debug4") Logging::max_level(utils::Logging::kDebug4);
    } catch (TCLAP::ArgException &e) {
        std::cerr << "ERROR: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }

    Tree ref_tree = DefaultTreeNewickReader().from_file(pathToReferenceTree);

    LOG_INFO << PrinterCompact().print(ref_tree);

    size_t m = countEvalTrees(pathToEvaluationTrees);
    if (m < (size_t(1) << 8)) doStuff<uint8_t>(pathToEvaluationTrees, m);
    else if (m < (size_t(1) << 16)) doStuff<uint16_t>(pathToEvaluationTrees, m);
    else if (m < (size_t(1) << 32)) doStuff<uint32_t>(pathToEvaluationTrees, m);
    else doStuff<uint64_t>(pathToEvaluationTrees, m);

    LOG_BOLD << "Done" << std::endl;

    return 0;
}
