#include <catch.hpp>

#include "genesis/genesis.hpp"
#include "genesis/tree/function/operators.hpp"
using namespace genesis;
using namespace genesis::tree;

#include "nni.hpp"

void test_tree_manipulation(
    std::string newickIn, std::string newickExpected, std::function<Tree(Tree)> manipulateTree) {
    Tree tree = DefaultTreeNewickReader().from_string(newickIn);
    Tree treeExpected = DefaultTreeNewickReader().from_string(newickExpected);
    Tree treeOut = manipulateTree(tree);
    REQUIRE(validate_topology(treeOut));
    auto writer = DefaultTreeNewickWriter();
    writer.write_values(false);
    std::string newickOut = writer.to_string(treeOut);
    REQUIRE(newickOut == newickExpected);
}

TEST_CASE("nni_a") {
    test_tree_manipulation("((A,B),C,D);", "((C,B),A,D);",
         [](Tree tree) { return nni_a(tree, tree.root_link().edge().index()); });

    test_tree_manipulation("(((A1,A2),B),C,D);", "((C,B),(A1,A2),D);",
                           [](Tree tree) { return nni_a(tree, tree.root_link().edge().index()); });
}

TEST_CASE("nni_b") {
    test_tree_manipulation("((A,B),C,D);", "((A,C),B,D);",
                           [](Tree tree) { return nni_b(tree, tree.root_link().edge().index()); });
    test_tree_manipulation("(((A1,A2),(B1,B2)),(C1,C2),(D1,D2));", "(((A1,A2),(C1,C2)),(B1,B2),(D1,D2));",
                           [](Tree tree) { return nni_b(tree, tree.root_link().edge().index()); });
}

TEST_CASE("nni_a_reverse") {
    test_tree_manipulation("((A,B),C,D);", "((A,B),C,D);",
         [](Tree tree) {
             Tree t = nni_a(tree, tree.root_link().edge().index());
             return nni_a(t, tree.root_link().edge().index()); });
}

TEST_CASE("nni_b_reverse") {
    test_tree_manipulation("((A,B),C,D);", "((A,B),C,D);",
         [](Tree tree) {
             Tree t = nni_b(tree, tree.root_link().edge().index());
             return nni_b(t, tree.root_link().edge().index()); });
}

TEST_CASE("LQIC after NNI") {
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);
    size_t e = 2;
    while (tree.edge_at(e).secondary_link().is_leaf()) ++e;

    SECTION("nni_a"){
        nni_a_with_lqic_update<uint64_t>(tree, e, qsc);
        std::vector<double> lqic1 = qsc.getLQICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> lqic2 = qsc.getLQICScores();
        REQUIRE(lqic1 == lqic2);
    }

    SECTION("nni_b"){
        nni_b_with_lqic_update<uint64_t>(tree, e, qsc);
        std::vector<double> lqic1 = qsc.getLQICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> lqic2 = qsc.getLQICScores();
        REQUIRE(lqic1 == lqic2);
    }
}
