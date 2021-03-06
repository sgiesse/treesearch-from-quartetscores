#include <catch.hpp>

#include "genesis/genesis.hpp"
#include "genesis/tree/function/operators.hpp"
using namespace genesis;
using namespace genesis::tree;

#include "nni.hpp"
#include "spr.hpp"
#include "../externals/generator/generator.hpp"
#include "starttree.hpp"

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
    test_tree_manipulation("((A,B),C,D);", "((A,C),B,D);",
         [](Tree tree) { return nni_a(tree, tree.root_link().edge().index()); });

    test_tree_manipulation("(((A1,A2),B),C,D);", "(((A1,A2),C),B,D);",
                           [](Tree tree) { return nni_a(tree, tree.root_link().edge().index()); });

    test_tree_manipulation("(D,E,(C,(A,B)));", "(D,E,(A,(C,B)));",
                           [](Tree tree) {
                               return nni_a(tree, find_node(tree, "A")->link().edge().primary_link().node().link().edge().index()); });

}

TEST_CASE("nni_b") {
    test_tree_manipulation("((A,B),C,D);", "((C,B),A,D);",
                           [](Tree tree) { return nni_b(tree, tree.root_link().edge().index()); });

    test_tree_manipulation("(((A1,A2),B),C,D);", "((C,B),(A1,A2),D);",
                           [](Tree tree) { return nni_b(tree, tree.root_link().edge().index()); });

    test_tree_manipulation("(D,E,(C,(A,B)));", "(D,E,(B,(A,C)));",
                           [](Tree tree) {
                               return nni_b(tree, find_node(tree, "A")->link().edge().primary_link().node().link().edge().index()); });
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

        for (size_t i = 0; i < 10; ++i) {
            e = (e+i) % tree.edge_count();
            while (tree.edge_at(e).secondary_link().is_leaf()) ++e;
            std::cout << "edge: " << e << std::endl;
            nni_a_with_lqic_update<uint64_t>(tree, e, qsc);
            std::vector<double> lqic1 = qsc.getLQICScores();
            qsc.recomputeScores(tree, false);
            std::vector<double> lqic2 = qsc.getLQICScores();
            //REQUIRE(lqic1 == lqic2);
            bool eq = true;
            for (size_t j = 0; j < lqic1.size(); ++j) {
                if (Approx(lqic1[j]) != lqic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }

    SECTION("nni_b"){
        nni_b_with_lqic_update<uint64_t>(tree, e, qsc);
        std::vector<double> lqic1 = qsc.getLQICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> lqic2 = qsc.getLQICScores();
        REQUIRE(lqic1 == lqic2);

        for (size_t i = 0; i < 10; ++i) {
            e = (e+i) % tree.edge_count();
            while (tree.edge_at(e).secondary_link().is_leaf()) ++e;
            std::cout << "edge: " << e << std::endl;
            nni_b_with_lqic_update<uint64_t>(tree, e, qsc);
            std::vector<double> lqic1 = qsc.getLQICScores();
            qsc.recomputeScores(tree, false);
            std::vector<double> lqic2 = qsc.getLQICScores();
            //REQUIRE(lqic1 == lqic2);
            bool eq = true;
            for (size_t j = 0; j < lqic1.size(); ++j) {
                if (Approx(lqic1[j]) != lqic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }
}

TEST_CASE("QPIC after NNI") {
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);
    size_t e = 2;
    while (tree.edge_at(e).secondary_link().is_leaf()) ++e;

    SECTION("nni_a"){
        nni_a_with_qpic_update<uint64_t>(tree, e, qsc);
        std::vector<double> qpic1 = qsc.getQPICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> qpic2 = qsc.getQPICScores();
        //REQUIRE(qpic1 == qpic2);
        bool eq = true;
        for (size_t j = 0; j < qpic1.size(); ++j) {
            if (Approx(qpic1[j]) != qpic2[j]) { eq = false; continue; }
        }
        REQUIRE(eq);

        for (size_t i = 0; i < 10; ++i) {
            e = (e+i) % tree.edge_count();
            while (tree.edge_at(e).secondary_link().is_leaf()) ++e;
            std::cout << "edge: " << e << std::endl;
            nni_a_with_qpic_update<uint64_t>(tree, e, qsc);
            std::vector<double> qpic1 = qsc.getQPICScores();
            qsc.recomputeScores(tree, false);
            std::vector<double> qpic2 = qsc.getQPICScores();
            //REQUIRE(qpic1 == qpic2);
            bool eq = true;
            for (size_t j = 0; j < qpic1.size(); ++j) {
                if (Approx(qpic1[j]) != qpic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }

    SECTION("nni_b"){
        std::cout << "edge: " << e << std::endl;
        nni_b_with_qpic_update<uint64_t>(tree, e, qsc);
        std::vector<double> qpic1 = qsc.getQPICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> qpic2 = qsc.getQPICScores();
        //REQUIRE(qpic1 == qpic2);
        bool eq = true;
        for (size_t j = 0; j < qpic1.size(); ++j) {
            if (Approx(qpic1[j]) != qpic2[j]) { eq = false; continue; }
        }
        REQUIRE(eq);

        for (size_t i = 0; i < 10; ++i) {
            e = (e+i) % tree.edge_count();
            while (tree.edge_at(e).secondary_link().is_leaf()) ++e;
            std::cout << "edge: " << e << std::endl;
            nni_b_with_qpic_update<uint64_t>(tree, e, qsc);
            std::vector<double> qpic1 = qsc.getQPICScores();
            qsc.recomputeScores(tree, false);
            std::vector<double> qpic2 = qsc.getQPICScores();
            //REQUIRE(qpic1 == qpic2);
            bool eq = true;
            for (size_t j = 0; j < qpic1.size(); ++j) {
                if (Approx(qpic1[j]) != qpic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }
}

TEST_CASE("EQPIC after NNI") {
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);
    size_t e = 2;
    while (tree.edge_at(e).secondary_link().is_leaf()) ++e;

    SECTION("nni_a"){
        nni_a_with_eqpic_update<uint64_t>(tree, e, qsc);
        std::vector<double> eqpic1 = qsc.getEQPICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> eqpic2 = qsc.getEQPICScores();
        //REQUIRE(eqpic1 == eqpic2);
        bool eq = true;
        for (size_t j = 0; j < eqpic1.size(); ++j) {
            if (Approx(eqpic1[j]) != eqpic2[j]) { eq = false; continue; }
        }
        REQUIRE(eq);

        for (size_t i = 0; i < 10; ++i) {
            e = (e+i) % tree.edge_count();
            while (tree.edge_at(e).secondary_link().is_leaf()) ++e;
            std::cout << "edge: " << e << std::endl;
            nni_a_with_eqpic_update<uint64_t>(tree, e, qsc);
            std::vector<double> eqpic1 = qsc.getEQPICScores();
            qsc.recomputeScores(tree, false);
            std::vector<double> eqpic2 = qsc.getEQPICScores();
            //REQUIRE(eqpic1 == eqpic2);
            bool eq = true;
            for (size_t j = 0; j < eqpic1.size(); ++j) {
                if (Approx(eqpic1[j]) != eqpic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }

    SECTION("nni_b"){
        std::cout << "edge: " << e << std::endl;
        nni_b_with_eqpic_update<uint64_t>(tree, e, qsc);
        std::vector<double> eqpic1 = qsc.getEQPICScores();
        qsc.recomputeScores(tree, false);
        std::vector<double> eqpic2 = qsc.getEQPICScores();
        //REQUIRE(eqpic1 == eqpic2);
        bool eq = true;
        for (size_t j = 0; j < eqpic1.size(); ++j) {
            if (Approx(eqpic1[j]) != eqpic2[j]) { eq = false; continue; }
        }
        REQUIRE(eq);

        for (size_t i = 0; i < 10; ++i) {
            e = (e+i) % tree.edge_count();
            while (tree.edge_at(e).secondary_link().is_leaf()) ++e;
            std::cout << "edge: " << e << std::endl;
            nni_b_with_eqpic_update<uint64_t>(tree, e, qsc);
            std::vector<double> eqpic1 = qsc.getEQPICScores();
            qsc.recomputeScores(tree, false);
            std::vector<double> eqpic2 = qsc.getEQPICScores();
            //REQUIRE(eqpic1 == eqpic2);
            bool eq = true;
            for (size_t j = 0; j < eqpic1.size(); ++j) {
                if (Approx(eqpic1[j]) != eqpic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }
}


TEST_CASE("NNI Generator") {
    std::string newickIn = "(((A1,A2),B),C,D);";
    Tree tree = DefaultTreeNewickReader().from_string(newickIn);

    std::vector<Tree> nnis = nni(tree);
    size_t c = 0;
    nni_generator genNNI(tree);
    //genNNI.tree = tree;
    for (Tree t; genNNI(t);) {
        REQUIRE(validate_topology(t));

        auto node_comparator = [] (TreeNode const& node_l,TreeNode const& node_r) {return (node_r.is_leaf() and node_l.is_leaf()) or (node_r.data<DefaultNodeData>().name == node_l.data<DefaultNodeData>().name); };
        auto edge_comparator = [] (TreeEdge const& edge_l,TreeEdge const& edge_r) {(void) edge_l; (void) edge_r; return true;};
        REQUIRE(genesis::tree::equal(t, nnis[c], node_comparator, edge_comparator));
        c++;
    }
    REQUIRE(c == nnis.size());
}

TEST_CASE("NNI Generator with LQIC Updates") {
    omp_set_num_threads(1);
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);
    QuartetScoreComputer<uint64_t> qsc2 = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);
    //REQUIRE(qsc.getLQICScores() == qsc2.getLQICScores());

    nni_generator_qsc<uint64_t> genNNI(tree, &qsc);
    int c = 0;
    for (Tree t; genNNI(t);) {
        std::cout << "Tree #:" << c++ << std::endl;
        qsc2.recomputeScores(t, false);
        //REQUIRE(qsc.getLQICScores() == qsc2.getLQICScores());

        auto lqic1 = qsc.getLQICScores();
        auto lqic2 = qsc2.getLQICScores();
        bool eq = true;
        for (size_t j = 0; j < lqic1.size(); ++j) {
            if (Approx(lqic1[j]) != lqic2[j]) { eq = false; continue; }
        }
        REQUIRE(eq);
    }
}


TEST_CASE("SPR LQIC update") {
    omp_set_num_threads(1);
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);

    for (size_t i = 0; i < tree.edge_count(); ++i) {
        for (size_t j = 0; j < tree.edge_count(); ++j) {
            if (!validSprMove(tree, i, j)) continue;
            qsc.recomputeScores(tree, false);
            Tree t(tree);
            spr(t, i, j);
            spr_lqic_update(t, i, j, qsc);
            std::vector<double> lqic1 = qsc.getLQICScores();
            qsc.recomputeScores(t, false);
            std::vector<double> lqic2 = qsc.getLQICScores();
            //REQUIRE(lqic1 == lqic2);
            bool eq = true;
            for (size_t j = 0; j < lqic1.size(); ++j) {
                if (Approx(lqic1[j]) != lqic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }
}


TEST_CASE("SPR EQPIC update") {
    omp_set_num_threads(1);
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);

    for (size_t i = 0; i < tree.edge_count(); ++i) {
        std::cout << i << "/" << tree.edge_count() << std::endl;
        for (size_t j = 0; j < tree.edge_count(); ++j) {
            if (!validSprMove(tree, i, j)) continue;
            qsc.recomputeScores(tree, false);
            Tree t(tree);
            spr(t, i, j);
            spr_eqpic_update(t, i, j, qsc);
            std::vector<double> eqpic1 = qsc.getEQPICScores();
            qsc.recomputeScores(t, false);
            std::vector<double> eqpic2 = qsc.getEQPICScores();
            //REQUIRE(eqpic1 == eqpic2);
            bool eq = true;
            for (size_t j = 0; j < eqpic1.size(); ++j) {
                if (Approx(eqpic1[j]) != eqpic2[j]) { eq = false; continue; }
            }
            REQUIRE(eq);
        }
    }
}

TEST_CASE("SPR reverse") {
    omp_set_num_threads(1);
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);

    for (size_t i = 0; i < tree.edge_count(); ++i) {
        for (size_t j = 0; j < tree.edge_count(); ++j) {
            if (validSprMove(tree, i, j)) {
                qsc.recomputeScores(tree, false);
                Tree t(tree);
                auto lqic1 = qsc.getLQICScores();
                spr(t, i, j);
                spr_lqic_update(t, i, j, qsc);
                spr(t, i, j);
                spr_lqic_update(t, i, j, qsc);
                auto lqic2 = qsc.getLQICScores();

                auto node_comparator = [] (TreeNode const& node_l,TreeNode const& node_r) {return node_r.data<DefaultNodeData>().name == node_l.data<DefaultNodeData>().name; };
                auto edge_comparator = [] (TreeEdge const& edge_l,TreeEdge const& edge_r) {(void) edge_l; (void) edge_r; return true;};
                REQUIRE(validate_topology(t));

                REQUIRE(genesis::tree::equal(tree, t, node_comparator, edge_comparator));
                REQUIRE(lqic1 == lqic2);
            }
        }
    }
}

TEST_CASE("SPR Generator") {
    omp_set_num_threads(1);
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);

    spr_generator_qsc<uint64_t> genSPR(tree, &qsc, false);
    for (Tree t; genSPR(t);) {
        REQUIRE(validate_topology(t));
        std::vector<double> lqic1 = qsc.getLQICScores();
        qsc.recomputeScores(t, false);
        std::vector<double> lqic2 = qsc.getLQICScores();
        //REQUIRE(lqic1 == lqic2);
        bool eq = true;
        for (size_t j = 0; j < lqic1.size(); ++j) {
            if (Approx(lqic1[j]) != lqic2[j]) { eq = false; continue; }
        }
        REQUIRE(eq);
    }
}


TEST_CASE("SPR") {
    test_tree_manipulation("((A,B),((C,(D,E)),((F,(G,H)),I)),(J,K));", "((A,B),(C,(((F,(G,H)),(D,E)),I)),(J,K));",
                           [](Tree tree) {
                               size_t i = find_node(tree, "D")->link().edge().primary_link().node().link().edge().index();
                               size_t j = find_node(tree, "F")->link().edge().primary_link().node().link().edge().index();
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               std::cout << i << " " << j << std::endl;
                               spr(tree, i, j);
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               return tree; });

    test_tree_manipulation("((A,B),((C,(D,E)),((F,(G,H)),I)),(J,K));", "(B,((C,(D,E)),((F,(G,H)),I)),(J,(A,K)));",
                           [](Tree tree) {
                               size_t i = find_node(tree, "A")->link().edge().index();
                               size_t j = find_node(tree, "K")->link().edge().index();
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               std::cout << i << " " << j << std::endl;
                               spr(tree, i, j);
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               return tree; });

    test_tree_manipulation("((A,B),((C,(D,E)),((F,(G,H)),I)),(J,K));", "((A,B),(F,(G,H)),((C,(D,E)),((J,K),I)));",
                           [](Tree tree) {
                               size_t i = find_node(tree, "A")->link().edge().primary_link().node().link().edge().index();
                               size_t j = find_node(tree, "F")->link().edge().primary_link().node().link().edge().index();
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               std::cout << i << " " << j << std::endl;
                               spr(tree, i, j);
                               std::cout << validate_topology(tree) << std::endl;
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               return tree; });

    test_tree_manipulation("((A,B),((C,(D,E)),((F,(G,H)),I)),(J,K));", "(((C,(D,E)),((A,B),I)),(F,(G,H)),(J,K));",
                           [](Tree tree) {
                               size_t i = find_node(tree, "J")->link().edge().primary_link().node().link().edge().index();
                               size_t j = find_node(tree, "F")->link().edge().primary_link().node().link().edge().index();
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               std::cout << i << " " << j << std::endl;
                               spr(tree, i, j);
                               std::cout << validate_topology(tree) << std::endl;
                               std::cout << PrinterCompact().print(tree, print_help) << std::endl;
                               return tree; });
}


TEST_CASE("valid SPR move") {
    Tree tree = DefaultTreeNewickReader().from_string("((A,B),((C,(D,E)),((F,(G,H)),I)),(J,K));");
    REQUIRE(validSprMove(tree, 16, 6) == true);
    REQUIRE(validSprMove(tree, 16, 2) == true);
    REQUIRE(validSprMove(tree, 16, 11) == true);
    REQUIRE(validSprMove(tree, 12, 6) == true);
    REQUIRE(validSprMove(tree, 12, 3) == true);

    REQUIRE(validSprMove(tree, 16, 3) == false);
    REQUIRE(validSprMove(tree, 16, 0) == false);
    REQUIRE(validSprMove(tree, 12, 11) == false);
    REQUIRE(validSprMove(tree, 12, 15) == false);
    REQUIRE(validSprMove(tree, 12, 14) == false);
    REQUIRE(validSprMove(tree, 12, 12) == false);
}


TEST_CASE("recomputeScores") {
    Tree tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    /*for (size_t i = 0; i < tree.edge_count(); ++i) {
        for (size_t j = 0; j < tree.edge_count(); ++j) {
            if (validSprMove(tree, i, j)) {
                spr(tree, i, j);
            }
        }
        }*/
    size_t m = countEvalTrees("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc = QuartetScoreComputer<uint64_t>(tree, "../tests/data/yeast_all.tre", m, true, true);
    std::vector<double> lqic1 = qsc.getLQICScores();

    Tree rand_tree = random_tree("../tests/data/yeast_all.tre");
    QuartetScoreComputer<uint64_t> qsc2 =
        QuartetScoreComputer<uint64_t>(rand_tree, "../tests/data/yeast_all.tre", m, true, true);
    tree = DefaultTreeNewickReader().from_file("../tests/data/yeast_reference.tre");
    /*for (size_t i = 0; i < tree.edge_count(); ++i) {
        for (size_t j = 0; j < tree.edge_count(); ++j) {
            if (validSprMove(tree, i, j)) {
                spr(tree, i, j);
            }
        }
        }*/
    qsc2.recomputeScores(tree, false);
    std::vector<double> lqic2 = qsc2.getLQICScores();
    REQUIRE(lqic1==lqic2);
    /*bool eq = true;
    for (size_t j = 0; j < lqic1.size(); ++j) {
        if (Approx(lqic1[j]) != lqic2[j]) { eq = false; continue; }
    }
    REQUIRE(eq);*/
}
