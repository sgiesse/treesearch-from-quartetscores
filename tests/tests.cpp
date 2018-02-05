#include <catch.hpp>

#include "genesis/genesis.hpp"
#include "genesis/tree/function/operators.hpp"
using namespace genesis;
using namespace genesis::tree;

#include "tree_operations.hpp"

TEST_CASE("nni_a") {
    std::string newickIn = "((A,B),C,D);";
    Tree tree = DefaultTreeNewickReader().from_string(newickIn);
    std::string newickExpected = "((C,B),A,D);";
    Tree tree2 = DefaultTreeNewickReader().from_string(newickExpected);

    size_t idx = tree.root_link().edge().index();
    tree = nni_a(tree, idx);

    REQUIRE(validate_topology(tree));

    auto writer = DefaultTreeNewickWriter();
    writer.write_values(false);
    std::string newickOut = writer.to_string(tree);
    REQUIRE(newickOut == newickExpected);
}

TEST_CASE("nni_b") {
    std::string newickIn = "((A,B),C,D);";
    Tree tree = DefaultTreeNewickReader().from_string(newickIn);
    std::string newickExpected = "((A,C),B,D);";
    Tree tree2 = DefaultTreeNewickReader().from_string(newickExpected);

    size_t idx = tree.root_link().edge().index();
    tree = nni_b(tree, idx);

    REQUIRE(validate_topology(tree));

    auto writer = DefaultTreeNewickWriter();
    writer.write_values(false);
    std::string newickOut = writer.to_string(tree);
    REQUIRE(newickOut == newickExpected);
}
