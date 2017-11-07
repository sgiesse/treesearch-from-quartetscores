#ifndef TREEOPERATIONS_HPP
#define TREEOPERATIONS_HPP

// --------- Forward Declarations
std::vector<Tree> nni(Tree& tree);
Tree nni_a(Tree& tree, int i);
Tree nni_b(Tree& tree, int i);
Tree make_random_nni_moves(Tree& tree, int n, std::uniform_int_distribution<int> distribution_edges, std::uniform_int_distribution<int> distribution_ab, std::mt19937 mt);

// -----------------------------


std::vector<Tree> nni(Tree& tree) {
    std::vector<Tree> trees;
    for(size_t i = 0; i < tree.edge_count(); ++i ) {
        auto const& edge = tree.edge_at( i );
        if (!(edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
            continue; //edge is no internode

        Tree tnew = nni_a(tree, i);
        //std::cout << "valid tnew1:  " << validate_topology(tnew) << std::endl;
        trees.push_back(tnew);

        Tree tnew2 = nni_b(tree, i);
        //std::cout << "valid tnew2:  " << validate_topology(tnew2) << std::endl;
        trees.push_back(tnew2);
    }
    return trees;
}

Tree nni_a(Tree& tree, int i) {
    Tree tnew(tree);

    if (tnew.edge_at(i).primary_link().next().edge().secondary_link().index() ==
        tnew.edge_at(i).primary_link().next().index()) {

        tnew.edge_at( i ).primary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next());
        tnew.edge_at(i).secondary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next().next());
        auto l1 = &tnew.edge_at( i ).primary_link().next().next().outer();
        tnew.edge_at( i ).primary_link().next().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().reset_outer(l1);

        tnew.edge_at(i).primary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().outer());

        tnew.edge_at(i).primary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().next().edge());
        tnew.edge_at(i).secondary_link().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().edge());
    } else {

        tnew.edge_at( i ).primary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next());
        tnew.edge_at(i).secondary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next());
        auto l1 = &tnew.edge_at( i ).primary_link().next().outer();
        tnew.edge_at( i ).primary_link().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().reset_outer(l1);

        tnew.edge_at(i).primary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().outer());

        tnew.edge_at(i).primary_link().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().edge());
        tnew.edge_at(i).secondary_link().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().edge());
    }

    return tnew;
}

Tree nni_b(Tree& tree, int i) {
    Tree tnew(tree);
    if (tnew.edge_at(i).primary_link().next().edge().secondary_link().index() ==
        tnew.edge_at(i).primary_link().next().index()) {

        tnew.edge_at( i ).primary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next().next());
        auto l2 = &tnew.edge_at( i ).primary_link().next().next().outer();
        tnew.edge_at( i ).primary_link().next().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().next().reset_outer(l2);

        tnew.edge_at(i).primary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().next().outer());

        tnew.edge_at(i).primary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().next().edge());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().next().edge());
    } else {
        tnew.edge_at( i ).primary_link().next().outer().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_outer(
            &tnew.edge_at( i ).primary_link().next());
        auto l2 = &tnew.edge_at( i ).primary_link().next().outer();
        tnew.edge_at( i ).primary_link().next().reset_outer(
            &tnew.edge_at( i ).secondary_link().next().next().outer());
        tnew.edge_at(i).secondary_link().next().next().reset_outer(l2);

        tnew.edge_at(i).primary_link().next().edge().reset_secondary_link(
            &tnew.edge_at(i).primary_link().next().outer());
        tnew.edge_at(i).secondary_link().next().next().edge().reset_secondary_link(
            &tnew.edge_at(i).secondary_link().next().next().outer());

        tnew.edge_at(i).primary_link().next().outer().reset_edge(
            &tnew.edge_at(i).primary_link().next().edge());
        tnew.edge_at(i).secondary_link().next().next().outer().reset_edge(
            &tnew.edge_at(i).secondary_link().next().next().edge());
    }
    return tnew;
}

Tree make_random_nni_moves(Tree& tree, int n, std::uniform_int_distribution<int> distribution_edges, std::uniform_int_distribution<int> distribution_ab, std::mt19937 mt) {
    Tree tnew = tree;
    for (int j = 0; j < n; ++j) {
        int i = distribution_edges(mt);
        while (!(tnew.edge_at(i).primary_link().node().is_inner() && tnew.edge_at(i).secondary_link().node().is_inner())) {
            i = distribution_edges(mt);
        }

        int ab = distribution_ab(mt);
        if (ab == 0)
            tnew = nni_a(tnew, i);
        else
            tnew = nni_b(tnew, i);
    }
    return tnew;
}

#endif
