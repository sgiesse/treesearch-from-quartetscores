#ifndef NNI_HPP
#define NNI_HPP

#include <algorithm>
#include <cmath>
#include <memory>
#include "genesis/tree/function/manipulation.hpp"
#include "utils.hpp"
#include "../externals/generator/generator.hpp"

// --------- Forward Declarations
std::vector<Tree> nni(Tree& tree);
Tree nni_a(Tree& tree, int i);
Tree nni_b(Tree& tree, int i);
void nni_b_inplace(Tree& tree, int i);
void nni_a_inplace(Tree& tree, int i);
Tree make_random_nni_moves(Tree& tree, int n);
// -----------------------------


GENERATOR(nni_generator) {
    size_t i;
    Tree tree;
    nni_generator(Tree t) { tree = t;}
    EMIT(Tree)
        for (i = 0; i < tree.edge_count(); i++){
            if (!(tree.edge_at(i).primary_link().node().is_inner() && tree.edge_at(i).secondary_link().node().is_inner()))
                continue; //edge is no internode

            nni_a_inplace(tree, i);
            YIELD(tree);
            nni_a_inplace(tree, i); // Reverse previous NNI move
            nni_b_inplace(tree, i);
            YIELD(tree);
            nni_b_inplace(tree, i);
        }
    STOP;
};

template<typename CINT>
GENERATOR(nni_generator_qsc) {
    size_t i;
    Tree tree;
    QuartetScoreComputer<CINT>* qsc;
    nni_generator_qsc(Tree t, QuartetScoreComputer<CINT>* _qsc) { tree = t; qsc = _qsc; }
    EMIT(Tree)
        for (i = 0; i < tree.edge_count(); i++){
            if (!(tree.edge_at(i).primary_link().node().is_inner() && tree.edge_at(i).secondary_link().node().is_inner()))
                continue; //edge is no internode

            nni_a_with_lqic_update(tree, i, *qsc);
            YIELD(tree);
            nni_a_with_lqic_update(tree, i, *qsc); // Reverse previous NNI move
            nni_b_with_lqic_update(tree, i, *qsc);
            YIELD(tree);
            nni_b_with_lqic_update(tree, i, *qsc);
        }
    STOP;
};


std::vector<Tree> nni(Tree& tree) {
    std::vector<Tree> trees;
    for(size_t i = 0; i < tree.edge_count(); ++i ) {
        auto const& edge = tree.edge_at( i );
        if (!(edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
            continue; //edge is no internode

        Tree tnew = nni_a(tree, i);
        trees.push_back(tnew);

        Tree tnew2 = nni_b(tree, i);
        trees.push_back(tnew2);

        if (!validate_topology(tnew) or !validate_topology(tnew2)) {
            throw std::runtime_error("NNI produced invalid topology!");
        }
    }
    return trees;
}

std::vector<Tree> nni_only_negative_lqic(Tree& tree, const std::vector<double>& lqic) {
    std::vector<Tree> trees;
    for(size_t i = 0; i < tree.edge_count(); ++i ) {
        auto const& edge = tree.edge_at( i );
        if (!(edge.primary_link().node().is_inner() && edge.secondary_link().node().is_inner()))
            continue; //edge is no internode
        if (lqic[i] > 0)
            continue;

        Tree tnew = nni_a(tree, i);
        trees.push_back(tnew);

        Tree tnew2 = nni_b(tree, i);
        trees.push_back(tnew2);

        if (!validate_topology(tnew) or !validate_topology(tnew2)) {
            throw std::runtime_error("NNI produced invalid topology!");
        }
    }
    return trees;
}

Tree nni_a(Tree& tree, int i) {
    Tree tnew(tree);
    nni_a_inplace(tnew, i);
    return tnew;
}

template<typename CINT>
void swap_LQIC(size_t i, size_t j, QuartetScoreComputer<CINT>& qsc) {
    double t = qsc.getLQICScores()[i];
    qsc.setLQIC(i, qsc.getLQICScores()[j]);
    qsc.setLQIC(j, t);
}

template<typename CINT>
void nni_a_with_lqic_update(Tree& tree, size_t e, QuartetScoreComputer<CINT>& qsc) {
    bool case1 =
        tree.edge_at(e).primary_link().next().edge().secondary_link().index() ==
        tree.edge_at(e).primary_link().next().index();

    nni_a_inplace(tree, e);
    qsc.recomputeLqicForEdge(tree, e);
    if (case1)
        swap_LQIC<CINT>(tree.edge_at(e).primary_link().next().next().edge().index(),
                        tree.edge_at(e).secondary_link().next().edge().index(), qsc);
    else
        swap_LQIC<CINT>(tree.edge_at(e).primary_link().next().edge().index(),
        tree.edge_at(e).secondary_link().next().next().edge().index(), qsc);
}


template<typename CINT>
void nni_b_with_lqic_update(Tree& tree, size_t e, QuartetScoreComputer<CINT>& qsc) {
    bool case1 =
        tree.edge_at(e).primary_link().next().edge().secondary_link().index() ==
        tree.edge_at(e).primary_link().next().index();

    nni_b_inplace(tree, e);
    qsc.recomputeLqicForEdge(tree, e);
    if (case1)
        swap_LQIC<CINT>(tree.edge_at(e).primary_link().next().next().edge().index(),
                        tree.edge_at(e).secondary_link().next().next().edge().index(), qsc);
    else
        swap_LQIC<CINT>(tree.edge_at(e).primary_link().next().edge().index(),
        tree.edge_at(e).secondary_link().next().edge().index(), qsc);

}

void nni_a_inplace(Tree& tree, int i) {
    if (tree.edge_at(i).primary_link().next().edge().secondary_link().index() ==
        tree.edge_at(i).primary_link().next().index()) {

        tree.edge_at( i ).primary_link().next().next().outer().reset_outer(
            &tree.edge_at( i ).secondary_link().next());
        tree.edge_at(i).secondary_link().next().outer().reset_outer(
            &tree.edge_at( i ).primary_link().next().next());
        auto l1 = &tree.edge_at( i ).primary_link().next().next().outer();
        tree.edge_at( i ).primary_link().next().next().reset_outer(
            &tree.edge_at( i ).secondary_link().next().outer());
        tree.edge_at(i).secondary_link().next().reset_outer(l1);

        tree.edge_at(i).primary_link().next().next().edge().reset_secondary_link(
            &tree.edge_at(i).primary_link().next().next().outer());
        tree.edge_at(i).secondary_link().next().edge().reset_secondary_link(
            &tree.edge_at(i).secondary_link().next().outer());

        tree.edge_at(i).primary_link().next().next().outer().reset_edge(
            &tree.edge_at(i).primary_link().next().next().edge());
        tree.edge_at(i).secondary_link().next().outer().reset_edge(
            &tree.edge_at(i).secondary_link().next().edge());
    } else {

        tree.edge_at( i ).primary_link().next().outer().reset_outer(
            &tree.edge_at( i ).secondary_link().next().next());
        tree.edge_at(i).secondary_link().next().next().outer().reset_outer(
            &tree.edge_at( i ).primary_link().next());
        auto l1 = &tree.edge_at( i ).primary_link().next().outer();
        tree.edge_at( i ).primary_link().next().reset_outer(
            &tree.edge_at( i ).secondary_link().next().next().outer());
        tree.edge_at(i).secondary_link().next().next().reset_outer(l1);

        tree.edge_at(i).primary_link().next().edge().reset_secondary_link(
            &tree.edge_at(i).primary_link().next().outer());
        tree.edge_at(i).secondary_link().next().next().edge().reset_secondary_link(
            &tree.edge_at(i).secondary_link().next().next().outer());

        tree.edge_at(i).primary_link().next().outer().reset_edge(
            &tree.edge_at(i).primary_link().next().edge());
        tree.edge_at(i).secondary_link().next().next().outer().reset_edge(
            &tree.edge_at(i).secondary_link().next().next().edge());
    }
}

Tree nni_b(Tree& tree, int i) {
    Tree tnew(tree);
    nni_b_inplace(tnew, i);
    return tnew;
}

void nni_b_inplace(Tree& tree, int i) {
    if (tree.edge_at(i).primary_link().next().edge().secondary_link().index() ==
        tree.edge_at(i).primary_link().next().index()) {

        tree.edge_at( i ).primary_link().next().next().outer().reset_outer(
            &tree.edge_at( i ).secondary_link().next().next());
        tree.edge_at(i).secondary_link().next().next().outer().reset_outer(
            &tree.edge_at( i ).primary_link().next().next());
        auto l2 = &tree.edge_at( i ).primary_link().next().next().outer();
        tree.edge_at( i ).primary_link().next().next().reset_outer(
            &tree.edge_at( i ).secondary_link().next().next().outer());
        tree.edge_at(i).secondary_link().next().next().reset_outer(l2);

        tree.edge_at(i).primary_link().next().next().edge().reset_secondary_link(
            &tree.edge_at(i).primary_link().next().next().outer());
        tree.edge_at(i).secondary_link().next().next().edge().reset_secondary_link(
            &tree.edge_at(i).secondary_link().next().next().outer());

        tree.edge_at(i).primary_link().next().next().outer().reset_edge(
            &tree.edge_at(i).primary_link().next().next().edge());
        tree.edge_at(i).secondary_link().next().next().outer().reset_edge(
            &tree.edge_at(i).secondary_link().next().next().edge());
    } else {
        tree.edge_at( i ).primary_link().next().outer().reset_outer(
            &tree.edge_at( i ).secondary_link().next());
        tree.edge_at(i).secondary_link().next().outer().reset_outer(
            &tree.edge_at( i ).primary_link().next());
        auto l2 = &tree.edge_at( i ).primary_link().next().outer();
        tree.edge_at( i ).primary_link().next().reset_outer(
            &tree.edge_at( i ).secondary_link().next().outer());
        tree.edge_at(i).secondary_link().next().reset_outer(l2);

        tree.edge_at(i).primary_link().next().edge().reset_secondary_link(
            &tree.edge_at(i).primary_link().next().outer());
        tree.edge_at(i).secondary_link().next().edge().reset_secondary_link(
            &tree.edge_at(i).secondary_link().next().outer());

        tree.edge_at(i).primary_link().next().outer().reset_edge(
            &tree.edge_at(i).primary_link().next().edge());
        tree.edge_at(i).secondary_link().next().outer().reset_edge(
            &tree.edge_at(i).secondary_link().next().edge());
    }
}

Tree make_random_nni_moves(Tree& tree, int n) {
    Tree tnew = tree;
    for (int j = 0; j < n; ++j) {
        int i = Random::get_rand_int(0, tree.edge_count()-1);
        while (!(tnew.edge_at(i).primary_link().node().is_inner() && tnew.edge_at(i).secondary_link().node().is_inner())) {
            i = Random::get_rand_int(0, tree.edge_count()-1);
        }

        int ab = Random::get_rand_int(0, 1);
        if (ab == 0)
            tnew = nni_a(tnew, i);
        else
            tnew = nni_b(tnew, i);
    }
    return tnew;
}

#endif
