#ifndef SPR_NNI
#define SPR_NNI

#include "tree_operations.hpp"
#include "../externals/generator/generator.hpp"

//-----------------------------------------------------
void spr(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx);
bool validSprMove(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx);
template<typename CINT> void spr_lqic_update(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx, QuartetScoreComputer<CINT>& qsc);
bool has_negative_lqic_on_spr_path(const Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx, const std::vector<double>& lqic);
//------------------------------------------------------


template<typename CINT>
GENERATOR(spr_generator_qsc) {
    size_t i;
    size_t j;
    std::vector<double> lqic;
    std::vector<size_t> invalidLQIC;
    Tree tree;
    QuartetScoreComputer<CINT>* qsc;
    bool restrict_by_lqic;
    spr_generator_qsc(Tree t, QuartetScoreComputer<CINT>* _qsc, bool _restrict_by_lqic) { tree = t; qsc = _qsc; restrict_by_lqic = _restrict_by_lqic; }

    EMIT(Tree)
        for (i = 0; i < tree.edge_count(); ++i) {
            for (j = 0; j < tree.edge_count(); ++j) {
                if (!validSprMove(tree, i, j)) continue;
                if (restrict_by_lqic and !has_negative_lqic_on_spr_path(tree, i, j, qsc->getLQICScores())) continue;

                spr(tree, i, j);
                spr_lqic_update(tree, i, j, *qsc);
                YIELD(tree);
                spr(tree, i, j);
                spr_lqic_update(tree, i, j, *qsc);
            }
        }
    STOP;
};

template<typename CINT>
void spr_lqic_update(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx, QuartetScoreComputer<CINT>& qsc) {
    std::vector<size_t> i1;
    std::vector<size_t> i2;
    std::vector<size_t> invalidLQIC;

    size_t pruneLinkIdx = tree.edge_at(pruneEdgeIdx).primary_link().index();
    size_t link_prune_no = tree.link_at(pruneLinkIdx).next().outer().index();
    size_t link_prune_nno = tree.link_at(pruneLinkIdx).next().next().outer().index();
    if (tree.link_at(pruneLinkIdx).node().is_root()) {
        i1.push_back(tree.link_at(link_prune_no).edge().index());
        i1.push_back(tree.link_at(link_prune_nno).edge().index());
    } else {
        size_t link_prune_sibling;
        if (link_prune_no == tree.link_at(pruneLinkIdx).node().link().edge().primary_link().index())
            link_prune_sibling = link_prune_nno;
        else link_prune_sibling = link_prune_no;
        size_t e = tree.link_at(link_prune_sibling).edge().index();
        while (!tree.edge_at(e).primary_link().node().is_root()) {
            i1.push_back(e);
            e = tree.edge_at(e).primary_link().node().link().edge().index();
        } i1.push_back(e);
    }

    size_t e = regraftEdgeIdx;
    while (!tree.edge_at(e).primary_link().node().is_root()) {
        i2.push_back(e);
        e = tree.edge_at(e).primary_link().node().link().edge().index();
    } i2.push_back(e);

    /* TODO
    if (!tree.link_at(pruneLinkIdx).node().is_root()) {
        //remove LCA -> Root from vector
        while (i1.size() > 0 and i2.size() > 0 and i1.back() == i2.back()) {
            i1.pop_back();
            i2.pop_back();
        }
    }*/

    invalidLQIC.reserve(i1.size() + i2.size() + 1);
    for (size_t i = 0; i < i1.size(); ++i) invalidLQIC.push_back(i1[i]);
    for (size_t i = 0; i < i2.size(); ++i) invalidLQIC.push_back(i2[i]);

    qsc.recomputeLqicForEdge(tree, invalidLQIC[0]);
    for (auto it = invalidLQIC.begin()+1; it != invalidLQIC.end(); ++it) qsc.recomputeLqicForEdge(*it);
}

bool validSprMove(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx) {
    if (pruneEdgeIdx >= tree.edge_count() or regraftEdgeIdx >= tree.edge_count()) return false;
    if (pruneEdgeIdx == regraftEdgeIdx) return false;
    if (regraftEdgeIdx == tree.edge_at(pruneEdgeIdx).primary_link().next().edge().index()) return false;
    if (regraftEdgeIdx == tree.edge_at(pruneEdgeIdx).primary_link().next().next().edge().index()) return false;

    for (auto it : eulertour(tree.edge_at(pruneEdgeIdx).primary_link())) {
        if (it.edge().index() == pruneEdgeIdx and it.link().index() == it.edge().secondary_link().index()) break;
        if (it.edge().index() == regraftEdgeIdx) return false;
    }

    return true;
}

bool has_negative_lqic_on_spr_path(const Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx, const std::vector<double>& lqic) {
    std::vector<size_t> i1;
    std::vector<size_t> i2;

    if (!tree.edge_at(pruneEdgeIdx).primary_link().node().is_root()) {
        size_t e = tree.edge_at(pruneEdgeIdx).primary_link().node().link().edge().index();
        while (!tree.edge_at(e).primary_link().node().is_root()) {
            i1.push_back(e);
            e = tree.edge_at(e).primary_link().node().link().edge().index();
        }
        i1.push_back(e);
    }

    if (!tree.edge_at(regraftEdgeIdx).primary_link().node().is_root()) {
        size_t e = tree.edge_at(regraftEdgeIdx).primary_link().node().link().edge().index();
        while (!tree.edge_at(e).primary_link().node().is_root()) {
            i2.push_back(e);
            e = tree.edge_at(e).primary_link().node().link().edge().index();
        }
        i2.push_back(e);
    }

    while (i1.size() > 0 and i2.size() > 0 and i1.back() == i2.back()) {
        i1.pop_back();
        i2.pop_back();
    }

    for (size_t i = 0; i < i1.size(); ++i) if (lqic[i1[i]] < 0) return true;
    for (size_t i = 0; i < i2.size(); ++i) if (lqic[i2[i]] < 0) return true;

    return false;
}

void spr(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx) {
    size_t pruneLinkIdx = tree.edge_at(pruneEdgeIdx).primary_link().index();
    LOG_DBG << "SPR(" << pruneEdgeIdx << " " << regraftEdgeIdx << ")" << std::endl;

    if (!tree.link_at(pruneLinkIdx).node().is_root()) {
        // Case 1: Prune Node is not root
        size_t link_prune_no = tree.link_at(pruneLinkIdx).next().outer().index();
        size_t link_prune_nno = tree.link_at(pruneLinkIdx).next().next().outer().index();
        size_t link_prune_parent, link_prune_sibling;
        if (link_prune_no == tree.link_at(pruneLinkIdx).node().link().edge().primary_link().index()) {
            link_prune_parent = link_prune_no;
            link_prune_sibling = link_prune_nno;
        } else {
            link_prune_parent = link_prune_nno;
            link_prune_sibling = link_prune_no;
        }

        size_t edge_prune_sibling = tree.link_at(link_prune_sibling).edge().index();
        size_t edge_prune_parent = tree.link_at(link_prune_parent).edge().index();
        size_t link_prune_parent_secondary = tree.edge_at(edge_prune_parent).secondary_link().index();
        size_t link_regraft_secondary = tree.edge_at(regraftEdgeIdx).secondary_link().index();

        reconnect_node_secondary(tree, edge_prune_parent, link_prune_sibling);
        reconnect_node_secondary(tree, regraftEdgeIdx, link_prune_parent_secondary);
        reconnect_node_secondary(tree, edge_prune_sibling, link_regraft_secondary);

        swap_edges(tree, regraftEdgeIdx, edge_prune_parent);
    }

    else {
        // Case 2: Prune Node is root
        size_t link_prune_no = tree.link_at(pruneLinkIdx).next().outer().index();
        size_t link_prune_nno = tree.link_at(pruneLinkIdx).next().next().index();
        size_t edge_prune_n = tree.link_at(link_prune_no).edge().index();
        size_t edge_prune_nn = tree.link_at(link_prune_nno).edge().index();
        size_t link_regraft_secondary = tree.edge_at(regraftEdgeIdx).secondary_link().index();

        bool regraftInNext = false;
        for(auto it : eulertour(tree.link_at(pruneLinkIdx).next().outer().next())) {
            if (it.link().index() == tree.link_at(pruneLinkIdx).next().next().index()) break;
            if (it.edge().index() == regraftEdgeIdx) regraftInNext = true;
        }

        size_t link_prune_n_sec = tree.edge_at(edge_prune_n).secondary_link().index();
        size_t link_prune_nn_sec = tree.edge_at(edge_prune_nn).secondary_link().index();
        if (regraftInNext) {
            reconnect_node_secondary(tree, edge_prune_n, link_regraft_secondary);
            reconnect_node_secondary(tree, edge_prune_nn, link_prune_n_sec);
            reconnect_node_secondary(tree, regraftEdgeIdx, link_prune_nn_sec);
        } else {
            reconnect_node_secondary(tree, edge_prune_nn, link_regraft_secondary);
            reconnect_node_secondary(tree, edge_prune_n, link_prune_nn_sec);
            reconnect_node_secondary(tree, regraftEdgeIdx, link_prune_n_sec);
        }
    }
}


#endif
