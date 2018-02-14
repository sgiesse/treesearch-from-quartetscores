#ifndef SPR_NNI
#define SPR_NNI

#include "tree_operations.hpp"

void new_spr(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx);
bool validSprMove(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx);

template<typename CINT>
void spr_lqic_update(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx, size_t lcaIdx, QuartetScoreComputer<CINT>& qsc) {

    for (auto it : path_set(
             tree.edge_at(pruneEdgeIdx).secondary_link().node(),
             tree.edge_at(regraftEdgeIdx).secondary_link().node(),
             tree.node_at(lcaIdx))) {
        //if (it.is_lca()) continue;
        //u_link = &it.link().outer();

        std::cout << it.edge().index() << " ";
        qsc.recomputeLqicForEdge(tree, it.edge().index());
    }
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

void new_spr(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx) {
    size_t pruneLinkIdx = tree.edge_at(pruneEdgeIdx).primary_link().index();
    size_t regraftLinkIdx = tree.edge_at(regraftEdgeIdx).primary_link().index();
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
        std::cout << link_prune_no << " " << link_prune_nno << " " << link_prune_parent << " " << link_prune_sibling << std::endl;

        size_t link_regraft_no = tree.link_at(regraftLinkIdx).next().outer().index();
        size_t link_regraft_nno = tree.link_at(regraftLinkIdx).next().next().outer().index();

        size_t edge_prune_sibling = tree.link_at(link_prune_sibling).edge().index();
        size_t edge_prune_parent = tree.link_at(link_prune_parent).edge().index();
        size_t link_prune_parent_secondary = tree.edge_at(edge_prune_parent).secondary_link().index();
        size_t link_regraft_secondary = tree.edge_at(regraftEdgeIdx).secondary_link().index();

        reconnect_node_secondary(tree, edge_prune_parent, link_prune_sibling);
        reconnect_node_secondary(tree, regraftEdgeIdx, link_prune_parent_secondary);
        reconnect_node_secondary(tree, edge_prune_sibling, link_regraft_secondary);
    }

    else {
        // Case 2: Prune Node is root

        size_t link_prune_no = tree.link_at(pruneLinkIdx).next().outer().index();
        size_t link_prune_nno = tree.link_at(pruneLinkIdx).next().next().index();
        size_t link_regraft_no = tree.link_at(regraftLinkIdx).next().outer().index();
        size_t link_regraft_nno = tree.link_at(regraftLinkIdx).next().next().outer().index();
        size_t edge_prune_n = tree.link_at(link_prune_no).edge().index();
        size_t edge_prune_nn = tree.link_at(link_prune_nno).edge().index();
        size_t link_regraft_secondary = tree.edge_at(regraftEdgeIdx).secondary_link().index();

        reconnect_node_primary(tree, edge_prune_nn, link_prune_no);
        reconnect_node_secondary(tree, regraftEdgeIdx, link_prune_nno);
        reconnect_node_secondary(tree, edge_prune_n, link_regraft_secondary);

        int c = 0;
        for(auto it : preorder(tree.root_node())) {
            if (it.link().index() != it.edge().secondary_link().index()
                and it.node().index() != tree.root_node().index()) {
                // revert edge direction
                size_t lp = it.edge().primary_link().index();
                size_t ls = it.edge().secondary_link().index();
                it.edge().reset_primary_link(&tree.link_at(ls));
                it.edge().reset_secondary_link(&tree.link_at(lp));
            }
            if (it.node().primary_link().index() != it.edge().primary_link().index()
                and it.node().index() != tree.root_node().index()) {
                it.node().reset_primary_link(&it.edge().secondary_link());
            }
            c++;
            if (c > (int)tree.link_count()) {
                LOG_ERR << "SPR done: preorder endless loop detected" << std::endl;
                validate_topology(tree);
                LOG_ERR << PrinterTable().print(tree);
                //return false;
                throw std::runtime_error("SPR not possible");
            }
        }

    }
}


#endif
