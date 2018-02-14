#ifndef TREEOPERATIONS_HPP
#define TREEOPERATIONS_HPP

#include <algorithm>
#include <cmath>
#include <memory>
#include "genesis/tree/function/manipulation.hpp"
#include "utils.hpp"

// --------- Forward Declarations
bool spr(Tree& tree, int i, int j);
std::vector<std::string> leafNames(const std::string &evalTreesPath);
bool verify_leaf_ids_match(Tree tree1, Tree tree2, bool verbose);

// -----------------------------

void reconnect_node_primary(Tree& tree, size_t edge, size_t new_primary_link) {
    std::cout << "reconnect_node_primary " << edge << " " << new_primary_link << std::endl;
    size_t secondary_link = tree.edge_at(edge).secondary_link().index();
    tree.edge_at(edge).reset_primary_link(&tree.link_at(new_primary_link));
    tree.link_at(new_primary_link).reset_edge(&tree.edge_at(edge));
    tree.link_at(new_primary_link).reset_outer(&tree.link_at(secondary_link));
    tree.link_at(secondary_link).reset_outer(&tree.link_at(new_primary_link));
}

void reconnect_node_secondary(Tree& tree, size_t edge, size_t new_secondary_link) {
    std::cout << "reconnect_node_secondary " << edge << " " << new_secondary_link << std::endl;
    size_t primary_link = tree.edge_at(edge).primary_link().index();
    tree.edge_at(edge).reset_secondary_link(&tree.link_at(new_secondary_link));
    tree.link_at(new_secondary_link).reset_edge(&tree.edge_at(edge));
    tree.link_at(new_secondary_link).reset_outer(&tree.link_at(primary_link));
    tree.link_at(primary_link).reset_outer(&tree.link_at(new_secondary_link));
}

bool spr(Tree& tree, size_t pruneEdgeIdx, size_t regraftEdgeIdx) {
    size_t pruneLinkIdx = tree.edge_at(pruneEdgeIdx).primary_link().index();
    size_t regraftLinkIdx = tree.edge_at(regraftEdgeIdx).primary_link().index();
    LOG_DBG << "SPR(" << pruneEdgeIdx << " " << regraftEdgeIdx << ")" << std::endl;
    size_t l0, l2, l3, l7, l10, l25;

    if (tree.edge_at(pruneEdgeIdx).primary_link().node().index() == tree.edge_at(regraftEdgeIdx).primary_link().node().index()) {
        LOG_WARN << "SPR not possible: 1\n";
        return false;
    }

    if (!tree.link_at(pruneLinkIdx).node().is_root()) {

        // Case 1: Prune Node is not root

        bool child = false;
        for(auto it : eulertour(tree.link_at(pruneLinkIdx))) {
            if (it.link().index() == tree.link_at(pruneLinkIdx).next().index()) break;
            if (it.edge().index() == regraftEdgeIdx) child = true;
        }
        bool parent = false;
        for(auto it : eulertour(tree.link_at(regraftLinkIdx))) {
            if (it.link().index() == tree.link_at(regraftLinkIdx).next().index()) break;
            if (it.edge().index() == pruneEdgeIdx) parent = true;
        }

        if (child or parent) {
            LOG_WARN << "SPR not possible: 2 ";
            if (child) { LOG_WARN << "child\n"; }
            if (parent) { LOG_WARN << "parent\n"; }
            return false;
        }

        size_t lps = tree.link_at(pruneLinkIdx).outer().index();
        size_t lrs = tree.link_at(regraftLinkIdx).outer().index();

        tree.edge_at(pruneEdgeIdx).reset_secondary_link(&tree.link_at(lrs));
        tree.edge_at(regraftEdgeIdx).reset_secondary_link(&tree.link_at(lps));

        tree.link_at(lrs).reset_edge(&tree.edge_at(pruneEdgeIdx));
        tree.link_at(lps).reset_edge(&tree.edge_at(regraftEdgeIdx));

        tree.link_at(pruneLinkIdx).reset_outer(&tree.link_at(lrs));
        tree.link_at(lrs).reset_outer(&tree.link_at(pruneLinkIdx));
        tree.link_at(regraftLinkIdx).reset_outer(&tree.link_at(lps));
        tree.link_at(lps).reset_outer(&tree.link_at(regraftLinkIdx));

        if (!validate_topology(tree)) {
            LOG_ERR << "SPR done: not valid" << std::endl;
            LOG_ERR << PrinterTable().print(tree);
            throw std::runtime_error("SPR not possible");
            return false;
        }
    }

    else {
        // Case 2: Prune Node is root
        bool regraftInNext = false;
        for(auto it : eulertour(tree.link_at(pruneLinkIdx).next().outer().next())) {
            if (it.link().index() == tree.link_at(pruneLinkIdx).next().next().index()) break;
            if (it.edge().index() == regraftEdgeIdx) regraftInNext = true;
        }
        bool regraftInNextNext = false;
        for(auto it : eulertour(tree.link_at(pruneLinkIdx).next().next().outer().next())) {
            if (it.link().index() == tree.link_at(pruneLinkIdx).index()) break;
            if (it.edge().index() == regraftEdgeIdx) regraftInNextNext = true;
        }
        if (!regraftInNext and !regraftInNextNext) {
            LOG_INFO << "SPR not possible: 3" << std::endl;
            return false;
            //throw std::runtime_error("SPR not possible");
        }

        if (regraftInNext) {
            l0 = tree.link_at(pruneLinkIdx).next().next().index();
            l10 = tree.link_at(regraftLinkIdx).outer().index();
            l2 = tree.link_at(pruneLinkIdx).next().index();
            l3 = tree.link_at(pruneLinkIdx).next().outer().index();
            l7 = regraftLinkIdx;
            l25 = tree.link_at(pruneLinkIdx).next().next().outer().index();
        } else {
            l0 = tree.link_at(pruneLinkIdx).next().index();
            l10 = tree.link_at(regraftLinkIdx).outer().index();
            l2 = tree.link_at(pruneLinkIdx).next().next().index();
            l3 = tree.link_at(pruneLinkIdx).next().next().outer().index();
            l7 = regraftLinkIdx;
            l25 = tree.link_at(pruneLinkIdx).next().outer().index();
        }

        size_t e0 = tree.link_at(l2).edge().index();
        size_t e3 = tree.link_at(l7).edge().index();
        size_t e12 = tree.link_at(l0).edge().index();

        tree.edge_at(e0).reset_secondary_link(&tree.link_at(l7));
        tree.edge_at(e3).reset_primary_link(&tree.link_at(l3));
        tree.edge_at(e3).reset_secondary_link(&tree.link_at(l25));
        tree.edge_at(e12).reset_secondary_link(&tree.link_at(l10));

        tree.link_at(l3).reset_edge(&tree.edge_at(e3));
        tree.link_at(l7).reset_edge(&tree.edge_at(e0));
        tree.link_at(l10).reset_edge(&tree.edge_at(e12));
        tree.link_at(l25).reset_edge(&tree.edge_at(e3));

        tree.link_at(l2).reset_outer(&tree.link_at(l7));
        tree.link_at(l7).reset_outer(&tree.link_at(l2));
        tree.link_at(l3).reset_outer(&tree.link_at(l25));
        tree.link_at(l25).reset_outer(&tree.link_at(l3));
        tree.link_at(l0).reset_outer(&tree.link_at(l10));
        tree.link_at(l10).reset_outer(&tree.link_at(l0));

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

        if (!validate_topology(tree)) {
            LOG_WARN << "SPR done: not valid" << std::endl;
            LOG_WARN << l0 << " "  << l2 << " "  << l3 << " "  << l7 << " "  << l10 << " "  << l25 << std::endl;
            LOG_WARN << PrinterTable().print(tree);
            throw std::runtime_error("SPR not possible");
            return false;
        }
    }

    LOG_DBG4 << "SPR done: success" << std::endl;

    return true;
}

std::vector<std::string> leafNames(const std::string &evalTreesPath) {
    std::set<std::string> leaf_names;
    utils::InputStream instream(utils::make_unique<utils::FileInputSource>(evalTreesPath));
    auto itTree = NewickInputIterator(instream, DefaultTreeNewickReader());
    while (itTree) {
        Tree const& tree = *itTree;
        for(auto const& node : tree.nodes()) {
            if (node->is_leaf()) {
                auto const& name = node->data<DefaultNodeData>().name;
                leaf_names.insert(name);
            }
        }
        ++itTree;
    }

    return std::vector<std::string>(leaf_names.begin(), leaf_names.end());
}

bool verify_leaf_ids_match(Tree tree1, Tree tree2, bool verbose = false) {
    std::map<std::string, size_t> leafToID1;
    std::map<std::string, size_t> leafToID2;

    for (size_t i = 0; i < tree1.node_count(); ++i) {
        if (tree1.node_at(i).is_leaf())
            leafToID1[tree1.node_at(i).data<DefaultNodeData>().name] = i;
    }
    for (size_t i = 0; i < tree2.node_count(); ++i) {
        if (tree2.node_at(i).is_leaf())
            leafToID2[tree2.node_at(i).data<DefaultNodeData>().name] = i;
    }

    bool ok = true;
    if (leafToID1.size() < leafToID2.size()) {
        for (auto it = leafToID1.begin(); it != leafToID1.end(); ++it) {
            if (leafToID2.find(it->first) == leafToID2.end() or
                leafToID1[it->first] != leafToID2[it->first]) {
                if (verbose) {
                    std::cout << "leafToID1[" << it->first << "] = " << leafToID1[it->first] << " != leafToID2[" << it->first << "] = "  << leafToID2[it->first] << std::endl;
                } else {
                    return false;
                }
            }
        }
    } else {
        for (auto it = leafToID2.begin(); it != leafToID2.end(); ++it) {
            if (leafToID1.find(it->first) == leafToID1.end() or
                leafToID1[it->first] != leafToID2[it->first]) {
                if (verbose) {
                    std::cout << "leafToID1[" << it->first << "] = " << leafToID1[it->first] << " != leafToID2[" << it->first << "] = "  << leafToID2[it->first] << std::endl;
                } else {
                    return false;
                }
            }
        }
    }

    return ok;
}

#endif
