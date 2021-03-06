#ifndef SPR_ITERATOR_HPP
#define SPR_ITERATOR_HPP

#include "TreeInformation.hpp"

#include "tree_operations.hpp"

class SPRiterator;

class SPRtree {
    friend class SPRiterator;
private:
    Tree tree;
    Tree current;
    std::vector<bool> spr_ok;
    size_t i;
    std::vector<double> lqic;
    void fill_spr_ok_for_i();
public:
    SPRtree(Tree& _tree, size_t i);
    SPRtree(const SPRtree& m);
    SPRtree();
    Tree get();
    void next_i();
    SPRiterator begin();
    SPRiterator end();
    void restrict_by_lqic(std::vector<double>& _lqic);
};

class SPRiterator {
private:
    SPRtree tree;
    size_t i;
    size_t j;
    void jump_to_next();
public:
    SPRiterator(SPRtree _tree);
    SPRiterator(SPRtree _tree, size_t _i, size_t _j);
    const SPRtree& operator*() const;
    SPRiterator* operator++();
    bool operator==(const SPRiterator& rhs) const;
    bool operator!=(const SPRiterator& rhs) const;
};

// --- implementation for SPRtree
SPRtree::SPRtree(Tree& _tree, size_t _i = 0) {
    tree = _tree;
    current = tree;
    i = _i;
    fill_spr_ok_for_i();
}

SPRtree::SPRtree(const SPRtree& m) {
    tree = m.tree;
    current = m.current;
    i = m.i;
    spr_ok = m.spr_ok;
    lqic = m.lqic;
}

SPRtree::SPRtree() { }

Tree SPRtree::get(){
    return current;
}

void SPRtree::next_i() {
    i++;
    LOG_DBG << "SPRtree::next_i() --> " << i;
    if (i < tree.edge_count()) fill_spr_ok_for_i();
}

void SPRtree::fill_spr_ok_for_i() {
    spr_ok = std::vector<bool> (tree.edge_count(), true);
    for (auto it : eulertour(tree.edge_at(i).primary_link())) {
        if (it.edge().index() == i and it.link().index() == it.edge().secondary_link().index()) break;
        spr_ok[it.edge().index()] = false;
    }
    size_t j = i;
    do {
        spr_ok[j] = false;
        j = tree.edge_at(j).primary_link().node().link().edge().index();
    } while (!tree.edge_at(j).primary_link().node().is_root());
    spr_ok[j] = false;

    spr_ok[tree.edge_at(i).primary_link().next().edge().index()] = false;
    spr_ok[tree.edge_at(i).primary_link().next().next().edge().index()] = false;

    if (lqic.size() > 0) {
        TreeInformation t_inf;
        t_inf.init(tree);

        for (size_t j = 0; j < lqic.size(); ++j) {
            size_t lca_idx = t_inf.lowestCommonAncestorIdx(
                tree.edge_at(i).primary_link().node().index(),
                tree.edge_at(j).primary_link().node().index(), tree.root_node().index());
            
            bool skip = true;
            if (!tree.edge_at(i).primary_link().node().is_root()) {
                if (tree.edge_at(i).primary_link().node().index() != lca_idx) {
                    size_t e = tree.edge_at(i).primary_link().node().link().edge().index();
                    while (tree.edge_at(e).secondary_link().node().index() != lca_idx
                           and e != tree.edge_at(i).primary_link().node().link().edge().index()) {
                        if (lqic[e] < 0) skip = false;
                        e = tree.edge_at(i).primary_link().node().link().edge().index();
                    }
                }
            }
            if (!tree.edge_at(j).primary_link().node().is_root()) {
                if (tree.edge_at(j).primary_link().node().index() != lca_idx) {
                    size_t e = tree.edge_at(j).primary_link().node().link().edge().index();
                    while (tree.edge_at(e).secondary_link().node().index() != lca_idx
                           and e != tree.edge_at(i).primary_link().node().link().edge().index()) {
                        if (lqic[e] < 0) skip = false;
                        e = tree.edge_at(i).primary_link().node().link().edge().index();
                    }
                }
            }
            if (skip) spr_ok[j] = false;
        }
    }
}

SPRiterator SPRtree::begin() {
    return SPRiterator(*this);
}

SPRiterator SPRtree::end() {
    return SPRiterator(*this, tree.edge_count(), 0);
}

void SPRtree::restrict_by_lqic(std::vector<double>& _lqic) {
    lqic = _lqic;
    fill_spr_ok_for_i();
}

// --- end implementation for SPRtree

// --- implementation for SPRiterator

SPRiterator::SPRiterator(SPRtree _tree) : SPRiterator::SPRiterator(_tree, 0, 0) { }

SPRiterator::SPRiterator(SPRtree _tree, size_t _i, size_t _j) {
    tree = _tree;
    i = _i;
    j = _j;
    if (!tree.spr_ok[j]) jump_to_next();
}

const SPRtree& SPRiterator::operator*() const {
    return tree;
}

void SPRiterator::jump_to_next() {
    while (true) {
        j++;
        if (j >= tree.tree.edge_count()) {
            i++; j = 0;
            if (i >= tree.tree.edge_count()) {
                //tree = nullptr;
                break;
            }
            tree.next_i();
        }
        if (tree.spr_ok[j]) {
            break;
        }
    }

}

SPRiterator* SPRiterator::operator++() {
    jump_to_next();
    //if (tree != nullptr) {
    if (i < tree.tree.edge_count()) {
        Tree t(tree.tree);
        spr(t, i, j);
        tree.current = t;
    }
    return new SPRiterator(tree, i, j);
}

bool SPRiterator::operator==(const SPRiterator& rhs) const {
    return (i == rhs.i and j == rhs.j);
}

bool SPRiterator::operator!=(const SPRiterator& rhs) const {
    if (i >= tree.tree.edge_count() and rhs.i >= rhs.tree.tree.edge_count())
        return false;
    return i != rhs.i or j != rhs.j;
}

// --- end implementation for SPRiterator

#endif
