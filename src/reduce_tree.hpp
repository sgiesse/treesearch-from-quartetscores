#ifndef REDUCE_TREE_HPP
#define REDUCE_TREE_HPP

//#include "treesearch.hpp"
#include "starttree.hpp"

void calculateAveragePairwiseDistance(std::vector<std::vector<double> >& D, Tree refTree, std::string pathToEvaluationTrees) {
	std::vector<Tree> evalTrees;

	size_t N = refTree.node_count();
	std::vector<std::vector<size_t> > occurInEval(N);
	D.resize(N);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			occurInEval[i].push_back(0);
			D[i].push_back(0);
		}
	}

	utils::InputStream instream(utils::make_unique<utils::FileInputSource>(pathToEvaluationTrees));
	auto itTree = NewickInputIterator(instream, DefaultTreeNewickReader());
	while (itTree) {
		evalTrees.push_back(*itTree);
		++itTree;
	}

	// Mapping node index <-> leaf name
	std::map<std::string, size_t> mapLeaves;
	std::map<size_t, std::string> revMapLeaves;
	for (size_t i = 0; i < refTree.node_count(); ++i) {
		if (refTree.node_at(i).is_leaf()){
			mapLeaves[refTree.node_at(i).data<DefaultNodeData>().name] = i;
			revMapLeaves[i] = refTree.node_at(i).data<DefaultNodeData>().name;
		}
	}

	// count how often pairs appear together
	for (Tree tree : evalTrees) {
		for (size_t i = 0; i < tree.node_count(); ++i) {
			if (! tree.node_at(i).is_leaf()) continue;
			for (size_t j = 0; j < tree.node_count(); ++j) {
				if (! tree.node_at(j).is_leaf()) continue;
				size_t u = mapLeaves[tree.node_at(i).data<DefaultNodeData>().name];
				size_t v = mapLeaves[tree.node_at(j).data<DefaultNodeData>().name];
				occurInEval[u][v]++;
				occurInEval[v][u]++;
			}
		}
	}

	for (Tree tree : evalTrees) {
		// TreeInformation (for LCA)
		TreeInformation tinf;
		tinf.init(tree);

		// add pairwise distances to D
		for (size_t i = 0; i < tree.node_count(); ++i) {
			if (! tree.node_at(i).is_leaf()) continue;
			for (size_t j = 0; j < tree.node_count(); ++j) {
				if (! tree.node_at(j).is_leaf()) continue;
				size_t u = mapLeaves[tree.node_at(i).data<DefaultNodeData>().name];
				size_t v = mapLeaves[tree.node_at(j).data<DefaultNodeData>().name];
				double d = tinf.distanceInEdges(i, j);
        if (d < 2 and i != j) throw std::runtime_error("Distance too small.");
				D[u][v] += d;
				D[v][u] += d;
			}
		}
	}
	for (size_t i = 0; i < D.size(); ++i) {
		for (size_t j = 0; j < D.size(); ++j) {
        if (occurInEval[i][j] > 0) 
            D[i][j] /= occurInEval[i][j];

        else D[i][j] = std::numeric_limits<double>::infinity();
        if (D[i][j] < 2 and i != j) throw std::runtime_error("Avg. Distance too small.");
		}
	}
}

std::vector<std::vector<std::string> > leaf_sets(std::string pathToEvaluationTrees) {
	Tree r_tree = random_tree(pathToEvaluationTrees);
    std::vector<std::vector<double> > D;
    calculateAveragePairwiseDistance(D, r_tree, pathToEvaluationTrees);

	std::vector<std::vector<size_t> > sets(r_tree.node_count());
	for (size_t i = 0; i < r_tree.node_count(); ++i) {
		if (r_tree.node_at(i).is_leaf()) sets[i].push_back(i);
	}
    const double D_MAX = 4;
    bool continueMergingLeaves = true;
    while (continueMergingLeaves) {
    	size_t best_i, best_j;
    	double max = D_MAX;

	    for (size_t i = 0; i < r_tree.node_count(); ++i) {
	        if (! r_tree.node_at(i).is_leaf()) continue;
	        if (sets[i].empty()) continue;
        	for (size_t j = i+1; j < r_tree.node_count(); ++j) {
	            if (! r_tree.node_at(j).is_leaf()) continue;
	            if (sets[j].empty()) continue;
	            if (D[i][j] < max) {
	            	max = D[i][j];
	            	best_i = i;
	            	best_j = j;
	            }
	        }
	    }

	    if (max < D_MAX) {
	    	for (size_t k = 0; k < r_tree.node_count(); ++k) {
	    		if (! r_tree.node_at(k).is_leaf()) continue;
	    		double sum = D[best_i][k];
	    		for (auto x : sets[best_j]) sum += D[x][k];
    			D[best_i][k] = sum / (sets[best_j].size() + 1);
	    	}
	    	for (auto x : sets[best_j]) sets[best_i].push_back(x);
    		sets[best_j].clear();
	    }
	    else continueMergingLeaves = false;
	}

	std::vector<std::vector<std::string> > leafSets;
	for (size_t i = 0; i < sets.size(); ++i) {
		if (sets[i].size() == 0) continue;
		leafSets.push_back(std::vector<std::string>());
		
		for (size_t j = 0; j < sets[i].size(); ++j) 
			leafSets[leafSets.size()-1].push_back(r_tree.node_at(sets[i][j]).data<DefaultNodeData>().name);
	}
	for (auto x : leafSets) {
		for (auto y : x) std::cout << y << " ";
		std::cout << std::endl << std::endl;
	}
	return leafSets;
}

/*std::pair<Tree, std::vector<std::vector<std::string> > > cluster_tree(const std::string & pathToEvaluationTrees) {
    std::vector<std::vector<std::string> > leafSets = leaf_sets(pathToEvaluationTrees);
    std::vector<std::string> leaves;
    for (auto x : leafSets) leaves.push_back(x[0]);
    std::shuffle(leaves.begin(), leaves.end(), Random::getMT());
    return std::pair<Tree, std::vector<std::vector<std::string> > >(random_tree_from_leaves(leaves), leafSets);
    }*/

Tree expanded_cluster_tree(Tree& clusterTree, std::vector<std::vector<std::string> >& leafSets) {
    for (size_t i = 0; i < clusterTree.node_count(); ++i) {
        if (clusterTree.node_at(i).is_leaf()) {
            size_t l_i = 0;
            while (l_i < leafSets.size() and leafSets[l_i][0] != clusterTree.node_at(i).data<DefaultNodeData>().name) ++l_i;
            if (l_i == leafSets.size()) continue; 
            //for (size_t j = 1; j < leafSets[l_i].size(); ++j) {
        	for (int j = leafSets[l_i].size()-1; j >= 1; --j) {
                add_new_node(clusterTree, clusterTree.node_at(i).link().edge()).secondary_link().node().data_cast<DefaultNodeData>()->name = leafSets[l_i][j];
        	}
        }
    }

    return clusterTree;
}


#endif
