#ifndef DISTANCE_SCORE_HPP
#define DISTANCE_SCORE_HPP

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
		for (size_t i = 0; i < refTree.node_count(); ++i) {
			if (! refTree.node_at(i).is_leaf()) continue;
			for (size_t j = 0; j < refTree.node_count(); ++j) {
				if (! refTree.node_at(j).is_leaf()) continue;
				size_t u = mapLeaves[refTree.node_at(i).data<DefaultNodeData>().name];
				size_t v = mapLeaves[refTree.node_at(j).data<DefaultNodeData>().name];
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
		}
	}
}


double distance_score4(std::vector<std::vector<double> >& D, Tree tree) {
	size_t N = tree.node_count();
	std::vector<std::vector<double> > Dt(N);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) Dt[i].push_back(0);
	}

	TreeInformation tinf;
	tinf.init(tree);
	for (size_t i = 0; i < tree.node_count(); ++i) {
		if (! tree.node_at(i).is_leaf()) continue;
		for (size_t j = 0; j < tree.node_count(); ++j) {
			if (! tree.node_at(j).is_leaf()) continue;
			double d = tinf.distanceInEdges(i, j);
			d = abs(d-D[i][j]);
			Dt[i][j] += d;
			Dt[j][i] += d;
		}	
	}

	double sum = 0;
	for (size_t i = 0; i < N; ++i) {
		// select and sum max 4 in D[i]
		if (N < 4) for (size_t j = 0; j < N; ++j) sum += Dt[i][j];
		else {
			std::vector<double> Dti(Dt[i].begin(), Dt[i].end());
			std::nth_element(Dti.begin(), Dti.begin()+N-4, Dti.end());
			for (auto it = Dti.begin()+N-4; it != Dti.end(); ++it) sum += *it;
		}
	}

	return sum;
}




#endif