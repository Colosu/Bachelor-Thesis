/*
 * Checkups.cpp
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#include <fst/fstlib.h>
#include "Graph.h"
#include "Checkups.h"
#include "Operations.h"

namespace fst {

Checkups::Checkups() {

}

Checkups::~Checkups() {

}

bool Checkups::are_equivalent(Graph* g1, Graph* g2) {

	if (Equal(*(g1->getTransducer()->Copy()), *(g2->getTransducer()->Copy()))) {
		cout << "are equal" << endl;
		return true;
	} else {
		return false;
	}
}

bool Checkups::is_valid(Graph* g) {

	if (g->getTransducer()->Properties(kIDeterministic, true) == kIDeterministic) {
		if (g->getTransducer()->Properties(kAccessible, true) == kAccessible) {
			if (g->getTransducer()->Properties(kCoAccessible, true) == kCoAccessible) {
				return true;
			}
		}
	}
	return false;
}

bool Checkups::is_validMutation(Graph* g) {

	if (g->getTransducer()->Properties(kIDeterministic, true) == kIDeterministic) {
		if (g->getTransducer()->Properties(kCoAccessible, true) == kCoAccessible) {
			return true;
		}
	}
	return false;
}

bool Checkups::has_FEP(Graph* g1, Graph* g2, int length, int &k) {

	StdMutableFst* orig = g1->getTransducer()->Copy();
	StdMutableFst* mut = g2->getTransducer()->Copy();
	k = 0;

	return FEPaux(orig, mut, orig->Start(), mut->Start(), "", "", 0, k, length);
}


bool Checkups::FEPaux(StdMutableFst* fsm1, StdMutableFst* fsm2, int qid1, int qid2, string output1, string output2, int iter, int &k, int length) {

	bool result = false;
	int it = iter;
	int q1 = qid1;
	int q2 = qid2;
	string out1 = output1;
	string out2 = output2;
	string o1 = out1;
	string o2 = out2;

	while (it <= length) {
		if (out1 != out2) {
			return false;
		} else if (k != 0) {
			result = true;
		}
		if (q1 != q2) {
			if (k == 0 || it < k) {
				k = it;
			}
			result = true;
		}

		if (it == length) {
			return result;
		}

		ArcIterator<StdMutableFst> arcIter1(*fsm1, q1);
		if (!arcIter1.Done()) {
			const StdArc &arc1 = arcIter1.Value();
			Matcher<StdMutableFst> matcher(*fsm2, MATCH_INPUT);
			matcher.SetState(q2);
			if (!matcher.Done()) {
				if (matcher.Find(arc1.ilabel)) {
					const StdArc &arc2 = matcher.Value();
					q1 = arc1.nextstate;
					q2 = arc2.nextstate;
					o1 = out1;
					o2 = out2;
					out1 = out1 + std::to_string(arc1.olabel);
					out2 = out2 + std::to_string(arc2.olabel);
					matcher.Next();
					for (; !matcher.Done(); matcher.Next()) {
						result = result || FEPaux(fsm1, fsm2, arc1.nextstate, arc2.nextstate, o1 + std::to_string(arc1.olabel), o2 + std::to_string(arc2.olabel), it + 1, k, length);
					}
				}
			}

			arcIter1.Next();
			for (; !arcIter1.Done(); arcIter1.Next()) {
				const StdArc &arc1 = arcIter1.Value();
				Matcher<StdMutableFst> matcher(*fsm2, MATCH_INPUT);
				matcher.SetState(q2);
				if (matcher.Find(arc1.ilabel)) {
					for (; !matcher.Done(); matcher.Next()) {
						const StdArc &arc2 = matcher.Value();
						result = result || FEPaux(fsm1, fsm2, arc1.nextstate, arc2.nextstate, o1 + std::to_string(arc1.olabel), o2 + std::to_string(arc2.olabel), it + 1, k, length);
					  }
				}
			}
		}
		it++;
	}


	return result;
}

} /* namespace std */
