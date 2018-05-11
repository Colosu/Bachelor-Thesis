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

	if (Equivalent(*(g1->getTransducer()->Copy()), *(g2->getTransducer()->Copy()))) {
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
		//if (g->getTransducer()->Properties(kCoAccessible, true) == kCoAccessible) {
			return true;
		//}
	}
	return false;
}

bool Checkups::is_validInput(Graph* g, string input, int length) {

	int qid = g->getTransducer()->Start();
	for (int i = 0; i < length; i++) {
		Matcher<StdMutableFst> matcher(*g->getTransducer(), MATCH_INPUT);
		matcher.SetState(qid);
		if (matcher.Find(input[i])) {
			const StdArc &arc1 = matcher.Value();
			qid = arc1.nextstate;
		} else {
			return false;
		}
	}
	return true;
}

double Checkups::PFEPState(Graph* g1, Graph* g2, int length) {

	StdMutableFst* orig = g1->getTransducer()->Copy();
	StdMutableFst* mut = g2->getTransducer()->Copy();
	int sucess = 0;
	int fails = 0;
	double result = 0;

	PFEPauxState(orig, mut, orig->Start(), mut->Start(), "", "", 0, sucess, fails, length, false);

	delete orig;
	delete mut;
	if(sucess != 0) {
		result = sucess / (double)(sucess + fails);
	}
	//if (result < 0.0001) {
	//	result = 0;
	//}
	return result;
}

bool Checkups::has_FEPState(Graph* g1, Graph* g2, int length, int &k) {

	StdMutableFst* orig = g1->getTransducer()->Copy();
	StdMutableFst* mut = g2->getTransducer()->Copy();
	k = 0;
	bool result = FEPauxState(orig, mut, orig->Start(), mut->Start(), "", "", 0, k, length);

	delete orig;
	delete mut;
	return result;
}

bool Checkups::has_FEPInput(Graph* g, string input1, string input2, int length) {

	StdMutableFst* fsm = g->getTransducer()->Copy();

	string output1 = "";
	string output2 = "";
	int qid = fsm->Start();

	for (int i = 0; i < length; i++) {
		Matcher<StdMutableFst> matcher(*fsm, MATCH_INPUT);
		matcher.SetState(qid);
		if (matcher.Find(input1[i])) {
			const StdArc &arc = matcher.Value();
			output1 += arc.olabel;
			qid = arc.nextstate;
		}
	}

	qid = fsm->Start();
	for (int i = 0; i < length; i++) {
		Matcher<StdMutableFst> matcher(*fsm, MATCH_INPUT);
		matcher.SetState(qid);
		if (matcher.Find(input2[i])) {
			const StdArc &arc = matcher.Value();
			output2 += arc.olabel;
			qid = arc.nextstate;
		}
	}

	delete fsm;
	if (output1 == output2) {
		return true;
	} else {
		return false;
	}
}


void Checkups::PFEPauxState(StdMutableFst* fsm1, StdMutableFst* fsm2, int qid1, int qid2, string output1, string output2, int iter, int &sucess, int &fails, int length, bool aux) {

	bool zero = true;
	if (output1 != output2) {
		fails++;
		return;
	}
	if (!aux && qid1 != qid2) {
		aux = true;
	}
	if (iter == length) {
		if(aux) {
			sucess++;
		}
		return;
	}

	for (ArcIterator<StdMutableFst> arcIter1(*fsm1, qid1); !arcIter1.Done(); arcIter1.Next()) {
		const StdArc &arc1 = arcIter1.Value();
		for (ArcIterator<StdMutableFst> arcIter2(*fsm2, qid2); !arcIter2.Done(); arcIter2.Next()) {
			const StdArc &arc2 = arcIter2.Value();
			if (arc1.ilabel == arc2.ilabel) {
				zero = false;
				PFEPauxState(fsm1, fsm2, arc1.nextstate, arc2.nextstate, output1 + (char)arc1.olabel, output2 + (char)arc2.olabel, iter + 1, sucess, fails, length, aux);
			}
		}
	}
	if(zero && aux) {
		sucess++;
	}
}


bool Checkups::FEPauxState(StdMutableFst* fsm1, StdMutableFst* fsm2, int qid1, int qid2, string output1, string output2, int iter, int &k, int length) {

	bool result = false;
	bool aux;
	int fails = 0;
	int sucess = 0;

	if (output1 != output2) {
		return false;
	} else if (k != 0) {
		result = true;
	} else if (qid1 != qid2) {
		k = iter;
		result = true;
	}

	if (iter == length) {
		return result;
	}

	for (ArcIterator<StdMutableFst> arcIter1(*fsm1, qid1); !arcIter1.Done(); arcIter1.Next()) {
		const StdArc &arc1 = arcIter1.Value();
		Matcher<StdMutableFst> matcher(*fsm2, MATCH_INPUT);
		matcher.SetState(qid2);
		if (matcher.Find(arc1.ilabel)) {
			for (; !matcher.Done(); matcher.Next()) {
				const StdArc &arc2 = matcher.Value();
				aux = FEPauxState(fsm1, fsm2, arc1.nextstate, arc2.nextstate, output1 + (char)arc1.olabel, output2 + (char)arc2.olabel, iter + 1, k, length);
				if (result && aux) {
					sucess++;
				} else if (result && !aux) {
					fails++;
				} else {
					result = result || aux;
				}
			}
		}
	}
	if (result && sucess > 0 && fails == 0) {
		result = true;
	} else {
		result = false;
	}

	return result;
}

} /* namespace std */
