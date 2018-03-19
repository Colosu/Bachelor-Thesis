/*
 * Mutations.cpp
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fst/fstlib.h>
#include "Graph.h"
#include "Mutations.h"

namespace fst {

Mutations::Mutations() {

	std::srand(time(NULL));
}

Mutations::~Mutations() {

}

Graph* Mutations::mutate(Graph* g, int length) {

	Graph* m = new Graph;
	StdMutableFst* mut = g->getTransducer()->Copy();
	StdArc arc;
	int count = 0;
	int state;

	while (count == 0) {
		state = rand() % mut->NumStates();
		for (ArcIterator<StdMutableFst> arcI(*mut, state); !arcI.Done(); arcI.Next()) {
			count++;
		}
	}

	MutableArcIterator<StdMutableFst>* arcIter = new MutableArcIterator<StdMutableFst>(mut, state);

	arcIter->Seek(rand() % count);

	arc = arcIter->Value();

	int obj = rand() % mut->NumStates();
	while (obj == arc.nextstate) {
		obj = rand() % mut->NumStates();
	}

	arc.nextstate = obj;
	arcIter->SetValue(arc);

	arcIter->~MutableArcIterator();

	m->setTransducer(mut);

	return m;
}

} /* namespace std */
