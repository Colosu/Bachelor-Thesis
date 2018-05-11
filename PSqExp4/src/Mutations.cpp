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
#include "Checkups.h"

namespace fst {

Mutations::Mutations() {

	std::srand(time(NULL));
}

Mutations::~Mutations() {

}

Graph* Mutations::mutateState(Graph* g, int length) {

	Graph* m = new Graph();
	//Checkups* check = new Checkups();
	StdMutableFst* mut = g->getTransducer()->Copy();
	StdArc arc;

	//do {

		int count = 0;
		int state;

		while (count == 0 || count == 1) {
			state = rand() % (min(length - 1, mut->NumStates() - 1));
			for (ArcIterator<StdMutableFst> arcI(*mut, state); !arcI.Done(); arcI.Next()) {
				count++;
			}
		}

		MutableArcIterator<StdMutableFst>* arcIter = new MutableArcIterator<StdMutableFst>(mut, state);

		arcIter->Seek(rand() % count);

		arc = arcIter->Value();

		int obj = rand() % mut->NumStates();
		while (obj == arc.nextstate || obj == state) {
			obj = rand() % mut->NumStates();
		}

		arc.nextstate = obj;
		arcIter->SetValue(arc);

		//delete arcIter;
	//} while(check->are_equivalent(g, new Graph(mut)));

	m->setTransducer(mut);

	return m;
}


string Mutations::mutateInput(string input1, int length) {

	string input = input1;
	int pos = rand()%length;
	int val = rand()%124;
	while (val < 33 || val == input1[pos]) {
		val = rand()%124;
	}
	input[pos] = val;

	return input;
}

} /* namespace std */
