/*
 * Graph.cpp
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#include "Graph.h"

namespace fst {

Graph::Graph() {
	transducer = NULL;
}

Graph::Graph(StdMutableFst* trans) {
	transducer = NULL;
	setTransducer(trans);
}

Graph::~Graph() {
	if (transducer != NULL) {
//		delete transducer;
	}
	transducer = NULL;
}

StdMutableFst* Graph::getTransducer() {
	return transducer;
}

void Graph::setTransducer(StdMutableFst* trans) {
	if (transducer != NULL) {
		delete transducer;
	}
	transducer = trans;
}

string Graph::getRandInput(int length) {

	string input = "";
	int qid = transducer->Start();
	int num = 0;
	int count = 0;
	for (int i = 0; i < length; i++) {
		ArcIterator<StdMutableFst> arcIter1(*transducer, qid);
		count = 0;
		for (; !arcIter1.Done(); arcIter1.Next()) {
			count++;
		}
		if(count != 0) {
			num = rand() % count;
			arcIter1.Seek(num);
			const StdArc &arc1 = arcIter1.Value();
			input = input + (char)arc1.ilabel;
			qid = arc1.nextstate;
		} else {
			i = length;
		}
	}
	return input;
}

} /* namespace std */
