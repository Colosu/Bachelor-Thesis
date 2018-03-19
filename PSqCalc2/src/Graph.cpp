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
	transducer->~Fst();
	transducer = NULL;
}

StdMutableFst* Graph::getTransducer() {
	return transducer;
}

void Graph::setTransducer(StdMutableFst* trans) {
	if (transducer != NULL) {
		transducer->~Fst();
	}
	transducer = trans;
}

} /* namespace std */
