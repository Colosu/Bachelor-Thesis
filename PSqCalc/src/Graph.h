/*
 * Graph.h
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <fst/fstlib.h>

namespace fst {

class Graph {
public:
	Graph();
	Graph(StdMutableFst* trans);
	~Graph();
	StdMutableFst* getTransducer();
	void setTransducer(StdMutableFst* trans);
	string getRandInput(int length);

private:
	StdMutableFst* transducer;

};

} /* namespace std */

#endif /* GRAPH_H_ */
