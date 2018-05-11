/*
 * Mutations.h
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#ifndef MUTATIONS_H_
#define MUTATIONS_H_

#include "Graph.h"

namespace fst {

class Mutations {
public:
	Mutations();
	~Mutations();
	Graph* mutateState(Graph* g, int length);
	string mutateInput(string input1, int length);
};

} /* namespace std */

#endif /* MUTATIONS_H_ */
