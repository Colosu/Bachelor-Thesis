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
	Graph* mutate(Graph* g, int length);
};

} /* namespace std */

#endif /* MUTATIONS_H_ */
