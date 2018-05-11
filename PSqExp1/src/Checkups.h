/*
 * Checkups.h
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#ifndef CHECKUPS_H_
#define CHECKUPS_H_

#include <fst/fstlib.h>
#include <string>
#include "Graph.h"

namespace fst {

class Checkups {
public:
	Checkups();
	~Checkups();
	bool are_equivalent(Graph* g1, Graph* g2);
	bool is_valid(Graph* g);
	bool is_validMutation(Graph* g);
	bool has_FEP(Graph* g1, Graph* g2, int length, int &k);

private:
	bool FEPaux(StdMutableFst* fsm1, StdMutableFst* fsm2, int qid1, int qid2, std::string output1, std::string output2, int iter, int &k, int length);
};

} /* namespace std */

#endif /* CHECKUPS_H_ */
