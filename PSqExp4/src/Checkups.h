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
	bool is_validInput(Graph* g, string input, int length);
	double PFEPState(Graph* g1, Graph* g2, int length);
	bool has_FEPState(Graph* g1, Graph* g2, int length, int &k);
	bool has_FEPInput(Graph* g, string input1, string input2, int length);

private:
	void PFEPauxState(StdMutableFst* fsm1, StdMutableFst* fsm2, int qid1, int qid2, std::string output1, std::string output2, int iter, int &sucess, int &fails, int length, bool aux);
	bool FEPauxState(StdMutableFst* fsm1, StdMutableFst* fsm2, int qid1, int qid2, std::string output1, std::string output2, int iter, int &k, int length);
};

} /* namespace std */

#endif /* CHECKUPS_H_ */
