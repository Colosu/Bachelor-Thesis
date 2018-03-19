/*
 * Operations.h
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include <semaphore.h>
#include <fst/fstlib.h>
#include "Graph.h"

namespace fst {

typedef struct {
	StdMutableFst* fsm;
	int qid;
	int iter;
	int length;
	int* inputs;
	std::map<string, int>* mapOtoI;
	string output;
	sem_t* sem;
} args;

void SqueezinessAux(void * argum);

class Operations {
public:
	Operations();
	~Operations();
	void minimization(Graph* g);
	StdMutableFst* product(Graph* g1, Graph* g2);
	void Squeeziness(Graph* g, int length, double Sq[], int I);
	void ProbSqueeziness(Graph* g, int length, double PSq[], int I);
	void ProbAndSqueeziness(Graph* g, int length, double Sq[], double PSq[], int I);

private:

};

} /* namespace std */

#endif /* OPERATIONS_H_ */
