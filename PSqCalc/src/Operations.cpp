/*
 * Operations.cpp
 *
 *  Created on: 29 jul. 2017
 *      Author: colosu
 */

#include <string>
#include <cmath>
#include <pthread.h>
#include <semaphore.h>
#include <fst/fstlib.h>
#include "Graph.h"
#include "Operations.h"

namespace fst {

Operations::Operations() {

}

Operations::~Operations() {

}

void Operations::minimization(Graph* g) {
	Minimize(g->getTransducer());
}

StdMutableFst* Operations::product(Graph* g1, Graph* g2) {
	StdMutableFst* a1 = g1->getTransducer()->Copy();
	StdMutableFst* a2 = g2->getTransducer()->Copy();
	//StdMutableFst* prod = Times(a1, a2);
	delete a1;
	delete a2;
	return NULL; //prod;
}

void Operations::Squeeziness(Graph* g, int length, double Sq[], int I) {


	if (length <= 0) {
		return;
	}

	int* inputs = new int;
	*inputs = 0;
	Sq[I] = 0;
	std::map<string, int>* mapOtoI = new std::map<string, int>;
	sem_t* sem = new sem_t;
	args* argum = new args;
	argum->fsm = g->getTransducer();
	argum->qid = g->getTransducer()->Start();
	argum->iter = 0;
	argum->length = length;
	argum->inputs = inputs;
	argum->mapOtoI = mapOtoI;
	argum->output = "";
	argum->sem = sem;

	sem_init(sem, 0, 1);

	SqueezinessAux(argum);

	sem_close(sem);

	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		Sq[I] += it->second * std::log2((double)it->second);
	}
	Sq[I] = Sq[I]/(*inputs);
	delete inputs;
	delete sem;
	delete argum;
}

void Operations::ProbSqueeziness(Graph* g, int length, double PSq[], int I) {


	if (length <= 0) {
		return;
	}

	double Sq = 0;
	int* inputs = new int;
	*inputs = 0;
	PSq[I] = 0;
	std::map<string, int>* mapOtoI = new std::map<string, int>;
	sem_t* sem = new sem_t;
	args* argum = new args;
	argum->fsm = g->getTransducer();
	argum->qid = g->getTransducer()->Start();
	argum->iter = 0;
	argum->length = length;
	argum->inputs = inputs;
	argum->mapOtoI = mapOtoI;
	argum->output = "";
	argum->sem = sem;

	sem_init(sem, 0, 1);

	SqueezinessAux(argum);

	sem_close(sem);

	long double max = 0;

	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		Sq += it->second * std::log2((double)it->second);
		if (it->second > max) {
			max = it->second;
		}
	}
	Sq = Sq/(*inputs);

	if (max > 1) {
		PSq[I] = Sq/std::log2(max);
	} else {
		PSq[I] = 0;
	}
	delete inputs;
	delete sem;
	delete argum;
}

void Operations::ProbAndSqueeziness(Graph* g, int length, double Sq[], double PSq[], int I) {


	if (length <= 0) {
		return;
	}

	int* inputs = new int;
	*inputs = 0;
	Sq[I] = 0;
	PSq[I] = 0;
	std::map<string, int>* mapOtoI = new std::map<string, int>;
	sem_t* sem = new sem_t;
	args* argum = new args;
	argum->fsm = g->getTransducer();
	argum->qid = g->getTransducer()->Start();
	argum->iter = 0;
	argum->length = length;
	argum->inputs = inputs;
	argum->mapOtoI = mapOtoI;
	argum->output = "";
	argum->sem = sem;
	//pthread_t th;

	sem_init(sem, 0, 1);

	//pthread_create(&th, NULL, SqueezinessAux, (void *)argum);
	//pthread_join(th, NULL);

	SqueezinessAux((void *)argum);

	sem_close(sem);

	long double max = 0;
	/*
	long double HIm = 0;
	long double HDom = std::log2((long double)(*inputs));
	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		HIm += (double)it->second/(*inputs) * std::log2((long double)(*inputs));
		if (it->second > max) {
			max = it->second;
		}
	}
	Sq[I] = HDom - HIm;
	*/

	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		Sq[I] += it->second * std::log2((long double)it->second);
		if (it->second > max) {
			max = it->second;
		}
	}
	Sq[I] = Sq[I]/(double)(*inputs);

	/*
	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		Sq[I] += std::pow(it->second / (double)(*inputs), 2) * std::log2((long double)(*inputs));
		if (it->second > max) {
			max = it->second;
		}
	}
	*/

	//if (Sq[I] < 0.0001) {
	//	Sq[I] = 0;
	//}

	if (max > 1 && Sq[I] != 0) {
		PSq[I] = Sq[I]/std::log2(max);
	} else {
		PSq[I] = 0;
	}
	delete inputs;
	delete sem;
}


void SqueezinessAux(void * argum) {

	args arg = *((args*)argum);
	args* arguments[10];
	ArcIterator<StdMutableFst> arcIter(*(arg.fsm), arg.qid);
	if (!arcIter.Done()) {
		if (arg.iter == arg.length) {
			sem_wait(arg.sem);
			*(arg.inputs) = *(arg.inputs) + 1;
			if (arg.mapOtoI->find(arg.output) != arg.mapOtoI->end()) {
				arg.mapOtoI->at(arg.output)++;
			} else {
				arg.mapOtoI->emplace(arg.output, 1);
			}
			sem_post(arg.sem);
		} else {
			for (int i = 0; i < 10; i++) {
				arguments[i] = new args();
				arguments[i]->fsm = arg.fsm;
				arguments[i]->qid = arg.qid;
				arguments[i]->iter = arg.iter + 1;
				arguments[i]->length = arg.length;
				arguments[i]->inputs = arg.inputs;
				arguments[i]->mapOtoI = arg.mapOtoI;
				arguments[i]->output = arg.output;
				arguments[i]->sem = arg.sem;
			}
			/*if (arg->iter < 7) {
				pthread_t th[10] = {0};
				int count = 0;
				for ( ; !arcIter.Done(); arcIter.Next()) {
					arguments[count]->qid = arcIter.Value().nextstate;
					arguments[count]->output = arg->output + to_string(arcIter.Value().olabel);
					pthread_create(&th[count], NULL, SqueezinessAux, (void *)arguments[count]);
					count++;
				}
				for (int i = 0; i < count; i++) {
					if(th[i] != 0) {
						pthread_join(th[i], NULL);
					}
				}
			} else {*/
				int count = 0;
				for ( ; !arcIter.Done(); arcIter.Next()) {
					arguments[count]->qid = arcIter.Value().nextstate;
					arguments[count]->output = arg.output + to_string(arcIter.Value().olabel);
					SqueezinessAux(arguments[count]);
					count++;
			//	}
			}
		}
	} else {
		if (arg.iter > 0) {
			sem_wait(arg.sem);
			*(arg.inputs) = *(arg.inputs) + 1;
			if (arg.mapOtoI->find(arg.output) != arg.mapOtoI->end()) {
				arg.mapOtoI->at(arg.output)++;
			} else {
				arg.mapOtoI->emplace(arg.output, 1);
			}
			sem_post(arg.sem);
		}
	}
	//delete arg;
	//pthread_exit(0);
}

} /* namespace std */
