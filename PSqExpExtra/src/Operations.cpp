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

	SqueezinessAux(argum);

	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		Sq[I] += it->second * std::log2((long double)it->second);
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

	SqueezinessAux(argum);

	long double max = 0;

	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		Sq += it->second * std::log2((long double)it->second);
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

void Operations::ProbAndSqueeziness(Graph* g, int length, double &Sq, double &PSq) {


	if (length <= 0) {
		return;
	}

	int* inputs = new int;
	*inputs = 0;
	double TSq = 0;
	std::map<string, int>* mapOtoI = new std::map<string, int>;
	args* argum = new args;
	argum->fsm = g->getTransducer();
	argum->qid = g->getTransducer()->Start();
	argum->iter = 0;
	argum->length = length;
	argum->inputs = inputs;
	argum->mapOtoI = mapOtoI;
	argum->output = "";

	SqueezinessAux((void *)argum);

	int max = 0;
	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		TSq += it->second * std::log2((long double)it->second);
		if (it->second > max) {
			max = it->second;
		}
	}

	Sq = TSq/(double)(*inputs);

	if (max > 1 && Sq != 0) {
		PSq = Sq/std::log2((long double)max);
	} else {
		PSq = 0;
	}
	delete inputs;
}

void Operations::ProbAndSqueezinessAndPColl(Graph* g, int length, double &Sq, double &PSq, double &PColl) {


	if (length <= 0) {
		return;
	}

	int* inputs = new int;
	*inputs = 0;
	double TSq = 0;
	double TPColl = 0;
	double test = 0;
	sem_t sem;
	sem_init(&sem, 0, 1);
	std::map<string, int>* mapOtoI = new std::map<string, int>;
	args* argum = new args;
	argum->fsm = g->getTransducer();
	argum->qid = g->getTransducer()->Start();
	argum->iter = 0;
	argum->length = length;
	argum->inputs = inputs;
	argum->mapOtoI = mapOtoI;
	argum->output = "";
	argum->sem = &sem;

	SqueezinessAux((void *)argum);


	sem_close(&sem);
	int max = 0;
	for (std::map<string, int>::iterator it = mapOtoI->begin(); it != mapOtoI->end(); it++) {
		TSq += it->second * std::log2((long double)it->second);
		TPColl += (double)(it->second) * (it->second - 1);
		test += it->second;
		if (it->second > max) {
			max = it->second;
		}
	}

	Sq = TSq/((double)*inputs);
	PColl = TPColl/((double)(*inputs) * (*inputs - 1));

	if (max > 1 && Sq != 0) {
		PSq = Sq/std::log2((long double)max);
	} else {
		PSq = 0;
	}
	delete inputs;
	delete mapOtoI;
}


void* SqueezinessAux(void * argum) {

	int N = 5;
	int M = 0;
	args arg = *((args*)argum);
	args* arguments[N];
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
			int count = 0;
			for ( ; !arcIter.Done(); arcIter.Next()) {
				arguments[count] = new args();
				arguments[count]->fsm = arg.fsm;
				arguments[count]->qid = arcIter.Value().nextstate;
				arguments[count]->iter = arg.iter + 1;
				arguments[count]->length = arg.length;
				arguments[count]->inputs = arg.inputs;
				arguments[count]->mapOtoI = arg.mapOtoI;
				arguments[count]->output = arg.output + to_string(arcIter.Value().olabel);
				arguments[count]->sem = arg.sem;
				count++;
			}
			if (arg.iter < M) {
				pthread_t th[N];
				for (int i = 0; i < N; i++){
					pthread_create(&th[i], NULL, SqueezinessAux, (void *)arguments[i]);
				}
				for (int i = 0; i < N; i++) {
					pthread_join(th[i], NULL);
				}
			} else {
				for (int i = 0; i < count; i++){
					SqueezinessAux(arguments[i]);
				}
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
	if(arg.iter > 0 && arg.iter < M + 1) {
		pthread_exit(0);
	}
	return new void*;
}

} /* namespace std */
