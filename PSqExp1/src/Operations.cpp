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
	a1->~Fst();
	a2->~Fst();
	return NULL; //prod;
}

void Operations::Squeeziness(Graph* g, int length, double Sq[]) {


	if (length <= 0) {
		return;
	}

	int* inputs = new int[length];
	for (int i = 0; i < length; i++) {
		inputs[i] = 0;
		Sq[i] = 0;
	}
	std::map<string, int>* mapOtoI = new std::map<string, int>[length];
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

	for (int i = 0; i < length; i++) {
		for (std::map<string, int>::iterator it = mapOtoI[i].begin(); it != mapOtoI[i].end(); it++) {
			Sq[i] += it->second * std::log2((double)it->second);
		}
		Sq[i] = Sq[i]/(*inputs);
	}
}

void Operations::ProbSqueeziness(Graph* g, int length, double PSq[]) {


	if (length <= 0) {
		return;
	}

	double* Sq = new double[length];
	int* inputs = new int[length];
	for (int i = 0; i < length; i++) {
		inputs[i] = 0;
		Sq[i] = 0;
		PSq[i] = 0;
	}
	std::map<string, int>* mapOtoI = new std::map<string, int>[length];
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
	for (int i = 0; i < length; i++) {
		max = 0;
		for (std::map<string, int>::iterator it = mapOtoI[i].begin(); it != mapOtoI[i].end(); it++) {
			Sq[i] += it->second * std::log2((double)it->second);
			if (it->second > max) {
				max = it->second;
			}
		}
		Sq[i] = Sq[i]/(*inputs);
		if (max > 1) {
			PSq[i] = Sq[i]/std::log2(max);
		} else {
			PSq[i] = Sq[i];
		}
	}
}

void Operations::ProbAndSqueeziness(Graph* g, int length, double Sq[], double PSq[]) {


	if (length <= 0) {
		return;
	}

	int* inputs = new int[length];
	for (int i = 0; i < length; i++) {
		inputs[i] = 0;
		Sq[i] = 0;
		PSq[i] = 0;
	}
	std::map<string, int>* mapOtoI = new std::map<string, int>[length];
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

	SqueezinessAux((void *)argum);

	sem_close(sem);

	long double max = 0;
	for (int i = 0; i < length; i++) {
		max = 0;
		for (std::map<string, int>::iterator it = mapOtoI[i].begin(); it != mapOtoI[i].end(); it++) {
			Sq[i] += it->second * std::log2((double)it->second);
			if (it->second > max) {
				max = it->second;
			}
		}
		Sq[i] = Sq[i]/(*inputs);
		if (max > 1) {
			PSq[i] = Sq[i]/std::log2(max);
		} else {
			PSq[i] = Sq[i];
		}
	}
}


void SqueezinessAux(void * argum) {

	args arg = *((args*)argum);
	args* arguments[10];
	ArcIterator<StdMutableFst> arcIter(*(arg.fsm), arg.qid);
	if (!arcIter.Done()) {
		if (arg.iter > 0) {
			sem_wait(arg.sem);
			arg.inputs[arg.iter-1] = arg.inputs[arg.iter-1] + 1;
			if (arg.mapOtoI[arg.iter-1].find(arg.output) != arg.mapOtoI[arg.iter-1].end()) {
				arg.mapOtoI[arg.iter-1].at(arg.output)++;
			} else {
				arg.mapOtoI[arg.iter-1].emplace(arg.output, 1);
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
			arg.inputs[arg.iter-1] = arg.inputs[arg.iter-1] + 1;
			if (arg.mapOtoI[arg.iter-1].find(arg.output) != arg.mapOtoI[arg.iter-1].end()) {
				arg.mapOtoI[arg.iter-1].at(arg.output)++;
			} else {
				arg.mapOtoI[arg.iter-1].emplace(arg.output, 1);
			}
			sem_post(arg.sem);
		}
		arg.iter = arg.length;
	}
	//delete arg;
	//pthread_exit(0);
}

} /* namespace std */
