/*
 * main.cpp
 *
 *  Created on: 19 sept. 2017
 *      Author: colosu
 */

#include <iostream>
#include <fstream>
#include <string>
#include <pthread.h>
#include <semaphore.h>
#include <gsl/gsl_statistics.h>
#include <fst/fst-decl.h>
#include <fst/fstlib.h>
#include "src/SqueezinessLib.h"
using namespace fst;

#define REP 100
#define LEN 25

Graph* G;
double FEP[LEN];
sem_t sem;

void * mutation(void * i);

int main() {

	//Initialization
	IOHandler* IOH = new IOHandler();
	Checkups* Checker = new Checkups();
	Operations* Ops = new Operations();

	std::string Ifile = "binary.fst";
	std::string Ofile = "Results.txt";

	std::ofstream OFile;

	OFile.open(Ofile);
	if (!OFile.is_open()) {
		std::cerr << "I can't create the output file." << std::endl;
		return 1;
	}

	OFile << "| Number of States | Pearson correlation PSq | Spearman correlation PSq | Pearson correlation Sq | Spearman correlation Sq |" << std::endl;

	double GSq[LEN];
	double GPSq[LEN];
	double PPSqcorr = 0;
	double PSqcorr = 0;
	double SPSqcorr = 0;
	double SSqcorr = 0;
	double aux[2*LEN];
	pthread_t th[REP];

	for (int J = 0; J < 50; J++) {
		sem_init(&sem, 0, 1);

		Ifile = "Tests/test" + to_string(J+1) + "/binary.fst";

		G = IOH->readGraph(Ifile);

		if (G == NULL) {
			return 1;
		}

		if (!Checker->is_valid(G)) {
			std::cerr << "Not valid graph." << std::endl;
			return 1;
		}

		Ops->minimization(G);

		Ops->ProbAndSqueeziness(G, LEN, GSq, GPSq);

		for (int j = 0; j < 2; j++) {

			for (int i = 0; i < LEN; i++) {
				FEP[i] = 0;
			}

			for (int i = 0; i < REP; i++) {
				pthread_create(&th[i], NULL, mutation, Checker);
			}

			for (int i = 0; i < REP; i++) {
				pthread_join(th[i], NULL);
			}

			PPSqcorr = gsl_stats_correlation(GPSq, 1, FEP, 1, LEN);
			PSqcorr = gsl_stats_correlation(GSq, 1, FEP, 1, LEN);
			SPSqcorr = gsl_stats_spearman(GPSq, 1, FEP, 1, LEN, aux);
			SSqcorr = gsl_stats_spearman(GSq, 1, FEP, 1, LEN, aux);
			OFile << G->getTransducer()->NumStates() << " & " << PPSqcorr << " & " << SPSqcorr << " & " << PSqcorr << " & " << SSqcorr << "\\\\" << std::endl;
			OFile << "\\hline" << std::endl;
			cout << "test " << to_string(J+1) << endl;
		}
		sem_close(&sem);
	}

	OFile.close();

	IOH->~IOHandler();
	Checker->~Checkups();
	Ops->~Operations();

	return 0;
}

void * mutation(void * i) {

	Mutations* Mutator = new Mutations();
	Checkups* Checker = (Checkups *) i;
	Graph* GM;
	GM = Mutator->mutate(G, LEN);

	int k = 0;
	int count = 0;
	while ((!Checker->is_validMutation(GM) || !Checker->has_FEP(G, GM, LEN, k)) && count < 500) {
		GM = Mutator->mutate(G, LEN);
		count++;
	}

	sem_wait(&sem);
	if (k > 0 && k <= LEN) {
		FEP[k]++;
	} else if (k != 0) {
		std::cerr << "k outside boundaries. Something went wrong." << std::endl;
	}
	sem_post(&sem);

	Mutator->~Mutations();
	pthread_exit(0);
}
