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
#define N 50
#define EX 1
#define INI 10

sem_t sem;

typedef struct {
	Graph* G;
	double* FEP;
	double* GSq;
	double* GPSq;
	int i;
} arguments;

void * mutationState(void * i);
bool isAlready(double fep, double* FEP, int I);
// void * mutationInput(void * i);

int main() {

	//Initialization
	IOHandler* IOH = new IOHandler();
	Checkups* Checker = new Checkups();
	//Operations* Ops = new Operations();

	std::string Ifile = "binary.fst";
	std::string Ofile = "Results.txt";

	std::ofstream OFile;

	OFile.open(Ofile);
	if (!OFile.is_open()) {
		std::cerr << "I can't create the output file." << std::endl;
		return 1;
	}

	OFile << "| Pearson correlation PSq | Spearman correlation PSq | Pearson correlation Sq | Spearman correlation Sq |" << std::endl;

	Graph* G;
	double FEP[N];
	double GSq[N];
	double GPSq[N];
	double PPSqcorr = 0;
	double PSqcorr = 0;
	double SPSqcorr = 0;
	double SSqcorr = 0;
	double aux[2*N];
	pthread_t th[N];
	string input1 = "";
	string input2 = "";

	sem_init(&sem, 0, 1);
	for (int J = INI; J < INI + EX; J++) {

		//for (int I = 0; I < N; I++) {

			Ifile = "Tests/test" + to_string(J+1) + "/binary.fst";

			G = IOH->readGraph(Ifile);

			if (G == NULL) {
				return 1;
			}

			if (!Checker->is_valid(G)) {
				std::cerr << "Not valid graph." << std::endl;
				return 1;
			}

			for (int i = 0; i < N; i++) {
				FEP[i] = 0;
				GSq[i] = 0;
				GPSq[i] = 0;
			}

			//Ops->minimization(G);

			//Ops->ProbAndSqueeziness(G, LEN, GSq, GPSq, I);

			/*
			for (int i = 0; i < REP; i++) {

				input1 = G->getRandInput(LEN);
				input2 = G->getRandInput(LEN);
				while (input1 == input2) {
					input2 = G->getRandInput(LEN);
				}

				if (Checker->has_FEPInput(G, input1, input2, LEN)) {
					FEP[I]++;
				}
			}
			*/

			for (int i = 0; i < N; i++) {
				arguments args;
				args.G = G;
				args.FEP = FEP;
				args.GPSq = GPSq;
				args.GSq = GSq;
				args.i = i;
				pthread_create(&th[i], NULL, mutationState, (void *)&args);
			}

			for (int i = 0; i < N; i++) {
				pthread_join(th[i], NULL);
			}

			//FEP[I] /= REP;

			cout << "test " << to_string(J+1) << endl;
		//}



		PPSqcorr = gsl_stats_correlation(GPSq, 1, FEP, 1, N);
		PSqcorr = gsl_stats_correlation(GSq, 1, FEP, 1, N);
		SPSqcorr = gsl_stats_spearman(GPSq, 1, FEP, 1, N, aux);
		SSqcorr = gsl_stats_spearman(GSq, 1, FEP, 1, N, aux);
		OFile << PPSqcorr << " & " << SPSqcorr << " & " << PSqcorr << " & " << SSqcorr << "\\\\" << std::endl;
		OFile << "\\hline" << std::endl;
		delete G;
	}

	sem_close(&sem);

	OFile.close();

	delete IOH;
	delete Checker;
	//Ops->~Operations();

	return 0;
}


void * mutationState(void * i) {

	Mutations* Mutator = new Mutations();
	Checkups* Checker = new Checkups;
	Operations* Ops = new Operations();
	Graph* GM;
	double fep = 0;
	arguments args = *((arguments *)i);
	GM = Mutator->mutateState(args.G, LEN);
	fep = Checker->PFEPState(args.G, GM, LEN);

	int count = 0;
	sem_wait(&sem);
	while (!Checker->is_validMutation(GM) || ((fep == 0 || fep == 1 /*|| isAlready(args.GSq[args.i], args.GSq, args.i)*/) && count < 500)) {
		GM->~Graph();
		GM = Mutator->mutateState(args.G, LEN);
		fep = Checker->PFEPState(args.G, GM, LEN);
		count++;
	}
	sem_post(&sem);

	Ops->ProbAndSqueeziness(GM, LEN, args.GSq, args.GPSq, args.i);
	args.FEP[args.i] = fep;

	delete Mutator;
	delete Checker;
	delete Ops;
	delete GM;
	pthread_exit(0);
}

bool isAlready(double fep, double* FEP, int I) {

	//sem_wait(&sem);
	for (int i = 0; i < N; i++) {
		if (FEP[i] == fep && i != I) {
			//sem_post(&sem);
			return true;
		}
	}
	//sem_post(&sem);
	return false;
}

/*
void * mutationInput(void * i) {

	Mutations* Mutator = new Mutations();
	Checkups* Checker = new Checkups;
	int I = *((int *)i);
	string input1 = "";
	string input2 = "";

	int count = 500;

	while (count >= 500) {
		input1 = G->getRandInput(LEN);
		count = 0;
		while (!Checker->is_validInput(G, input2, LEN) && count < 500) {
			input2 = Mutator->mutateInput(input1, LEN);
			count++;
		}
	}


	if (Checker->has_FEPInput(G, input1, input2, LEN)) {
		sem_wait(&sem);
		FEP[I]++;
		sem_post(&sem);
	}

	delete Mutator;
	delete Checker;
	pthread_exit(0);
}
*/
