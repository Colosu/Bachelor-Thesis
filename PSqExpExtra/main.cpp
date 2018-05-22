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
#include <gsl/gsl_blas.h>
#include <fst/fst-decl.h>
#include <fst/fstlib.h>
#include "src/SqueezinessLib.h"

using namespace fst;

#define REP 2
#define LEN 7
#define N 10
#define EX 17
#define INI 483

sem_t sem;

typedef struct {
	Mutations* Mutator;
	Checkups* Checker;
	Operations* Ops;
	Graph* G;
	double* FEP;
	double* GSq;
	double* GPSq;
	int i;
} arguments;

void * mutationState(void * i);
bool isAlready(double fep, double* FEP, int I);
bool lt3(double* FEP);
// void * mutationInput(void * i);

int main(int argc, char * argv[]) {

	//Initialization
	IOHandler* IOH = new IOHandler();
	Mutations* Mutator = new Mutations();
	Checkups* Checker = new Checkups();
	Operations* Ops = new Operations();

	std::string Ifile = "binary.fst";
	std::string FEPfile = "fep.txt";
	std::string Sqfile = "Sq.txt";
	std::string PSqfile = "PSq.txt";
	std::string Ofile = "Results.txt";

	std::ofstream OFile;
	std::ofstream FEPFile;
	std::ofstream PSqFile;
	std::ofstream SqFile;

	OFile.open(Ofile);
	if (!OFile.is_open()) {
		std::cerr << "I can't create the output file." << std::endl;
		return 1;
	}

	OFile << "| Machine | PColl deviation | PSq deviations | Sq deviation |" << std::endl;

	Graph* G;
	double FEP[N];
	double GSq[N];
	double GPSq[N];
	double FEPDesv = 0;
	double SqDesv = 0;
	double PSqDesv = 0;
	//double aux[2*N];
	pthread_t th[N];
	string input1 = "";
	string input2 = "";

	sem_init(&sem, 0, 1);
	for (int J = INI; J < INI + EX; J++) {
		for (int K = 0; K < REP; K++) {
			for (int i = 0; i < N; i++) {
				FEP[i] = 0;
				GSq[i] = 0;
				GPSq[i] = 0;
			}


			Ifile = "./Tests/test" + to_string((J+1)) + "/binary.fst";
			//Ifile = "./Tests/Phone/binary.fst";
			//Ifile = "./War of the Worlds/binary.fst";
			FEPfile = "./FEP/fep" + to_string(J+1) + ".txt";
			Sqfile = "./Sq/Sq" + to_string(J+1) + ".txt";
			PSqfile = "./PSq/PSq" + to_string(J+1) + ".txt";

			G = IOH->readGraph(Ifile);

			if (G == NULL) {
				return 1;
			}

			if (!Checker->is_valid(G)) {
				std::cerr << "Not valid graph." << std::endl;
				return 1;
			}

			try {
				arguments args[N];
				for (int i = 0; i < N; i++) {
					args[i].Mutator = Mutator;
					args[i].Checker = Checker;
					args[i].Ops = Ops;
					args[i].G = G;
					args[i].FEP = FEP;
					args[i].GPSq = GPSq;
					args[i].GSq = GSq;
					args[i].i = i;
				}
				for (int i = 0; i < N; i++) {
					pthread_create(&th[i], NULL, mutationState, (void *)&args[i]);
				}

				for (int i = 0; i < N; i++) {
					pthread_join(th[i], NULL);
				}

			} catch (exception &e) {
				cerr << "Exception: " << e.what() << endl;
			}

			cout << "test " << to_string(J+1) << endl;


			FEPFile.open(FEPfile);
			if (!FEPFile.is_open()) {
				std::cerr << "I can't create the fep output file." << std::endl;
				return 1;
			}
			SqFile.open(Sqfile);
			if (!SqFile.is_open()) {
				std::cerr << "I can't create the Sq output file." << std::endl;
				return 1;
			}
			PSqFile.open(PSqfile);
			if (!PSqFile.is_open()) {
				std::cerr << "I can't create the PSq output file." << std::endl;
				return 1;
			}

			for (int i = 0; i < N; i++){
				FEPFile << FEP[i] << endl;
				SqFile << GSq[i] << endl;
				PSqFile << GPSq[i] << endl;
			}

			FEPFile.close();
			SqFile.close();
			PSqFile.close();

			gsl_vector* FEPv1 = gsl_vector_alloc(N);
			gsl_vector* FEPv2 = gsl_vector_alloc(N);
			double FEPMean = 0;
			for (int i = 0; i < N; i++) {
				gsl_vector_set(FEPv1, i, FEP[i]);
				FEPMean += FEP[i];
			}
			FEPMean /= N;
			for (int i = 0; i < N; i++) {
				gsl_vector_set(FEPv2, i, FEPMean);
			}
			gsl_vector_sub(FEPv1, FEPv2);
			FEPDesv = gsl_blas_dnrm2(FEPv1) / FEPMean;

			gsl_vector* GSqv1 = gsl_vector_alloc(N);
			gsl_vector* GSqv2 = gsl_vector_alloc(N);
			double GSqMean = 0;
			for (int i = 0; i < N; i++) {
				gsl_vector_set(GSqv1, i, GSq[i]);
				GSqMean += GSq[i];
			}
			GSqMean /= N;
			for (int i = 0; i < N; i++) {
				gsl_vector_set(GSqv2, i, GSqMean);
			}
			gsl_vector_sub(GSqv1, GSqv2);
			SqDesv = gsl_blas_dnrm2(GSqv1) / GSqMean;

			gsl_vector* GPSqv1 = gsl_vector_alloc(N);
			gsl_vector* GPSqv2 = gsl_vector_alloc(N);
			double GPSqMean = 0;
			for (int i = 0; i < N; i++) {
				gsl_vector_set(GPSqv1, i, GPSq[i]);
				GPSqMean += GPSq[i];
			}
			GPSqMean /= N;
			for (int i = 0; i < N; i++) {
				gsl_vector_set(GPSqv2, i, GPSqMean);
			}
			gsl_vector_sub(GPSqv1, GPSqv2);
			PSqDesv = gsl_blas_dnrm2(GPSqv1) / GPSqMean;

			OFile << "M_{" << J+1 << "}^" << K+1 << " & " << FEPDesv << " & " << PSqDesv << " & " << SqDesv << "\\\\" << std::endl;
			OFile << "\\hline" << std::endl;

			delete G;

		}
	}

	sem_close(&sem);

	OFile.close();

	delete IOH;
	delete Mutator;
	delete Checker;
	delete Ops;

	return 0;
}


void * mutationState(void * i) {

	Graph* GM;
	arguments args = *((arguments *)i);

	double SqG = 0;
	double PSqG = 0;
	double PCollG = 0;
	GM = args.Mutator->mutateState(args.G, LEN);;

	while (!args.Checker->is_validMutation(GM)) {
		delete GM;
		GM = args.Mutator->mutateState(args.G, LEN);
	}

	args.Ops->ProbAndSqueezinessAndPColl(GM, LEN, SqG, PSqG, PCollG);

	sem_wait(&sem);
	args.GSq[args.i] = SqG;
	args.GPSq[args.i] = PSqG;
	args.FEP[args.i] = PCollG;
	sem_post(&sem);

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

bool lt3(double* FEP) {

	double fep[N];

	for (int i = 0; i < N; i++) {
		fep[i] = -1;
	}

	int I = 0;
	int j;
	for (int i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (FEP[i] == fep[j]) {
				break;
			}
		}
		if (j == N) {
			fep[I] = FEP[i];
			I++;
		}
	}

	if (I < 3)
		return true;
	else
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
