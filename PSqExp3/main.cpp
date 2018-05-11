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

#define REP 2
#define LEN 20
#define N 10
#define EX 50
#define INI 0

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

	OFile << "| Machine | #states | Pearson correlation PSq | Spearman correlation PSq | Pearson correlation Sq | Spearman correlation Sq |" << std::endl;

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
		for (int I = 0; I < REP; I++) {

			Ifile = "./Tests/test" + to_string(J+1) + "/binary.fst";
			//Ifile = "./Tests/Phone/binary.fst";
			//Ifile = "./War of the Worlds/binary.fst";
			FEPfile = "./Tests/test" + to_string(J+1) + "/fep" + to_string(I) + ".txt";
			Sqfile = "./Tests/test" + to_string(J+1) + "/Sq" + to_string(I) + ".txt";
			PSqfile = "./Tests/test" + to_string(J+1) + "/PSq" + to_string(I) + ".txt";

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

			try {
				while (lt3(FEP)) {
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

			PPSqcorr = gsl_stats_correlation(GPSq, 1, FEP, 1, N);
			PSqcorr = gsl_stats_correlation(GSq, 1, FEP, 1, N);
			SPSqcorr = gsl_stats_spearman(GPSq, 1, FEP, 1, N, aux);
			SSqcorr = gsl_stats_spearman(GSq, 1, FEP, 1, N, aux);
			OFile << "M" /* _{ << J+1 << }*/ " & " << G->getTransducer()->NumStates() << " & " << PPSqcorr << " & " << SPSqcorr << " & " << PSqcorr << " & " << SSqcorr << "\\\\" << std::endl;
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
	double Sq = 0;
	double PSq = 0;
	double fep = 0;
	GM = args.Mutator->mutateState(args.G, LEN);

	while (!args.Checker->is_validMutation(GM)) {
		delete GM;
		GM = args.Mutator->mutateState(args.G, LEN);
	}

	args.Ops->ProbAndSqueeziness(GM, LEN, Sq, PSq);
	fep = args.Checker->PFEPState(args.G, GM, LEN);
	sem_wait(&sem);
	args.FEP[args.i] = fep;
	args.GSq[args.i] = Sq;
	args.GPSq[args.i] = PSq;
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

	if (I < 3) {
		//cout << "lt3" << endl;
		return true;
	} else {
		return false;
	}
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
