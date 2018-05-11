/*
 * main.cpp
 *
 *  Created on: 8 sept. 2017
 *      Author: colosu
 */

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <gsl/gsl_statistics.h>
using namespace std;

#define REP 200

/*
double PColl(unsigned int d, unsigned int div[]);
double Sq(unsigned int d, unsigned int div[]);
double PSq(unsigned int d, unsigned int div[], double s);
double DRR(unsigned int d, unsigned int div[]);
double Pearson(double pc[100], double ps[100]);
double Spearman(double pc[100], double ps[100]);
unsigned int* setSubdivision(unsigned int d, unsigned int m);
*/

int main() {

	srand(time(NULL));
	string Ifile = "test.txt";
	string Ofile = "PearsonResults.txt";
	string Ofile2 = "SpearmanResults.txt";

	ifstream IFile;
	ofstream OFile;
	ofstream OFile2;

	IFile.open(Ifile);
	if (!IFile.is_open()) {
		cerr << "I can't open the input file." << endl;
		return 1;
	}
	OFile.open(Ofile);
	if (!OFile.is_open()) {
		cerr << "I can't create the output file." << endl;
		return 1;
	}
	OFile2.open(Ofile2);
	if (!OFile2.is_open()) {
		cerr << "I can't create the output file." << endl;
		return 1;
	}
	unsigned int d = 0;
	unsigned int m = 0;
	double pc[REP];
	double ps[REP];
	double s[REP];
	double dr[REP];
	double Pps = 0;
	double Ps = 0;
	double Pdr = 0;
	double Sps = 0;
	double Ss = 0;
	double Sdr = 0;
	unsigned int sum = 0;
	unsigned int max = 0;
	unsigned int count = 0;
	double aux[2*REP];

	OFile << "| Lenght | Maximun size | PSq | Sq | DRR |" << endl;
	OFile2 << "| Lenght | Maximun size | PSq | Sq | DRR |" << endl;
	while (IFile.peek() != EOF) {
		IFile >> d >> m;
		for (int I = 0; I < 2; I++) {
			for (int i = 0; i < REP; i++) {
				sum = 0;
				unsigned int* div = new unsigned int[d];

				/*
				for (unsigned int j = 0; j < d; j++) {
					div[j] = 0;
				}
				*/

				unsigned int n = 0;
				unsigned int j = 0;
				while (sum < d && j < d) {
					n = rand() % m;
					while (n <= 0) {
						n = rand() % m;
					}
					sum += n;
					div[j] = n;
					j++;
				}

				if (sum > d) {
					div[j-1] = div[j-1] - (sum - d);
				}

				div[j] = 0;

				/*
				unsigned int* div;
				div = setSubdivision(d, m);
				pc[i] = PColl(d, div);
				s[i] = Sq(d, div);
				ps[i] = PSq(d, div, s[i]);
				dr[i] = DRR(d, div);
				*/

				pc[i] = 0;
				s[i] = 0;
				ps[i] = 0;
				dr[i] = 0;
				max = 0;
				count = 0;
				int k = 0;
				while (div[k] != 0) {
					pc[i] += ((double)div[k]*(div[k]-1))/(d*(d-1));
					s[i] += (double)div[k]/d * log2(div[k]);
					if (div[k] > max) {
						max = div[k];
					}
					count++;
					k++;
				}
				ps[i] = s[i]/log2(max);
				dr[i] = (double)d/count;
				delete div;
			}
			Pps = gsl_stats_correlation(pc, 1, ps, 1, REP);
			Ps = gsl_stats_correlation(pc, 1, s, 1, REP);
			Pdr = gsl_stats_correlation(pc, 1, dr, 1, REP);
			Sps = gsl_stats_spearman(pc, 1, ps, 1, REP, aux);
			Ss = gsl_stats_spearman(pc, 1, s, 1, REP, aux);
			Sdr = gsl_stats_spearman(pc, 1, dr, 1, REP, aux);
			OFile << d << " & " << m << " & " << Pps << " & " << Ps << " & " << Pdr << "\\\\" << endl;
			OFile << "\\hline" << endl;
			OFile2 << d << " & " << m << " & " << Sps << " & " << Ss << " & " << Sdr << "\\\\" << endl;
			OFile2 << "\\hline" << endl;
		}
		if(IFile.peek() == '\n') {
			IFile.ignore(1);
		}
	}
	IFile.close();
	OFile.close();
	OFile2.close();
	return 0;
}

/*
double PColl(unsigned int d, unsigned int div[]) {
	int i = 0;
	double sum = 0;
	while (div[i] != 0) {
		sum += ((double)div[i]*(div[i]-1))/(d*(d-1));
		i++;
	}
	return sum;
}

double Sq(unsigned int d, unsigned int div[]) {
	int i = 0;
	double sum = 0;
	while (div[i] != 0) {
		sum += pow((double)div[i]/d, 2) * log2(d);
		i++;
	}
	return sum;
}

double PSq(unsigned int d, unsigned int div[], double s) {
	int i = 0;
	unsigned int max = 0;
	while (div[i] != 0) {
		if (div[i] > max) {
			max = div[i];
		}
		i++;
	}
	return s/log2(max);
}

double DRR(unsigned int d, unsigned int div[]) {
	int i = 0;
	unsigned int sum = 0;
	while (div[i] != 0) {
		sum++;
		i++;
	}
	return (double)d/sum;
}

double Pearson(double pc[100], double ps[100]) {

	return gsl_stats_correlation(pc, 1, ps, 1, REP);
}

double Spearman(double pc[100], double ps[100]) {

	double aux[2*REP];
	return gsl_stats_spearman(pc, 1, ps, 1, REP, aux);
}

unsigned int* setSubdivision(unsigned int d, unsigned int m) {
	unsigned int sum = 0;
	unsigned int* div = new unsigned int[d];
	for (unsigned int j = 0; j < d; j++) {
		div[j] = 0;
	}
	unsigned int n = 0;
	unsigned int j = 0;
	while (sum < d && j < d) {
		n = rand() % m;
		while (n <= 0) {
			n = rand() % m;
		}
		sum += n;
		div[j] = n;
		j++;
	}

	if (sum > d) {
		div[j-1] = div[j-1] - (sum - d);
	}

	return div;
}
*/
