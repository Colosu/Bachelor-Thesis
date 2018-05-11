/*
 * main.cpp
 *
 *  Created on: 25 sept. 2017
 *      Author: colosu
 */


#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <fst/fstlib.h>
using namespace std;

#define NREP 500
#define MAX_STATES 25
#define MIN_STATES 25
#define MAX_TRANSITIONS 5 //PER STATE
#define MIN_TRANSITIONS 5 //PER STATE
#define NCHARS 5
#define MIN_CHAR 65


int main() {

	srand(time(NULL));

	string command;
	string ofile;
	ofstream* OFile;
	int nstates = 0;
	int ntransitions = 0;
	int nextstate = 0;
	int input = 0;
	int output = 0;
	fst::StdMutableFst* fst = NULL;
	int* inputs = new int [NCHARS];

	for (int I = 0; I < NREP; I++) {

		command = "mkdir test" + to_string(I+1);
		system(command.c_str());
		command = "cp ascii.syms test" + to_string(I+1) + "/isyms.syms";
		system(command.c_str());
		command = "cp ascii.syms test" + to_string(I+1) + "/osyms.syms";
		system(command.c_str());

		ofile = "test" + to_string(I+1) + "/graph.fst";

		while (fst == NULL || fst->Properties(fst::kIDeterministic, true) != fst::kIDeterministic || fst->Properties(fst::kAccessible, true) != fst::kAccessible || fst->Properties(fst::kCoAccessible, true) != fst::kCoAccessible) {

			OFile = new ofstream();
			OFile->open(ofile, ofstream::trunc);
			if (!OFile->is_open() || !OFile->good()) {
				std::cerr << "I can't create the output file." << std::endl;
				return 1;
			}

			nstates = 0;
			nstates = MIN_STATES + (rand() % (MAX_STATES - MIN_STATES + 1));

			for (int i = 0; i < nstates; i++) {

				ntransitions = 0;
				if (i != nstates - 1) {
					ntransitions = MIN_TRANSITIONS + (rand() % (MAX_TRANSITIONS - MIN_TRANSITIONS + 1));
				}

				if (ntransitions == 0) {
					*OFile << to_string(i) << endl;
				} else {

					for (int j = 0; j < NCHARS; j++) {
						inputs[j] = 0;
					}

					input = 0;
					output = 0;
					while (input == 0 || inputs[input - MIN_CHAR] != 0) {
						input = MIN_CHAR + (rand() % NCHARS);
					}
					inputs[input - MIN_CHAR] = 1;
					output = MIN_CHAR + (rand() % NCHARS);

					*OFile << to_string(i) << " " << to_string(i + 1) << " " << (char)input << " " << (char)output << endl;

					for (int j = 1; j < ntransitions; j++) {

						nextstate = 0;
						input = 0;
						output = 0;
						while (nextstate == i || nextstate == 0) {
							nextstate = rand() % nstates;
						}
						while (input == 0 || inputs[input - MIN_CHAR] != 0) {
							input = MIN_CHAR + (rand() % NCHARS);
						}
						inputs[input - MIN_CHAR] = 1;
						output = MIN_CHAR + (rand() % NCHARS);
						*OFile << to_string(i) << " " << to_string(nextstate) << " " << (char)input << " " << (char)output << endl;
					}
				}
			}

			OFile->close();
			OFile->~basic_ofstream();

			command = "fstcompile --isymbols=test" + to_string(I+1) + "/isyms.syms --osymbols=test" + to_string(I+1) + "/osyms.syms --keep_isymbols --keep_osymbols test" + to_string(I+1) + "/graph.fst test" + to_string(I+1) + "/binary.fst";
			system(command.c_str());

			command = "test" + to_string(I+1) + "/binary.fst";
			fst = fst->Read(command);

			command = "fstdraw --isymbols=test" + to_string(I+1) + "/isyms.syms --osymbols=test" + to_string(I+1) + "/osyms.syms -portrait test" + to_string(I+1) + "/binary.fst | dot -Tps >test" + to_string(I+1) + "/binary.ps";
			system(command.c_str());
		}
		fst->~Fst();
		fst = NULL;
	}
	delete inputs;

	return 0;
}
