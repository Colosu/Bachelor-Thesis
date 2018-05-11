/*
 * IOHandler.cpp
 *
 *  Created on: 30 jul. 2017
 *      Author: colosu
 */

#include <iostream>
#include <string>
#include <fst/fstlib.h>
#include "Graph.h"
#include "IOHandler.h"

namespace fst {

IOHandler::IOHandler() {

}

IOHandler::~IOHandler() {

}

Graph* IOHandler::readGraph(std::string file) {

	Graph* g = new Graph();
	g->setTransducer(StdMutableFst::Read(file));
	if (g->getTransducer() == NULL) {
		std::cerr << "Unable to load the automaton." << std::endl;
		return NULL;
	}
	return g;

}

} /* namespace std */
