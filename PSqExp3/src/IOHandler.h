/*
 * IOHandler.h
 *
 *  Created on: 30 jul. 2017
 *      Author: colosu
 */

#ifndef IOHANDLER_H_
#define IOHANDLER_H_

#include <string>
#include "Graph.h"

namespace fst {

class IOHandler {
public:
	IOHandler();
	~IOHandler();
	Graph* readGraph(std::string file);

private:

};

} /* namespace std */

#endif /* IOHANDLER_H_ */
