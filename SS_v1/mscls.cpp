/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              mscls.cpp
/***********************************************/

#include "mscls.h"
#include <time.h>


int main(int argc, char **argv)
{
	clock_t start, end;
	start = clock();

	multisort ml;
	//	ml.free_vals_automatic=false;

	// obtain param from command line
	ml.getParam(argc, argv);
	ml.exec();

	end = clock();
	cerr << "All the procedures are done in " << (double)(end - start) / CLOCKS_PER_SEC << " sec.\n";
	return (0);
}

