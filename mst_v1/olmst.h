/***********************************************/
//     Copyright (c) 2010, Kana Shimizu
//         All rights reserved.
//              onmst.h
/***********************************************/

//#include "..\mscls\\mscls.h"
#ifdef SINGLE_CORE
#include "..\mscls-1.0-fix\\mscls.h"
#else
#include "..\mslidesort\\parallelslidesort.h"
#endif

#include <time.h>
#include <list>

#define MIN_CLS 1

#define TYPE_SIGNED_INDEX int

class Node{	
public:
	char* id;
	TYPE_SIGNED_INDEX parent;
	double edge;
	TYPE_SIGNED_INDEX member; //only valid for root
	TYPE_SIGNED_INDEX rank; //only valid for root
	TYPE_SIGNED_INDEX num_of_sod;
	TYPE_SIGNED_INDEX *sod;
	TYPE_SIGNED_INDEX sod_ptr;
	TYPE_SIGNED_INDEX degree;
	Node(){
		parent = (int)NULL;
		edge = (int)NULL;
		rank = (int)NULL;
	}
};

class MSTree{
public:
	static int mkMST_main(int argc, char** argv);
	static int OutPair(const char* head1, const char* head2,TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size);
	static int mkForest(const char* head1, const char* head2,TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size);
	static int VOID(const char* head1, const char* head2,TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size);
	static int linkChildren();
	static int linkChildrenAsForest();
	static int outputRoot();
	static int outputLargeTreeForWalrus(char *ofname);
	static int outputLargeTreeForText(char *ofname);
	static int outputLargeTreeForSIF(char *ofname);
	static int outputVisualResult(char *ofname);
};

