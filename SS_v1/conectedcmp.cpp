
/***********************************************/
//     Copyright (c) 2010-2012, Kana Shimizu
//         All rights reserved.
//              test_libslidesort.cpp
/***********************************************/

#include<stdio.h>
#include "..\\mscls-1.0-fix\\mscls.h"

class Node{
public:
	Node* next;
	Node* end;
	TYPE_INDEX seqID;
};


Node *cc_list;
Node **cc_table;

int degree(const char* fasta_head1, const char* fasta_head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size)
{
	Node* seq1 = cc_table[seqid1];
	Node* seq2 = cc_table[seqid2];

	cerr<<seqid1<<" - "<<seqid2<<"\n";
	if( seq1 != seq2){
		cerr<<"!!\n";
		seq2->end->next = seq1;
		seq2->end = seq1->end;
		while( seq1->next != seq1 ){
			seq1 = seq1->next;
			cc_table[ seq1->seqID ] = seq2;
		}
		cc_table[ seq1->seqID ] = seq2;
	}

	for(int i=0;i<8;i++){
		cerr<<cc_table[i]->seqID<<"\n";
	}
  return(0);
}

int main(int argc, char **argv)
{
 
  multisort ml;  // create object 
  ml.free_vals_automatic=false;

  //set callback function
  ml.setResGetFuncPtr(degree);

  // if obtain param from command line
  ml.getParam(argc, argv);

  cc_list = (Node *)malloc(sizeof(Node)*ml.showNumOfSeqs());
  cc_table = (Node **)malloc(sizeof(Node*)*ml.showNumOfSeqs());

  for(int i=0;i<ml.showNumOfSeqs();i++){
	  cc_list[i].seqID = i;
	  cc_list[i].next = &cc_list[i];
	  cc_list[i].end = &cc_list[i];
	  cc_table[i] = &cc_list[i];
  }

  ml.exec();

  cerr<<"num of similar pairs = "<<ml.showNumOfPairs()<<"\n";
  cerr<<"num of sequences = "<<ml.showNumOfSeqs()<<"\n";

  for(int i=0;i<ml.showNumOfSeqs();i++){
	cerr<<i<<" : "<<cc_table[i]->seqID<<"\n";
  }

  return (0);
}
