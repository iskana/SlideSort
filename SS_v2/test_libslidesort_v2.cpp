/***********************************************/
//     Copyright (c) 2010-2012, Kana Shimizu
//         All rights reserved.
//              test_libslidesort.cpp
/***********************************************/

#include<stdio.h>
#include "parallelslidesort.h"

int *res_mscls;
int cnt_p;
int num_of_seq_of_all_input;

/*******
Similar pairs are obtained through callback function.
Callback function should be implemented as follows.

int FunctionName(const char*, const char*, TYPE_INDEX, TYPE_INDEX, char*, char*, double);

Each argument shows following information.

- Fasta header of 1st sequence
- Fasta header of 2nd sequence
- ID of 1st sequence
- ID of 2st sequence
- Alignment of 1st sequence
- Alignment of 2st sequence
- Distance of a pair
- Aligned size of a pair
*********/

int degree(const char* fasta_head1, const char* fasta_head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size)
{
  res_mscls[seqid1]++;
  res_mscls[seqid2]++;

  /**
  //if dist=0, NULL is returned.
  cnt_p++;
  if(cnt_p%1000==0){
    cerr<<cnt_p<<" pairs have been found\n";
  }
  **/
  return(0);
}

int degree2(const char* fasta_head1, const char* fasta_head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size)
{
  res_mscls[seqid1/2]++;
  res_mscls[seqid2/2]++;

  /**
  //if dist=0, NULL is returned.
  cnt_p++;
  if(cnt_p%1000==0){
    cerr<<cnt_p<<" pairs have been found\n";
  }
  **/
  return(0);
}

int main(int argc, char **argv)
{
 
/**	
	In default setting, slidesort output all pairs to STDOUT.
	I recommend to use "-l" option which disables outputting pairs, if you run slidesort on big data.
**/

  cnt_p=0;
  //multisort ml;  // create object
  parallelslidesort ml;  // create object

  ml.free_vals_automatic=false;
/**
  In default setting, original sequence and related information is released after exec().
  If you would like to keep these data on memory,
  
  set free_vals_automatic=false.
  
  The data is kept on memory while the object (in this case "ml") is alive.
  Please take care of memory size when you use this option.

  Those information are accessible by following methods.

  string showFastaHeader(TYPE_INDEX seqid); // show fasta header of a sequence of seqid
  string showSequence(TYPE_INDEX seqid); // show a sequence of seqid
  TYPE_LONG showNumOfPairs(); // show total number of pairs
  TYPE_INDEX showNumOfSeqs(); // show total number of sequences. If you run slidesort with default setting, sequences including unknown character (e.g. n, N, Z, X) are excluded. If using -v, number of sequence is doubled for including reverse complement sequences.
  static TYPE_INDEX showRevCompID(TYPE_INDEX seqid); // If you use -v option, this function shows ID of Reverse complement sequence of seqid.

**/

  // if parameters are obtained from command line
  if(ml.getParam(argc, argv)<0) return (1);

  // -o ml.co.outputfile set to false
  ml.co.outputfile = false;

  //set callback function
  if (ml.co.isRevComp) {
    ml.setResGetFuncPtr(degree2);
  } else {
    ml.setResGetFuncPtr(degree);
  }

  /**
  If you would like to set parameters by the code, 
  please use following options.

  cmlOptions inco;  // An object for setting options.

  inco.distance = 3; // set distance
  inco.charType = DNA; //set char type
  inco.exclude_unknown_character=true; // exclude sequences with unknown character (e.g. n, N, Z, X) or not.
  inco.key_size = 0; // key size is automatically chose if sat 0.

  //set distance type
  inco.distanceType = EDIT_DISTANCE; // set EDIT_DISTANCE for edit distance, or HAMMING_DISTANCE for hamming distance
  // set input
  inco.inputFileName = "input.fst"; // set input filename

  // set output
  inco.outputfile=true; // set true if outputting file. set false if not.
  inco.outputFileName="-"; // set output file name. set "-" if outputting result in stdout.
  inco.outputaln = true; // set true if outputting alignment for each pair, set false if not.

  // set gap cost
  inco.gap_limit=-1;  // gap width
  inco.gap_cost = 1.0; // gap cost
  inco.gap_open = 0; // gap open cost

  //Options for comparing two datasets
  inco.isMapMode = false; //set true if comparing two datasets. set false if not.
  inco.fst_head = 0; // set position of a top sequence of 1st dataset. Positions of the 1st sequence, 2nd sequence, 3rd sequence, ... are referred to as 0, 1, 2, .. .  
  inco.fst_size = 0; //size of first dataset. if 0, automatic detect
  inco.snd_head = 100; // Set position of a top sequence of 2nd dataset.
  inco.snd_size = 0; //set size of second dataset. if 0, automatic detect.
  inco.exclude_unknown_character=true;

  //Options for searching on subset of an input sequences
  inco.isPartialMode=false; // set true if part of datasets are used. set false if not.
  inco.fst_head = 30; // set position of a top sequence.
  inco.fst_size = 0; // set size of dataset.

  //Options for reverse complement search
  inco.isRevComp=false;// set true if searching with reverse complement
  inco.isOutputBothStrand=false;// set true if outputting a same pair twice if both strand make similar pairs (dist(A,B)<=d and dist(A,B')<=d. B' is reverse complement of B).

  ml.getParam(inco);
  **/

  num_of_seq_of_all_input = ml.showNumOfSeqs();
  if (ml.co.isRevComp) {
     num_of_seq_of_all_input /= 2;
  }

  res_mscls = (int *)malloc(sizeof(TYPE_INDEX)*num_of_seq_of_all_input);
  memset(res_mscls,0x00,sizeof(TYPE_INDEX)*num_of_seq_of_all_input);
  if(ml.exec() < 0) return (1);

  cerr<<"num of similar pairs = "<<ml.showNumOfPairs()<<"\n";
  cerr<<"num of sequences = "<<ml.showNumOfSeqs()<<"\n";

  if (ml.co.isRevComp) {
    for(int i=0;i<num_of_seq_of_all_input;i++){
      cout<<ml.showFastaHeader(i*2)<<": "<<res_mscls[i]<<"\n";
    }
  } else {
    for(int i=0;i<num_of_seq_of_all_input;i++){
      cout<<ml.showFastaHeader(i)<<": "<<res_mscls[i]<<"\n";
    }
  }

  free(res_mscls);

  return (0);
}
