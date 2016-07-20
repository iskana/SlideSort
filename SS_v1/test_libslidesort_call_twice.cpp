#include<stdio.h>
#include"mscls.h"

int *res_mscls;
int cnt_p;

/*  An example of callback function: int funcname(const char*, const char*, int, int, char*, char*, double)  */
int degree(const char* fasta_head1, const char* fasta_head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size){
  res_mscls[seqid1]++;
  res_mscls[seqid2]++;

  //if dist=0, NULL is returned.
  cnt_p++;
  if(cnt_p%1000==0){
    cerr<<cnt_p<<" pairs have been found\n";
  }
  return(0);
}

int outPairs(const char* fasta_head1, const char* fasta_head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size){
	cerr<<fasta_head2<<"-"<<fasta_head2<<"\n";
	return(0);
}

int main(int argc, char **argv)
{
  cnt_p=0;
  /**
     slidesort -t E -i inputfile -d distance -l
     -l option <-   do not outputfile.
     please see, " slidesort -help" for further info.
   **/

  multisort ml;
  //set callback function
  ml.setResGetFuncPtr(degree);
  // if obtain param from command line
  ml.getParam(argc, argv);
  res_mscls = (int *)malloc(sizeof(int)*ml.sq.seqName.size());
  memset(res_mscls,0x00,sizeof(int)*ml.sq.seqName.size());
  ml.exec();
  for(int i=0;i<ml.sq.seqName.size();i++){
    cerr<<ml.sq.seqName[i]<<":"<<res_mscls[i]<<"\n";
  }
  free(res_mscls);
  
  cerr<<"input a next filename\n";
  char fname[256];
  fgets(fname,sizeof(fname),stdin);
  fname[strlen(fname)-1]='\0';
  
  ml.setResGetFuncPtr(outPairs);
  ml.co.inputFileName=fname;
  ml.getParam(ml.co);
  ml.exec();
  
  return (0);
}

