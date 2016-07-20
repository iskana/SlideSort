/***********************************************/
//     Copyright (c) 2010-2012, Kana Shimizu
//         All rights reserved.
//              mscls.cpp
/***********************************************/

#include "parallelslidesort.h"
#include <time.h>

//#define TEST_PARALLEL

FILE *outfp;
cmlOptions mco;
void openOutputFile();
void closeOutputFile();

int getSimilarPairs(const char* fasta_head1, const char* fasta_head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1,
char* aln2, double dist,int aln_size)
{
    // distがゼロの場合は、0(小数点無し)で表示。それ以外は、小数点以下2桁で表示。
	if(mco.outputfile){
		if (dist<0.005) fprintf(outfp,"%s,%s,%d\n",fasta_head1,fasta_head2,(int)dist);
		else            fprintf(outfp,"%s,%s,%.2f\n",fasta_head1,fasta_head2,dist);
	}
    if(mco.outputaln){
         fprintf(outfp, "%s\n", aln1);
         fprintf(outfp, "%s\n", aln2);
    }
    return(0);
}

int main(int argc, char **argv)
{
//    clock_t start, end;
//    start=clock();

    parallelslidesort ml;
    ml.free_vals_automatic=false;
#ifdef CALL_BACK_FUNC_TEST
    ml.setResGetFuncPtr(test_degree);
#endif
/**
//  /** set param manually
    cmlOptions inco;
    inco.charType=DNA; //char type
    inco.distanceType=EDIT_DISTANCE;
    inco.gap_limit=5;  //gap width
    inco.distance=5; // distance
    inco.key_size = 0; //key size is automatically chose if sat 0
                    //      inco.distanceType=HAMMING_DISTANCE;
    inco.isMapMode = false;
                            //inco.fst_head = 0;
                            //inco.fst_size = 10;
                            //inco.snd_head = 100;
                            //inco.snd_size = 100;
    inco.exclude_unknown_character=true;
    inco.isPartialMode=false;
    inco.inputFileName="test5.txt";
    inco.outputfile=false;
    inco.outputFileName="out.txt";

    ml.getParam(inco);
//      **/

    // obtain param from command line
    ml.getParam(argc, argv);
#ifdef CALL_BACK_FUNC_TEST
    res_mscls = (int *)malloc(sizeof(TYPE_INDEX)*ml.sq.num_of_valid_seq);
    memset(res_mscls,0x00,sizeof(TYPE_INDEX)*ml.sq.num_of_valid_seq);
#endif

    mco.outputFileName = ml.co.outputFileName;
    mco.outputaln = ml.co.outputaln;
	mco.outputfile = ml.co.outputfile;

    ml.co.outputfile = false;
    ml.setResGetFuncPtr(getSimilarPairs);
    openOutputFile();
#ifdef TEST_PARALLEL
    ml.exec2();
#else
	ml.exec();
#endif
    closeOutputFile();

#ifdef CALL_BACK_FUNC_TEST
    for(int i=0;i<ml.sq.num_of_valid_seq;i++){
//      cerr<<ml.sq.seqName[i]<<":"<<res_mscls[i]<<"\n";
    }
    free(res_mscls);
#endif

//    cerr<<ml.num_of_comparison<<" pairs were compared\n";
//    cerr<<ml.num_of_similar_pairs<<" pairs were found\n";

    /**
    for(int i=0;i<ml.sq.num_of_valid_seq;i++){
        cerr<<ml.showFastaHeader(i)<<"\n";
        cerr<<ml.showSequence(i)<<"\n";
    }
    **/
//    end=clock();
//    cerr<<"All the procedures are done in "<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";
    return (0);
}

/**************************************************
 *
 * slidesort結果出力ファイルポインタオープン関数.
 *
 * Argument:
 *
 * Return:
 *
 *************************************************/
void openOutputFile(){
    if(mco.outputFileName=="-"){
        outfp = stdout;
    }else{
        outfp = fopen(mco.outputFileName.c_str(),"w");
    }
}

/**************************************************
 *
 * slidesort結果出力ファイルポインタクローズ関数.
 *
 * Argument:
 *
 * Return:
 *
 *************************************************/
void closeOutputFile(){
    if(mco.outputFileName=="-"){
    }else{
        fclose(outfp);
    }
}

