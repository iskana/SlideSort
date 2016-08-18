/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              mscls.h
/***********************************************/

#include "common.h"
#include "param.h"

#ifndef VCPP
#include<stdint.h>
#endif

class cmlOptions
{
public:
	int inputFileType;
	int distanceType;
	int distance;
	int allowedgap;
	int gap_limit;
	double gap_cost;
	double gap_open;
	int ptr_of_second_dataset;
	int charType;
	int key_size;
	bool exclude_unknown_character;
	string inputFileName;
	string outputFileName;
	string usage;
	bool outputaln;
	bool outputfile;
	bool isRevComp;
	bool isOutputBothStrand;

	//for mapping mode
	int fst_head;
	int snd_head;
	int fst_size;
	int snd_size;
	bool isMapMode;

	bool isPartialMode;
	bool isInputFromStdin;
	int tmp_file_size;

	//for parallel SS
	bool isSortOrgSeq;
	bool isDevSort; // ソートを分割するか
	int startBlc; //ソート対象とする開始ブロック
	bool isOutputIdenticalPairs; // execが呼ばれるごとにidentpairsが報告されると、重複するため。

	cmlOptions(){
		inputFileType=FORMAT_FASTA;
		distanceType=HAMMING_DISTANCE;
		isInputFromStdin=false;
		usage="mscls";
		tmp_file_size = TMP_FILE_SIZE;
		isSortOrgSeq = true;
		isDevSort = false;
		isOutputIdenticalPairs=true;
		startBlc=0;
	}
	int setCmlOptions(int &argc, char **&argv);
};


class varChar
{
public:
	TYPE_INDEX size;
	TYPE_CHAR *val;
	varChar(TYPE_INDEX n){
		val=(TYPE_CHAR *)calloc(n,sizeof(TYPE_CHAR));
		if(val==NULL){cerr<<"calloc error (varChar)\n";}
		size=n;
	}	
	~varChar(){ 
		if(val)
			free(val); 
	}
};

class varInt
{
public:
	TYPE_INDEX size;
	TYPE_LABEL *val;
	varInt(TYPE_INDEX n){
		val=(TYPE_LABEL *)calloc(n,sizeof(TYPE_LABEL));
		if(val==NULL){cerr<<"calloc error (varInt)\n";}
		size=n;
	}	
	~varInt(){ 
		if(val)
			free(val); 
	}
};

class varLongInt
{
public:
	TYPE_INDEX size;
	TYPE_INDEX *val;
	varLongInt(TYPE_INDEX n){
		val=(TYPE_INDEX *)calloc(n,sizeof(TYPE_INDEX));
		if(val==NULL){cerr<<"calloc error (varLongInt)\n";}
		size=n;
	}	
	~varLongInt(){ 
		if(val)
			free(val); 
	}
};

class charTable
{
public:
	int *toInt;
	char *toChar;
	int num_of_character;
	int unknown_character;
	int overlap_character;
	int lim_wild_card;
	bool isUseWildCard;
	void setCharTable(cmlOptions co);

	charTable(){
		toInt=NULL;
		toChar=NULL;
	}

	//no destruct for copying in getparam
	void freeTables();
};


//データ headは常に数+1
class seq
{
//入力配列から、冗長性を取り除いたもの。
	void freetoOrgSeqIndex();
public:
	TYPE_INDEX num_of_seq_of_all_input; //元データの数　N含有リードも含む
	TYPE_INDEX num_of_seq; //non-redundant num for access of ss internal index
	TYPE_INDEX num_of_valid_seq; // redundant num of access of org order  
	TYPE_LONG seq_length;
	int max_seq_length;
	int min_seq_length;
	TYPE_CHARACTER *nr_seq;   //non redundant sequence
	TYPE_INDEX *head;

	varInt** toOrgSeqIndex; // 配列数は32bit(の本数)を超えない前提。

	TYPE_INDEX *revCompIdx;
	int *num_of_unknown_char_matchs;
	string *seqName;
	bool mem_aloc_flag;

	int *revCompSeqDist;

	//for map
	char *org_seq_mapID;
	char *nr_seq_mapID;

	void freeSeq();
	void outputSeqInfo();

	void makeRevCompIdx();
	void makeRevCompSeqDist();

	// for libSS
	TYPE_INDEX *toSSIndex;
	int makeMapOfOrgIdxtoSSIdx();
#ifdef INDEXING_INPUT_ORDER
	vector<TYPE_INDEX> toInputOrder;
#endif

	seq(){
		seq_length=0;
		num_of_seq=0;
		num_of_valid_seq = 0;
		mem_aloc_flag=false;
		org_seq_mapID=NULL;
		nr_seq_mapID=NULL;
		nr_seq=NULL;
		head=NULL;
		revCompIdx=NULL;
		toSSIndex=NULL;
		toOrgSeqIndex=NULL;
		seqName=NULL;
		revCompSeqDist=NULL;
	}
	~seq(){
		if(mem_aloc_flag){
			if(nr_seq)
				free(nr_seq);
			if(head){
				free(head);
			}
			freetoOrgSeqIndex();
			if(num_of_unknown_char_matchs){
				free(num_of_unknown_char_matchs);
			}
			if(nr_seq_mapID){
				free(nr_seq_mapID);
			}
			if(org_seq_mapID){
				free(org_seq_mapID);
			}
			if(seqName)
				delete(seqName);
			if(revCompIdx)
				free(revCompIdx);
			if(toSSIndex)
				free(toSSIndex);
			if(revCompSeqDist)
				free(revCompSeqDist);
		}
	}
};

class box
{
//入力配列を、同一の長さに切った部分配列。
//入力配列が全て同一の長さで、部分一致を探さない場合はseqと同一。

public:
	TYPE_INDEX num_of_box;
	int box_length;

	// 1 box 当たりの・・
	// boxのサイズがblockの数で割り切れない場合は、block長を調整。long-blockの lshort-block順番
	int num_of_blocks;
	int num_of_long_blocks;
	int short_block_length;
	int long_block_length;
	int key_size;
	TYPE_INDEX *head;
	int *block_offset;

	// for edit distance set of non_indel box and indel boxes
	int box_set_size;
	int box_set_center;

	void outputBoxInfo();
	void freeBox();

	box(){
		head=NULL;
		block_offset=NULL;
	}
	~box(){
		if(head){
			free(head);
		}
		if(block_offset){
			free(block_offset);
		}
	}
};

class bitstream
{
public:
	int num_of_stream;
	TYPE_INDEX stream_length;
	TYPE_INDEX bit_block_length;
	int num_of_regbox_perseq;
	int num_of_bb_per_reg;
	int num_of_bb_per_seq;

	BIT_BLOCK **stream;
	varInt **brmap; //br[block]:reg1,...regn

	bitstream(int size, TYPE_INDEX length);
	void free();
	~bitstream(){
		free();
	}

};

// for call-back function
typedef int (*GETRESFUNC)(const char *seqid1, const char *seqid2, TYPE_INDEX index_of_seq1, TYPE_INDEX index_of_seq2, char* aln1, char* aln2, double dist, int aln_size);

class multisort
{
	//normal SACHICA
	void setBox_nonOV_nonIC_HD();
	void setBox_nonOV_nonIC_ED();

	void resetOrder(int pos);
	void resetIsRegionTop(int pos);

	void findSimilarHDpairs();
	void findSimilarHDpairsHP();
	void findSimilarHDpairsHP128();
	void findSimilarEDpairs();

	// for test only
	void findSimilarHDpairsHP128_for_kmer_test();

	void findSimilarHDRegion(int ptr);
	void findSimilarRegion(int ptr);

	void freeBucket();
	void freeVals();
	void freeIsRegionTop();
	void freeOrder();
	void freeIsOffsetZero();
	void freeOrderCpy();
	void freeIsOffsetZeroCpy();

	//for SSE
	void setBitstream();
	void setBitstream128();
	void freeBlockMask();
	void freeLowerBlockMask();

	int calcHammingDist128(TYPE_INDEX id1, TYPE_INDEX id2);
	int calcHammingDist(TYPE_INDEX id1, TYPE_INDEX id2);


	FILE *outfp;
	//for Sorting
	int char_size; //setBox
	int word_size[2]; //setBox
	int num_of_bucket[2]; //blocksort_ED
	TYPE_INDEX **bucket; //blocksort_ED
	int **bucketID; //blocksort_ED
	int num_of_words[2];//blocksort_ED
	int l_word_size[2];//blocksort_ED

	bool is_set_param;

	varLongInt **order; // for edit distance
	varChar **is_offset_zero;//setBox
	varChar **is_region_top;//setBox
	varChar **is_offset_zero_cpy; //メモリを食う。（計算量削減のため。） //setBox
	varLongInt **order_cpy; //メモリを食う。（計算量削減のため。） //setBox
	vector<int> key_stack; 

	//main routine
	void setBox();
	int blockSort(int next);
	int blockSort_ED(void);
	void openOutputFile();
	void closeOutputFile();
	void outputIdenticalPairs();
	double calcEditDistance(TYPE_INDEX seqid1, TYPE_INDEX seqid2, TYPE_CHARACTER *seq1,  TYPE_CHARACTER *seq2);
	int calcHDistanceFromED(TYPE_INDEX seqid1, TYPE_INDEX seqid2, TYPE_CHARACTER *seq1,  TYPE_CHARACTER *seq2);

	//virtual function to wrap call-back function
	virtual int execCallbackFunc(const char *seqid1, const char *seqid2, TYPE_INDEX index_of_seq1, TYPE_INDEX index_of_seq2, char* aln1, char* aln2, double dist, int aln_size);
	/**
	// implimating following function
	int multisort::execCallbackFunc(const char *seqid1, const char *seqid2, TYPE_INDEX index_of_seq1, TYPE_INDEX index_of_seq2, char* aln1, char* aln2, double dist, int aln_size){
	 
	}
	**/

public:
	seq sq;
	box bx;
	cmlOptions co;
	charTable ct;

	//for SSE
	bitstream *bs; //setBox
	BIT_BLOCK *block_mask; //setBox
	BIT_BLOCK *lower_block_mask; //setBox
	int dup_chk_thr;
	
	TYPE_LONG num_of_similar_pairs;
	TYPE_LONG num_of_comparison;
	TYPE_LONG sum_of_region;

	bool free_vals_automatic;

	int setResGetFuncPtr(GETRESFUNC fptr);
	int getParam(int argc, char **argv);
	int getParam(cmlOptions co);
	int exec();

	//func pointer
	GETRESFUNC grfptr;

	// for lib
	virtual TYPE_LONG showNumOfPairs();
	virtual TYPE_INDEX showNumOfSeqs();
	virtual string showSequence(TYPE_INDEX id);
	virtual string showFastaHeader(TYPE_INDEX id);
	static TYPE_INDEX showRevCompID(TYPE_INDEX id);
#ifdef INDEXING_INPUT_ORDER
	virtual TYPE_INDEX showInputOrder(TYPE_INDEX id);
#endif
	//	getResult(char* id1, char* id2, char *out_seq1, char *out_seq2,double distance,int size_of_alignment);
//	virtual int getResult(int a, int b);

	//for parallel SS
	int debug();
	varLongInt *resetordertmplt;
	void getCmlOptions(cmlOptions &org);
	void attachCharTable(charTable &in_ct);
	void detachCharTable();
	void attachSeq(seq &in_sq);
	void detachSeq();
	void attachBox(box &in_bx);
	void detachBox();
	void attachBitstream(bitstream *in_bs, BIT_BLOCK *in_block_mask, BIT_BLOCK *in_lower_block_mask);
	void detachBitstream();
	void setSortRange(TYPE_INDEX fst_head, TYPE_INDEX fst_size);
	void setSortRange(TYPE_INDEX fst_head, TYPE_INDEX fst_size, TYPE_INDEX snd_head, TYPE_INDEX snd_size);
	void initMapRanges();
	void setMapRanges(TYPE_INDEX fst_head, TYPE_INDEX fst_size, int map_id);

	multisort(){
		is_set_param=false;
		free_vals_automatic=true;

		bs = NULL;
		block_mask=NULL;
		lower_block_mask=NULL;

		is_region_top=NULL;
		is_offset_zero=NULL;
		is_offset_zero_cpy=NULL;
		order=NULL;
		order_cpy=NULL;

		num_of_similar_pairs=0;
		num_of_comparison=0;

		bucket=NULL;
		bucketID=NULL;

		resetordertmplt=NULL;

		grfptr=NULL;
	}
	~multisort(){
		if(free_vals_automatic==false)
			freeVals();
	}
};


// ファイルの読み込み、冗長性の除去、Nの扱いなど。
class preprocess
{
	cmlOptions co;
	TYPE_CHARACTER *org_seq;//exclude_unknown_characterの場合は、有効な配列のみ。
	TYPE_INDEX *sort, *work, *head;
	TYPE_INDEX num_of_seq_of_all_input; //元データの数　N含有リードも含む
	TYPE_INDEX num_of_seq;
	TYPE_LONG seq_length; // total seq length
	int max_seq_length;
	int min_seq_length;
	charTable ct;

	int readOrgFastaSeq();
	int readOrgFastqSeq();
	int readOrgFastaSeqViaStdin();
	int readOrgFastqSeqViaStdin();

	void sortOrgSeq();
	void getSeq(seq &s);
	void cpOrgseq(seq &s);
	void readRevComp();

	// 論文用。配る時はムダ。
	void findKmarPair();

public:
	vector<string> seqName;
#ifdef INDEXING_INPUT_ORDER
	vector<TYPE_INDEX> toInputOrder;
#endif
	preprocess(int argc, char **argv);
	preprocess(cmlOptions inco);
	void setParams();
	void getParams(seq &s, cmlOptions &co, charTable &ct);
	void mkNonredundantData(char *outputfile);

	//for parallel SS
	void keepOrgSeq(); // sortorgseqの代わり。ソートしない。

	~preprocess(){
		if(org_seq){
			free(org_seq);
		}
		if(sort){
			free(sort);
		}
		if(work){
			free(work);
		}
		if(head){
			free(head);
		}
	}
};

