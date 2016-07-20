/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              multisort.cpp
/***********************************************/

#include"mscls.h"
#include<time.h>

int multisort::execCallbackFunc(const char *seqid1, const char *seqid2, TYPE_INDEX index_of_seq1, TYPE_INDEX index_of_seq2, char* aln1, char* aln2, double dist, int aln_size){
	return( grfptr(seqid1, seqid2, index_of_seq1, index_of_seq2, aln1, aln2, dist, aln_size) );
}

int multisort::getParam(cmlOptions inco){
	preprocess *p = new preprocess(inco);
	p->setParams();
	p->getParams(sq, co, ct);
	setBox();
	sq.outputSeqInfo();
	bx.outputBoxInfo();
	if(!free_vals_automatic && sq.toSSIndex==NULL){
		sq.makeMapOfOrgIdxtoSSIdx();
	}
	cerr<<"get param finish\n";
	is_set_param=true;
	delete(p);
	return(0);
}

int multisort::getParam(int argc, char **argv){
	preprocess *p = new preprocess(argc, argv);
	p->setParams();
	p->getParams(sq, co, ct);
	setBox();
	sq.outputSeqInfo();
	bx.outputBoxInfo();
	if(!free_vals_automatic && sq.toSSIndex==NULL){
		sq.makeMapOfOrgIdxtoSSIdx();
	}
	cerr<<"get param finish\n";
	is_set_param=true;
	delete(p);
	return(0);
}

int multisort::exec(){
//	clock_t start, end;
//	start=clock();

#ifdef CALL_BACK_FUNC
	if(grfptr==NULL){
		cerr<<"ERROR: set callbackfunction.\n";
		exit(1);		
	}
#endif

	if(is_set_param==false){
		cerr<<"ERROR: set parameters.\n";
		exit(1);
	}
	if(co.outputfile){
		openOutputFile();
	}
	if(co.isDevSort && co.startBlc>0){
		co.isOutputIdenticalPairs=false;
	}
	if(co.isOutputIdenticalPairs)
		outputIdenticalPairs();
	if(co.distanceType==EDIT_DISTANCE){
		blockSort_ED();
	}else{
		if(co.isDevSort){
			blockSort(co.startBlc);
		}else{
			blockSort(0);
		}
	}
	if(co.outputfile){
		closeOutputFile();
	}
	//	end=clock();
	cerr<<num_of_comparison<<" pairs were compared\n";
	cerr<<num_of_similar_pairs<<" pairs were found\n";
	//	cerr<<"All the sliding sort procedures are done in "<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";
	freeBucket();
	if(free_vals_automatic)
		freeVals();
	return(0);
}


void multisort::openOutputFile(){
	if(co.outputFileName=="-"){
		outfp = stdout;
	}else{
		outfp = fopen(co.outputFileName.c_str(),"w");
	}
}

void multisort::closeOutputFile(){
	if(co.outputFileName=="-"){
	}else{
		fclose(outfp);
	}
}

void multisort::setBox()
{
	clock_t start, end;
	if(co.distanceType == HAMMING_DISTANCE){
		if(sq.max_seq_length!=sq.min_seq_length){
			cerr<<"ERROR: Input sequences are not same length\n       min_length="<<sq.min_seq_length<<"  max_length="<<sq.max_seq_length<<"\n";
			cerr<<"       Use overlap/inclusion option: -overlap, -inclusion\n";
			if(sq.max_seq_length-sq.min_seq_length < MAX_RECOMMENDED_EDIT_DISTANCE){
				cerr<<"       or use edit distance option: -distance_type E -d "<<sq.max_seq_length-sq.min_seq_length<<"\n";				
			}
			exit(1);
		}
		setBox_nonOV_nonIC_HD();
		// for SSE
		cerr<<"SET BITSTREAM\n";
		start=clock();
		if(sq.max_seq_length<SIZE_OF_REGISTER){
			setBitstream128();
		}else{
			setBitstream();
		}
		end=clock();
		cerr<<" done in ";
		cerr<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";
		dup_chk_thr = (bx.num_of_blocks + bx.key_size)/2; 
	}else{
		if((sq.max_seq_length-co.distance) > sq.min_seq_length){
			cerr<<"ERROR: The length of input sequences shoud range from minimum length of sequence to minimum length of sequence + edit_distance\n       min_length="<<sq.min_seq_length<<"  max_length="<<sq.max_seq_length<<"\n";
			cerr<<"       Use longer edit distance\n";
			exit(1);
		}
		setBox_nonOV_nonIC_ED();
	}

	//set region info
	order = (varLongInt **)malloc(sizeof(varLongInt *)*bx.key_size);
	is_region_top = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
	for(int i=0;i<bx.key_size;i++){
		order[i] = NULL;
		is_region_top[i] = NULL;
	}
	if(co.distanceType==HAMMING_DISTANCE){
		for(int i=0;i<bx.key_size;i++){
			varLongInt *od = new varLongInt(bx.num_of_box);
			varChar *irt = new varChar(bx.num_of_box);
			order[i] = od;
			is_region_top[i] = irt;
		}
	}else{
		is_offset_zero = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
		order_cpy  = (varLongInt **)malloc(sizeof(varLongInt *)*bx.key_size);
		is_offset_zero_cpy = (varChar **)malloc(sizeof(varChar *)*bx.key_size);

		for(int i=0;i<bx.key_size;i++){
			is_offset_zero[i]=NULL;
			order_cpy[i]=NULL;
			is_offset_zero_cpy[i]=NULL;
		}
	}

}

void multisort::resetIsRegionTop(int pos)
{
	if(pos==0){
		memset(is_region_top[0]->val,0x00,sizeof(TYPE_CHAR)*bx.num_of_box);
		is_region_top[0]->val[0]=REGION_TOP;
	}else{
		memcpy(is_region_top[pos]->val, is_region_top[pos-1]->val, sizeof(TYPE_CHAR)*bx.num_of_box); 
	}
}

void multisort::resetOrder(int pos)
{
	if(pos==0){
		if(resetordertmplt==NULL){
			for(TYPE_INDEX i=0;i<bx.num_of_box;i++)
				order[pos]->val[i] =i;
		}else{
			memcpy(order[pos]->val,resetordertmplt->val,sizeof(TYPE_INDEX)*order[pos]->size);
		}
	}else{
		for(TYPE_INDEX i=0;i<bx.num_of_box;i++)
			order[pos]->val[i] =order[pos-1]->val[i];		
	}
}

void multisort::freeIsRegionTop()
{
	if(is_region_top){
		for(int i=0;i<bx.key_size;i++){
			if(is_region_top[i]) delete(is_region_top[i]);
		}
		free(is_region_top);
		is_region_top=NULL;
	}
}

void multisort::freeOrder(){
	if(order){
		for(int i=0;i<bx.key_size;i++) { 
			if(order[i]) delete(order[i]);
		}
		free(order);
		order=NULL;
	}
}

void multisort::freeOrderCpy(){
	if(order_cpy){
		for(int i=0;i<bx.key_size;i++) { 
			if(order_cpy[i]) delete(order_cpy[i]);
		}
		free(order_cpy);
		order_cpy=NULL;
	}
}

void multisort::freeIsOffsetZero(){
	if(is_offset_zero){
		for(int i=0;i<bx.key_size;i++){
			if(is_offset_zero[i]) delete(is_offset_zero[i]);
		}
		free(is_offset_zero);
		is_offset_zero=NULL;
	}
}

void multisort::freeIsOffsetZeroCpy(){
	if(is_offset_zero_cpy){
		for(int i=0;i<bx.key_size;i++){
			if(is_offset_zero_cpy[i]) delete(is_offset_zero_cpy[i]);
		}
		free(is_offset_zero_cpy);
		is_offset_zero_cpy=NULL;
	}
}

void multisort::freeVals(){
	sq.freeSeq();
	bx.freeBox();
	ct.freeTables();
	num_of_similar_pairs=0;
	num_of_comparison=0;

	if(block_mask){
		freeBlockMask();
	}
	if(lower_block_mask){
		freeLowerBlockMask();
	}
	if(bs){
		delete(bs);
		bs=NULL;
	}
	if(is_region_top){
		freeIsRegionTop();
	}
	if(order){
		freeOrder();
	}
	if(order_cpy){
		freeOrderCpy();
	}
	if(is_offset_zero){
		freeIsOffsetZero();
	}
	if(is_offset_zero_cpy){
		freeIsOffsetZeroCpy();
	}
}

void multisort::freeBucket(){
	if(bucket){
		free(bucket[LONG_WORD]);
		free(bucket[SHORT_WORD]);
		free(bucket);
		bucket=NULL;
	}
	if(bucketID){
		free(bucketID[LONG_WORD]);
		free(bucketID[SHORT_WORD]);
		free(bucketID);
		bucketID=NULL;
	}
}


void multisort::outputIdenticalPairs(){
	TYPE_CHARACTER *rev_seq1, *rev_seq2;
	rev_seq1 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq1,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	TYPE_INDEX otmp;
	int num_of_unknown_char;
	int dlen=0;
	TYPE_INDEX seq_start=0;
	TYPE_INDEX i=0;
	for(TYPE_INDEX ptr=0 ;ptr<bx.num_of_box; ptr+=bx.box_set_size){
		if(resetordertmplt==NULL){
			i = ptr/bx.box_set_size;
		}else{
			i = resetordertmplt->val[ptr]/bx.box_set_size;	
		}
		num_of_unknown_char = sq.num_of_unknown_char_matchs[i];
		if(num_of_unknown_char>co.distance){
			goto skipNonPairs;
		}
		otmp=sq.head[i] + bx.box_set_center;
		seq_start=0;
		dlen=0;
		while(sq.nr_seq[otmp]>=ct.num_of_character){
			dlen++;
			otmp++;
		}
		seq_start = otmp;
		otmp = sq.head[i] + bx.box_set_center + sq.max_seq_length;
		while(sq.nr_seq[otmp]>=ct.num_of_character){
			otmp--;
			dlen++;
		}
		dlen--;

#ifdef CALL_BACK_FUNC
		for(otmp=0;otmp<sq.max_seq_length-dlen;otmp++){
			rev_seq1[otmp] = ct.toChar[sq.nr_seq[seq_start+otmp]];
		}
		rev_seq1[sq.max_seq_length-dlen]='\0';
#endif
		TYPE_INDEX ret_id1, ret_id2;
		for(TYPE_INDEX j=0;j<sq.toOrgSeqIndex[i]->size;j++){
			for(TYPE_INDEX k=j+1;k<sq.toOrgSeqIndex[i]->size;k++){
				ret_id1 = sq.toOrgSeqIndex[i]->val[j];
				ret_id2 = sq.toOrgSeqIndex[i]->val[k];
				if(co.isRevComp){
					if(min(ret_id1, ret_id2)%2==1 || ret_id1/2 == ret_id2/2){
						goto revDifChk;
					}
					if(co.isOutputBothStrand==false){
						TYPE_INDEX rev_id = sq.toSSIndex[ (TYPE_INDEX)(ret_id2/2)*2 + abs(ret_id2%2-1) ];
						if( rev_id==i && ret_id2%2==1){
							goto revDifChk;					
						}
					}
				}

				if(co.isMapMode){
					if(sq.org_seq_mapID[ret_id1] != sq.org_seq_mapID[ret_id2]){
						num_of_similar_pairs++;
#ifdef CALL_BACK_FUNC
						if(co.outputaln){
							CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, rev_seq1,rev_seq1,num_of_unknown_char,sq.max_seq_length-dlen);
						}else{
							CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, NULL,NULL,num_of_unknown_char,NULL);
						}
#endif 
#ifdef FILE_OUTPUT_PAIRS
						if(co.outputfile){
							fprintf(outfp,"%s,%s,%d\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),num_of_unknown_char);
							if(co.outputaln){
								int otmp=0;
								for(otmp=0; otmp<sq.max_seq_length-dlen; otmp++){
									fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + seq_start ] ]);
								}
								fprintf(outfp,"\n");
								for(otmp=0; otmp<sq.max_seq_length-dlen; otmp++){
									fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + seq_start ] ]);
								}
								fprintf(outfp,"\n");
							}
						}
#endif
						//cout<<sq.seqName[ sq.toOrgSeqIndex[i]->val[j] ]<<"\t";
						//cout<<sq.seqName[ sq.toOrgSeqIndex[i]->val[k] ]<<"\t"<<"0\n";	
					}
				}else{
						num_of_similar_pairs++;
#ifdef CALL_BACK_FUNC
						if(co.outputaln){
							CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, rev_seq1,rev_seq1,num_of_unknown_char,sq.max_seq_length-dlen);
						}else{
							CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, NULL,NULL,num_of_unknown_char,NULL);
						}
#endif 
#ifdef FILE_OUTPUT_PAIRS
						if(co.outputfile){
							fprintf(outfp,"%s,%s,%d\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),num_of_unknown_char);
							if(co.outputaln){
								for(otmp=0; otmp<sq.max_seq_length-dlen; otmp++){
									fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + seq_start ] ]);
								}
								fprintf(outfp,"\n");
								for(otmp=0; otmp<sq.max_seq_length-dlen; otmp++){
									fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + seq_start ] ]);
								}
								fprintf(outfp,"\n");
							}
						}
#endif
					//cout<<sq.seqName[ sq.toOrgSeqIndex[i]->val[j] ]<<"\t";
					//cout<<sq.seqName[ sq.toOrgSeqIndex[i]->val[k] ]<<"\t"<<"0\n";
				}
revDifChk:;
			}
		}
skipNonPairs:;
	}

	free(rev_seq1);
}

int multisort::setResGetFuncPtr(GETRESFUNC infptr){
	grfptr = infptr;
	return(0);
}

string multisort::showFastaHeader(TYPE_INDEX id){
	string fh = sq.seqName[id];
	return(fh);
}

TYPE_INDEX multisort::showNumOfSeqs(){
	return(sq.num_of_valid_seq);
}

TYPE_LONG multisort::showNumOfPairs(){
	return(num_of_similar_pairs);
}

#ifdef INDEXING_INPUT_ORDER
TYPE_INDEX multisort::showInputOrder(TYPE_INDEX id){
	if(co.exclude_unknown_character){
		return(sq.toInputOrder[id]);
	}else{
		return(id);
	}
}
#endif

string multisort::showSequence(TYPE_INDEX id){
	string ret;
	int cnt=0;
	TYPE_INDEX head = sq.head[ sq.toSSIndex[id] ] + bx.box_set_center;
	char* tmp = (char *)malloc(sizeof(char)*(sq.max_seq_length+1));
	while(sq.nr_seq[head] != ct.overlap_character && cnt<sq.max_seq_length){
		if(sq.nr_seq[head] != ct.lim_wild_card)
			tmp[cnt++] = ct.toChar[ sq.nr_seq[head] ];
		head++;
	}
	tmp[cnt]='\0';
	ret = tmp;
	free(tmp);
	return(ret);
}

TYPE_INDEX multisort::showRevCompID(TYPE_INDEX id){
	if(id%2==0){
		id++;
	}else{
		id--;
	}
	return(id);
}


void multisort::attachSeq(seq &in_sq){
	sq = in_sq;
}

void multisort::detachSeq(){
	sq.mem_aloc_flag = false;
	/**
	sq.head = NULL;
	sq.nr_seq = NULL;
	sq.toOrgSeqIndex = NULL;
	sq.revCompIdx = NULL;
	sq.num_of_unknown_char_matchs = NULL;
	sq.seqName = NULL;
	sq.revCompSeqDist = NULL;
	sq.org_seq_mapID = NULL;
	sq.nr_seq_mapID = NULL;
	sq.toSSIndex = NULL;
	**/
}

void multisort::attachBox(box &in_bx){
	bx = in_bx;
}

void multisort::detachBox(){
	bx.head = NULL;
	bx.block_offset = NULL;

}

void multisort::attachBitstream(bitstream *in_bs,  BIT_BLOCK *in_block_mask, BIT_BLOCK *in_lower_block_mask){
	bs = in_bs;
	block_mask = in_block_mask;
	lower_block_mask = in_lower_block_mask;
}

void multisort::detachBitstream(){
	bs = NULL;
	block_mask = NULL;
	lower_block_mask = NULL;
}

void multisort::initMapRanges(){
	if(sq.nr_seq_mapID){
		cerr<<"WARNING: org_seq_map_ID will be reallocated.\n";
	}
	if(sq.org_seq_mapID){
		cerr<<"WARNING: nr_seq_map_ID will be reallocated.\n";
	}
	sq.org_seq_mapID = (char*)malloc(sizeof(char)*sq.num_of_valid_seq);
	sq.nr_seq_mapID = (char*)malloc(sizeof(char)*sq.num_of_seq);
}

void multisort::setMapRanges(TYPE_INDEX fst_head, TYPE_INDEX fst_size, int map_id){
	if(sq.org_seq_mapID==NULL || sq.nr_seq_mapID==NULL){
		cerr<<"ERROR: seq_map_ID is not allocated.\n";		
	}
	for(TYPE_INDEX i=fst_head; i<fst_head+fst_size; i++){
		sq.nr_seq_mapID[i] = map_id;
		for(int j=0;j<sq.toOrgSeqIndex[i]->size;j++){
			sq.org_seq_mapID[ sq.toOrgSeqIndex[i]->val[j] ] = map_id;
		}
	}
}

void multisort::setSortRange(TYPE_INDEX fst_head, TYPE_INDEX fst_size, TYPE_INDEX snd_head, TYPE_INDEX snd_size){
	TYPE_INDEX fst_start_box, snd_start_box, num_of_box_of_fst;
	//set region info
	order = (varLongInt **)malloc(sizeof(varLongInt *)*bx.key_size);
	is_region_top = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
	for(int i=0;i<bx.key_size;i++){
		order[i] = NULL;
		is_region_top[i] = NULL;
	}
	if(co.distanceType==HAMMING_DISTANCE){
		bx.num_of_box = fst_size+snd_size;
		fst_start_box = fst_head;
		snd_start_box = snd_head;
		num_of_box_of_fst = fst_size;
		for(int i=0;i<bx.key_size;i++){
			varLongInt *od = new varLongInt(bx.num_of_box);
			varChar *irt = new varChar(bx.num_of_box);
			order[i] = od;
			is_region_top[i] = irt;
		}
	}else{
		bx.num_of_box = (fst_size+snd_size)*bx.box_set_size;
		num_of_box_of_fst = fst_size*bx.box_set_size;
		fst_start_box = fst_head*bx.box_set_size;
		snd_start_box = snd_head*bx.box_set_size;
		is_offset_zero = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
		order_cpy  = (varLongInt **)malloc(sizeof(varLongInt *)*bx.key_size);
		is_offset_zero_cpy = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
		for(int i=0;i<bx.key_size;i++){
			is_offset_zero[i]=NULL;
			order_cpy[i]=NULL;
			is_offset_zero_cpy[i]=NULL;
		}
	}
	resetordertmplt = new varLongInt(bx.num_of_box);
	TYPE_INDEX i=0, j=0;
	for(;i<num_of_box_of_fst;i++){
		resetordertmplt->val[i] = fst_start_box + i;
	}
	for(;i<bx.num_of_box;i++,j++){
		resetordertmplt->val[i] = snd_start_box + j;
	}

	if(sq.nr_seq != NULL && ct.toChar != NULL && ct.toInt != NULL){
		is_set_param = true;	
	}

	if(co.distanceType==HAMMING_DISTANCE){
		if(bs == NULL || block_mask == NULL ){
			is_set_param = false;
		}
	}
}

void multisort::setSortRange(TYPE_INDEX fst_head, TYPE_INDEX fst_size){
	TYPE_INDEX start_box;
	//set region info
	order = (varLongInt **)malloc(sizeof(varLongInt *)*bx.key_size);
	is_region_top = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
	for(int i=0;i<bx.key_size;i++){
		order[i] = NULL;
		is_region_top[i] = NULL;
	}
	if(co.distanceType==HAMMING_DISTANCE){
		bx.num_of_box = fst_size;
		start_box = fst_head;
		for(int i=0;i<bx.key_size;i++){
			varLongInt *od = new varLongInt(bx.num_of_box);
			varChar *irt = new varChar(bx.num_of_box);
			order[i] = od;
			is_region_top[i] = irt;
		}
	}else{
		bx.num_of_box = fst_size*bx.box_set_size;
		start_box = fst_head*bx.box_set_size;
		is_offset_zero = (varChar **)malloc(sizeof(varChar *)*bx.key_size);
		order_cpy  = (varLongInt **)malloc(sizeof(varLongInt *)*bx.key_size);
		is_offset_zero_cpy = (varChar **)malloc(sizeof(varChar *)*bx.key_size);

		for(int i=0;i<bx.key_size;i++){
			is_offset_zero[i]=NULL;
			order_cpy[i]=NULL;
			is_offset_zero_cpy[i]=NULL;
		}
	}
	resetordertmplt = new varLongInt(bx.num_of_box);
	for(TYPE_INDEX i=0;i<bx.num_of_box;i++){
		resetordertmplt->val[i] = start_box + i;
	}

	if(sq.nr_seq != NULL && ct.toChar != NULL && ct.toInt != NULL){
		is_set_param = true;	
	}

	if(co.distanceType==HAMMING_DISTANCE){
		if(bs == NULL || block_mask == NULL ){
			is_set_param = false;			
		}
	}
}

void multisort::getCmlOptions(cmlOptions &in_co){
	co = in_co;
}

void multisort::attachCharTable(charTable &in_ct){
	ct = in_ct;
}

void multisort::detachCharTable(){
	ct.toInt = NULL;
	ct.toChar = NULL;
}


