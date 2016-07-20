/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              files.cpp
/***********************************************/

#include"mscls.h"
#include<map>
#include<time.h>

preprocess::preprocess(int argc, char **argv)
{
	org_seq=NULL;
	sort=NULL;
	work=NULL;
	head=NULL;
	co.setCmlOptions(argc, argv);
}

preprocess::preprocess(cmlOptions inco)
{
	org_seq=NULL;
	sort=NULL;
	work=NULL;
	head=NULL;
	co=inco;
}

void preprocess::setParams(){
	clock_t start, end;
	org_seq=NULL;sort=NULL;work=NULL;head=NULL;
	ct.setCharTable(co);
	cerr<<"read file ...";
	start=clock();
	if(co.isInputFromStdin){
		FILE* fp = stdin;
		char ch = getc(fp);
		if(co.inputFileType==FORMAT_FASTA){
			switch(ch){
			case '@':
				cerr<<"ERROR: Format of the input file is incorrect.\n";
				cerr<<"WARNING: RUN SLIDESORT WITH FASTQ\n";
				readOrgFastqSeqViaStdin();
				break;
			case '>':
				readOrgFastaSeqViaStdin();
				break;
			default:
				cerr<<"ERROR: Format of the input file is incorrect.\n";
				exit(1);
				break;
			}
		}else if(co.inputFileType==FORMAT_FASTQ){		
			switch(ch){
			case '@':
				readOrgFastqSeqViaStdin();
				break;
			case '>':
				cerr<<"ERROR: Format of the input file is incorrect.\n";
				cerr<<"WARNING: RUN SLIDESORT WITH FASTA\n";
				readOrgFastaSeqViaStdin();
				break;
			default:
				cerr<<"ERROR: Format of the input file is incorrect.\n";
				exit(1);
				break;
			}
		}else{
			cerr<<"ERROR: Set input file type.\n";
			exit(1);
		}		
	}else{
		if(co.inputFileType==FORMAT_FASTA){
			if(readOrgFastaSeq()==FORMAT_FASTA_TO_FASTQ)
				readOrgFastqSeq();			

		}else if(co.inputFileType==FORMAT_FASTQ){
			if(readOrgFastqSeq()==FORMAT_FASTQ_TO_FASTA)
				readOrgFastaSeq();			
		}else{
			cerr<<"ERROR: Set input file type.\n";
			exit(1);
		}	
	}



	if(co.isRevComp){
		readRevComp();
	}
	end=clock();
	cerr<<" done in ";
	cerr<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";
	// set search width.
#ifndef INTERNAL_GAP_OPEN
	if(max_seq_length == max_seq_length){
		co.gap_limit = (co.distance - 2*co.gap_open)/co.gap_cost;
	}else{
		co.gap_limit = (co.distance - co.gap_open)/co.gap_cost;
	}
#endif
	co.allowedgap = co.gap_limit/2;
	ct.isUseWildCard=false;
	if(max_seq_length-min_seq_length != 0){
		co.allowedgap += co.gap_limit%2;
		ct.isUseWildCard=true;
		//debg
		//			cerr<<"-----------"<<ct.num_of_character<<"\n";
	}else if(co.allowedgap==0){
		cerr<<"Switch to Hamming Distance Mode.\n";
		co.distanceType=HAMMING_DISTANCE;
	}
	//		cerr<<"slide-w:"<<co.allowedgap<<"\n";
	if(co.charType==INPUT_INT){
		if(ct.num_of_character>MIN_CHARACTER_FOR_HPMODE){
			co.distanceType = EDIT_DISTANCE;
		}
		co.outputaln=false;
	}

	if(co.key_size<=0){
		cerr<<"warning: key size is automatically determined.\n";
		if(co.distanceType==EDIT_DISTANCE){
			co.key_size = max_seq_length/10 - co.distance;
			if(co.key_size<=1){
				co.key_size = 2;
			}
		}else{
			if(num_of_seq < HAMMING_KEY_DEC_SEQ_NUM){
				co.key_size = 2;
			}else{
				co.key_size = 3;			
			}
		}
	}
	if(co.isSortOrgSeq){
		cerr<<"remove redundancy ... ";
		start=clock();
		sortOrgSeq();
		cerr<<"done in ";
		end=clock();
		cerr<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";
	}else{
		cerr<<"cpy orgseq ... ";
		keepOrgSeq();
		cerr<<"done.\n";				
	}
	/* RONBUN */
	//	findKmarPair();
	//	exit(1);
}

void preprocess::getParams(seq &s, cmlOptions &in_co, charTable &in_ct)
{
	in_co = co;
	in_ct = ct;
	getSeq(s);
	if(co.isRevComp){
		s.makeRevCompIdx();
		if(co.isOutputBothStrand==false){
			s.makeRevCompSeqDist();
			s.makeMapOfOrgIdxtoSSIdx();
		}
	}
	//   cpOrgseq(s);
	/**
	for(int i=0;i<s.seq_length;i++){
	cerr<<s.nr_seq[i];
	}
	**/
}


int preprocess::readOrgFastaSeq()
{
	FILE *fp;

	fp = fopen(co.inputFileName.c_str(),"r");
	if(fp==NULL){
		fprintf(stderr,"ERROR: Could not find %s\n",co.inputFileName.c_str());
		exit(1);
	}

	char ch;
	max_seq_length=0;
	min_seq_length=MIN_SEQ_LENGTH_DEFAULT;
	seq_length=0;
	num_of_seq=0;
	num_of_seq_of_all_input=0;
	TYPE_INDEX ltmp=0;
	bool lflag=false;

	TYPE_INDEX ttlSeq=0;
	int mflag=false;
	TYPE_INDEX fst_size_tmp=0, snd_size_tmp=0;

	map<int,int> char_map;
	int map_int=0, map_int_pos=0;
	char map_int_str[SIZE_OF_TYPE_CHARACTER];
	int num_of_unknown_chars=0;

	// check if file type is collect
	if((ch=getc(fp))!='>'){
		cerr<<"ERROR: Format of the input file is incorrect.\n";
		fclose(fp);
		if(ch=='@'){
			cerr<<"WARNING: RUN SLIDESORT WITH FASTQ\n";			
			return(FORMAT_FASTA_TO_FASTQ);
		}
		exit(1);
	}
	fseek(fp,0,SEEK_SET);

	while((ch=getc(fp))!=EOF){
		switch(ch){
		case '>':
			num_of_seq_of_all_input++;
			if(!co.exclude_unknown_character){
				if(num_of_unknown_chars<=co.distance)
					lflag=false;
			}		
			if(lflag){
#ifdef INDEXING_INPUT_ORDER
				toInputOrder.pop_back();
#endif
				num_of_seq--;
				seq_length -= ltmp;
			}else{
				if(max_seq_length<ltmp) max_seq_length = ltmp;
				if(min_seq_length>ltmp && ltmp>0) min_seq_length = ltmp;
			}
			lflag=false;
			ltmp=0;
			num_of_unknown_chars = 0;

			//Clear the line with <
			while((ch=getc(fp))!='\n');

			if(co.isMapMode){
				ttlSeq++;
				if( ttlSeq > co.snd_head && ((co.snd_size==0) || ( ttlSeq <= co.snd_size + co.snd_head ) ) ){
					mflag=false;
				}else if( ttlSeq > co.fst_head && ((co.fst_size==0) || ( ttlSeq <= co.fst_size + co.fst_head ) ) ){
					mflag=false;
				}else{
					mflag=true;
					break;
				}
			}

			if(co.isPartialMode){
				ttlSeq++;
				if( ttlSeq > co.fst_head && ((co.fst_size==0) || ( ttlSeq <= co.fst_size + co.fst_head ) ) ){
					mflag=false;
				}else{
					mflag=true;
					break;
				}
			}
#ifdef INDEXING_INPUT_ORDER
			toInputOrder.push_back(num_of_seq_of_all_input-1);
#endif
			num_of_seq++;

			break;
		case ' ':
			break;
		case '\n':
			break;
		case '\r':
			break;
		case ',':
			if(co.charType==INPUT_INT){
				map_int_str[map_int_pos]='\0';
				int tmp_i = atoi(map_int_str); 
				if( char_map.insert( map<int, int>::value_type( tmp_i,map_int ) ).second ){
					map_int++;
				}
				map_int_pos=0;
				ltmp++;
				seq_length++;
			}
			break;
		default:
			if( (co.isMapMode || co.isPartialMode) && mflag==true){
				break;
			}
			if(co.charType==INPUT_INT){
				map_int_str[map_int_pos++]=ch;
				break;
			}
			if(ct.toInt[ch]==ct.unknown_character){
				lflag=true;
				num_of_unknown_chars++;
			}
			ltmp++;
			seq_length++;
			break;
		}
	}
	if(!co.exclude_unknown_character){
		if(num_of_unknown_chars<=co.distance)
			lflag=false;
	}

	if(lflag){
		seq_length -= ltmp;
		num_of_seq--;
#ifdef INDEXING_INPUT_ORDER
		toInputOrder.pop_back();
#endif
	}else{
		if(max_seq_length<ltmp) max_seq_length = ltmp;
		if(min_seq_length>ltmp && ltmp>0) min_seq_length = ltmp;
	}

	if(!num_of_seq){
		cerr<<"\nERROR: no valid input sequence.\n       use -u option if sequence including unknown characters.\n";
		exit(1);		
	}

	if(co.distanceType == HAMMING_DISTANCE){
		if(max_seq_length!=min_seq_length){
			cerr<<"ERROR: Input sequences are not same length\n       min_length="<<min_seq_length<<"  max_length="<<max_seq_length<<"\n";
			if(max_seq_length-min_seq_length < MAX_RECOMMENDED_EDIT_DISTANCE){
				cerr<<"       or use edit distance option: -t E -d "<<max_seq_length-min_seq_length<<"\n";				
			}
			exit(1);
		}
	}else{
		if((max_seq_length-co.distance) > min_seq_length){
			cerr<<"ERROR: The length of input sequences shoud range from minimum length of sequence to minimum length of sequence + edit_distance\n       min_length="<<min_seq_length<<"  max_length="<<max_seq_length<<"\n";
			cerr<<"       Use longer edit distance\n";
			exit(1);
		}
	}

	long long size_lim=1;
	for(int i=0;i<sizeof(long)*8-1;i++){
		size_lim = size_lim*2;
	}
	if(co.isRevComp){
		size_lim /=2;
	}
	if(seq_length>(size_lim-1)){
		cerr<<"\nERROR: Too large data. Use LL package.\n";
		exit(1);
	}

	fseek(fp,0,SEEK_SET);

	// cpy data from file.
	org_seq = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*seq_length);
	if(org_seq==NULL){
		cerr<<"ERROR: cannot allocate memory for org_seq\n";
		exit(1);
	}
	//+1:head[i]+j<head[i+1]のため。
	head = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*(num_of_seq+1));
	if(head==NULL){
		cerr<<"ERROR: cannot allocate memory for head\n";
		exit(1);
	}

	ltmp=0;
	lflag=false;
	TYPE_INDEX cur=0;
	TYPE_INDEX seqID=0;
	char tid[MAX_FASTA_ID_LENGTH];
	int i;
	num_of_unknown_chars = 0;

	ttlSeq=0;
	mflag=false;

	while(1){
		ch = getc(fp);
		if(ch == EOF){
			break;
		}
		switch(ch){
		case '>':
			if(!co.exclude_unknown_character){
				if(num_of_unknown_chars<=co.distance)
					lflag=false;
			}
			if(lflag){
				cur -=ltmp;
				seqID--;
				seqName.pop_back();
				if(co.isMapMode){
					if(ttlSeq>=co.snd_head){
						snd_size_tmp--;
					}else{
						fst_size_tmp--;
					}
				}
				if(co.isPartialMode){
					fst_size_tmp--;
				}
			}
			lflag=false;
			ltmp=0;
			num_of_unknown_chars=0;

			if(co.isMapMode){
				ttlSeq++;
				if( ttlSeq > co.snd_head && ((co.snd_size==0) || ( ttlSeq <= co.snd_size + co.snd_head ) ) ){
					snd_size_tmp++;
					mflag=false;
				}else if( ttlSeq > co.fst_head && ((co.fst_size==0) || ( ttlSeq <= co.fst_size + co.fst_head ) ) ){
					fst_size_tmp++;
					mflag=false;
				}else{
					mflag=true;
					break;
				}
			}
			if(co.isPartialMode){
				ttlSeq++;
				if( ttlSeq > co.fst_head && ((co.fst_size==0) || ( ttlSeq <= co.fst_size + co.fst_head ) ) ){
					fst_size_tmp++;
					mflag=false;
				}else{
					mflag=true;
					break;
				}
			}

			i=0;
			while(1){
				ch = getc(fp);
				if(ch == '\n'){
					break;
				}
				if(i<MAX_FASTA_ID_LENGTH-1) tid[i++]=ch;
			}
			tid[i]='\0';
			head[seqID++] = cur;
			seqName.push_back(tid);
			break;
		case ' ':
			break;
		case '\n':
			break;
		case '\r':
			break;
		case ',':
			if(co.charType==INPUT_INT){
				map_int_str[map_int_pos]='\0';
				int tmp_i = atoi(map_int_str);
				org_seq[cur++]=char_map[tmp_i];
				map_int_pos=0;
				ltmp++;
				seq_length++;
			}
			break;
		default:
			if((co.isMapMode || co.isPartialMode) && mflag==true){
				break;
			}
			if(co.charType==INPUT_INT){
				map_int_str[map_int_pos++]=ch;
				break;
			}
			if(cur<seq_length){
				org_seq[cur++]=ct.toInt[ch];
			}else{
				cur++;
			}
			if(ct.toInt[ch]==ct.unknown_character){
				lflag=true;
				num_of_unknown_chars++;
			}
			ltmp++;
			break;		
		}
	}
	if(!co.exclude_unknown_character){
		if(num_of_unknown_chars<=co.distance){
			lflag=false;
		}
	}
	if(lflag){
		cur -=ltmp;
		seqID--;
		seqName.pop_back();
		if(co.isMapMode){
			if(seqID>=co.snd_head){
				snd_size_tmp--;
			}else{
				fst_size_tmp--;
			}
		}
		if(co.isPartialMode){
			fst_size_tmp--;			
		}
	}
	if(co.isMapMode||co.isPartialMode){
		co.fst_size = fst_size_tmp;
		co.snd_size = snd_size_tmp;
	}
	head[seqID] = cur;

	fclose(fp);

	/***
	char set for input_int
	***/
	if(co.charType==INPUT_INT){
		ct.num_of_character = map_int;
		ct.overlap_character =map_int;
		ct.lim_wild_card = map_int+1;
		//cerr<<ct.overlap_character<<"\n";
		//cerr<<ct.num_of_character<<"\n";
	}

	/**
	for(int i=0;i<num_of_seq;i++){
	cerr<<seqName[i]<<" "<<head[i]<<"\n";
	}
	**/
	return(0);
}

void preprocess::readRevComp(){
	/**
	seq1, seq2, seq3 
	-> 	seq1, seq1-rev, seq2, seq2-rev, seq3, seq3-rev,
	**/

	TYPE_CHARACTER *revtable = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*ct.num_of_character);
	revtable[ct.num_of_character-1] = ct.unknown_character;

	revtable[DNA_A] = REV_DNA_A;
	revtable[DNA_T] = REV_DNA_T;
	revtable[DNA_G] = REV_DNA_G;
	revtable[DNA_C] = REV_DNA_C;

	TYPE_LONG double_seq_length = seq_length*2;
	TYPE_INDEX double_num_of_seq = num_of_seq*2;
	TYPE_CHARACTER *double_org_seq = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*double_seq_length);
	if(double_org_seq==NULL){
		cerr<<"ERROR: cannot allocate memory for d_org_seq\n";
		exit(1);
	}
	TYPE_INDEX *double_head = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*(double_num_of_seq+1));
	if(double_head==NULL){
		cerr<<"ERROR: cannot allocate memory for d_head\n";
		exit(1);
	}

	vector<string> tmp = seqName;
	seqName.clear();

	for(TYPE_INDEX cnt=0;cnt<num_of_seq;cnt++){
		seqName.push_back(tmp[cnt]);
		seqName.push_back(tmp[cnt].append("-rev"));
	}
	tmp.clear();

#ifdef INDEXING_INPUT_ORDER
	vector<TYPE_INDEX> tmpd = toInputOrder;
	toInputOrder.clear();
	for(int i=0;i<tmpd.size();i++){
		toInputOrder.push_back(tmpd[i]*2);
		toInputOrder.push_back(tmpd[i]*2+1);
	}
	tmpd.clear();
#endif

	int ptr = head[0];
	for(int cnt=0;cnt<num_of_seq;cnt++){
		int len = head[cnt+1]-head[cnt];
		for(int i=0;i<len;i++){
			*(double_org_seq++) = *(org_seq++);
		}
		org_seq--;
		for(int i=0;i<len;i++){
			*(double_org_seq++) = revtable[ *(org_seq--) ];
		}
		org_seq++;
		org_seq += len;
		double_head[cnt*2] = ptr;
		ptr += len;
		double_head[cnt*2+1] = ptr;
		ptr += len;
	}
	co.fst_size *= 2;
	co.snd_size *= 2;

	double_head[double_num_of_seq] = double_seq_length;
	org_seq -= seq_length;
	double_org_seq -= double_seq_length;

	num_of_seq = double_num_of_seq;
	seq_length = double_seq_length;
	TYPE_INDEX *thp;
	thp = head;
	head = double_head;
	free(thp);

	TYPE_CHARACTER *tcp;
	tcp = org_seq;
	org_seq = double_org_seq;
	free(tcp);
	//	for(int i=0;i<seq_length;i++) cerr<<(int)org_seq[i]<<" ";
}

/**

else{
for(int cnt=0;cnt<num_of_seq;cnt++){
seqName.push_back(tmp[cnt]);
}
for(int cnt=0;cnt<num_of_seq;cnt++){
seqName.push_back(tmp[cnt].append("-rev"));
}
tmp.clear();
int ptr = head[0];
for(int cnt=0;cnt<num_of_seq;cnt++){
int len = head[cnt+1]-head[cnt];
for(int i=0;i<len;i++){
*(double_org_seq++) = *(org_seq++);
}
double_head[cnt] = ptr;
ptr += len;
}
org_seq -= (seq_length+1);
for(int cnt=0;cnt<num_of_seq;cnt++){
int len = head[cnt+1]-head[cnt];
org_seq += len;
for(int i=0;i<len;i++){
*(double_org_seq++) = revtable[ *(org_seq--) ];
}
org_seq += len;
double_head[cnt+num_of_seq] = ptr;
ptr += len;
}
org_seq++;
co.fst_size = num_of_seq;
co.snd_size = num_of_seq;
}

**/

void preprocess::cpOrgseq(seq &sq)
{
	sq.num_of_seq = num_of_seq;
	sq.head = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*(sq.num_of_seq+1));

	if(co.distanceType==HAMMING_DISTANCE){
		sq.seq_length = seq_length;
	}else{
		sq.seq_length = sq.num_of_seq*(max_seq_length+co.allowedgap) + co.allowedgap;
	}
	sq.nr_seq = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*sq.seq_length);
	sq.mem_aloc_flag=true;
	if(co.distanceType==HAMMING_DISTANCE){
		memcpy(sq.nr_seq,org_seq,sizeof(TYPE_CHARACTER)*seq_length);
		memcpy(sq.head,head,sizeof(TYPE_INDEX)*(sq.num_of_seq+1));
	}else{
		TYPE_INDEX i,j,iCnt=0,h;
		for(i=0;i<num_of_seq;i++){
			sq.head[i] = i*(max_seq_length+co.allowedgap);
			for(j=0;j<co.allowedgap;j++){
				sq.nr_seq[iCnt++] = ct.overlap_character;
			}
			h = (max_seq_length-(head[i+1]-head[i]))/2;
			for(j=0;j<h;j++){
				sq.nr_seq[iCnt++] = ct.lim_wild_card;
			}
			for(j=0;j<head[i+1]-head[i];j++){
				sq.nr_seq[iCnt++] = org_seq[ j+head[i] ];
			}
			for(j=0;j<max_seq_length-(head[i+1]-head[i])-h;j++){
				sq.nr_seq[iCnt++] = ct.lim_wild_card;
			}
		}
		for(j=0;j<co.allowedgap;j++){
			sq.nr_seq[iCnt++] = ct.overlap_character;
		}
	}
	sq.max_seq_length = max_seq_length;
	sq.min_seq_length = min_seq_length;
	sq.seqName = (string*)malloc(sizeof(string)*num_of_seq);
	for(TYPE_INDEX i=0;i<num_of_seq;i++){
		sq.seqName[i] = seqName[i];
		varInt *vTmp = new varInt(1);
		vTmp->val[0] = i;
//		sq.toOrgSeqIndex.push_back(vTmp);
	}

	if(co.isMapMode){
		TYPE_INDEX i;
		sq.org_seq_mapID = (char*)malloc(sizeof(char)*num_of_seq);
		for(i=0;i<co.fst_size;i++){
			sq.org_seq_mapID[i] = FST_DATASET;
		}
		for(;i<num_of_seq;i++){
			sq.org_seq_mapID[i] = SND_DATASET;
		}
		sq.nr_seq_mapID = sq.org_seq_mapID;
	}
}
#if 0
void preprocess::sortOrgSeq()
{
	sort = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*(num_of_seq));
	work = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*(num_of_seq));

	if(sort==NULL){
		cerr<<"ERROR: cannot allocate memory for sort\n";
		exit(1);
	}
	if(work==NULL){
		cerr<<"ERROR: cannot allocate memory for work\n";
		exit(1);
	}

	// +1は桁用
	int num_of_bucket = ct.num_of_character+1;
	int *bucket = (int *)malloc(sizeof(int)*num_of_bucket);
	if(bucket==NULL){
		cerr<<"ERROR: cannot allocate memory for bucket\n";
		exit(1);
	}

	TYPE_INDEX i,j;
	for(i=0;i<num_of_seq;i++){
		sort[i]=i;
	}

	for(j=0;j<max_seq_length;j++){
		for(i=0;i<num_of_bucket;i++) bucket[i]=0;
		for(i=0;i<num_of_seq;i++){
			if(head[i]+j<head[i+1]){
				bucket[ org_seq[ head[i]+j ] ]++;
			}else{
				bucket[num_of_bucket-1]++;
			}
		}
		for(i=1;i<num_of_bucket;i++) bucket[i] += bucket[i-1];
		// minus check
		//		for(i=num_of_seq-1;i>=0;i--){
		for(i=num_of_seq-1;;i--){
			if(head[ sort[i] ]+j<head[ sort[i]+1 ]){
				work[ --bucket[ org_seq[ head[ sort[i] ]+j ] ] ] = sort[i];
			}else{
				work[ --bucket[num_of_bucket-1] ] = sort[i];
			}
			if(i==0) break;
		}
		for(i=0;i<num_of_seq;i++) sort[i] = work[i];

		/** for debug 
		for(i=0;i<num_of_seq;i++) fprintf(stderr," %d -",sort[i]);	
		fprintf(stderr,"\n");
		**/
	}
	free(bucket);
}
#endif


void preprocess::mkNonredundantData(char *outputfile)
{
	if(sort==NULL){
		cerr<<"not yet sorted..\n";
		exit(1);
	}
	FILE *ofp =fopen(outputfile,"w");
	seq sq;
	getSeq(sq);

	TYPE_INDEX cnt=0;
	for(TYPE_INDEX i=0;i<sq.seq_length;i++){
		if(sq.head[cnt]==i){
			if(cnt!=0) putc('\n',ofp);
			putc('>',ofp);
			fputs(sq.seqName[ sq.toOrgSeqIndex[cnt++]->val[0] ].c_str(),ofp);
			putc('\n',ofp);
		}
		putc(ct.toChar[ sq.nr_seq[i] ],ofp);
	}
	fclose(ofp);
}
