#include"mscls.h"
#include<map>
#include<time.h>


int preprocess::readOrgFastqSeqViaStdin()
{
	FILE *fp;
	char *tmp_stdin=NULL;
	char *org_tmp_stdin;
	long long tmp_stdin_ptr=0;
	fp = stdin;
	tmp_stdin = (char *)malloc(sizeof(char)*co.tmp_file_size);		

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

	int map_int=0, map_int_pos=0;
	char map_int_str[SIZE_OF_TYPE_CHARACTER];
	int num_of_unknown_chars=0;

	int dnflag=false;


	// check if file type is collect
	ch='@';
	while(1){
		tmp_stdin[tmp_stdin_ptr]=ch;
		tmp_stdin_ptr++;
		if(tmp_stdin_ptr > co.tmp_file_size){
			int alpha;
			if(tmp_stdin_ptr > 1000000000){
				alpha = tmp_stdin_ptr/20;
			}else{
				alpha = tmp_stdin_ptr;
			}
			if(realloc(tmp_stdin,sizeof(char)*(tmp_stdin_ptr + alpha))==NULL){
				cerr<<"can not allocate tmp_file\n";
				exit(1);
			}
		}

		switch(ch){
		case '@':
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

			//Clear the line with @ 
			while((ch=getc(fp))!='\n'){
				tmp_stdin[tmp_stdin_ptr]=ch;
				tmp_stdin_ptr++;
				if(tmp_stdin_ptr > co.tmp_file_size){
					int alpha;
					if(tmp_stdin_ptr > 1000000000){
						alpha = tmp_stdin_ptr/20;
					}else{
						alpha = tmp_stdin_ptr;
					}
					//？？？　tmp_stdin = (char *)realloc としなくて良い？？
					if(realloc(tmp_stdin,sizeof(char)*(tmp_stdin_ptr + alpha))==NULL){
						cerr<<"can not allocate tmp_file\n";
						exit(1);
					}
				}			
			}
			tmp_stdin[tmp_stdin_ptr]=ch;
			tmp_stdin_ptr++;
			
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
			//	map or partial モードで範囲外の時はインクリメントされない。 
			num_of_seq++;

			break;
		case ' ':
			break;
		case '\n':
		    dnflag=false;
			while(1){
				ch=getc(fp);
				tmp_stdin[tmp_stdin_ptr]=ch;
				tmp_stdin_ptr++;
				if(tmp_stdin_ptr > co.tmp_file_size){
					int alpha;
					if(tmp_stdin_ptr > 1000000000){
						alpha = tmp_stdin_ptr/20;
					}else{
						alpha = tmp_stdin_ptr;
					}
					// tmp_stdin
					if(realloc(tmp_stdin,sizeof(char)*(tmp_stdin_ptr + alpha))==NULL){
						cerr<<"can not allocate tmp_file\n";
						exit(1);
					}
				}
				if(ch=='\n'){
					if(dnflag){
						break;
					}else{
						dnflag=true;
					}
				}
			}
			break;
		case '\r':
			break;
		default:
			if( (co.isMapMode || co.isPartialMode) && mflag==true){
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

		ch=getc(fp);
		if(ch==EOF){
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

		tmp_stdin[tmp_stdin_ptr]=EOF;

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

	org_tmp_stdin = tmp_stdin;
	while(1){
		ch = *(tmp_stdin++);
		if(ch == EOF){
			break;
		}
		switch(ch){
		case '@':
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
					while(1){
						ch = *(tmp_stdin++);
						if(ch == '\n'){
							break;
						}
					}	
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
					while(1){
						ch = *(tmp_stdin++);
						if(ch == '\n'){
							break;
						}
					}
					break;
				}
			}

			i=0;
			while(1){
				ch = *(tmp_stdin++);
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
			while(1){
				ch = *(tmp_stdin++);
				if(ch == '\n'){
					break;
				}
			}
			while(1){
				ch = *(tmp_stdin++);
				if(ch == '\n'){
					break;
				}
			}
			break;
		case '\r':
			break;
		default:
			if((co.isMapMode || co.isPartialMode) && mflag==true){
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
	free(org_tmp_stdin);

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


int preprocess::readOrgFastaSeqViaStdin()
{
	FILE *fp;
	char *tmp_stdin=NULL;
	char *org_tmp_stdin;
	long long tmp_stdin_ptr=0;

	fp = stdin;
	tmp_stdin = (char *)malloc(sizeof(char)*co.tmp_file_size);		


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

	ch='>';
	while(1){
		tmp_stdin[tmp_stdin_ptr]=ch;
		tmp_stdin_ptr++;
		if(tmp_stdin_ptr > co.tmp_file_size){
			int alpha;
			if(tmp_stdin_ptr > 1000000000){
				alpha = tmp_stdin_ptr/20;
			}else{
				alpha = tmp_stdin_ptr;
			}
			if(realloc(tmp_stdin,sizeof(char)*(tmp_stdin_ptr + alpha))==NULL){
				cerr<<"can not allocate tmp_file\n";
				exit(1);
			}
		}

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
			while((ch=getc(fp))!='\n'){
				tmp_stdin[tmp_stdin_ptr]=ch;
				tmp_stdin_ptr++;
				if(tmp_stdin_ptr > co.tmp_file_size){
					int alpha;
					if(tmp_stdin_ptr > 1000000000){
						alpha = tmp_stdin_ptr/20;
					}else{
						alpha = tmp_stdin_ptr;
					}
					if(realloc(tmp_stdin,sizeof(char)*(tmp_stdin_ptr + alpha))==NULL){
						cerr<<"can not allocate tmp_file\n";
						exit(1);
					}
				}			
			}
			tmp_stdin[tmp_stdin_ptr]=ch;
			tmp_stdin_ptr++;

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
		ch=getc(fp);
		if(ch==EOF){
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

	tmp_stdin[tmp_stdin_ptr]=EOF;

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

	org_tmp_stdin = tmp_stdin;
	while(1){
		ch = *(tmp_stdin++);
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
				ch = *(tmp_stdin++);
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