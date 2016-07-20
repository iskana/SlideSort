#include"mscls.h"
#include<map>
#include<set>

//#define DEBUG
#ifdef DEBUG
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

	int i,j;
	for(i=0;i<num_of_seq;i++){
		sort[i]=i;
	}

	for(j=0;j<max_seq_length;j++){
		memset(bucket,0x00,sizeof(int)*num_of_bucket);
		int tmp;
		for(i=0;i<num_of_seq;i++){
			tmp = org_seq[ head[i]+j ];
				bucket[ tmp ]++;
		}
		for(i=1;i<num_of_bucket;i++) bucket[i] += bucket[i-1];
		for(i=num_of_seq-1;;i--){
			tmp = org_seq[ head[ sort[i] ]+j ];
			work[ --bucket[ tmp ] ] = sort[i];
			if(i==0) break;
		}
		for(i=0;i<num_of_seq;i++) sort[i] = work[i];

		/** for debug 
		for(i=0;i<num_of_seq;i++) fprintf(stderr," %d -",sort[i]);	
		fprintf(stderr,"\n");
		**/
	}
	free(bucket);
	exit(1);
}

#else

void preprocess::keepOrgSeq()
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

	TYPE_INDEX i;
	for(i=0;i<num_of_seq;i++){
		sort[i]=i;
	}
}

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
	int char_size = ct.num_of_character;
	if(ct.isUseWildCard){
		char_size++;
	}
	int word_size = log((double)num_of_seq)/log((double)char_size)-1;
	if(word_size<=0) word_size=1;
	if(min_seq_length<word_size) word_size = min_seq_length;
	int num_of_bucket = pow(char_size,word_size);
	TYPE_INDEX *bucket = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*num_of_bucket);
	int *bucketID = (int *)malloc(sizeof(int)*num_of_bucket);
	int bid_ptr=0;

	int num = ((min_seq_length%word_size)+(max_seq_length-min_seq_length))/word_size + min_seq_length/word_size;
	int l_word_size = ((min_seq_length%word_size)+(max_seq_length-min_seq_length))%word_size;

//	cerr<<"ws:"<<word_size<<" num:"<<num<<"  lws:"<<l_word_size<<"\n";
//	cerr<<"mid:"<<min_seq_length/word_size<<"\n";
//	cerr<<"cs:"<<char_size<<"\n";

	if(bucket==NULL){
		cerr<<"ERROR: cannot allocate memory for bucket\n";
		exit(1);
	}

	TYPE_INDEX i,j;
	for(i=0;i<num_of_seq;i++){
		sort[i]=i;
	}
	TYPE_INDEX tmp=0;
	for(j=0;j< min_seq_length/word_size;j++){
		memset(bucket,0x00,sizeof(TYPE_INDEX)*num_of_bucket);
		bid_ptr=0;
		for(i=0;i<num_of_seq;i++){
			tmp=0;
			for(int cnt=0;cnt<word_size;cnt++){
				tmp = org_seq[ head[i]+cnt+j*word_size ] + tmp*char_size;
			}
			if(bucket[ tmp ]==0){
				bucketID[bid_ptr]=tmp;
				bid_ptr++;
			}
			bucket[ tmp ]++;
		}
		for(i=1;i<bid_ptr;i++) bucket[bucketID[i]] += bucket[bucketID[i-1]];
		for(i=num_of_seq-1;;i--){
			tmp=0;
			for(int cnt=0;cnt<word_size;cnt++){
				tmp = org_seq[ head[ sort[i] ]+ cnt + j*word_size ] + tmp*char_size;
			}
			work[ --bucket[ tmp ] ] = sort[i];
			if(i==0) break;
		}
		for(i=0;i<num_of_seq;i++) sort[i] = work[i];
	//	/** for debug 
//		for(i=0;i<num_of_seq;i++) fprintf(stderr," %d -",sort[i]);	
//		fprintf(stderr,"\n");
	//	**/
	}

	for(;j<num;j++){
		memset(bucket,0x00,sizeof(TYPE_INDEX)*num_of_bucket);
		bid_ptr=0;
		for(i=0;i<num_of_seq;i++){
			tmp=0;
			for(int cnt=0;cnt<word_size;cnt++){
				if( head[i]+cnt+j*word_size<head[i+1] ){
					tmp = org_seq[ head[i]+cnt+j*word_size ] + tmp*char_size;
				}else{
					tmp = ct.overlap_character + tmp*char_size;
				}
			}
			if(bucket[ tmp ]==0){
				bucketID[bid_ptr]=tmp;
				bid_ptr++;
			}
			bucket[ tmp ]++;
		}
		for(i=1;i<bid_ptr;i++) bucket[bucketID[i]] += bucket[bucketID[i-1]];
		for(i=num_of_seq-1;;i--){
			tmp=0;
			for(int cnt=0;cnt<word_size;cnt++){
				if( head[ sort[i] ]+cnt+j*word_size<head[ sort[i]+1 ] ){
					tmp = org_seq[ head[ sort[i] ]+ cnt + j*word_size ] + tmp*char_size;
				}else{
					tmp = ct.overlap_character + tmp*char_size;					
				}
			}
			work[ --bucket[ tmp ] ] = sort[i];
			if(i==0) break;
		}
		for(i=0;i<num_of_seq;i++) sort[i] = work[i];

		//	/** for debug 
//		cerr<<j<<"   ";
//		for(i=0;i<num_of_seq;i++) fprintf(stderr," %d -",sort[i]);	
//		fprintf(stderr,"\n");
		//	**/
	}
	if(l_word_size>0){
		memset(bucket,0x00,sizeof(TYPE_INDEX)*num_of_bucket);
		bid_ptr=0;
		for(i=0;i<num_of_seq;i++){
			tmp=0;
			for(int cnt=0;cnt<l_word_size;cnt++){
				if( head[i]+cnt+j*word_size<head[i+1] ){
					tmp = org_seq[ head[i]+cnt+j*word_size ] + tmp*char_size;
				}else{
					tmp = ct.overlap_character + tmp*char_size;
				}
			}
			if(bucket[ tmp ]==0){
				bucketID[bid_ptr]=tmp;
				bid_ptr++;
			}
			bucket[ tmp ]++;
		}
		for(i=1;i<bid_ptr;i++) bucket[bucketID[i]] += bucket[bucketID[i-1]];
		for(i=num_of_seq-1;;i--){
			tmp=0;
			for(int cnt=0;cnt<l_word_size;cnt++){
				if( head[ sort[i] ]+cnt+j*word_size<head[ sort[i]+1] ){
					tmp = org_seq[ head[ sort[i] ]+ cnt + j*word_size ] + tmp*char_size;
				}else{
					tmp = ct.overlap_character + tmp*char_size;					
				}
			}
			work[ --bucket[ tmp ] ] = sort[i];
			if(i==0) break;
		}
		for(i=0;i<num_of_seq;i++) sort[i] = work[i];

		//	/** for debug 
//		for(i=0;i<num_of_seq;i++) fprintf(stderr," %d -",sort[i]);	
//		fprintf(stderr,"\n");
		//	**/
	}

	free(bucket);
//	exit(1);
}
#endif

void preprocess::getSeq(seq &sq)
{
	//count sequence length
	TYPE_INDEX i,j,rTmp=0;
	int flag, len1, len2;

	vector<varInt*> tmpidx;

	// for edit distance
	int box_size;

	box_size = max_seq_length;
	sq.num_of_seq_of_all_input = num_of_seq_of_all_input;

	if(co.isSortOrgSeq){
		sq.num_of_seq++;
		sq.seq_length += head[ sort[0]+1 ] - head[ sort[0] ];
		rTmp++;
		for(i=1;i<num_of_seq;i++){
			flag=0;
			len1 = head[ sort[i]+1 ] - head[ sort[i] ];
			len2 = head[ sort[i-1]+1 ] - head[ sort[i-1] ];
			if(len1 == len2){
				for(j=0;j<len1;j++)
					if(org_seq[ head[sort[i]] + j ] != org_seq[ head[sort[i-1]] + j ]) {flag=1; break;}
			}else{
				flag=1;
			}
			if(flag){
				sq.seq_length += len1;
				sq.num_of_seq++;
				varInt *vTmp = new varInt(rTmp);
				tmpidx.push_back(vTmp);
				rTmp=0;
			}
			rTmp++;
		}
		varInt *vTmp = new varInt(rTmp);
		tmpidx.push_back(vTmp);
		sq.toOrgSeqIndex = (varInt**)malloc(sizeof(varInt*)*tmpidx.size());

		for(i=0;i<tmpidx.size();i++){
			sq.toOrgSeqIndex[i] = tmpidx[i];
		}
	}else{
		sq.num_of_seq = num_of_seq;
		sq.seq_length = seq_length;
		sq.toOrgSeqIndex = (varInt**)malloc(sizeof(varInt*)*num_of_seq);
		for(i=0;i<num_of_seq;i++){
			sq.toOrgSeqIndex[i] = new varInt(1);
		}
	}

	// map mode
	int fstFlag=false;
	int sndFlag=false;
	if(co.isMapMode){
		sq.org_seq_mapID = (char*)malloc(sizeof(char)*num_of_seq);
		sq.nr_seq_mapID = (char*)malloc(sizeof(char)*sq.num_of_seq);
		for(i=0;i<co.fst_size;i++){
			sq.org_seq_mapID[i] = FST_DATASET;
		}
		for(;i<num_of_seq;i++){
			sq.org_seq_mapID[i] = SND_DATASET;
		}
	}

	// copy org_seq to nr_seq
	sq.head = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*(sq.num_of_seq+1));
	sq.mem_aloc_flag=true;
	if(co.distanceType==EDIT_DISTANCE){
		sq.seq_length = sq.num_of_seq*(box_size+co.allowedgap) + co.allowedgap;
	}
	sq.nr_seq = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*sq.seq_length);

	TYPE_INDEX sCnt=0, iCnt=0;
	int h=0;
	rTmp=0;

	sq.toOrgSeqIndex[iCnt]->val[rTmp++]=sort[0];
	sq.head[iCnt]=sCnt;
	if(co.distanceType==EDIT_DISTANCE){
		for(j=0;j<co.allowedgap;j++){
			sq.nr_seq[sCnt++] = ct.overlap_character;
		}
		h = (max_seq_length-(head[sort[0]+1]-head[sort[0]]))/2;

		for(j=0;j<h;j++){
			sq.nr_seq[sCnt++] = ct.lim_wild_card;
		}
	}
	for(j=head[ sort[0] ];j<head[ sort[0]+1 ];j++){
		sq.nr_seq[sCnt++]=org_seq[ j ];
	}
	if(co.distanceType==EDIT_DISTANCE){
		for(j=0;j<max_seq_length-(head[sort[0]+1]-head[sort[0]])-h;j++){
			sq.nr_seq[sCnt++] = ct.lim_wild_card;
		}
	}
	iCnt++;

	if(co.isMapMode){
		if(sq.org_seq_mapID[sort[0]]==FST_DATASET){
			fstFlag=true;		
		}else{
			sndFlag=true;		
		}
	}

	// for match of unknown characters.
	int num_unknown_char;
	sq.num_of_unknown_char_matchs = (int *)malloc(sizeof(int)*sq.num_of_seq);
	memset(sq.num_of_unknown_char_matchs,0x00,sizeof(int)*sq.num_of_seq);
	if(co.exclude_unknown_character){
		sq.num_of_unknown_char_matchs[0] = 0;
	}else{
		num_unknown_char=0;
		for(j=head[ sort[0] ];j<head[ sort[0]+1 ];j++){
			if(org_seq[ j ]==ct.unknown_character)
				num_unknown_char++;				
		}
		sq.num_of_unknown_char_matchs[0] = num_unknown_char;
	}
	
	for(i=1;i<num_of_seq;i++){
		flag=0;
		len1 = head[ sort[i]+1 ] - head[ sort[i] ];
		len2 = head[ sort[i-1]+1 ] - head[ sort[i-1] ];
		if(len1 == len2 && co.isSortOrgSeq==true){
			for(j=0;j<len1;j++)
				if(org_seq[ head[sort[i]] + j ] != org_seq[ head[sort[i-1]] + j ]) {flag=1; break;}
		}else{
			flag=1;
		}
		if(flag){
			rTmp=0;
			sq.head[iCnt]=sCnt;
			if(co.distanceType==EDIT_DISTANCE){
				for(j=0;j<co.allowedgap;j++){
					sq.nr_seq[sCnt++] = ct.overlap_character;
				}
				h = (max_seq_length-len1)/2;
				for(j=0;j<h;j++){
					sq.nr_seq[sCnt++] = ct.lim_wild_card;
				}
			}
			for(j=0;j<len1;j++){
				sq.nr_seq[sCnt++] = org_seq[ head[sort[i]] + j ];
			}
			if(co.distanceType==EDIT_DISTANCE){
				for(j=0;j<max_seq_length-len1-h;j++){
					sq.nr_seq[sCnt++] = ct.lim_wild_card;
				}
			}
			if(co.isMapMode){
				if(fstFlag==true){
					sq.nr_seq_mapID[iCnt-1] = FST_DATASET;
				}
				if(sndFlag==true){
					sq.nr_seq_mapID[iCnt-1] = SND_DATASET;
				}
				if( fstFlag==true && sndFlag==true){
					sq.nr_seq_mapID[iCnt-1] = BOTH_DATASET;
				}
				fstFlag=false;
				sndFlag=false;
			}
			if(co.exclude_unknown_character){
				sq.num_of_unknown_char_matchs[iCnt] = 0;
			}else{
				num_unknown_char=0;
				for(j=0;j<len1;j++){
					if(org_seq[ head[sort[i]] + j ]==ct.unknown_character)
						num_unknown_char++;				
				}
				sq.num_of_unknown_char_matchs[iCnt] = num_unknown_char;
			}
			iCnt++;
		}
		if(co.isMapMode){
			if(sq.org_seq_mapID[sort[i]]==FST_DATASET){
				fstFlag=true;
			}else{
				sndFlag=true;		
			}
		}
		sq.toOrgSeqIndex[iCnt-1]->val[rTmp++]=sort[i];

	}

	if(co.isMapMode){
		if(fstFlag==true){
			sq.nr_seq_mapID[iCnt-1] = FST_DATASET;
		}
		if(sndFlag==true){
			sq.nr_seq_mapID[iCnt-1] = SND_DATASET;
		}
		if( fstFlag==true && sndFlag==true){
			sq.nr_seq_mapID[iCnt-1] = BOTH_DATASET;
		}
	}

	//20100315 確かめてない。
	while(sCnt<sq.seq_length)
		sq.nr_seq[sCnt++] = ct.overlap_character; 
//		sq.nr_seq[sq.seq_length-1] = ct.overlap_character; 
	sq.head[sq.num_of_seq]=sq.seq_length;

	sq.seqName = new string[seqName.size()];
	sq.num_of_valid_seq = seqName.size();
	for(i=0;i<seqName.size();i++){
		sq.seqName[i] = seqName[i];
	}

#ifdef INDEXING_INPUT_ORDER
	for(i=0;i<toInputOrder.size();i++){
		sq.toInputOrder.push_back(toInputOrder[i]);
	}
#endif

	sq.max_seq_length = max_seq_length;
	sq.min_seq_length = min_seq_length;

	/** for debug**/
#if 0
	for(i=0;i<seqName.size();i++){
		cerr<<seqName[i]<<"\n";
		if(co.isMapMode)
			cerr<<"\t"<<(int)sq.org_seq_mapID[i]<<"\n";
	}
	cerr<<"--\n";
	for(i=0;i<sq.num_of_seq;i++){
		for(j=0;j<sq.toOrgSeqIndex[i]->size;j++){
			cerr<< seqName[sq.toOrgSeqIndex[i]->val[j]]<<"-";
		}
		cerr<<"\n";
		if(co.isMapMode){
			cerr<<(int)sq.nr_seq_mapID[i]<<"  -  ";
			for(j=0;j<sq.toOrgSeqIndex[i]->size;j++){
				cerr<< (int)sq.org_seq_mapID[ sq.toOrgSeqIndex[i]->val[j] ]<<"-";
			}
			cerr<<"\n";
		}
	}
#endif
	/**/
}

