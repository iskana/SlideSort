/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              hamming_distance.cpp
/***********************************************/

#include"mscls.h"
#include <time.h>

void multisort::setBox_nonOV_nonIC_HD()
{
	bx.num_of_box = sq.num_of_seq;
	bx.box_length = sq.max_seq_length;

	// error when free 
	//	bx.head = sq.head;

	bx.key_size = co.key_size;

	//for output
	bx.box_set_center=0;
	bx.box_set_size=1;

	//block数の決定
	bx.num_of_blocks = co.distance + bx.key_size;
	bx.short_block_length = bx.box_length/bx.num_of_blocks;
	if(bx.box_length%bx.num_of_blocks==0){
		bx.long_block_length = bx.short_block_length;
		bx.num_of_long_blocks = bx.num_of_blocks;
	}else{
		bx.long_block_length = bx.short_block_length + 1;
		bx.num_of_long_blocks = bx.box_length - bx.short_block_length*bx.num_of_blocks;
	}

	/**
	char_size = ct.num_of_character;
	word_size[SHORT_WORD] = log((double)sq.num_of_seq)/log((double)char_size)-3;
	if(word_size[SHORT_WORD]<=0) word_size[SHORT_WORD]=1;
	word_size[LONG_WORD] = word_size[SHORT_WORD];
	if(bx.short_block_length<word_size[SHORT_WORD]) word_size[SHORT_WORD] = bx.short_block_length;
	if(bx.long_block_length<word_size[LONG_WORD]) word_size[LONG_WORD] = bx.long_block_length;

	num_of_bucket[SHORT_WORD] = pow(char_size,word_size[SHORT_WORD]);
	num_of_bucket[LONG_WORD] = pow(char_size,word_size[LONG_WORD]);

	bucket = (TYPE_INDEX **)malloc(sizeof(TYPE_INDEX*)*2);
	bucket[SHORT_WORD] = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*num_of_bucket[SHORT_WORD]);
	bucket[LONG_WORD] = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*num_of_bucket[LONG_WORD]);

	if(bucket[SHORT_WORD]==NULL || bucket[LONG_WORD] == NULL){
		cerr<<"ERROR: can not allocate memory for bucket\n";
		exit(1);
	}

	bucketID = (int **)malloc(sizeof(int*)*2);
	bucketID[SHORT_WORD] = (int *)malloc(sizeof(int)*num_of_bucket[SHORT_WORD]);
	bucketID[LONG_WORD] = (int *)malloc(sizeof(int)*num_of_bucket[LONG_WORD]);
	if(bucketID[SHORT_WORD]==NULL || bucketID[LONG_WORD] == NULL){
		cerr<<"ERROR: can not allocate memory for bucketID\n";
		exit(1);
	}

	num_of_words[SHORT_WORD] = bx.short_block_length/word_size[SHORT_WORD];
	num_of_words[LONG_WORD] = bx.long_block_length/word_size[LONG_WORD];

	l_word_size[SHORT_WORD] = bx.short_block_length%word_size[SHORT_WORD]; 
	l_word_size[LONG_WORD] = bx.long_block_length%word_size[LONG_WORD];
	**/
}

void multisort::findSimilarHDpairs(){
	if(sq.max_seq_length<SIZE_OF_REGISTER){
		findSimilarHDpairsHP128();
//		FOR TEST ONLY
//		findSimilarHDpairsHP128_for_kmer_test();
	}else{
		findSimilarHDpairsHP();
	}
}

#if 0
void multisort::findSimilarpairs(){

	int cHead=0,nHead=0;
	int ptr = bx.key_size -1;

	int flag = NON_REDUNDANT, id1, id2, cPtr, hd;

	for(int cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
		if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
			nHead = cnt;
			for(int i=cHead;i<nHead;i++){
				for(int j=i+1;j<nHead;j++){
					//128bitなので。
					id1 = order[ptr]->val[i];
					id2 = order[ptr]->val[j];;
					cPtr = ptr-1;

					for(int key = key_stack[ptr]-1 ;key>=0 ;key--){
//						cerr<<key<<" , "<<key_stack[cPtr]<<"\n";
						//CHK duplication
						flag = NON_REDUNDANT;
						for(int bit_pos=0;bit_pos<bs->num_of_stream;bit_pos++){
							BIT_BLOCK* sPtr1 = &(bs->stream[bit_pos][ id1 ]);
							BIT_BLOCK* sPtr2 = &(bs->stream[bit_pos][ id2 ]);
							if(key == key_stack[cPtr] ){
								if(cPtr>0)	cPtr--;
								//								cerr<<"break\n";
								break;
							}else{
								seq1 = _mm_load_si128((__m128i *)sPtr1);
								seq2 = _mm_load_si128((__m128i *)sPtr2);
								mask = _mm_load_si128((__m128i *)&block_mask[key*2]);
								mSeq1 = _mm_and_si128(seq1, mask);
								mSeq2 = _mm_and_si128(seq2, mask);
								hdc = _mm_cmpeq_epi32(mSeq1, mSeq2);
								_mm_store_si128((__m128i *)res,hdc);
								if(res[0]==0xffffffffffffffff && res[1]==0xffffffffffffffff){
									if(flag==ONE_BIT_REDUNDANT){
										flag = REDUNDANT;
										break;
									}
									flag = ONE_BIT_REDUNDANT;
								}else{
									break;
								}
								//for debug
//								if(key_stack[0]==0 && key_stack[1]==3 && key==1 && bit_pos==1 && order[ptr]->val[i] == 2 && order[ptr]->val[j] == 3)
//									order[100]->val[j]=100;
							}
						}
						if(flag==REDUNDANT){
//							cerr<<"redundant\n";
							break;
						}
					}

					//check hamming distance
					if(flag!=REDUNDANT){
						num_of_comparison++;
						for(int bit_pos=0;bit_pos<bs->num_of_stream;bit_pos++){
							BIT_BLOCK* sPtr1 = &(bs->stream[bit_pos][ id1 ]);
							BIT_BLOCK* sPtr2 = &(bs->stream[bit_pos][ id2 ]);
							seq1 = _mm_load_si128((__m128i *)sPtr1);
							seq2 = _mm_load_si128((__m128i *)sPtr2);
							xor_res = _mm_xor_si128(seq1,seq2);
							if(bit_pos==0){
								hdn = xor_res;
							}else{
								hdn = _mm_or_si128(xor_res,hdc);
							}
							hdc = hdn;
						}
						_mm_store_si128((__m128i *)res,hdc);
						//N以内をさがす。hd<Nとしたいので。
						hd = -1;

						while(res[0]!=0 && hd<co.distance){
							res[0] &=(res[0]-1);
							hd++;
						}
						while(res[1]!=0 && hd<co.distance){
							res[1] &=(res[1]-1);
							hd++;
						}
						if(hd<co.distance){
							num_of_similar_pairs++;
						}
//						cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<" , "<<sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[j] ]->val[0] ]<<"--";
//						cerr<<"hamming distance:"<<hd+1<<"\n";
					}else{
//						cerr<<flag<<"\n";
//						cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<" , "<<sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[j] ]->val[0] ]<<"\n";
					}

				}
			}
			cHead = cnt;
		}
	}
	if(res){
		_mm_free(res);
	}
}
#endif

int multisort::blockSort(int next)
{
	if(key_stack.size()==0 && (next==0 || (co.isDevSort && next == co.startBlc))){
		char_size = ct.num_of_character;
		word_size[SHORT_WORD] = log((double)sq.num_of_seq)/log((double)char_size)-3;
		if(word_size[SHORT_WORD]<=0) word_size[SHORT_WORD]=1;
		word_size[LONG_WORD] = word_size[SHORT_WORD];
		if(bx.short_block_length<word_size[SHORT_WORD]) word_size[SHORT_WORD] = bx.short_block_length;
		if(bx.long_block_length<word_size[LONG_WORD]) word_size[LONG_WORD] = bx.long_block_length;

		num_of_bucket[SHORT_WORD] = pow(char_size,word_size[SHORT_WORD]);
		num_of_bucket[LONG_WORD] = pow(char_size,word_size[LONG_WORD]);

		bucket = (TYPE_INDEX **)malloc(sizeof(TYPE_INDEX*)*2);
		bucket[SHORT_WORD] = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*num_of_bucket[SHORT_WORD]);
		bucket[LONG_WORD] = (TYPE_INDEX *)malloc(sizeof(TYPE_INDEX)*num_of_bucket[LONG_WORD]);

		if(bucket[SHORT_WORD]==NULL || bucket[LONG_WORD] == NULL){
			cerr<<"ERROR: can not allocate memory for bucket\n";
			exit(1);
		}

		bucketID = (int **)malloc(sizeof(int*)*2);
		bucketID[SHORT_WORD] = (int *)malloc(sizeof(int)*num_of_bucket[SHORT_WORD]);
		bucketID[LONG_WORD] = (int *)malloc(sizeof(int)*num_of_bucket[LONG_WORD]);
		if(bucketID[SHORT_WORD]==NULL || bucketID[LONG_WORD] == NULL){
			cerr<<"ERROR: can not allocate memory for bucketID\n";
			exit(1);
		}

		num_of_words[SHORT_WORD] = bx.short_block_length/word_size[SHORT_WORD];
		num_of_words[LONG_WORD] = bx.long_block_length/word_size[LONG_WORD];

		l_word_size[SHORT_WORD] = bx.short_block_length%word_size[SHORT_WORD]; 
		l_word_size[LONG_WORD] = bx.long_block_length%word_size[LONG_WORD];

	}
	if(key_stack.size()==bx.key_size){
		cerr<<"Finish sorting with keys .. ";
		for(int j=0;j<key_stack.size();j++) cerr<<key_stack[j]<<",";
		cerr<<"\n";

		//find similar pair
		findSimilarHDpairs();

		int c=bx.key_size-1;
		int t=0;
		while(c>=0 && key_stack[c]==bx.num_of_blocks - bx.key_size + c){
			key_stack.pop_back();
			c--;
		}

		if(c<0) return(0);

		if(co.isDevSort && c==0 && key_stack[0]==co.startBlc){
			return(0);
		}

		t=key_stack[c]+1;
		key_stack.pop_back();

		blockSort(t);
		return(0);
	}
	key_stack.push_back(next);

	findSimilarRegion(key_stack.size()-1);

	//for debug
	/**
	cerr<<"ptr: "<<ptr<<"\n";
	cerr<<"size :"<<order[ptr]->size<<"\n";
	for(int i=0;i<bx.num_of_box;i++){
		cerr<<"         val:"<<is_region_top[0]->val[i]<<"  "<<is_region_top[1]->val[i]<<"  "<<is_region_top[2]->val[i]<<" \n";
	}
	**/
	blockSort(next+1);
	return(0);
}

#ifdef HPSORT
void multisort::findSimilarRegion(int ptr)
{

//	clock_t start, end;
//	start=clock();

	TYPE_INDEX *work = (TYPE_INDEX *)calloc(sizeof(TYPE_INDEX),bx.num_of_box);

	int block_length;
	int block_type;
	int tmp_word;
	int bid_ptr;

	//ptr-1のクラスタ状況を引き継ぐ。
	resetOrder(ptr);
	resetIsRegionTop(ptr);

	TYPE_INDEX cHead=0, nHead=0;
	if(key_stack[ptr]<bx.num_of_long_blocks){
		block_length = bx.long_block_length;
		block_type = LONG_WORD;
	}else{
		block_length = bx.short_block_length;
		block_type = SHORT_WORD;
	}

	int offset;//blockの先頭
	if(key_stack[ptr]<bx.num_of_long_blocks){
		offset = key_stack[ptr]*bx.long_block_length;
	}else{
		offset = bx.num_of_long_blocks*bx.long_block_length + (key_stack[ptr] - bx.num_of_long_blocks)*bx.short_block_length;
	}

	int col;
	tmp_word=0;
	for(col=0;col<num_of_words[block_type];col++){
		cHead=0;nHead=0;
		for(TYPE_INDEX cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
			if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
				if(cnt-cHead==1){
				}else{
					nHead = cnt;
					memset(bucket[block_type],0x00,sizeof(TYPE_INDEX)*num_of_bucket[block_type]);
					bid_ptr=0;
					//クラスの頻度の数え上げ
					for(TYPE_INDEX i=cHead;i<nHead;i++){
						tmp_word=0;
						for(int cls_cnt=0;cls_cnt<word_size[block_type];cls_cnt++){
							tmp_word = sq.nr_seq[ sq.head[ order[ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
						}
						if(bucket[block_type][tmp_word]==0){
							bucketID[block_type][bid_ptr] = tmp_word;
							bid_ptr++;
						}
						bucket[block_type][ tmp_word ]++;
					}
					for(int i=1;i<bid_ptr;i++){
						bucket[block_type][ bucketID[block_type][i] ] += bucket[block_type][  bucketID[block_type][i-1] ];
					}

					//新しいクラスの追加
					for(int i=0;i<bid_ptr;i++){
						if(bucket[block_type][ bucketID[block_type][i] ] + cHead < bx.num_of_box){
							is_region_top[ptr]->val[ bucket[block_type][ bucketID[block_type][i] ] + cHead ]=REGION_TOP;
						}
					}
					//edit distance用はsortingしない。
					//for(int i=bucket[ct.overlap_character-1];i<bucket[ct.overlap_character];i++){
					//	is_region_top[ptr]->val[ i + cHead ]=REGION_TOP;
					//	cerr<<cHead<<","<< i + cHead<<"\n";
					//}
					//ソート
					for(TYPE_INDEX i=nHead-1;;i--){
						tmp_word=0;
						for(int cls_cnt=0;cls_cnt<word_size[block_type];cls_cnt++){
							tmp_word = sq.nr_seq[ sq.head[ order[ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
						}
						work[ --bucket[block_type][tmp_word] + cHead ] = order[ptr]->val[i];
						if(i==cHead) break;
					}
				}
				cHead = cnt;
			}
		}
		for(TYPE_INDEX i=0;i<bx.num_of_box;i++)
			order[ptr]->val[i]=work[i];

		/**	for debug	**/
		//		cerr<<"--- ptr_pos="<<ptr<<"  key="<<key_stack[ptr]<<"  "<<offset<<"---\n";
		//		for(int i=0;i<bx.num_of_box;i++){

		//							cerr<< order[ptr]->val[i]<<" , ";
		//							cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<"\n";
		//							cerr<<is_region_top[ptr]->val[i]<<" , "<<work[i]<<" , ";
		//							cerr<<order[ptr]->val[i]<<"\n";
		//		}
		//**/

	}
//	cerr<<"test";

	cHead=0;nHead=0;
	for(TYPE_INDEX cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
		if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
			if(cnt-cHead==1){
			}else{
				nHead = cnt;
				memset(bucket[block_type],0x00,sizeof(TYPE_INDEX)*num_of_bucket[block_type]);
				bid_ptr=0;
				//クラスの頻度の数え上げ
				for(TYPE_INDEX i=cHead;i<nHead;i++){
					tmp_word=0;
					for(int cls_cnt=0;cls_cnt<l_word_size[block_type];cls_cnt++){
						tmp_word = sq.nr_seq[ sq.head[ order[ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
					}
					if(bucket[block_type][tmp_word]==0){
						bucketID[block_type][bid_ptr] = tmp_word;
						bid_ptr++;
					}
					bucket[block_type][ tmp_word ]++;
				}
				for(int i=1;i<bid_ptr;i++){
					bucket[block_type][ bucketID[block_type][i] ] += bucket[block_type][  bucketID[block_type][i-1] ];
				}

				//新しいクラスの追加
				for(int i=0;i<bid_ptr;i++){
					if(bucket[block_type][ bucketID[block_type][i] ] + cHead < bx.num_of_box){
						is_region_top[ptr]->val[ bucket[block_type][ bucketID[block_type][i] ] + cHead ]=REGION_TOP;
					}
				}
				//edit distance用はsortingしない。
				//for(int i=bucket[ct.overlap_character-1];i<bucket[ct.overlap_character];i++){
				//	is_region_top[ptr]->val[ i + cHead ]=REGION_TOP;
				//	cerr<<cHead<<","<< i + cHead<<"\n";
				//}
				//ソート
				for(TYPE_INDEX i=nHead-1;;i--){
					tmp_word=0;
					for(int cls_cnt=0;cls_cnt<l_word_size[block_type];cls_cnt++){
						tmp_word = sq.nr_seq[ sq.head[ order[ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
					}
					work[ --bucket[block_type][tmp_word] + cHead ] = order[ptr]->val[i];
					if(i==cHead) break;
				}
			}
			cHead = cnt;
		}
	}
	for(TYPE_INDEX i=0;i<bx.num_of_box;i++)
		order[ptr]->val[i]=work[i];

	if(work){
		free(work);
	}

//	end=clock();
//	cerr<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";

}
#endif
#ifndef HPSORT
void multisort::findSimilarRegion(int ptr)
{
	clock_t start, end;
	start=clock();
    
	int bucket_size;
	TYPE_INDEX *work = (int *)calloc(sizeof(int),bx.num_of_box);
	if(co.distanceType==HAMMING_DISTANCE){
		bucket_size = ct.num_of_character;
	}else{
		bucket_size = ct.num_of_character + 1;
	}
	int *bucket = (int *)calloc(sizeof(int),bucket_size);
	if(bucket==NULL){
		cerr<<"ERROR: can not allocate memory for bucket\n";
		exit(1);
	}

	//ptr-1のクラスタ状況を引き継ぐ。
	resetOrder(ptr);
	resetIsRegionTop(ptr);

	TYPE_INDEX cHead=0, nHead=0;
	int block_length;
	if(key_stack[ptr]<bx.num_of_long_blocks){
		block_length = bx.long_block_length;
	}else{
		block_length = bx.short_block_length;
	}

	int offset;//blockの先頭
	if(key_stack[ptr]<bx.num_of_long_blocks){
		offset = key_stack[ptr]*bx.long_block_length;
	}else{
		offset = bx.num_of_long_blocks*bx.long_block_length + (key_stack[ptr] - bx.num_of_long_blocks)*bx.short_block_length;
	}
	
	for(int col=0;col<block_length;col++){
		cHead=0;nHead=0;
		for(TYPE_INDEX cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
			if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
				if(cnt-cHead==1){
				}else{
					nHead = cnt;
					memset(bucket,0x00,sizeof(int)*bucket_size);
					//クラスの頻度の数え上げ
					for(TYPE_INDEX i=cHead;i<nHead;i++){
						bucket[ sq.nr_seq[ sq.head[ order[ptr]->val[i] ] + offset + col ] ]++;
					}
					for(int i=1;i<bucket_size;i++){
						bucket[i] += bucket[i-1];
					}
					//新しいクラスの追加
					for(int i=0;i<ct.num_of_character;i++){
						if(bucket[i] + cHead<bx.num_of_box)
							is_region_top[ptr]->val[ bucket[i] + cHead ]=REGION_TOP;
					}
					//edit distance用はsortingしない。
					for(int i=bucket[ct.overlap_character-1];i<bucket[ct.overlap_character];i++){
						is_region_top[ptr]->val[ i + cHead ]=REGION_TOP;
						cerr<<cHead<<","<< i + cHead<<"\n";
					}
					//ソート
					for(TYPE_INDEX i=nHead-1;i>=cHead;i--){
						work[ --bucket[ sq.nr_seq[ sq.head[ order[ptr]->val[i] ] + offset + col ]  ] + cHead ] = order[ptr]->val[i];						
					}
				}
				cHead = cnt;
			}
		}
		for(TYPE_INDEX i=0;i<bx.num_of_box;i++)
			order[ptr]->val[i]=work[i];
/**	for debug	**/
//		cerr<<"--- ptr_pos="<<ptr<<"  key="<<key_stack[ptr]<<"  "<<offset<<"---\n";
//		for(int i=0;i<bx.num_of_box;i++){
//			cerr<<is_region_top[ptr]->val[i]<<" , "<<work[i]<<" , ";
//			cerr<<order[ptr]->val[i]<<"\n";
//		}
//**/
	}
//	cerr<<"test";
	if(bucket){
		free(bucket);
	}
	if(work){
		free(work);
	}

	end=clock();
	cerr<<(double)(end-start)/CLOCKS_PER_SEC<<" sec.\n";

}
#endif

