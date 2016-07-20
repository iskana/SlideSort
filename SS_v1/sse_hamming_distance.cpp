/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              sse_hamming_distance.cpp
/***********************************************/

#include"mscls.h"
#include <emmintrin.h>
#ifdef SSE4
#include <nmmintrin.h>
#endif

void bit_stream_debug(BIT_BLOCK bits);


int multisort::calcHammingDist128(TYPE_INDEX id1, TYPE_INDEX id2){

	__m128i seq1, seq2, chd;
	int hdr;
	TYPE_INDEX cPtr;
	BIT_BLOCK *res = (BIT_BLOCK *)_mm_malloc(sizeof(BIT_BLOCK)*2,ALIGNED_BORDER);
	if(res==NULL){
		exit(1);
	}
	//128bitなので。
	BIT_BLOCK *bits = bs->stream[0];

	memset(res,0x00,sizeof(BIT_BLOCK)*2);
	__m128i zero = _mm_set_epi32(0x00000000,0x00000000,0x00000000,0x00000000);
	__m128i one = _mm_set_epi32(0x00000000, 0x00000001,0x00000000, 0x00000001);

	unsigned short s_flg_1;

	seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
	seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
	chd = _mm_xor_si128(seq1,seq2);

	for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
		bits = bs->stream[bit_pos];
		seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
		seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
		seq1 = _mm_xor_si128(seq1,seq2);
		chd = _mm_or_si128(seq1, chd);
	}
	hdr = -1; // hd++と合わせる
#ifdef SSE4
	_mm_store_si128((__m128i *)res,chd);
	hdr += _mm_popcnt_u64(res[0]);
	hdr += _mm_popcnt_u64(res[1]);
#else
	while(hdr<co.distance){
		seq1 = _mm_cmpeq_epi8(chd,zero);
		s_flg_1 = _mm_movemask_epi8 (seq1);
		if(s_flg_1==0xffff){
			break;
		}
		hdr++;
		if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
			hdr++;
		}
		seq1 = _mm_sub_epi64(chd,one);
		chd = _mm_and_si128(chd,seq1);
	}
#endif
	//for count matched unknown_chars
	if(!co.exclude_unknown_character){
		seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
		seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
		chd = _mm_and_si128(seq1,seq2);
		for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
			bits = bs->stream[bit_pos];
			seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
			seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
			seq1 = _mm_and_si128(seq1,seq2);
			chd = _mm_and_si128(seq1, chd);
		}
#ifdef SSE4
		_mm_store_si128((__m128i *)res,chd);
		hdr += _mm_popcnt_u64(res[0]);
		hdr += _mm_popcnt_u64(res[1]);
#else
		while(hdr<co.distance){
			seq1 = _mm_cmpeq_epi8(chd,zero);
			s_flg_1 = _mm_movemask_epi8 (seq1);
			if(s_flg_1==0xffff){
				break;
			}
			hdr++;
			if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
				hdr++;
			}
			seq1 = _mm_sub_epi64(chd,one);
			chd = _mm_and_si128(chd,seq1);
		}
#endif
	}
	hdr++; // -1と合わせる
	if(res){
		_mm_free(res);
	}
	return(hdr);
}


void multisort::findSimilarHDpairsHP128()
{

#ifdef CALL_BACK_FUNC
	TYPE_CHARACTER *rev_seq1, *rev_seq2;
	rev_seq1 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq1,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	rev_seq2 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq2,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
#endif

	TYPE_INDEX cHead=0,nHead=0;
	int ptr = bx.key_size -1;

	__m128i seq1, seq2, chd;
	int flag = NON_REDUNDANT, hd, hdr;
	TYPE_INDEX id1, id2, cPtr, org_seq_id1, org_seq_id2;
	BIT_BLOCK *res = (BIT_BLOCK *)_mm_malloc(sizeof(BIT_BLOCK)*2,ALIGNED_BORDER);
	if(res==NULL){
		exit(1);
	}

	memset(res,0x00,sizeof(BIT_BLOCK)*2);
	__m128i zero = _mm_set_epi32(0x00000000,0x00000000,0x00000000,0x00000000);
	__m128i one = _mm_set_epi32(0x00000000, 0x00000001,0x00000000, 0x00000001);

	unsigned short s_flg_1;

	for(TYPE_INDEX cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
		if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
			nHead = cnt;
			for(TYPE_INDEX i=cHead;i<nHead;i++){
				org_seq_id1 = order[ptr]->val[i];
				id1 = 2*org_seq_id1; 				//128bitなので。
				for(TYPE_INDEX j=i+1;j<nHead;j++){
					flag = NON_REDUNDANT;
					org_seq_id2 = order[ptr]->val[j];
					id2 = 2*org_seq_id2;
#ifndef VCPP
					BIT_BLOCK *bits = bs->stream[0];
#endif

					if(co.isMapMode){
						if(sq.nr_seq_mapID[ org_seq_id1 ]==sq.nr_seq_mapID[ org_seq_id2 ]){
							if(sq.nr_seq_mapID[ org_seq_id1 ]!=BOTH_DATASET){
								goto mapDifChkHd128;
							}
						}
					}
					if(co.isRevComp){
						if( min(sq.revCompIdx[org_seq_id1] , sq.revCompIdx[org_seq_id2])%2 == 1 ){
							goto mapDifChkHd128;
						}
					}

#ifdef VCPP
					//128bitなので。
					BIT_BLOCK *bits = bs->stream[0];
#endif
					seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
					seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
					chd = _mm_xor_si128(seq1,seq2);

					for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
						bits = bs->stream[bit_pos];
						seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
						seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
						seq1 = _mm_xor_si128(seq1,seq2);
						chd = _mm_or_si128(seq1, chd);
					}

					cPtr = ptr-1;
					for(int key = key_stack[ptr]-1 ;key>=0 ;key--){
						if(cPtr>=0 && key == key_stack[cPtr]){
							if(cPtr>0)	cPtr--;
						}else{
							seq1 = _mm_load_si128((__m128i *)&block_mask[key*2]);
							seq1 = _mm_and_si128(chd, seq1);

							seq1 = _mm_cmpeq_epi8(seq1,zero);
							s_flg_1 = _mm_movemask_epi8 (seq1);
							if(s_flg_1==0xffff){
								flag = REDUNDANT;
								break;
							}
							/**
							s_flg_1 = _mm_movemask_epi8 (mhd);
							seq1 = _mm_add_epi8(mhd,fmsk);
							s_flg_2 = _mm_movemask_epi8 (seq1);
							if(s_flg_1==0 && s_flg_2==0){
							flag = REDUNDANT;
							break;
							}
							**/
						}
					}

					//check hamming distance
					if(flag==NON_REDUNDANT){
						num_of_comparison++;
						hd = -1; // hd++と合わせる
#ifdef SSE4
						_mm_store_si128((__m128i *)res,chd);
						hd += _mm_popcnt_u64(res[0]);
						hd += _mm_popcnt_u64(res[1]);
#else
						while(hd<co.distance){
							seq1 = _mm_cmpeq_epi8(chd,zero);
							s_flg_1 = _mm_movemask_epi8 (seq1);
							if(s_flg_1==0xffff){
								break;
							}
							hd++;
							if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
								hd++;
							}
							seq1 = _mm_sub_epi64(chd,one);
							chd = _mm_and_si128(chd,seq1);
						}
#endif
						if(hd<co.distance){
							//for count matched unknown_chars
							if(!co.exclude_unknown_character){
								seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
								seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
								chd = _mm_and_si128(seq1,seq2);
								for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
									bits = bs->stream[bit_pos];
									seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
									seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
									seq1 = _mm_and_si128(seq1,seq2);
									chd = _mm_and_si128(seq1, chd);
								}
#ifdef SSE4
								_mm_store_si128((__m128i *)res,chd);
								hd += _mm_popcnt_u64(res[0]);
								hd += _mm_popcnt_u64(res[1]);
#else
								while(hd<co.distance){
									seq1 = _mm_cmpeq_epi8(chd,zero);
									s_flg_1 = _mm_movemask_epi8 (seq1);
									if(s_flg_1==0xffff){
										break;
									}
									hd++;
									if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
										hd++;
									}
									seq1 = _mm_sub_epi64(chd,one);
									chd = _mm_and_si128(chd,seq1);
								}
#endif
								if(hd>=co.distance){
									goto mapDifChkHd128;
								}
							}


#ifdef CALL_BACK_FUNC
							if(co.outputaln){
								for(int otmp=0; otmp<sq.max_seq_length; otmp++){
									*(rev_seq1++) = ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id1 ] ] ];
									*(rev_seq2++) = ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id2 ] ] ];
								}
								*(rev_seq1)='\0';
								*(rev_seq2)='\0';
								rev_seq1 -= sq.max_seq_length;
								rev_seq2 -= sq.max_seq_length;
							}
#endif
							hd++; // -1と合わせる
							TYPE_INDEX ret_id1, ret_id2;
							for(int t_id1=0;t_id1<sq.toOrgSeqIndex[ org_seq_id1 ]->size;t_id1++){
								ret_id1 = sq.toOrgSeqIndex[ org_seq_id1 ]->val[t_id1];
								for(int t_id2=0;t_id2<sq.toOrgSeqIndex[ org_seq_id2 ]->size;t_id2++){
									ret_id2 = sq.toOrgSeqIndex[ org_seq_id2 ]->val[t_id2];

									if(co.isRevComp){
										if( sq.toOrgSeqIndex[ org_seq_id1 ]->val[0]/2 == sq.toOrgSeqIndex[ org_seq_id2 ]->val[0]/2){
											if(t_id1>=t_id2){
												goto revDifChk;													
											}
										}
										if(co.isOutputBothStrand==false){
											if(abs(sq.revCompSeqDist[org_seq_id2]-hd)<=co.distance){
												TYPE_INDEX rev_id = sq.toSSIndex[ (TYPE_INDEX)(ret_id2/2)*2 + abs(ret_id2%2-1) ];

												/** calc hamming dist**/
												hdr = calcHammingDist128(2*org_seq_id1, 2*rev_id);

												if(hd==hdr){
													if(ret_id2%2==1){
														goto revDifChk;
													}
												}else if(hd>hdr){
													goto revDifChk;													
												}
											}
										}
									}

									if(co.isMapMode){
										if(sq.org_seq_mapID[ ret_id1 ] != sq.org_seq_mapID[ ret_id2 ]){
											num_of_similar_pairs++;
#ifdef CALL_BACK_FUNC
											if(co.outputaln){
												CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,rev_seq1,rev_seq2, hd,sq.max_seq_length);
											}else{
												CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,NULL,NULL,hd,NULL);												
											}
#endif

#ifdef FILE_OUTPUT_PAIRS
											if(co.outputfile){
												fprintf(outfp,"%s,%s,%d\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),hd);
												if(co.outputaln){
													for(int otmp=0; otmp<sq.max_seq_length; otmp++){
														fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id1] ] ]);
													}
													fprintf(outfp,"\n");
													for(int otmp=0; otmp<sq.max_seq_length; otmp++){
														fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id2]  ] ]);
													}
													fprintf(outfp,"\n");
												}
											}
#endif
										}												
									}else{
#ifdef CALL_BACK_FUNC
										if(co.outputaln){
											CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,rev_seq1,rev_seq2, hd,sq.max_seq_length);
										}else{
											CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,NULL,NULL,hd,NULL);												
										}
#endif
#ifdef FILE_OUTPUT_PAIRS
										if(co.outputfile){
											fprintf(outfp,"%s,%s,%d\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),hd);
											if(co.outputaln){
												for(int otmp=0; otmp<sq.max_seq_length; otmp++){
													fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id1] ] ]);
												}
												fprintf(outfp,"\n");
												for(int otmp=0; otmp<sq.max_seq_length; otmp++){
													fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id2]  ] ]);
												}
												fprintf(outfp,"\n");
											}
										}
#endif									
										num_of_similar_pairs++;
									}
revDifChk:;

								}
							}
						}					
						//						cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<" , "<<sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[j] ]->val[0] ]<<" ,";
						//						cerr<<"HD:"<<hd+1<<"\n";
					}else{
						//						cerr<<flag<<"\n";
						//						cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<" , "<<sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[j] ]->val[0] ]<<"\n";
					}
mapDifChkHd128:;
				}
			}
			cHead = cnt;
		}
	}
	if(res){
		_mm_free(res);
	}
}


bitstream::bitstream(int size, TYPE_INDEX length)
{	
	num_of_stream = size;
	stream_length = length;
	TYPE_INDEX bit_block_len = (length + BIT_BLOCK_LENGTH -1 )/BIT_BLOCK_LENGTH;
	bit_block_length=bit_block_len;
	stream = (BIT_BLOCK **)_mm_malloc(sizeof(BIT_BLOCK *)*size, ALIGNED_BORDER);
	if(stream == NULL){
		cerr<<"ERROR: _mm_malloc for stream\n";
		exit(1);
	}
	for(int i=0;i<size;i++){
		stream[i] = (BIT_BLOCK *)_mm_malloc(sizeof(BIT_BLOCK)*bit_block_len, ALIGNED_BORDER);
		if(stream[i] == NULL){
			cerr<<"ERROR: _mm_malloc for stream["<<i<<"]\n";
			exit(1);
		}
		memset(stream[i],0x00,sizeof(BIT_BLOCK)*bit_block_len);
	}
}

void bitstream::free()
{
	for(int i=0;i<num_of_stream;i++){
		if(stream[i])
			_mm_free(stream[i]);
	}
	if(stream)
		_mm_free(stream);
}

void multisort::freeBlockMask()
{
	if(block_mask){
		_mm_free(block_mask);
		block_mask=NULL;
	}
}

void multisort::freeLowerBlockMask()
{
	if(lower_block_mask){
		_mm_free(lower_block_mask);
		lower_block_mask=NULL;
	}
}

void multisort::setBitstream128()
{
	int bits=0;
	while(ct.num_of_character>pow(2,bits)){
		bits++;
	}
	bs = new bitstream(bits,sq.num_of_seq*SIZE_OF_REGISTER);
	// for unknown_character
	int all_one=0;
	for(int i=0;i<bits;i++){
		all_one += pow(2,i);
	}

	BIT_BLOCK bTmp=0;
	BIT_BLOCK sTmp;
	int offset;

	for(int ptr=0;ptr<bs->num_of_stream;ptr++){
		BIT_BLOCK mask = pow(2,ptr); //ないと思うけど・・32 bitを超える場合は無理。
		int block_pos=0;
		bTmp=0;
		for(TYPE_INDEX i=0;i<sq.num_of_seq;i++){
			bTmp=0;
			for(int j=0;j<sq.max_seq_length;j++){
				if(j%BIT_BLOCK_LENGTH==0 && j!=0){
					bs->stream[ptr][block_pos++]=bTmp;
					bTmp=0;
				}
				sTmp = sq.nr_seq[ sq.head[i] + j];
				if(!co.exclude_unknown_character){
					if(sTmp == ct.unknown_character){
						sTmp = all_one;
					}
				}
				offset = j%BIT_BLOCK_LENGTH-ptr;
				if(offset>=0){
					bTmp |= ( (mask & sTmp) << offset  );
				}else{
					bTmp |= ( (mask & sTmp) >> abs(offset) );
				}
			}
			//			cerr<<bTmp<<" ";
			bs->stream[ptr][block_pos++]=bTmp;
			if(sq.max_seq_length<64){
				bs->stream[ptr][block_pos++]=0x0000000000000000;
			}
		}
	}
	/** for debug
	for(int i=0;i<sq.num_of_seq*2;i+=2){
	cerr<<i<<"\n";
	bit_stream_debug(bs->stream[0][i]);
	bit_stream_debug(bs->stream[1][i]);
	cerr<<"\n";
	bit_stream_debug(bs->stream[0][i+1]);
	bit_stream_debug(bs->stream[1][i+1]);
	cerr<<"---\n";
	}
	**/
	//set block_mask
	block_mask = (BIT_BLOCK*)_mm_malloc(sizeof(BIT_BLOCK)*bx.num_of_blocks*2, ALIGNED_BORDER);
	lower_block_mask = (BIT_BLOCK*)_mm_malloc(sizeof(BIT_BLOCK)*bx.num_of_blocks*2, ALIGNED_BORDER);

	memset(block_mask,0x00,sizeof(BIT_BLOCK)*bx.num_of_blocks*2);
	memset(lower_block_mask,0x00,sizeof(BIT_BLOCK)*bx.num_of_blocks*2);

	BIT_BLOCK longB = 0xffffffffffffffff;
	BIT_BLOCK shortB = 0xffffffffffffffff;

	longB = longB >>(BIT_BLOCK_LENGTH-bx.long_block_length);
	shortB = shortB >>(BIT_BLOCK_LENGTH-bx.short_block_length);

	// for each block
	offset=0;
	for(int i=0;i<bx.num_of_blocks*2;i+=2){
		if(i<bx.num_of_long_blocks*2){
			bTmp = longB;
		}else{
			bTmp = shortB;			
		}
		if(offset!=64)
			block_mask[i] = bTmp << offset;
		if( BIT_BLOCK_LENGTH-offset >0 ){
			if(BIT_BLOCK_LENGTH - offset!=64)
				block_mask[i+1] = bTmp >> (BIT_BLOCK_LENGTH - offset);
		}else{
			if(offset-BIT_BLOCK_LENGTH!=64)
				block_mask[i+1] = bTmp << (offset - BIT_BLOCK_LENGTH);			
		}
		if(i<bx.num_of_long_blocks*2){
			offset += bx.long_block_length;
		}else{
			offset += bx.short_block_length;
		}
		//		bit_stream_debug(block_mask[i]);
		//		bit_stream_debug(block_mask[i+1]);
		//		cerr<<"--\n";
	}

	// for lower block mask
	for(int i=0;i<bx.num_of_blocks*2;i+=2){
		for(int j=i; j<bx.num_of_blocks*2;j+=2){
			lower_block_mask[i] |= block_mask[j];
			lower_block_mask[i+1] |= block_mask[j+1];
		}
		//		bit_stream_debug(lower_block_mask[i]);
		//		bit_stream_debug(lower_block_mask[i+1]);
	}

}






/** FOR TEST ONLY**/

void multisort::findSimilarHDpairsHP128_for_kmer_test()
{

#ifdef CALL_BACK_FUNC
	TYPE_CHARACTER *rev_seq1, *rev_seq2;
	rev_seq1 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq1,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	rev_seq2 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq2,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
#endif

	TYPE_INDEX cHead=0,nHead=0;
	int ptr = bx.key_size -1;

	__m128i seq1, seq2, chd, chd1, chd2;
	int flag = NON_REDUNDANT, hd, hdr;
	TYPE_INDEX id1, id2, cPtr, org_seq_id1, org_seq_id2;
	BIT_BLOCK *res = (BIT_BLOCK *)_mm_malloc(sizeof(BIT_BLOCK)*2,ALIGNED_BORDER);
	if(res==NULL){
		exit(1);
	}

	memset(res,0x00,sizeof(BIT_BLOCK)*2);
	__m128i zero = _mm_set_epi32(0x00000000,0x00000000,0x00000000,0x00000000);
	__m128i one = _mm_set_epi32(0x00000000, 0x00000001,0x00000000, 0x00000001);

	unsigned short s_flg_1;

	for(TYPE_INDEX cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
		if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
			nHead = cnt;
			for(TYPE_INDEX i=cHead;i<nHead;i++){
				org_seq_id1 = order[ptr]->val[i];
				id1 = 2*org_seq_id1; 				//128bitなので。
				for(TYPE_INDEX j=i+1;j<nHead;j++){
					flag = NON_REDUNDANT;
					org_seq_id2 = order[ptr]->val[j];
					id2 = 2*org_seq_id2;
#ifndef VCPP
					BIT_BLOCK *bits = bs->stream[0];
#endif

					if(co.isMapMode){
						if(sq.nr_seq_mapID[ org_seq_id1 ]==sq.nr_seq_mapID[ org_seq_id2 ]){
							if(sq.nr_seq_mapID[ org_seq_id1 ]!=BOTH_DATASET){
								goto mapDifChkHd128;
							}
						}
					}
					if(co.isRevComp){
						if( min(sq.revCompIdx[org_seq_id1] , sq.revCompIdx[org_seq_id2])%2 == 1 ){
							goto mapDifChkHd128;
						}
					}

#ifdef VCPP
					//128bitなので。
					BIT_BLOCK *bits = bs->stream[0];
#endif
					seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
					seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
					chd = _mm_xor_si128(seq1,seq2);
					chd1 = chd;

					for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
						bits = bs->stream[bit_pos];
						seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
						seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
						seq1 = _mm_xor_si128(seq1,seq2);
						chd = _mm_or_si128(seq1, chd);
					}
					chd2 = seq1;
					cPtr = ptr-1;
					for(int key = key_stack[ptr]-1 ;key>=0 ;key--){
						if(cPtr>=0 && key == key_stack[cPtr]){
							if(cPtr>0)	cPtr--;
						}else{
							seq1 = _mm_load_si128((__m128i *)&block_mask[key*2]);
							seq1 = _mm_and_si128(chd, seq1);

							seq1 = _mm_cmpeq_epi8(seq1,zero);
							s_flg_1 = _mm_movemask_epi8 (seq1);
							if(s_flg_1==0xffff){
								flag = REDUNDANT;
								break;
							}
							/**
							s_flg_1 = _mm_movemask_epi8 (mhd);
							seq1 = _mm_add_epi8(mhd,fmsk);
							s_flg_2 = _mm_movemask_epi8 (seq1);
							if(s_flg_1==0 && s_flg_2==0){
							flag = REDUNDANT;
							break;
							}
							**/
						}
					}

					//check hamming distance
					if(flag==NON_REDUNDANT){
						num_of_comparison++;
						hd = -1; // hd++と合わせる
#ifdef SSE4
						_mm_store_si128((__m128i *)res,chd);
						hd += _mm_popcnt_u64(res[0]);
						hd += _mm_popcnt_u64(res[1]);

						chd = _mm_and_si128(chd1,chd2);
						_mm_store_si128((__m128i *)res,chd);
						hd += _mm_popcnt_u64(res[0]);
						hd += _mm_popcnt_u64(res[1]);
#else
						while(hd<co.distance){
							seq1 = _mm_cmpeq_epi8(chd,zero);
							s_flg_1 = _mm_movemask_epi8 (seq1);
							if(s_flg_1==0xffff){
								break;
							}
							hd++;
							if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
								hd++;
							}
							seq1 = _mm_sub_epi64(chd,one);
							chd = _mm_and_si128(chd,seq1);
						}
#endif
						if(hd<co.distance){
							//for count matched unknown_chars
							if(!co.exclude_unknown_character){
								seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
								seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
								chd = _mm_and_si128(seq1,seq2);
								for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
									bits = bs->stream[bit_pos];
									seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
									seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
									seq1 = _mm_and_si128(seq1,seq2);
									chd = _mm_and_si128(seq1, chd);
								}
#ifdef SSE4
								_mm_store_si128((__m128i *)res,chd);
								hd += _mm_popcnt_u64(res[0]);
								hd += _mm_popcnt_u64(res[1]);
#else
								while(hd<co.distance){
									seq1 = _mm_cmpeq_epi8(chd,zero);
									s_flg_1 = _mm_movemask_epi8 (seq1);
									if(s_flg_1==0xffff){
										break;
									}
									hd++;
									if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
										hd++;
									}
									seq1 = _mm_sub_epi64(chd,one);
									chd = _mm_and_si128(chd,seq1);
								}
#endif
								if(hd>=co.distance){
									goto mapDifChkHd128;
								}
							}


#ifdef CALL_BACK_FUNC
							if(co.outputaln){
								for(int otmp=0; otmp<sq.max_seq_length; otmp++){
									*(rev_seq1++) = ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id1 ] ] ];
									*(rev_seq2++) = ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id2 ] ] ];
								}
								*(rev_seq1)='\0';
								*(rev_seq2)='\0';
								rev_seq1 -= sq.max_seq_length;
								rev_seq2 -= sq.max_seq_length;
							}
#endif
							hd++; // -1と合わせる
							TYPE_INDEX ret_id1, ret_id2;
							for(int t_id1=0;t_id1<sq.toOrgSeqIndex[ org_seq_id1 ]->size;t_id1++){
								ret_id1 = sq.toOrgSeqIndex[ org_seq_id1 ]->val[t_id1];
								for(int t_id2=0;t_id2<sq.toOrgSeqIndex[ org_seq_id2 ]->size;t_id2++){
									ret_id2 = sq.toOrgSeqIndex[ org_seq_id2 ]->val[t_id2];

									if(co.isRevComp){
										if( sq.toOrgSeqIndex[ org_seq_id1 ]->val[0]/2 == sq.toOrgSeqIndex[ org_seq_id2 ]->val[0]/2){
											if(t_id1>=t_id2){
												goto revDifChk;													
											}
										}
										if(co.isOutputBothStrand==false){
											if(abs(sq.revCompSeqDist[org_seq_id2]-hd)<=co.distance){
												TYPE_INDEX rev_id = sq.toSSIndex[ (TYPE_INDEX)(ret_id2/2)*2 + abs(ret_id2%2-1) ];

												/** calc hamming dist**/
												hdr = calcHammingDist128(2*org_seq_id1, 2*rev_id);

												if(hd==hdr){
													if(ret_id2%2==1){
														goto revDifChk;
													}
												}else if(hd>hdr){
													goto revDifChk;													
												}
											}
										}
									}

									if(co.isMapMode){
										if(sq.org_seq_mapID[ ret_id1 ] != sq.org_seq_mapID[ ret_id2 ]){
											num_of_similar_pairs++;
#ifdef CALL_BACK_FUNC
											if(co.outputaln){
												CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,rev_seq1,rev_seq2, hd,sq.max_seq_length);
											}else{
												CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,NULL,NULL,hd,NULL);												
											}
#endif

#ifdef FILE_OUTPUT_PAIRS
											if(co.outputfile){
												fprintf(outfp,"%s,%s,%d\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),hd);
												if(co.outputaln){
													for(int otmp=0; otmp<sq.max_seq_length; otmp++){
														fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id1] ] ]);
													}
													fprintf(outfp,"\n");
													for(int otmp=0; otmp<sq.max_seq_length; otmp++){
														fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id2]  ] ]);
													}
													fprintf(outfp,"\n");
												}
											}
#endif
										}												
									}else{
#ifdef CALL_BACK_FUNC
										if(co.outputaln){
											CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,rev_seq1,rev_seq2, hd,sq.max_seq_length);
										}else{
											CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2,NULL,NULL,hd,NULL);												
										}
#endif
#ifdef FILE_OUTPUT_PAIRS
										if(co.outputfile){
											fprintf(outfp,"%s,%s,%d\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),hd);
											if(co.outputaln){
												for(int otmp=0; otmp<sq.max_seq_length; otmp++){
													fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id1] ] ]);
												}
												fprintf(outfp,"\n");
												for(int otmp=0; otmp<sq.max_seq_length; otmp++){
													fprintf(outfp,"%c",ct.toChar[ sq.nr_seq[ otmp + sq.head[org_seq_id2]  ] ]);
												}
												fprintf(outfp,"\n");
											}
										}
#endif									
										num_of_similar_pairs++;
									}
revDifChk:;

								}
							}
						}					
						//						cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<" , "<<sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[j] ]->val[0] ]<<" ,";
						//						cerr<<"HD:"<<hd+1<<"\n";
					}else{
						//						cerr<<flag<<"\n";
						//						cerr<< sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[i] ]->val[0] ]<<" , "<<sq.seqName[ sq.toOrgSeqIndex[ order[ptr]->val[j] ]->val[0] ]<<"\n";
					}
mapDifChkHd128:;
				}
			}
			cHead = cnt;
		}
	}
	if(res){
		_mm_free(res);
	}
}

