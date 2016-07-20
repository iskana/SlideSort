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

void bit_stream_debug(BIT_BLOCK bits)
{

	BIT_BLOCK bTmp= bits;
	vector<int> out;

	for(int i=0;i<BIT_BLOCK_LENGTH;i++){
		if(bTmp==0){
			out.push_back(0);
		}else{
			out.push_back(bTmp%2);
		}
		bTmp=bTmp/2;
	}
	for(int i=0;i<out.size();i++){
		cerr<<out[i];
	}
	cerr<<"\n";
}

int multisort::calcHammingDist(TYPE_INDEX id1, TYPE_INDEX id2){

	__m128i seq1, seq2;
	int hd;
	__m128i *chd = (__m128i *)_mm_malloc(sizeof(__m128i)*bs->num_of_regbox_perseq,ALIGNED_BORDER);
	__m128i zero = _mm_set_epi32(0x00000000,0x00000000,0x00000000,0x00000000);
	__m128i one = _mm_set_epi32(0x00000000, 0x00000001,0x00000000, 0x00000001);

#ifdef SSE4
	BIT_BLOCK res[2];
#endif

	unsigned short s_flg_1;
	BIT_BLOCK *bits = bs->stream[0];

	seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
	seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
	chd[0] = _mm_xor_si128(seq1,seq2);
	for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
		seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
		seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
		chd[reg_pos] = _mm_xor_si128(seq1,seq2);
	}
	for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
		bits = bs->stream[bit_pos];
		// for reg_pos=0
		seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
		seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
		seq1 = _mm_xor_si128(seq1,seq2);
		chd[0] = _mm_or_si128(seq1, chd[0]);
		for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
			seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
			seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
			seq1 = _mm_xor_si128(seq1,seq2);
			chd[reg_pos] = _mm_or_si128(seq1, chd[reg_pos]);

		}
	}

	hd = -1;
	for(int reg_pos=0;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
#ifdef SSE4
		_mm_store_si128((__m128i *)res,chd[reg_pos]);
		hd += _mm_popcnt_u64(res[0]);
		hd += _mm_popcnt_u64(res[1]);
		if(hd>=co.distance){
			break;
		}
#else
		while(hd<co.distance){
			seq1 = _mm_cmpeq_epi8(chd[reg_pos],zero);
			s_flg_1 = _mm_movemask_epi8 (seq1);
			if(s_flg_1==0xffff){
				break;
			}
			hd++;
			if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
				hd++;
			}
			seq1 = _mm_sub_epi64(chd[reg_pos],one);
			chd[reg_pos] = _mm_and_si128(chd[reg_pos],seq1);
		}
#endif
	}
	//for count matched unknown_chars
	if(!co.exclude_unknown_character){
		seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
		seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
		chd[0] = _mm_and_si128(seq1,seq2);
		for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
			seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
			seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
			chd[reg_pos] = _mm_and_si128(seq1,seq2);
		}
		for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
			bits = bs->stream[bit_pos];
			// for reg_pos=0
			seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
			seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
			seq1 = _mm_and_si128(seq1,seq2);
			chd[0] = _mm_and_si128(seq1, chd[0]);
			for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
				seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
				seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
				seq1 = _mm_and_si128(seq1,seq2);
				chd[reg_pos] = _mm_and_si128(seq1, chd[reg_pos]);
			}
		}
		for(int reg_pos=0;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
#ifdef SSE4
			_mm_store_si128((__m128i *)res,chd[reg_pos]);
			hd += _mm_popcnt_u64(res[0]);
			hd += _mm_popcnt_u64(res[1]);
#else
			while(hd<co.distance){
				seq1 = _mm_cmpeq_epi8(chd[reg_pos],zero);
				s_flg_1 = _mm_movemask_epi8 (seq1);
				if(s_flg_1==0xffff){
					break;
				}
				hd++;
				if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
					hd++;
				}
				seq1 = _mm_sub_epi64(chd[reg_pos],one);
				chd[reg_pos] = _mm_and_si128(chd[reg_pos],seq1);
			}
#endif
		}
	}
	hd++;
	_mm_free(chd);
	return(hd);
}


void multisort::findSimilarHDpairsHP()
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

	__m128i seq1, seq2;
	int flag = NON_REDUNDANT, hd, hdr;
	TYPE_INDEX id1, id2, cPtr, org_seq_id1, org_seq_id2;
	__m128i *chd = (__m128i *)_mm_malloc(sizeof(__m128i)*bs->num_of_regbox_perseq,ALIGNED_BORDER);

	__m128i zero = _mm_set_epi32(0x00000000,0x00000000,0x00000000,0x00000000);
	__m128i one = _mm_set_epi32(0x00000000, 0x00000001,0x00000000, 0x00000001);

#ifdef SSE4
	BIT_BLOCK res[2];
#endif

	/**
	__m128i lo_one = _mm_set_epi32(0x00000000, 0x00000000,0x00000000, 0x00000001);
	__m128i up_one = _mm_set_epi32(0x00000000, 0x00000001,0x00000000, 0x00000000);
	__m128i up_f = _mm_set_epi32(0xffffffff, 0xffffffff,0x00000000, 0x00000000);
	__m128i lo_f = _mm_set_epi32(0x00000000, 0x00000000,0xffffffff, 0xffffffff);
	**/
	unsigned short s_flg_1;

	//	BIT_BLOCK *dbg = (BIT_BLOCK *)_mm_malloc(sizeof(BIT_BLOCK)*2,ALIGNED_BORDER);

	for(TYPE_INDEX cnt=1;cnt<=is_region_top[ptr]->size;cnt++){
		if(cnt==is_region_top[ptr]->size || is_region_top[ptr]->val[cnt]==REGION_TOP){
			nHead = cnt;
			for(TYPE_INDEX i=cHead;i<nHead;i++){
				org_seq_id1 = order[ptr]->val[i];
				id1 = (org_seq_id1)*bs->num_of_bb_per_seq;
				for(TYPE_INDEX j=i+1;j<nHead;j++){
					flag = NON_REDUNDANT;
					org_seq_id2 = order[ptr]->val[j];
					id2 = (org_seq_id2)*bs->num_of_bb_per_seq;

#ifndef VCPP
					BIT_BLOCK *bits = bs->stream[0];
#endif

					if(co.isMapMode){
						if(sq.nr_seq_mapID[ org_seq_id1 ]==sq.nr_seq_mapID[ org_seq_id2 ]){
							if(sq.nr_seq_mapID[ org_seq_id1 ]!=BOTH_DATASET){
								goto mapDifChkHd;
							}
						}
					}
					if(co.isRevComp){
						if( min(sq.revCompIdx[org_seq_id1] , sq.revCompIdx[org_seq_id2])%2 == 1 ){
							goto mapDifChkHd;
						}
					}
#ifdef VCPP
					BIT_BLOCK *bits = bs->stream[0];
#endif
					seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
					seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
					chd[0] = _mm_xor_si128(seq1,seq2);
					for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
						seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
						seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
						chd[reg_pos] = _mm_xor_si128(seq1,seq2);
					}
					for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
						bits = bs->stream[bit_pos];
						// for reg_pos=0
						seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
						seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
						seq1 = _mm_xor_si128(seq1,seq2);
						chd[0] = _mm_or_si128(seq1, chd[0]);
						for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
							seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
							seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
							seq1 = _mm_xor_si128(seq1,seq2);
							chd[reg_pos] = _mm_or_si128(seq1, chd[reg_pos]);

						}
					}

					//	for(int reg_pos=0;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
					//		_mm_store_si128((__m128i *)dbg,chd[reg_pos]);
					//		bit_stream_debug(dbg[0]);
					//		bit_stream_debug(dbg[1]);
					//	}

					//duplication check
					cPtr = ptr-1;
					for(int key = key_stack[ptr]-1 ;key>=0 ;key--){
						if(cPtr>=0 && key == key_stack[cPtr]){
							if(cPtr>0)	cPtr--;
						}else{
							flag = REDUNDANT;
							for(int reg_pos = bs->brmap[key]->val[0] ;reg_pos < bs->brmap[key]->val[1]+1;reg_pos++){									
								seq1 = _mm_load_si128((__m128i *)&block_mask[key*bs->num_of_bb_per_seq + reg_pos*bs->num_of_bb_per_reg]);
								seq1 = _mm_and_si128(chd[reg_pos], seq1);

								seq1 = _mm_cmpeq_epi8(seq1,zero);
								s_flg_1 = _mm_movemask_epi8 (seq1);

								if(s_flg_1!=0xffff){
									flag = NON_REDUNDANT;
									break;
								}
							}
							if(flag==REDUNDANT){
								break;
							}
						}
					}

					//check hamming distance

					if(flag==NON_REDUNDANT){
						num_of_comparison++;
						hd = -1;
						for(int reg_pos=0;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
#ifdef SSE4
							_mm_store_si128((__m128i *)res,chd[reg_pos]);
							hd += _mm_popcnt_u64(res[0]);
							hd += _mm_popcnt_u64(res[1]);
							if(hd>=co.distance){
								break;
							}
#else
							while(hd<co.distance){
								seq1 = _mm_cmpeq_epi8(chd[reg_pos],zero);
								s_flg_1 = _mm_movemask_epi8 (seq1);
								if(s_flg_1==0xffff){
									break;
								}
								hd++;
								if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
									hd++;
								}
								seq1 = _mm_sub_epi64(chd[reg_pos],one);
								chd[reg_pos] = _mm_and_si128(chd[reg_pos],seq1);
							}
#endif
						}
						if(hd<co.distance){
							//for count matched unknown_chars
							if(!co.exclude_unknown_character){
								seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
								seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
								chd[0] = _mm_and_si128(seq1,seq2);
								for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
									seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
									seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
									chd[reg_pos] = _mm_and_si128(seq1,seq2);
								}
								for(int bit_pos=1;bit_pos<bs->num_of_stream;bit_pos++){
									bits = bs->stream[bit_pos];
									// for reg_pos=0
									seq1 = _mm_load_si128((__m128i *)&(bits[ id1 ]));
									seq2 = _mm_load_si128((__m128i *)&(bits[ id2 ]));
									seq1 = _mm_and_si128(seq1,seq2);
									chd[0] = _mm_and_si128(seq1, chd[0]);
									for(int reg_pos=1;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
										seq1 = _mm_load_si128((__m128i *)&(bits[ id1+reg_pos*bs->num_of_bb_per_reg ]));
										seq2 = _mm_load_si128((__m128i *)&(bits[ id2+reg_pos*bs->num_of_bb_per_reg ]));
										seq1 = _mm_and_si128(seq1,seq2);
										chd[reg_pos] = _mm_and_si128(seq1, chd[reg_pos]);
									}
								}
								for(int reg_pos=0;reg_pos<bs->num_of_regbox_perseq;reg_pos++){
#ifdef SSE4
									_mm_store_si128((__m128i *)res,chd[reg_pos]);
									hd += _mm_popcnt_u64(res[0]);
									hd += _mm_popcnt_u64(res[1]);
#else
									while(hd<co.distance){
										seq1 = _mm_cmpeq_epi8(chd[reg_pos],zero);
										s_flg_1 = _mm_movemask_epi8 (seq1);
										if(s_flg_1==0xffff){
											break;
										}
										hd++;
										if((s_flg_1&0x00ff) != 0x00ff && (s_flg_1&0xff00) != 0xff00){
											hd++;
										}
										seq1 = _mm_sub_epi64(chd[reg_pos],one);
										chd[reg_pos] = _mm_and_si128(chd[reg_pos],seq1);
									}
#endif
									if(hd>=co.distance){
										goto mapDifChkHd;
									}
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
							hd++;
							TYPE_INDEX ret_id1, ret_id2;
							for(TYPE_INDEX t_id1=0;t_id1<sq.toOrgSeqIndex[ org_seq_id1 ]->size;t_id1++){
								ret_id1 = sq.toOrgSeqIndex[ org_seq_id1 ]->val[t_id1];
								for(TYPE_INDEX t_id2=0;t_id2<sq.toOrgSeqIndex[ org_seq_id2 ]->size;t_id2++){
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
												hdr = calcHammingDist((org_seq_id1)*bs->num_of_bb_per_seq, (rev_id)*bs->num_of_bb_per_seq);
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
mapDifChkHd:;
				}
			}
			cHead = cnt;
		}
	}

	_mm_free(chd);
}

void multisort::setBitstream()
{
	int bits=0;
	while(ct.num_of_character>pow(2,bits)){
		bits++;
	}
	int num_rb_ps=sq.max_seq_length/SIZE_OF_REGISTER;
	if(sq.max_seq_length%SIZE_OF_REGISTER>0){
		num_rb_ps++;
	}
	// for unknown_character
	int all_one=0;
	for(int i=0;i<bits;i++){
		all_one += pow(2,i);
	}
	bs = new bitstream(bits,sq.num_of_seq*SIZE_OF_REGISTER*num_rb_ps);
	bs->num_of_regbox_perseq = num_rb_ps;
	bs->num_of_bb_per_reg = SIZE_OF_REGISTER/BIT_BLOCK_LENGTH;
	bs->num_of_bb_per_seq = bs->num_of_bb_per_reg*bs->num_of_regbox_perseq;

	BIT_BLOCK bTmp=0;
	BIT_BLOCK sTmp;
	int offset;
	// stream[num]: bb1, bb2, ... , bbn
	//              ------ bbps -------
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
			while(block_pos%bs->num_of_bb_per_reg!=0){
				bs->stream[ptr][block_pos++] = BIT_BLOCK_0;		
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
	block_mask = (BIT_BLOCK*)_mm_malloc(sizeof(BIT_BLOCK)*bx.num_of_blocks*bs->num_of_bb_per_seq, ALIGNED_BORDER);
	memset(block_mask,0x00,sizeof(BIT_BLOCK)*bx.num_of_blocks*bs->num_of_bb_per_seq);

	int nlbb = bx.long_block_length/BIT_BLOCK_LENGTH;
	if(bx.long_block_length%BIT_BLOCK_LENGTH>0){
		nlbb++;
	}
	int nsbb;
	if(bx.long_block_length!=bx.short_block_length){
		nsbb = bx.short_block_length/BIT_BLOCK_LENGTH;
		if(bx.short_block_length%BIT_BLOCK_LENGTH>0){
			nsbb++;
		}
	}else{
		nsbb = nlbb;
	}

	BIT_BLOCK longB;
	BIT_BLOCK shortB;

	int olbb = BIT_BLOCK_LENGTH*nlbb - bx.long_block_length;
	int osbb = BIT_BLOCK_LENGTH*nsbb - bx.short_block_length; 

	// set 0xff.. to the region of each block
	int bptr = 0;
	int head = 0;
	//	cerr<<bx.long_block_length<<"   "<<nlbb<<"\n";
	for(int i=0;i<bx.num_of_blocks;i++){
		if(i<bx.num_of_long_blocks){
			for(int j=0;j<nlbb;j++){
				if(j==0){
					if(head<64)
						block_mask[bptr] = BIT_BLOCK_1<<head;					
					//					cerr<<bptr<<"    j:"<<j<<"   "<<head<<" shifted. \n";
				}else{
					block_mask[bptr] = BIT_BLOCK_1;
				}
				if(j==nlbb-1){
					if(olbb-head>=0){
						//						cerr<<olbb-head<<" <-\n";
						block_mask[bptr] &= BIT_BLOCK_1>>(olbb-head);
						head = BIT_BLOCK_LENGTH-(olbb-head);
						bptr += bs->num_of_bb_per_seq;
						//						cerr<<head<<" next->\n";
						break;
					}else{
						bptr++;
						//						cerr<<BIT_BLOCK_LENGTH+(olbb-head)<<" <- after+\n";
						block_mask[bptr] = BIT_BLOCK_1>>BIT_BLOCK_LENGTH+(olbb-head);
						head = head-olbb;
						bptr += bs->num_of_bb_per_seq;
						//						cerr<<head<<" next->\n";
						break;
					}
				}
				//				cerr<<"+";
				bptr++;
			}
		}else{
			for(int j=0;j<nsbb;j++){
				if(j==0){
					block_mask[bptr] = BIT_BLOCK_1<<head;
					//					cerr<<bptr<<"    j:"<<j<<"   "<<head<<" shifted. \n";
				}else{
					block_mask[bptr] = BIT_BLOCK_1;
				}
				if(j==nsbb-1){
					if(osbb-head>=0){
						//						cerr<<osbb-head<<" <-\n";
						block_mask[bptr] &= BIT_BLOCK_1>>(osbb-head);						
						head = BIT_BLOCK_LENGTH-(osbb-head);
						bptr += bs->num_of_bb_per_seq;
						//						cerr<<head<<" next->\n";
						break;
					}else{
						bptr++;
						//						cerr<<BIT_BLOCK_LENGTH+(osbb-head)<<" <- after+\n";
						block_mask[bptr] = BIT_BLOCK_1>>BIT_BLOCK_LENGTH+(osbb-head);
						head = head-osbb;
						bptr += bs->num_of_bb_per_seq;
						//						cerr<<head<<" next->\n";
						break;
					}
				}
				//				cerr<<"+";
				bptr++;
			}		
		}
	}

	bs->brmap = (varInt **)malloc(sizeof(varInt *)*bx.num_of_blocks);
	int bstart=0, bstop=0;
	for(int i=0;i<bx.num_of_blocks;i++){
		if(i<bx.num_of_long_blocks)	{
			bstart = i*bx.long_block_length/SIZE_OF_REGISTER;
			bstop = ((i+1)*bx.long_block_length-1)/SIZE_OF_REGISTER;
		}else{
			bstart = (bx.num_of_long_blocks*bx.long_block_length + (i-bx.num_of_long_blocks)*bx.short_block_length)/SIZE_OF_REGISTER;
			bstop = (bx.num_of_long_blocks*bx.long_block_length + ((i+1-bx.num_of_long_blocks)*bx.short_block_length)-1)/SIZE_OF_REGISTER;
		}

		bs->brmap[i] = new varInt(2);
		bs->brmap[i]->val[0] = bstart;
		bs->brmap[i]->val[1] = bstop;
		//		cerr<<bstart<<" to "<<bstop<<"\n";
	}

	//debug
	/**
	for(int i=0;i<bx.num_of_blocks;i++){
	for(int j=0;j<bs->num_of_bb_per_seq;j++){
	bit_stream_debug(block_mask[j+i*bs->num_of_bb_per_seq]);
	}
	cerr<<"\n-------\n";
	}
	cerr<<"DEBUG\n";
	**/
}


