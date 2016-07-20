/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              edit_distance2.cpp
/***********************************************/

#include"mscls.h"

#if 0
//#ifdef HPSORT
int multisort::blockSort_ED()
{
	/** print
	for(int i=0;i<sq.seq_length;i++){
		cerr<<ct.toChar[sq.nr_seq[i]];
	}
**/

	TYPE_CHARACTER *out_seq1, *out_seq2;
	out_seq1 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
	out_seq2 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
	memset(out_seq1,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
	memset(out_seq2,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
#ifdef CALL_BACK_FUNC
	TYPE_CHARACTER *rev_seq1, *rev_seq2;
	rev_seq1 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	rev_seq2 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq1,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq2,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
#endif

	// for sort
//	int bucket_size;
	TYPE_INDEX *work_order, *work_is_offset_zero;
	char_size = ct.num_of_character + 1;

	word_size[SHORT_WORD] = log((double)sq.num_of_seq)/log((double)char_size)-2;
	if(word_size[SHORT_WORD]<=0) word_size[SHORT_WORD]=1;
	word_size[LONG_WORD] = word_size[SHORT_WORD];
	if(bx.short_block_length<word_size[SHORT_WORD]) word_size[SHORT_WORD] = bx.short_block_length;
	if(bx.long_block_length<word_size[LONG_WORD]) word_size[LONG_WORD] = bx.long_block_length;

	num_of_bucket[SHORT_WORD] = pow(char_size,word_size[SHORT_WORD]);
	num_of_bucket[LONG_WORD] = pow(char_size,word_size[LONG_WORD]);

	bucket = (int **)malloc(sizeof(int*)*2);
	bucket[SHORT_WORD] = (int *)malloc(sizeof(int)*num_of_bucket[SHORT_WORD]);
	bucket[LONG_WORD] = (int *)malloc(sizeof(int)*num_of_bucket[LONG_WORD]);

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

	int bid_ptr;
	int block_length;
	int block_type;
	int tmp_word;


	//for 
	int *block_stack = (int *)calloc(bx.key_size,sizeof(int));
	int stack_ptr = 0;
	int cHead = 0, pHead=0, cnt, prev, next_block, copy_flag;
	int region_len = bx.num_of_box;
	int min_zero_id_flag;

	copy_flag=0;
	key_stack.push_back(0);
	order[0] = new varInt(bx.num_of_box);
	is_region_top[0] = new varInt(bx.num_of_box);
	is_offset_zero[0] = new varInt(bx.num_of_box);
	resetOrder(0);
	resetIsRegionTop(0);
	for(int i=0;i<bx.num_of_box;i++){
		if(i%bx.box_set_size == bx.box_set_center){
			is_offset_zero[0]->val[i]=OFFSET_ZERO;
		}
	}
	//cpy
	order_cpy[0] = new varInt(bx.num_of_box);
	is_offset_zero_cpy[0] = new varInt(bx.num_of_box);
	memcpy(order_cpy[0]->val, order[0]->val, sizeof(TYPE_INDEX)*order[0]->size);
	memcpy(is_offset_zero_cpy[0]->val, is_offset_zero[0]->val, sizeof(TYPE_INDEX)*is_offset_zero[0]->size);

	cerr<<"blockSort of ED\n";
//	cerr<<bx.box_set_center<<"\n";
	while(stack_ptr >= 0){
		if(stack_ptr == bx.key_size){
			stack_ptr--;
			cHead = block_stack[stack_ptr];
			key_stack.pop_back();

			region_len /= bx.box_set_size;
			TYPE_INDEX *pair_list = (TYPE_INDEX *)calloc(sizeof(TYPE_INDEX),region_len);
			TYPE_INDEX *ioz_list = (TYPE_INDEX *)calloc(sizeof(TYPE_INDEX),region_len);
			int cpy=0;

			for(int cnt=pHead;cnt<cHead;cnt++){
				if(is_region_top[stack_ptr]->val[cnt]>=COPY_TO_NEXT_LIST){
					pair_list[cpy] = order[stack_ptr]->val[cnt];
					if(is_region_top[stack_ptr]->val[cnt]==OFFSET_ZERO){
						//offsetzeroのものは必ずしも、copy_to_next_LISTになってないことに注意。
						ioz_list[cpy] = OFFSET_ZERO;
						pair_list[cpy] = (order[stack_ptr]->val[cnt]/bx.box_set_size)*bx.box_set_size+bx.box_set_center;
					}
					cpy++;
				}
			}
			/** print
			cerr<<"pair candidate\n";
			for(int cnt=pHead;cnt<cHead;cnt++){
				cerr<<order[stack_ptr]->val[cnt]<<" , ";
			}
			cerr<<"\n";
			for(int cnt=pHead;cnt<cHead;cnt++){
				cerr<<is_offset_zero[stack_ptr]->val[cnt]<<" , ";
			}
			cerr<<"\n";

			for(int cnt=0;cnt<region_len;cnt++){
				cerr<<pair_list[cnt]<<" ("<<pair_list[cnt]/bx.box_set_size<<") , ";
			}
			cerr<<"\n";
			for(int cnt=0;cnt<region_len;cnt++){
				cerr<<ioz_list[cnt]<<",\t";
			}
			cerr<<"\n";
**/

			/**
			cerr<<"keys = ";
			for(int i=0;i<key_stack.size();i++){
				cerr<<key_stack[i]<<",";
			}
			cerr<<"\n";
**/
//			cerr<<"calc candidate\n";
			int small_id;
			int big_id;
			int org_seq_id1;
			int org_seq_id2;

			for(int box1=0;box1<region_len-1;box1++){
				for(int box2=box1+1;box2<region_len;box2++){
					org_seq_id1 = (pair_list[box1]/bx.box_set_size);
					org_seq_id2 = (pair_list[box2]/bx.box_set_size);

					//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[0] ]<<"\t";
					//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[0] ]<<"\n";

					if(co.isMapMode){
						if(sq.nr_seq_mapID[org_seq_id1]==sq.nr_seq_mapID[org_seq_id2]){
							if(sq.nr_seq_mapID[org_seq_id1]!=BOTH_DATASET){
								goto mapDifChk;
							}
						}
					}

					if(org_seq_id1 != org_seq_id2){
						if(pair_list[box1]<pair_list[box2]){
							small_id = box1;
							big_id = box2;
						}else{
							small_id = box2;
							big_id = box1;
						}
						if(ioz_list[ small_id ]==OFFSET_ZERO){
							int flag = NON_REDUNDANT;
							/** print
							cerr<<pair_list[box1]<<" ("<<pair_list[box1]/bx.box_set_size<<") , ";
							cerr<<pair_list[box2]<<" ("<<pair_list[box2]/bx.box_set_size<<")\n";
							**/
							//duplication check
							int key_ptr = stack_ptr;
							int key = key_stack[ key_ptr ] - 1;
							key_ptr--;
							/**
							cerr<<"key="<<key<<"\n";
							cerr<<"key_stack"<<key_stack[key_ptr]<<"\n";
							cerr<<"key_ptr"<<key_ptr<<"\n";
							**/
							while(key>=0){
								if(key==key_stack[key_ptr] && key_ptr>=0){
//									cerr<<"key=_stack"<<key_stack[key_ptr]<<"\n";
									key_ptr--;
								}else{
									int small_seq_id = (pair_list[small_id]/bx.box_set_size)*bx.box_set_size + bx.box_set_center;
									int big_seq_id = (pair_list[big_id]/bx.box_set_size)*bx.box_set_size;
									for(int slide=0;slide<bx.box_set_size;slide++){
										int block_len;
										if(key<bx.num_of_long_blocks){
											block_len = bx.long_block_length;
										}else{
											block_len = bx.short_block_length;
										}
										int i;
										int small_c,big_c;
										for(i=0;i<block_len;i++){
											small_c = sq.nr_seq[ bx.head[ small_seq_id ] + bx.block_offset[key] + i ];
											big_c = sq.nr_seq[ bx.head[ big_seq_id ] + bx.block_offset[key]+slide + i ];
											if(small_c == big_c){
											//	cerr<< ct.toChar[ sq.nr_seq[ bx.head[ small_seq_id ] + bx.block_offset[key] + i ] ];
											}else{
												break;
											}
										}
										if(i==block_len){
											flag=REDUNDANT;
											/** print
											cerr<<"key "<<key<<"\n";
											cerr<<"slide "<<slide<<"\n";
											cerr<<bx.block_offset[key]<<"\n";
											cerr<<bx.head[ small_seq_id ] + bx.block_offset[key]<<" , "<<bx.head[ big_seq_id ] + bx.block_offset[key]+slide<<"\n";
											**/
											break;
										}
									}
									if(flag==REDUNDANT){
										break;
									}
								}
								key--;
							}

							if(flag==NON_REDUNDANT){
								num_of_comparison++;
								int ed;
								if((ed=calcEditDistance(org_seq_id1, org_seq_id2, out_seq1, out_seq2))<=co.distance){
									/** Print IDs**/
									for(int id1=0;id1<sq.toOrgSeqIndex[ org_seq_id1 ]->size;id1++){
										for(int id2=0;id2<sq.toOrgSeqIndex[ org_seq_id2 ]->size;id2++){
											if(co.isMapMode){
												if(sq.org_seq_mapID[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1] ] != sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2]){
													num_of_similar_pairs++;
#ifdef CALL_BACK_FUNC
													if(co.outputaln){													
														for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
															*(rev_seq1++) = ct.toChar[ out_seq1[otmp] ];
														}
														*(rev_seq1)='\0';
														rev_seq1 -= sq.max_seq_length+co.distance;
														for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
															*(rev_seq2++) = ct.toChar[ out_seq2[otmp] ];
														}
														*(rev_seq2)='\0';
														rev_seq2 -= sq.max_seq_length+co.distance;

														CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1], sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2],rev_seq1,rev_seq2,ed);
													}else{
														CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1], sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2],NULL,NULL,ed);													
													}
#endif 
#ifdef FILE_OUTPUT_PAIRS
													if(co.outputfile){
														fprintf(outfp,"%s,%s,%d\n",sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1] ].c_str(),sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2] ].c_str(),ed);
														if(co.outputaln){
															for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
																fprintf(outfp,"%c",ct.toChar[ out_seq1[otmp] ]);
															}
															fprintf(outfp,"\n");
															for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
																fprintf(outfp,"%c",ct.toChar[ out_seq2[otmp] ]);
															}
															fprintf(outfp,"\n");
														}
													}
#endif
													//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1] ]<<"\t";
													//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2] ]<<"\t"<<ed<<"\n";
												}												
											}else{
												num_of_similar_pairs++;
#ifdef FILE_OUTPUT_PAIRS
												if(co.outputfile){
													fprintf(outfp,"%s,%s,%d\n",sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1] ].c_str(),sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2] ].c_str(),ed);
													if(co.outputaln){
														for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
															fprintf(outfp,"%c",ct.toChar[ out_seq1[otmp] ]);
														}
														fprintf(outfp,"\n");
														for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
															fprintf(outfp,"%c",ct.toChar[ out_seq2[otmp] ]);
														}
														fprintf(outfp,"\n");
													}
												}
#endif
#ifdef CALL_BACK_FUNC
												if(co.outputaln){													
													for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
														*(rev_seq1++) = ct.toChar[ out_seq1[otmp] ];
													}
													*(rev_seq1)='\0';
													rev_seq1 -= sq.max_seq_length+co.distance;
													for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
														*(rev_seq2++) = ct.toChar[ out_seq2[otmp] ];
													}
													*(rev_seq2)='\0';
													rev_seq2 -= sq.max_seq_length+co.distance;
													CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1], sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2],rev_seq1,rev_seq2,ed);
												}else{
													CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1], sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2],NULL,NULL,ed);
												}

#endif 
												//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1] ]<<"\t";
												//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2] ]<<"\t"<<ed<<"\n";
											}
										}
									}
//									cerr<<"=ED=\n";
/**
									cerr<<pair_list[small_id]<<" ("<<pair_list[small_id]/bx.box_set_size<<" : "<<pair_list[small_id]-(pair_list[small_id]/bx.box_set_size)*bx.box_set_size-bx.box_set_center<<") ,";
									cerr<<pair_list[big_id]<<" ("<<pair_list[big_id]/bx.box_set_size<<" : "<<pair_list[big_id]-(pair_list[big_id]/bx.box_set_size)*bx.box_set_size-bx.box_set_center<<")\n";

									for(int i=0;i<bx.box_length;i++){										
										cerr<<ct.toChar[ sq.nr_seq[ bx.head[ (pair_list[small_id]/bx.box_set_size)*bx.box_set_size + bx.box_set_center ] +i ] ];
									}
									cerr<<"\n";
									for(int i=0;i<bx.box_length;i++){
										cerr<<ct.toChar[ sq.nr_seq[ bx.head[ (pair_list[big_id]/bx.box_set_size)*bx.box_set_size + bx.box_set_center ] +i] ];
									}
									cerr<<"\n--\n";
									**/
								}
							}
						}
					}
mapDifChk:;
				}
			}
//			cerr<<num_of_comparison<<"\n";

			free(pair_list);
			free(ioz_list);

		}
		if(cHead==0){
			if(copy_flag){
				order[stack_ptr] = new varInt(region_len);
				// 原因不明↓
				is_region_top[stack_ptr] = new varInt(region_len+1);
				is_region_top[stack_ptr]->size = region_len;
				is_offset_zero[stack_ptr] = new varInt(region_len);
				order_cpy[stack_ptr] = new varInt(region_len);
				is_offset_zero_cpy[stack_ptr] = new varInt(region_len);

				//copy
				int cnt=pHead;
				int cpy=0;
				while(1){
					if(is_region_top[stack_ptr-1]->val[cnt] >=COPY_TO_NEXT_LIST){
					//	cerr<<"select "<<order[stack_ptr-1]->val[cnt]<<"\n";
						if(is_region_top[stack_ptr-1]->val[cnt] == OFFSET_ZERO){
							is_offset_zero[stack_ptr]->val[cpy + bx.box_set_center] = OFFSET_ZERO;
						}
						for(int i=0;i<bx.box_set_size;i++){
							order[stack_ptr]->val[cpy++]
							= (order[stack_ptr-1]->val[cnt]/bx.box_set_size)*bx.box_set_size + i;	
						}
					}
					if(cpy>=region_len){
						break;
					}
					cnt++;
				}
				memcpy(order_cpy[stack_ptr]->val, order[stack_ptr]->val, sizeof(TYPE_INDEX)*order[stack_ptr]->size);
				memcpy(is_offset_zero_cpy[stack_ptr]->val, is_offset_zero[stack_ptr]->val, sizeof(TYPE_INDEX)*is_offset_zero[stack_ptr]->size);
//				cerr<<"malloc:"<<stack_ptr<<"\n";
			}else{
				//orderとiozの修正の必要あり。
				memcpy(order[stack_ptr]->val, order_cpy[stack_ptr]->val, sizeof(TYPE_INDEX)*order_cpy[stack_ptr]->size);
				memcpy(is_offset_zero[stack_ptr]->val, is_offset_zero_cpy[stack_ptr]->val, sizeof(TYPE_INDEX)*is_offset_zero_cpy[stack_ptr]->size);
//				cerr<<"cpy:"<<stack_ptr<<"\n";				
			}

			/****   SORT *****/
			cerr<<"sort "<<key_stack[stack_ptr]<<"...\n";
			memset(is_region_top[stack_ptr]->val,0x00,sizeof(TYPE_INDEX)*is_region_top[stack_ptr]->size);
			is_region_top[stack_ptr]->val[0] = REGION_TOP;

			work_order = (TYPE_INDEX *)calloc(order[stack_ptr]->size,sizeof(TYPE_INDEX));
			work_is_offset_zero = (TYPE_INDEX *)calloc(is_offset_zero[stack_ptr]->size,sizeof(TYPE_INDEX));

			if(key_stack[stack_ptr]<bx.num_of_long_blocks){
				block_length = bx.long_block_length;
				block_type=LONG_WORD;
			}else{
				block_length = bx.short_block_length;
				block_type=SHORT_WORD;
			}
/**
			int offset;//blockの先頭
			if(key_stack[stack_ptr]<bx.num_of_long_blocks){
				offset = key_stack[stack_ptr]*bx.long_block_length;
			}else{
				offset = bx.num_of_long_blocks*bx.long_block_length + (key_stack[stack_ptr] - bx.num_of_long_blocks)*bx.short_block_length;
			}
**/

			int offset = bx.block_offset[ key_stack[stack_ptr] ];
			int col;
			int cur, nxt;
			for(col=0;col<num_of_words[block_type];col++){
				cur=0;nxt=0;
				for(int cnt=1;cnt<=is_region_top[stack_ptr]->size;cnt++){
					if(cnt==is_region_top[stack_ptr]->size || is_region_top[stack_ptr]->val[cnt]==REGION_TOP){
						if(cnt-cur==1){
						}else{
							nxt = cnt;
							memset(bucket[block_type],0x00,sizeof(int)*num_of_bucket[block_type]);
							bid_ptr=0;
							//クラスの頻度の数え上げ
							for(int i=cur;i<nxt;i++){
								tmp_word=0;
								for(int cls_cnt=0;cls_cnt<word_size[block_type];cls_cnt++){
									tmp_word = sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
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

							//新しいクラスの追加 //
							for(int i=0;i<bid_ptr;i++){
								//↓　bucket[i] + cur < is_region_top[stack_ptr]->size?? 確かめ
								if(bucket[block_type][ bucketID[block_type][i] ] + cur < bx.num_of_box){
									is_region_top[stack_ptr]->val[ bucket[block_type][ bucketID[block_type][i] ] + cur ]=REGION_TOP;
								}
							}
							//edit distance用はsortingしない。
							//for(int i=bucket[ct.overlap_character-1];i<bucket[ct.overlap_character];i++){
							//	is_region_top[stack_ptr]->val[ i + cur ]=REGION_TOP;
							//	cerr<<cur<<","<< i + cur<<"\n";
							//							}
							//ソート
							for(int i=nxt-1;i>=cur;i--){
								tmp_word=0;
								for(int cls_cnt=0;cls_cnt<word_size[block_type];cls_cnt++){
									tmp_word = sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
								}
								work_order[ --bucket[block_type][tmp_word] + cur ] = order[stack_ptr]->val[i];	
								work_is_offset_zero[ bucket[block_type][tmp_word] + cur ] = is_offset_zero[stack_ptr]->val[i];
							}
						}
						cur = cnt;
					}
				}
				for(int i=0;i<order[stack_ptr]->size;i++){
					order[stack_ptr]->val[i] = work_order[i];
					is_offset_zero[stack_ptr]->val[i] = work_is_offset_zero[i];
				}
			}

			cur=0;nxt=0;
			for(int cnt=1;cnt<=is_region_top[stack_ptr]->size;cnt++){
				if(cnt==is_region_top[stack_ptr]->size || is_region_top[stack_ptr]->val[cnt]==REGION_TOP){
					if(cnt-cur==1){
					}else{
						nxt = cnt;
						memset(bucket[block_type],0x00,sizeof(int)*num_of_bucket[block_type]);
						bid_ptr=0;
						//クラスの頻度の数え上げ
						for(int i=cur;i<nxt;i++){
							tmp_word=0;
							for(int cls_cnt=0;cls_cnt<l_word_size[block_type];cls_cnt++){
								tmp_word = sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
							}
							if(bucket[block_type][tmp_word]==0){
								bucketID[block_type][bid_ptr]=tmp_word;
								bid_ptr++;
							}
							bucket[block_type][ tmp_word ]++;
						}
						for(int i=1;i<bid_ptr;i++){
							bucket[block_type][ bucketID[block_type][i] ] += bucket[block_type][  bucketID[block_type][i-1] ];
						}
						//新しいクラスの追加 //
						for(int i=0;i<bid_ptr;i++){
							//↓　bucket[i] + cur < is_region_top[stack_ptr]->size?? 確かめ
							if(bucket[block_type][ bucketID[block_type][i] ] + cur < bx.num_of_box){
								is_region_top[stack_ptr]->val[ bucket[block_type][ bucketID[block_type][i] ] + cur ]=REGION_TOP;
							}
						}
						//edit distance用はsortingしない。
						//for(int i=bucket[ct.overlap_character-1];i<bucket[ct.overlap_character];i++){
						//	is_region_top[stack_ptr]->val[ i + cur ]=REGION_TOP;
						//	cerr<<cur<<","<< i + cur<<"\n";
						//							}
						//ソート
						for(int i=nxt-1;i>=cur;i--){
							tmp_word=0;
							for(int cls_cnt=0;cls_cnt<l_word_size[block_type];cls_cnt++){
								tmp_word = sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
							}
							work_order[ --bucket[block_type][tmp_word] + cur ] = order[stack_ptr]->val[i];	
							work_is_offset_zero[ bucket[block_type][tmp_word] + cur ] = is_offset_zero[stack_ptr]->val[i];
						}
					}
					cur = cnt;
				}
			}
			for(int i=0;i<order[stack_ptr]->size;i++){
				order[stack_ptr]->val[i] = work_order[i];
				is_offset_zero[stack_ptr]->val[i] = work_is_offset_zero[i];
			}

			if(work_order) free(work_order);
			if(work_is_offset_zero) free(work_is_offset_zero);

			/**print
			for(int i=0;i<order[stack_ptr]->size;i++){
			cerr<<order[stack_ptr]->val[i]<<", ";
			}
			cerr<<"\n";
			for(int i=0;i<order[stack_ptr]->size;i++){
			cerr<<is_region_top[stack_ptr]->val[i]<<", ";
			}
			cerr<<"\n";
			for(int i=0;i<order[stack_ptr]->size;i++){
			cerr<<is_offset_zero[stack_ptr]->val[i]<<", ";
			}
			cerr<<"\n\n";
			print**/
			//			cerr<<"done\n";
		}

		if(cHead<is_region_top[stack_ptr]->size){
			while(1){
				prev = cHead;
				cnt = cHead+1;
				region_len=0;
				//orderは必ず昇順にソートされていることを利用。
				min_zero_id_flag = false;

				if(order[stack_ptr]->val[cHead]%bx.box_set_size == bx.box_set_center){
					if(is_offset_zero[stack_ptr]->val[cHead] == OFFSET_ZERO){
						is_region_top[stack_ptr]->val[cHead] = OFFSET_ZERO;
						min_zero_id_flag = true;
						region_len++;
					}
				}
				while( cnt < is_region_top[stack_ptr]->size && is_region_top[stack_ptr]->val[cnt]!=REGION_TOP ){
					if( min_zero_id_flag && order[stack_ptr]->val[cnt]/bx.box_set_size != order[stack_ptr]->val[cnt-1]/bx.box_set_size){
						is_region_top[stack_ptr]->val[cnt] = COPY_TO_NEXT_LIST;
						prev = cnt;
						region_len++;
					}
					if(order[stack_ptr]->val[cnt]%bx.box_set_size == bx.box_set_center){
						//コピーするindexにzeroをつける。
						if(is_offset_zero[stack_ptr]->val[cnt] == OFFSET_ZERO){
							if(min_zero_id_flag){
								is_region_top[stack_ptr]->val[prev] = OFFSET_ZERO;								
							}else{
								is_region_top[stack_ptr]->val[cnt] = OFFSET_ZERO;
								region_len++;
								min_zero_id_flag = true;
							}
						}
					}
					cnt++;
				}
//				cerr<<region_len<<"/"<<cnt-cHead<<"\n";
				if(region_len<2){ // offset_zeroが一個もない場合は、これ以上ソートしても意味ない。
					cHead=cnt;
					//				if(num_of_offset_zero==0)
					//					cerr<<".";
					if(cnt >= is_region_top[stack_ptr]->size){
						break;
					}
				}else{
					break;
				}
			}
		}
		pHead = cHead;
/** print
		for(int i=0;i<order[stack_ptr]->size;i++){
			cerr<<order[stack_ptr]->val[i]<<", ";
		}
		cerr<<"\n";
		for(int i=0;i<order[stack_ptr]->size;i++){
			cerr<<is_region_top[stack_ptr]->val[i]<<", ";
		}
		cerr<<"\n";
		for(int i=0;i<order[stack_ptr]->size;i++){
			cerr<<is_offset_zero[stack_ptr]->val[i]<<", ";
		}
		cerr<<"\n\n";
**/
		if(cHead>=is_region_top[stack_ptr]->size){
			next_block = key_stack[stack_ptr] +1;
			if(next_block <= co.distance + stack_ptr){
				key_stack.pop_back();
				key_stack.push_back(next_block);
				block_stack[stack_ptr] =0;
				cHead=0;
				copy_flag = 0;
//				cerr<<"move to next block "<<next_block<<"\n";
			}else{
//				cerr<<"free "<<stack_ptr<<"\n";
				if(order[stack_ptr]) delete(order[stack_ptr]);
//				cerr<<"order ";
				if(is_region_top[stack_ptr]) delete(is_region_top[stack_ptr]);
//				cerr<<" region";
				if(is_offset_zero[stack_ptr]) delete(is_offset_zero[stack_ptr]);
//				cerr<<" zero ";
				if(order_cpy[stack_ptr]) delete(order_cpy[stack_ptr]);
//				cerr<<" order c ";
				if(is_offset_zero_cpy[stack_ptr]) delete(is_offset_zero_cpy[stack_ptr]);
//				cerr<<" zero c\n ";
				key_stack.pop_back();
				block_stack[stack_ptr] = 0;
				stack_ptr--;
				if(stack_ptr>=0){
					cHead = block_stack[stack_ptr];
				}
			}
		}else{		
//			cerr<<"key_stack "<<key_stack[stack_ptr]<<"\n";
			block_stack[stack_ptr] = cnt;
			region_len *= bx.box_set_size;
//			cerr<<"region-len "<<region_len<<"\n";
			next_block = key_stack[stack_ptr] +1;
			key_stack.push_back(next_block);
			stack_ptr++;
			cHead=0;
			copy_flag = 1;
//			cerr<<"goto next\n";
//			cerr<<"next_ptr "<<cnt<<"\n";
//			cerr<<"from "<<pHead<<" to "<<cnt<<"\n";
		}
	}
	cerr<<"end of ED\n";
	// free 対策
	for(int i=0;i<bx.key_size;i++){
		order[i] = NULL;
		is_region_top[i] = NULL;
	}
	return(0);
}
#endif

