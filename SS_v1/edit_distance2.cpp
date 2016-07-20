/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              edit_distance2.cpp
/***********************************************/

#include"mscls.h"

#ifdef HPSORT

int multisort::debug(){
	return(0);
}

int debug_flg=0;

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

	TYPE_CHARACTER *out_seq1_rc, *out_seq2_rc;
	if(co.isOutputBothStrand==false){													
		out_seq1_rc = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
		out_seq2_rc = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
		memset(out_seq1_rc,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
		memset(out_seq2_rc,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance));
	}

#ifdef CALL_BACK_FUNC
	TYPE_CHARACTER *rev_seq1, *rev_seq2;
	rev_seq1 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	rev_seq2 = (TYPE_CHARACTER *)malloc(sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq1,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
	memset(rev_seq2,0x00,sizeof(TYPE_CHARACTER)*(sq.max_seq_length+co.distance + 1));
#endif

	// for sort
	TYPE_INDEX *work_order;
	TYPE_CHAR *work_is_offset_zero;
	char_size = ct.num_of_character + 1;
	if(ct.isUseWildCard){
		char_size++;
	}
	word_size[SHORT_WORD] = log((double)sq.num_of_seq)/log((double)char_size)-2;
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

	int bid_ptr;
	int block_length;
	int block_type;
	int tmp_word;

//for small sort
	int sm_bucket_size = char_size;
	int *sm_bucket = (int *)calloc(sizeof(int),sm_bucket_size);
	if(bucket==NULL){
		cerr<<"ERROR: can not allocate memory for sm_bucket\n";
		exit(1);
	}

	//for 
	int *block_stack = (int *)calloc(bx.key_size,sizeof(int));
	int stack_ptr = 0;
	TYPE_INDEX cHead = 0, pHead=0, cnt, prev;
	TYPE_INDEX region_len = bx.num_of_box;
	int next_block, copy_flag;	
	int min_zero_id_flag;

	copy_flag=0;
	if(co.isDevSort){
		key_stack.push_back(co.startBlc);
	}else{
		key_stack.push_back(0);
	}
	order[0] = new varLongInt(bx.num_of_box);
	is_region_top[0] = new varChar(bx.num_of_box);
	is_offset_zero[0] = new varChar(bx.num_of_box);
	resetOrder(0);
	resetIsRegionTop(0);
	for(TYPE_INDEX i=0;i<bx.num_of_box;i++){
		if(i%bx.box_set_size == bx.box_set_center){
			is_offset_zero[0]->val[i]=OFFSET_ZERO;
		}
	}
	//cpy
	if(!co.isDevSort){
		order_cpy[0] = new varLongInt(bx.num_of_box);
		is_offset_zero_cpy[0] = new varChar(bx.num_of_box);
		memcpy(order_cpy[0]->val, order[0]->val, sizeof(TYPE_INDEX)*order[0]->size);
		memcpy(is_offset_zero_cpy[0]->val, is_offset_zero[0]->val, sizeof(TYPE_CHAR)*is_offset_zero[0]->size);
	}

	cerr<<"Starting SlideSort\n";
	sum_of_region=0;
//	cerr<<bx.box_set_center<<"\n";
	while(stack_ptr >= 0){
		if(stack_ptr == bx.key_size){
			stack_ptr--;
			cHead = block_stack[stack_ptr];
			key_stack.pop_back();

			region_len /= bx.box_set_size;
			TYPE_INDEX *pair_list = (TYPE_INDEX *)calloc(sizeof(TYPE_INDEX),region_len);
			TYPE_CHAR *ioz_list = (TYPE_CHAR *)calloc(sizeof(TYPE_CHAR),region_len);
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
			TYPE_INDEX small_id, big_id, org_seq_id1, org_seq_id2;

			for(TYPE_INDEX box1=0;box1<region_len-1;box1++){
				for(TYPE_INDEX box2=box1+1;box2<region_len;box2++){
					org_seq_id1 = (pair_list[box1]/bx.box_set_size);
					org_seq_id2 = (pair_list[box2]/bx.box_set_size);

//					cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[0] ]<<"\t";
//					cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[0] ]<<"\n";

					if(co.isMapMode){
						if(sq.nr_seq_mapID[org_seq_id1]==sq.nr_seq_mapID[org_seq_id2] && sq.nr_seq_mapID[org_seq_id1]!=BOTH_DATASET){
							goto mapDifChk;
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
//							cerr<<sq.revCompIdx[org_seq_id1]<<" "<<sq.revCompIdx[org_seq_id2]<<"\n";
							if(co.isRevComp){
								if( min(sq.revCompIdx[org_seq_id1] , sq.revCompIdx[org_seq_id2])%2 == 1 ){
									goto mapDifChk;
								}
							}

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
								if(key_ptr>=0 && key==key_stack[key_ptr]){
//									cerr<<"key=_stack"<<key_stack[key_ptr]<<"\n";
									key_ptr--;
								}else{
									TYPE_INDEX small_seq_id = (pair_list[small_id]/bx.box_set_size)*bx.box_set_size + bx.box_set_center;
									TYPE_INDEX big_seq_id = (pair_list[big_id]/bx.box_set_size)*bx.box_set_size;
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
								double ed, edr;
								TYPE_INDEX ret_id1, ret_id2;
								if((ed=calcEditDistance(org_seq_id1, org_seq_id2, out_seq1, out_seq2))<=co.distance){
									/** Print IDs**/
									for(int id1=0;id1<sq.toOrgSeqIndex[ org_seq_id1 ]->size;id1++){
										ret_id1 = sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1];
										for(int id2=0;id2<sq.toOrgSeqIndex[ org_seq_id2 ]->size;id2++){
											ret_id2 = sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2];

											if(co.isRevComp){
												if( sq.toOrgSeqIndex[ org_seq_id1 ]->val[0]/2 == sq.toOrgSeqIndex[ org_seq_id2 ]->val[0]/2){
													if(id1>=id2){
														goto revDifChk;													
													}
												}
												if( org_seq_id1 ==  sq.toSSIndex[ (ret_id1/2)*2 + abs(ret_id1%2-1) ] && ret_id1%2==1){
														goto revDifChk;													
												}
												if(co.isOutputBothStrand==false){
													if(fabs((double)sq.revCompSeqDist[org_seq_id2]-ed)<=(double)co.distance){
														TYPE_INDEX rev_id = sq.toSSIndex[ (TYPE_INDEX)(ret_id2/2)*2 + abs(ret_id2%2-1) ];
														edr = calcEditDistance(org_seq_id1, rev_id, out_seq1_rc, out_seq2_rc);
														if(ed==edr){
															if(ret_id2%2==1){
																goto revDifChk;
															}
														}else if(ed>edr){
															cerr<<sq.seqName[  (ret_id2/2)*2 + abs(ret_id2%2-1)  ]<<"  "<<rev_id<<"  "<<ed<<" "<<edr<<" E2\n";
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
														//modify 20101208
														int top_margin=0;
														for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
															*(rev_seq1++) = ct.toChar[ out_seq1[otmp] ];
															*(rev_seq2++) = ct.toChar[ out_seq2[otmp] ];
															if(out_seq1[otmp] == out_seq2[otmp] && ct.toChar[ out_seq1[otmp] ] == '-'){
																top_margin++;
															}
														}
														*(rev_seq1)='\0';
														*(rev_seq2)='\0';
														rev_seq1 -= sq.max_seq_length+co.distance;
														rev_seq2 -= sq.max_seq_length+co.distance;
														CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, rev_seq1+top_margin, rev_seq2+top_margin, ed, sq.max_seq_length+co.distance-top_margin);

/**														for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
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
														**/
										//modify 20101208

													}else{
														CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1,ret_id2, NULL,NULL,ed,NULL);													
													}
#endif 
#ifdef FILE_OUTPUT_PAIRS
													if(co.outputfile){
														fprintf(outfp,"%s,%s,%.2f\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),ed);
														if(co.outputaln){
															int top_pos = sq.max_seq_length+co.distance-1;
															while(out_seq1[top_pos]==out_seq2[top_pos] && ct.toChar[ out_seq1[top_pos] ]=='-'){
																top_pos--;
															}
															for(int otmp=top_pos;otmp>=0;otmp--){
																fprintf(outfp,"%c",ct.toChar[ out_seq1[otmp] ]);
															}
															fprintf(outfp,"\n");
															for(int otmp=top_pos;otmp>=0;otmp--){
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
													fprintf(outfp,"%s,%s,%.2f\n",sq.seqName[ ret_id1 ].c_str(),sq.seqName[ ret_id2 ].c_str(),ed);
													if(co.outputaln){
														int top_pos = sq.max_seq_length+co.distance-1;
														while(out_seq1[top_pos]==out_seq2[top_pos] && ct.toChar[ out_seq1[top_pos] ]=='-'){
															top_pos--;
														}
														for(int otmp=top_pos;otmp>=0;otmp--){
															fprintf(outfp,"%c",ct.toChar[ out_seq1[otmp] ]);
														}
														fprintf(outfp,"\n");
														for(int otmp=top_pos;otmp>=0;otmp--){
															fprintf(outfp,"%c",ct.toChar[ out_seq2[otmp] ]);
														}
														fprintf(outfp,"\n");
													}
												}
#endif
#ifdef CALL_BACK_FUNC
												if(co.outputaln){	
													int top_margin=0;
													for(int otmp=sq.max_seq_length+co.distance-1;otmp>=0;otmp--){
														*(rev_seq1++) = ct.toChar[ out_seq1[otmp] ];
														*(rev_seq2++) = ct.toChar[ out_seq2[otmp] ];
														if(out_seq1[otmp] == out_seq2[otmp] && ct.toChar[ out_seq1[otmp] ] == '-'){
															top_margin++;
														}
													}
													*(rev_seq1)='\0';
													*(rev_seq2)='\0';
													rev_seq1 -= sq.max_seq_length+co.distance;
													rev_seq2 -= sq.max_seq_length+co.distance;
													CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, rev_seq1+top_margin, rev_seq2+top_margin, ed, sq.max_seq_length+co.distance-top_margin);
/**
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
													**/
												}else{
													CALL_BACK_FUN_PTR_OR_WRAP_FUNC(sq.seqName[ ret_id1 ].c_str(), sq.seqName[ ret_id2 ].c_str(),ret_id1, ret_id2, NULL,NULL,ed,NULL);
												}
#endif 
												//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id1 ]->val[id1] ]<<"\t";
												//cout<<sq.seqName[ sq.toOrgSeqIndex[ org_seq_id2 ]->val[id2] ]<<"\t"<<ed<<"\n";
											}

revDifChk:;
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
				order[stack_ptr] = new varLongInt(region_len);
				// 原因不明↓
				is_region_top[stack_ptr] = new varChar(region_len+1);
				is_region_top[stack_ptr]->size = region_len;
				is_offset_zero[stack_ptr] = new varChar(region_len);
				order_cpy[stack_ptr] = new varLongInt(region_len);
				is_offset_zero_cpy[stack_ptr] = new varChar(region_len);

				//copy
				TYPE_INDEX cnt=pHead;
				TYPE_INDEX cpy=0;
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
				memcpy(is_offset_zero_cpy[stack_ptr]->val, is_offset_zero[stack_ptr]->val, sizeof(TYPE_CHAR)*is_offset_zero[stack_ptr]->size);
//				cerr<<"malloc:"<<stack_ptr<<"\n";
			}else{
				// 異なるブロックでソート
				if(stack_ptr==0){
					if(co.isDevSort){
						if(key_stack[stack_ptr]!=co.startBlc)
							goto END_OF_BLOCKSORT_ED;
					}else{
						//orderとiozの修正の必要あり。
						memcpy(order[stack_ptr]->val, order_cpy[stack_ptr]->val, sizeof(TYPE_INDEX)*order_cpy[stack_ptr]->size);
						memcpy(is_offset_zero[stack_ptr]->val, is_offset_zero_cpy[stack_ptr]->val, sizeof(TYPE_CHAR)*is_offset_zero_cpy[stack_ptr]->size);
						//				cerr<<"cpy:"<<stack_ptr<<"\n";
					}
					cerr<<"Block combinations starting from "<<key_stack[stack_ptr]<<".\n";
				}else{
					//orderとiozの修正の必要あり。
					memcpy(order[stack_ptr]->val, order_cpy[stack_ptr]->val, sizeof(TYPE_INDEX)*order_cpy[stack_ptr]->size);
					memcpy(is_offset_zero[stack_ptr]->val, is_offset_zero_cpy[stack_ptr]->val, sizeof(TYPE_CHAR)*is_offset_zero_cpy[stack_ptr]->size);
					//				cerr<<"cpy:"<<stack_ptr<<"\n";
				}
			}

			/****   SORT *****/
//			cerr<<"sort "<<key_stack[stack_ptr]<<"...\n";

			memset(is_region_top[stack_ptr]->val,0x00,sizeof(TYPE_CHAR)*is_region_top[stack_ptr]->size);
			is_region_top[stack_ptr]->val[0] = REGION_TOP;

			work_order = (TYPE_INDEX *)calloc(order[stack_ptr]->size,sizeof(TYPE_INDEX));
			work_is_offset_zero = (TYPE_CHAR *)calloc(is_offset_zero[stack_ptr]->size,sizeof(TYPE_CHAR));

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
			TYPE_INDEX cur, nxt;
//			cerr<<word_size[block_type]<<"\n";
//			cerr<<num_of_words[block_type]<<"\n";
			if(is_region_top[stack_ptr]->size > num_of_bucket[block_type]){
				for(col=0;col<num_of_words[block_type];col++){
					cur=0;nxt=0;
					TYPE_INDEX num_ec=0;
					for(TYPE_INDEX cnt=1;cnt<=is_region_top[stack_ptr]->size;cnt++){
						if(cnt==is_region_top[stack_ptr]->size || is_region_top[stack_ptr]->val[cnt]==REGION_TOP){
							if(cnt-cur==1){
							}else{
								nxt = cnt;
								memset(bucket[block_type],0x00,sizeof(TYPE_INDEX)*num_of_bucket[block_type]);
								bid_ptr=0;
								//クラスの頻度の数え上げ
								for(TYPE_INDEX i=cur;i<nxt;i++){
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
								num_ec++;

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
//								for(TYPE_INDEX i=nxt-1;i>=cur;i--){
								for(TYPE_INDEX i=nxt-1;;i--){
									tmp_word=0;
									for(int cls_cnt=0;cls_cnt<word_size[block_type];cls_cnt++){
										tmp_word = sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col*word_size[block_type] + cls_cnt ] + tmp_word*char_size;
									}
									work_order[ --bucket[block_type][tmp_word] + cur ] = order[stack_ptr]->val[i];	
									work_is_offset_zero[ bucket[block_type][tmp_word] + cur ] = is_offset_zero[stack_ptr]->val[i];
									if(i==cur) break;
								}
							}
							cur = cnt;
						}
					}
					for(TYPE_INDEX i=0;i<order[stack_ptr]->size;i++){
						order[stack_ptr]->val[i] = work_order[i];
						is_offset_zero[stack_ptr]->val[i] = work_is_offset_zero[i];
					}
					/**ようけんとう**/
					if(num_ec > sq.num_of_seq/1000){
						col++;
						break;
					}
				}
				col = col*word_size[block_type];
			}else{
				col=0;
			}
			for(;col<block_length;col++){
				cur=0;nxt=0;
				for(int cnt=1;cnt<=is_region_top[stack_ptr]->size;cnt++){
					if(cnt==is_region_top[stack_ptr]->size || is_region_top[stack_ptr]->val[cnt]==REGION_TOP){
						if(cnt-cur==1){
						}else{
							nxt = cnt;
							memset(sm_bucket,0x00,sizeof(int)*sm_bucket_size);
							//クラスの頻度の数え上げ
							for(int i=cur;i<nxt;i++){
								sm_bucket[ sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col ] ]++;
							}
							for(int i=1;i<sm_bucket_size;i++){
								sm_bucket[i] += sm_bucket[i-1];
							}
							//新しいクラスの追加 //
							for(int i=0;i<char_size-1;i++){
								//↓　sm_bucket[i] + cur < is_region_top[stack_ptr]->size?? 確かめ
								if(sm_bucket[i] + cur < bx.num_of_box){
									is_region_top[stack_ptr]->val[ sm_bucket[i] + cur ]=REGION_TOP;
								}
							}
							//edit distance用はsortingしない。
							//for(int i=sm_bucket[ct.overlap_character-1];i<sm_bucket[ct.overlap_character];i++){
							//	is_region_top[stack_ptr]->val[ i + cur ]=REGION_TOP;
							//	cerr<<cur<<","<< i + cur<<"\n";
							//							}
							//ソート
//							for(TYPE_INDEX i=nxt-1;i>=cur;i--){
							for(TYPE_INDEX i=nxt-1;;i--){
								work_order[ --sm_bucket[ sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col ]  ] + cur ] = order[stack_ptr]->val[i];	
								work_is_offset_zero[ sm_bucket[ sq.nr_seq[ bx.head[ order[stack_ptr]->val[i] ] + offset + col ]  ] + cur ] = is_offset_zero[stack_ptr]->val[i];
								if(i==cur) break;
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

		//コピー元の指定
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
//			sum_of_region += region_len;
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

END_OF_BLOCKSORT_ED:

	cerr<<"Finish\n";
//	cerr<<"sum_of_region:"<<sum_of_region<<"\n";
	free(sm_bucket);
	free(out_seq1);
	free(out_seq2);
	if(co.isOutputBothStrand==false){
		free(out_seq1_rc);
		free(out_seq2_rc);
	}

	// free 対策
	for(int i=0;i<bx.key_size;i++){
		is_region_top[i] = NULL;
		order[i] = NULL;
		order_cpy[i] = NULL;
		is_offset_zero_cpy[i] = NULL;
		is_offset_zero[i] = NULL;
	}
	return(0);
}
#endif

