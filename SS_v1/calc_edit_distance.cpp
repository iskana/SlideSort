/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              sse_edit_distance.cpp
/***********************************************/

#include"mscls.h"
#include <emmintrin.h>

#define MATCH 0
#define MISMATCH 1
#define UP 2
#define LEFT 3
#define TERM 4


//#define PRINT_ALN
//#define DEBUG


int multisort::calcHDistanceFromED(TYPE_INDEX seq1, TYPE_INDEX seq2, TYPE_CHARACTER *out_seq1, TYPE_CHARACTER *out_seq2){
	int hd=0;
	TYPE_INDEX head1 = sq.head[seq1] + bx.box_set_center; 
	TYPE_INDEX head2 = sq.head[seq2] + bx.box_set_center;

	for(int i=0;i<sq.max_seq_length;i++){
		if(sq.nr_seq[head1+i]!=sq.nr_seq[head2+i]) hd++;
		if(hd>co.distance){
			break;
		}
	}
	return(hd);
}


double multisort::calcEditDistance(TYPE_INDEX seq1, TYPE_INDEX seq2, TYPE_CHARACTER *out_seq1, TYPE_CHARACTER *out_seq2){
	double ed=0;
#ifndef DEBUG
	int len  = sq.max_seq_length+1;
	int len2 = co.allowedgap*2+1;
	TYPE_INDEX head1 = sq.head[seq1] + bx.box_set_center; 
	TYPE_INDEX head2 = sq.head[seq2] + bx.box_set_center;
	int top_margin=0;
	int skip_margin=0;
	int last=0,wflag=0;
	int lg1=0,lg2=0;

	while(sq.nr_seq[head1+sq.max_seq_length-1-lg1]==ct.lim_wild_card){
		lg1++;
	}
	while(sq.nr_seq[head2+sq.max_seq_length-1-lg2]==ct.lim_wild_card){
		lg2++;
	}
	if(abs(lg1-lg2)>co.gap_limit){
		return(co.distance+1);
	}

	while(sq.nr_seq[head1]==ct.lim_wild_card && sq.nr_seq[head2]==ct.lim_wild_card){
		head1++;
		len--;
		head2++;
	}
	while(sq.nr_seq[head1+skip_margin]==ct.lim_wild_card){
		skip_margin++;
	}
	while(sq.nr_seq[head2+top_margin]==ct.lim_wild_card){
		top_margin++;
	}

	if(abs(lg1-lg2)+skip_margin+top_margin>co.gap_limit){
		return(co.distance+1);
	}

	double *tmp = (double *)calloc(sizeof(double),(len)*(len2));
	int *ptr = (int *)calloc(sizeof(int),(len)*(len2));

	int i=0,j=0,cur;
	double min;
	double val;

	int chk=0;
#ifndef INTERNAL_GAP_OPEN
	int is_head=1;
#endif
	for(i=0;i<len;i++){
		if(i<skip_margin+1){
#ifdef PRINT_ALN
			cerr<<"-  ";
#endif
			for(j=0;j<co.allowedgap+ top_margin - skip_margin;j++){
				tmp[(i)*len2+j] = bx.box_length*co.gap_cost;
				ptr[(i)*len2+j] = TERM; 
#ifdef PRINT_ALN
				cerr<<tmp[(i)*len2+j]<<"\t";
#endif
			}
			for(;j<len2;j++){
				tmp[(i)*len2+j] = (j-co.allowedgap - top_margin + skip_margin)*co.gap_cost;
#ifndef INTERNAL_GAP_OPEN
				if( (j-co.allowedgap - top_margin + skip_margin)>0 )
					tmp[(i)*len2+j] += co.gap_open;	
#endif
				ptr[(i)*len2+j] = LEFT; 
#ifdef PRINT_ALN
				cerr<<tmp[(i)*len2+j]<<"\t";
#endif
			}
			ptr[(i)*len2] = TERM;
#ifdef PRINT_ALN
			cerr<<"\n";
#endif
		}else{
			if( (sq.nr_seq[ head1 + i -1] == ct.overlap_character) || ( sq.nr_seq[ head1 + i -1] == ct.lim_wild_card )){
//				cerr<<"break by w\n";
				break;
			}
			if(chk==len2){
//				cerr<<"break by chk\n";
				free(tmp);
				free(ptr);
				return(co.distance+1);
			}
			chk=0;

#ifdef PRINT_ALN
//		cerr<<i<<"-";
		cerr<<ct.toChar[sq.nr_seq[ head1 + i -1 ]]<<"  ";
#endif
			for(cur=0,j=i-co.allowedgap;j<i+co.allowedgap+1;j++,cur++){
				wflag=0;
				if(i==skip_margin+1){
					wflag=1;
				}
				if( (sq.nr_seq[ head2 + j] == ct.lim_wild_card || sq.nr_seq[ head2 + j] == ct.overlap_character) && sq.nr_seq[ head2 + j - 1] !=  ct.lim_wild_card && sq.nr_seq[ head2 + j - 1] != ct.overlap_character ){
				//wの一つ手前。
					last = cur;
					wflag=1;
				}
				if(j==0){
					tmp[ (i)*len2 + cur ] = tmp[ (i-1)*len2 + cur + 1 ]+co.gap_cost;

#ifndef INTERNAL_GAP_OPEN
					if(is_head){
						tmp[ (i)*len2 + cur ] +=  co.gap_open;
						is_head=0;
					}
#endif

#ifdef PRINT_ALN
					cerr<<ct.toChar[sq.nr_seq[ head2 + j -1 ]];
					cerr<<tmp[(i)*len2 + cur]<<"\t";
#endif
				}else{
#ifdef PRINT_ALN
					cerr<<ct.toChar[sq.nr_seq[ head2 + j -1 ]];
#endif
					if(sq.nr_seq[ head1 + i -1] == sq.nr_seq[ head2 + j -1] && sq.nr_seq[ head1 + i -1]< ct.unknown_character){
						min = tmp[(i-1)*len2 + cur ];
						ptr[(i)*len2 + cur] = MATCH;
					}else{
						min = tmp[(i-1)*len2 + cur ] +1;
						ptr[(i)*len2 + cur] = MISMATCH;
					}
					if(cur<(co.allowedgap)*2){
						val = tmp[(i-1)*len2 + cur +1 ] + co.gap_cost;
						if(ptr[(i-1)*len2 + cur +1]!=UP){
							val += co.gap_open;
#ifdef INTERNAL_GAP_OPEN
							if(wflag==1){
								val -= co.gap_open;
							}
#endif
						}
						if( val<= min){
							min = val;
							ptr[(i)*len2 + cur] = UP;
						}
					}
					if(cur>0){					
						val = tmp[(i)*len2 + cur-1]+co.gap_cost;
						if(ptr[(i)*len2 + cur-1] != LEFT){
							val += co.gap_open;
#ifdef INTERNAL_GAP_OPEN
							if((sq.nr_seq[ head1 + i] == ct.overlap_character) || ( sq.nr_seq[ head1 + i] == ct.lim_wild_card )){
								val-=co.gap_open;
							}
#endif
						}
						if(val <= min){
							min = val;
							ptr[(i)*len2 + cur] = LEFT;
						}
					}
					tmp[(i)*len2 + cur] = min;

					if(min>co.distance){
						chk++;
					}
#ifdef PRINT_ALN
//					cerr<<ptr[(i)*len2 + cur]<<"\t";
					cerr<<tmp[(i)*len2 + cur]<<"\t";
#endif
				}
			}
#ifdef PRINT_ALN
			cerr<<"\n";
#endif
		}
	}

	ed=tmp[(i-1)*len2 + last];

#ifdef PRINT_ALN
	cerr<<"\ned="<<ed<<"\n";
	cerr<<last<<"\n";
#endif

#endif


	//trace back
	if(co.outputaln && ed<=co.distance){
		i--;
		int c = ptr[i*len2 + last];
		int cnt;
//		cerr<<i<<"   "<<last<<"\n";
		for(cnt=0;cnt<sq.max_seq_length+co.distance;cnt++){
			switch(c){
				case MATCH:
					out_seq1[cnt] = sq.nr_seq[ head1 + i -1 ];
					out_seq2[cnt] = sq.nr_seq[ head2 + i -1 + last - co.allowedgap ];
//					cerr<<ct.toChar[ out_seq1[cnt] ]<<"-"<<ct.toChar[ out_seq2[cnt] ]<<"\n";
					i--;
					break;
				case MISMATCH:
					out_seq1[cnt] = sq.nr_seq[ head1 + i -1 ];
					out_seq2[cnt] = sq.nr_seq[ head2 + i -1 + last - co.allowedgap ];
//					cerr<<ct.toChar[ out_seq1[cnt] ]<<"/"<<ct.toChar[ out_seq2[cnt] ]<<"\n";
					i--;
					break;
				case UP:
					out_seq1[cnt] = sq.nr_seq[ head1 + i -1 ];
					out_seq2[cnt] = ct.lim_wild_card+1;
//					cerr<<ct.toChar[ out_seq1[cnt] ]<<"|"<<ct.toChar[ out_seq2[cnt] ]<<"\n";
					i--; last++;
					break;
				case LEFT:
					out_seq1[cnt] = ct.lim_wild_card+1;
					out_seq2[cnt] = sq.nr_seq[ head2 + i - 1 + last - co.allowedgap ];
//					cerr<<ct.toChar[ out_seq1[cnt] ]<<"<"<<ct.toChar[ out_seq2[cnt] ]<<"\n";
					last--;
					break;
				case TERM:
					out_seq1[cnt] = ct.lim_wild_card+1;
					out_seq2[cnt] = ct.lim_wild_card+1;
//					cerr<<ct.toChar[ out_seq1[cnt] ]<<"@"<<ct.toChar[ out_seq2[cnt] ]<<"\n";
					i--;
					goto end_trace;
					break;
			}
			c = ptr[i*len2 + last];
		}
end_trace:;
		for(;cnt<sq.max_seq_length+co.distance;cnt++){
			out_seq1[cnt] = ct.lim_wild_card+1;
			out_seq2[cnt] = ct.lim_wild_card+1;
		}
		for(cnt=sq.max_seq_length+co.distance-1;cnt>=0;cnt--){
			if(out_seq1[cnt] == ct.lim_wild_card || out_seq1[cnt] == ct.overlap_character){
				out_seq1[cnt] = ct.lim_wild_card+1;
			}
#ifdef PRINT_ALN
			cerr<<ct.toChar[ out_seq1[cnt] ];
#endif			
		}
#ifdef PRINT_ALN
		cerr<<"\n";
#endif
		for(cnt=sq.max_seq_length+co.distance-1;cnt>=0;cnt--){
			if(out_seq2[cnt] == ct.lim_wild_card || out_seq2[cnt] == ct.overlap_character){
				out_seq2[cnt] = ct.lim_wild_card+1;
			}
#ifdef PRINT_ALN
			cerr<<ct.toChar[ out_seq2[cnt] ];
#endif
		}
#ifdef PRINT_ALN
		cerr<<"\n";
#endif
	}

	free(tmp);
	free(ptr);
	return(ed);
}


#if 0
int multisort::calcEditDistance(int seq1, int seq2){
	int ed=0;
#ifndef DEBUG
	int dist=0;
	int len=sq.max_seq_length+1;
	int len2 = len;
	int *tmp = (int *)calloc(sizeof(int),(len)*(len2));
	int *ptr = (int *)calloc(sizeof(int),(len)*(len2));
	int i=0,j=0,cur;
	int max;
	int head1 = sq.head[seq1] + bx.box_set_center; 
	int head2 = sq.head[seq2] + bx.box_set_center;

	for(i=1;i<len;i++){
		if(sq.nr_seq[ head1 + i -1] == ct.overlap_character){
			break;
		}
//			cout<<ct.toChar[sq.nr_seq[ sq.head[seq1] + bx.box_set_center + i ]];

		for(j=1,cur=1;j<i+co.distance;j++,cur++){
			if(sq.nr_seq[ head2 + j -1] == ct.overlap_character){
				break;
			}
			max=0;
			if(sq.nr_seq[ head1 + i -1] == sq.nr_seq[ head2 + j -1] ){
				max = tmp[(i-1)*len2 + cur-1 ] + 1;
				ptr[(i)*len2 + cur] = MATCH;
			}else{
				max = tmp[(i-1)*len2 + cur-1 ];
				ptr[(i)*len2 + cur] = MISMATCH;
			}
			if(tmp[(i-1)*len2 + cur] > max ){
				max = tmp[(i-1)*len2 + cur];
				ptr[(i)*len2 + cur] = UP;
			}
			if(tmp[(i)*len2 + cur-1] > max){
				max = tmp[(i)*len2 + cur-1];
				ptr[(i)*len2 + cur] = LEFT;
			}
			tmp[(i)*len2 + cur] = max;
//			cerr<<tmp[(i)*len + j]<<",";
		}
//		cerr<<"\n";
	}
	ed=0;
	dist=0;
	int p=(i-1)*len2 + cur-1;
	while(p>0){
		switch(ptr[p]){
			case MATCH:
				p -= len2;
				p -=1;
				break;
			case MISMATCH:
				p -= len2;
				p -=1;
				ed++;
				break;
			case UP:
				p -=len2;
				ed++;
				break;
			case LEFT:
				ed++;
				p -=1;
				break;
		}
	}
//	cerr<<"ed="<<ed<<"\n";

	free(tmp);
	free(ptr);
#endif
	return(ed);
}
#endif

