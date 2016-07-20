/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              seq.cpp
/***********************************************/

#include"mscls.h"

void seq::freeSeq(){
	if(mem_aloc_flag){
		if(nr_seq){
			free(nr_seq);
			nr_seq = NULL;
		}
		if(head){
			free(head);
			head = NULL;
		}
		if(toOrgSeqIndex)
			freetoOrgSeqIndex();
		if(seqName){
			delete [] (seqName);
			seqName=NULL;
		}
		if(revCompIdx){
			free(revCompIdx);
			revCompIdx=NULL;
		}
		if(toSSIndex){
			free(toSSIndex);
			toSSIndex=NULL;
		}
		if(num_of_unknown_char_matchs){
			free(num_of_unknown_char_matchs);
			num_of_unknown_char_matchs=NULL;
		}
		if(nr_seq_mapID){
			free(nr_seq_mapID);
			nr_seq_mapID=NULL;
		}
		if(org_seq_mapID){
			free(org_seq_mapID);
			org_seq_mapID=NULL;
		}
		if(revCompSeqDist){
			free(revCompSeqDist);
			revCompSeqDist=NULL;
		}
	}
	mem_aloc_flag=false;
	num_of_seq_of_all_input=0;
	seq_length=0;
	num_of_seq=0;
	num_of_valid_seq=0;
	max_seq_length=0;
	min_seq_length=0;
}

void seq::freetoOrgSeqIndex()
{
	for(TYPE_INDEX i=0;i<num_of_seq;i++){
		if(toOrgSeqIndex[i])
			delete(toOrgSeqIndex[i]);
	}
	if(toOrgSeqIndex)
		free(toOrgSeqIndex);
}

int seq::makeMapOfOrgIdxtoSSIdx()
{
	toSSIndex = (TYPE_INDEX*)malloc(sizeof(TYPE_INDEX)*num_of_valid_seq);
	memset(toSSIndex,0x00,sizeof(TYPE_INDEX)*num_of_valid_seq);
	for(TYPE_INDEX i=0;i<num_of_seq;i++){
		for(TYPE_INDEX j=0;j<toOrgSeqIndex[i]->size;j++){
			toSSIndex[ toOrgSeqIndex[i]->val[j] ] = i;
		}
	}
	return(0);
}

void seq::outputSeqInfo()
{
	cerr<<"--------SEQ INFO--------\n";
	cerr<<"num_of_original/valid_and_non_redundant seq (with revComp if -v): "<<num_of_seq_of_all_input<<"/"<<num_of_seq<<"\n";
	cerr<<"non_redundant_seq_length: "<<seq_length<<"\n";
	cerr<<"max_seq_length: "<<max_seq_length<<"\n";
	cerr<<"min_seq_length: "<<min_seq_length<<"\n";

	/** for debug
	for(int i=0;i<seq_length;i++)
		cerr<<nr_seq[i];
	cerr<<"\n";
	**/
}

void seq::makeRevCompIdx(){
	/**
	TYPE_INDEX min;
	for(int i=0;i<num_of_seq;i++){
		min = toOrgSeqIndex[i]->val[0];
		for(int j=1; j<toOrgSeqIndex[i]->size;j++){
			if(min>toOrgSeqIndex[i]->val[j]){
				min = toOrgSeqIndex[i]->val[j];
			}
		}
		revCompIdx.push_back(min);
	}

	for(int i=0;i<num_of_seq;i++){
		cerr<<revCompIdx[i]<<"\n";
	}
	**/
	revCompIdx = (TYPE_INDEX*)malloc(sizeof(TYPE_INDEX)*num_of_seq);
	for(TYPE_INDEX i=0;i<num_of_seq;i++){
		revCompIdx[i] = toOrgSeqIndex[i]->val[0];
	}
}

void seq::makeRevCompSeqDist(){
	int dna[4]={0,0,0,0};
	int match=0;
	TYPE_CHARACTER *tmp_seq = nr_seq;
	revCompSeqDist = (int*)malloc(sizeof(int)*num_of_seq);
	for(int i=0;i<num_of_seq;i++){
		for(int cnt=0;cnt<head[i+1]-head[i];cnt++){
			if(*tmp_seq<4)
				dna[ *tmp_seq ]++;
			tmp_seq++;
		}
		match = min(dna[DNA_A],dna[REV_DNA_A]) +  min(dna[DNA_T],dna[REV_DNA_T]) +  min(dna[DNA_G],dna[REV_DNA_G]) + min(dna[DNA_C],dna[REV_DNA_C]);
		revCompSeqDist[i]=dna[DNA_A]+dna[DNA_T]+dna[DNA_G]+dna[DNA_C]-match;
		dna[DNA_A] = 0;
		dna[DNA_T] = 0;
		dna[DNA_G] = 0;
		dna[DNA_C] = 0;
	}
}

