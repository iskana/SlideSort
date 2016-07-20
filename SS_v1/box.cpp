/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              box.cpp
/***********************************************/

#include"mscls.h"


void box::freeBox(){
	if(head){
		free(head);
		head = NULL;
	}
	if(block_offset){
		free(block_offset);
		block_offset=NULL;
	}
}

void box::outputBoxInfo()
{
	
	cerr<<"--------BOX INFO--------\n";
	cerr<<"num_of_box: "<<num_of_box<<"\n";
	cerr<<"box_length: "<<box_length<<"\n";

	cerr<<"num_of_blocks: "<<num_of_blocks<<"\n";
	cerr<<"num_of_long_blocks: "<<num_of_long_blocks<<"\n";
	cerr<<"short_block_length: "<<short_block_length<<"\n";
	cerr<<"long_block_length: "<<long_block_length<<"\n";
	cerr<<"key_size: "<<key_size<<"\n";

}