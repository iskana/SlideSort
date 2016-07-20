/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              chartable.cpp
/***********************************************/

#include"mscls.h"

void charTable::freeTables(){
	if(toInt){
		free(toInt);
		toInt=NULL;
	}
	if(toChar){
		free(toChar);
		toChar=NULL;
	}
}

void charTable::setCharTable(cmlOptions co)
{
	switch(co.charType){
		case DNA:
			num_of_character=4;
			unknown_character=4;
			if(!co.exclude_unknown_character){
				overlap_character = 5;
				lim_wild_card = 6;
				num_of_character++;
			}else{
				overlap_character = 4;
				lim_wild_card = 5;
			}
			toInt=(int*)malloc(sizeof(int)*128);
			for(int i=0;i<128;i++) toInt[i] = unknown_character;

			//FOR DNA SETTING
			toInt['a']=DNA_A;
			toInt['t']=DNA_T;
			toInt['g']=DNA_G;
			toInt['c']=DNA_C;
			toInt['A']=DNA_A;
			toInt['T']=DNA_T;
			toInt['G']=DNA_G;
			toInt['C']=DNA_C;

//			toChar = (char*)malloc(sizeof(char)*(num_of_character+1));
//20110112
			toChar = (char*)malloc(sizeof(char)*(num_of_character+3));
			toChar[DNA_A] = 'A';
			toChar[DNA_T] = 'T';
			toChar[DNA_G] = 'G';
			toChar[DNA_C] = 'C';
			toChar[lim_wild_card] = 'w';
			toChar[unknown_character] = 'n';
			toChar[overlap_character] = 'o';
			toChar[lim_wild_card+1] = '-';
			break;

		case PROTEIN:
			num_of_character=20;
			unknown_character=20;
			if(!co.exclude_unknown_character){
				overlap_character = 21;
				lim_wild_card = 22;
				num_of_character++;
			}else{
				overlap_character = 20;
				lim_wild_card = 21;
			}
			toInt=(int*)malloc(sizeof(int)*128);
			for(int i=0;i<128;i++) toInt[i] = unknown_character;

			toInt['A']=PROTEIN_A;
			toInt['C']=PROTEIN_C;
			toInt['D']=PROTEIN_D;
			toInt['E']=PROTEIN_E;
			toInt['F']=PROTEIN_F;
			toInt['G']=PROTEIN_G;
			toInt['H']=PROTEIN_H;
			toInt['I']=PROTEIN_I;
			toInt['K']=PROTEIN_K;
			toInt['L']=PROTEIN_L;
			toInt['M']=PROTEIN_M;
			toInt['N']=PROTEIN_N;
			toInt['P']=PROTEIN_P;
			toInt['Q']=PROTEIN_Q;
			toInt['R']=PROTEIN_R;
			toInt['S']=PROTEIN_S;
			toInt['T']=PROTEIN_T;
			toInt['V']=PROTEIN_V;
			toInt['W']=PROTEIN_W;
			toInt['Y']=PROTEIN_Y;

			toInt['a']=PROTEIN_A;
			toInt['c']=PROTEIN_C;
			toInt['d']=PROTEIN_D;
			toInt['e']=PROTEIN_E;
			toInt['f']=PROTEIN_F;
			toInt['g']=PROTEIN_G;
			toInt['h']=PROTEIN_H;
			toInt['i']=PROTEIN_I;
			toInt['k']=PROTEIN_K;
			toInt['l']=PROTEIN_L;
			toInt['m']=PROTEIN_M;
			toInt['n']=PROTEIN_N;
			toInt['p']=PROTEIN_P;
			toInt['q']=PROTEIN_Q;
			toInt['r']=PROTEIN_R;
			toInt['s']=PROTEIN_S;
			toInt['t']=PROTEIN_T;
			toInt['v']=PROTEIN_V;
			toInt['w']=PROTEIN_W;
			toInt['y']=PROTEIN_Y;

			toChar = (char*)malloc(sizeof(char)*(num_of_character+3));
			toChar[PROTEIN_A] = 'A';
			toChar[PROTEIN_C] = 'C';
			toChar[PROTEIN_D] = 'D';
			toChar[PROTEIN_E] = 'E';
			toChar[PROTEIN_F] = 'F';
			toChar[PROTEIN_G] = 'G';
			toChar[PROTEIN_H] = 'H';
			toChar[PROTEIN_I] = 'I';
			toChar[PROTEIN_K] = 'K';
			toChar[PROTEIN_L] = 'L';
			toChar[PROTEIN_M] = 'M';
			toChar[PROTEIN_N] = 'N';
			toChar[PROTEIN_P] = 'P';
			toChar[PROTEIN_Q] = 'Q';
			toChar[PROTEIN_R] = 'R';
			toChar[PROTEIN_S] = 'S';
			toChar[PROTEIN_T] = 'T';
			toChar[PROTEIN_V] = 'V';
			toChar[PROTEIN_W] = 'W';
			toChar[PROTEIN_Y] = 'Y';
			toChar[unknown_character] = 'z';

			toChar[lim_wild_card] = 'w';
			toChar[unknown_character] = 'n';
			toChar[overlap_character] = 'o';
			toChar[lim_wild_card+1] = '-';

			break;

		case INPUT_INT:
			break;

		default:
			break;
	}

}

