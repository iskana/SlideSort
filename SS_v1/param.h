/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              param.h
/***********************************************/

#define FORMAT_FASTA 0x00
#define FORMAT_FASTQ 0x01

#define FORMAT_FASTA_TO_FASTQ 0x01
#define FORMAT_FASTQ_TO_FASTA 0x02


#define HAMMING_DISTANCE 0x00
#define EDIT_DISTANCE 0x01

#define DNA 0x00
#define PROTEIN 0x01
#define INPUT_INT 0x02

#define DNA_A 0
#define DNA_T 1
#define DNA_G 2
#define DNA_C 3
#define GAP 100

#define REV_DNA_A 1
#define REV_DNA_T 0
#define REV_DNA_G 3
#define REV_DNA_C 2

#define MAX_FASTA_ID_LENGTH 128

#define MAX_RECOMMENDED_EDIT_DISTANCE 10

#define MIN_CHARACTER_FOR_HPMODE 7

#define REGION_TOP 1
#define REGION_SKIP 2

#define COPY_TO_NEXT_LIST 3
//COPY_TO_NEXT_LISTよりも大きい値に。
#define OFFSET_ZERO 4

#define REDUNDANT 0
#define NON_REDUNDANT 1

#define SHORT_WORD 0
#define LONG_WORD 1

//for SSE
#ifdef VCPP
#define BIT_BLOCK unsigned __int64 
#else
#define BIT_BLOCK uint64_t 
#endif
// (sizeof(unsigned __int64)*8) とすると、マイナスの値が表現できない。
#define SIZE_OF_REGISTER 128
#define BIT_BLOCK_LENGTH 64
#define ALIGNED_BORDER 16
#define BIT_BLOCK_0 0x0000000000000000
#define BIT_BLOCK_1 0xffffffffffffffff



// type define
#define TYPE_CHARACTER char //used for holding character
#define SIZE_OF_TYPE_CHARACTER 64
#define TYPE_CHAR char
#define TYPE_INDEX long long
#define TYPE_LABEL int // for offset
#define TYPE_LONG unsigned long long
/* unsinged int を使ってる時は、flagにマイナスを使わないこと！！*/

#define MIN_SEQ_LENGTH_DEFAULT 1000000000

#define FST_DATASET 0
#define SND_DATASET 1
#define BOTH_DATASET 2

#define SET_REGION_START 10

#define PROTEIN_A 0
#define PROTEIN_C 1
#define PROTEIN_D 2
#define PROTEIN_E 3
#define PROTEIN_F 4
#define PROTEIN_G 5
#define PROTEIN_H 6
#define PROTEIN_I 7
#define PROTEIN_K 8
#define PROTEIN_L 9
#define PROTEIN_M 10
#define PROTEIN_N 11
#define PROTEIN_P 12
#define PROTEIN_Q 13
#define PROTEIN_R 14
#define PROTEIN_S 15
#define PROTEIN_T 16
#define PROTEIN_V 17
#define PROTEIN_W 18
#define PROTEIN_Y 19

#define HAMMING_KEY_DEC_SEQ_NUM 18000000

#define min(x,y) (x<y ? x:y)

#define TMP_FILE_SIZE 100000000
