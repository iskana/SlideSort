/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//             commoh.h
/***********************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<math.h>
#include<memory.h>

using namespace std;

#define VCPP
#define FILE_OUTPUT_PAIRS
#define HPSORT
#define NEWSSE
//#define SSE4
//#define CALL_BACK_FUNC
//#define INTERNAL_GAP_OPEN
#define CALL_BACK_FUNC_BY_VIRTUAL_FUNC

/******
#define INDEXING_INPUT_ORDER 
******/

#ifdef CALL_BACK_FUNC_BY_VIRTUAL_FUNC
#define CALL_BACK_FUN_PTR_OR_WRAP_FUNC execCallbackFunc
#else
#define CALL_BACK_FUN_PTR_OR_WRAP_FUNC grfptr
#endif

#ifdef VCPP
#define pow(x,y) (TYPE_LABEL)pow((double)x,(double)y)
#endif

#define SINGLE_CORE