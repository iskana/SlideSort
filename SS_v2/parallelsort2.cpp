#include <stdio.h>
#include <iostream>
#include <time.h>
#include "mscls.h"

int parallelslidesort::exec2()
{
    time_t start, end;
    start = time(NULL);

    char    msg[512];

    memset(msg, 0, sizeof(msg));

    // -mt または -mp オプションが設定されている場合は、parallelslidesort を実行する。
    if (isMTOptSet||isMPOptSet) {

		initMapRanges();
		setMapRanges(0,sq.num_of_seq,SET_REGION_START);

        omp_set_num_threads(mt); // マルチスレッド数（オプションで指定）を設定する
        omp_init_lock(&myLock); // ロックの生成

        #pragma omp parallel for
		for (int i = 0; i < co.distance+1; ++i) // Vector bc の要素数分ループ
        {
//            printf("Thread ID:%d, index:%d\n", omp_get_thread_num(), i);
//           fprintf(stderr, "parallelsort::exec ; Thread ID:%d, index:%d\n", omp_get_thread_num(), i);
            pssExecutor psse;

            psse.setResGetFuncPtr(grfptr);
            psse.getCmlOptions(co);

            psse.attachSeq(sq);
            psse.attachBox(bx);
            psse.attachCharTable(ct);
            psse.attachBitstream(bs, block_mask, lower_block_mask);
            psse.free_vals_automatic=false;

            // -o オプションを無効にする
//            psse.co.outputfile    = false;
            psse.co.isMapMode = false;

			psse.setSortRange(0,sq.num_of_seq);
            psse.setMyLock(&myLock);
            psse.exec();

            psse.detachBitstream();
            psse.detachCharTable();
            psse.detachBox();
            psse.detachSeq();
            // スレッド毎の類似配列ペア数を類似配列ペアの合計値に加算
            num_of_similar_pairs += psse.num_of_similar_pairs;
        }
    }
    // -mt オプションが設定されていない場合はシングルプロセスで slidesort を実行する。
    else {
        // slidesort のメイン処理を実行する関数を呼び出す。
        ss_main();
    }

    end = time(NULL);
    cerr<<"All the procedures are done in "<<difftime(end, start)<<" sec.\n";

    return(0);
}