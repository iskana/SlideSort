#include "parallelslidesort.h"

/**************************************************
 *
 * pssExecutor class
 *
 * multisort実行クラス.
 *
 *  Created on: 2012/12/05
 *
 *************************************************/


/**************************************************
 *
 * ロック設定関数.
 *
 * ロックを設定する.
 *
 * Argument:
 *   *myLock : I : ロック
 *
 * Return:
 *
 *************************************************/
void pssExecutor::setMyLock(omp_lock_t *myLock) {
    myLock_ptr = myLock;
}

/**************************************************
 *
 * コールバック関数実行関数.
 *
 * コールバック関数を呼び出す.
 * OpenMPの定義を呼び出し、排他制御を行う.
 * 開発時は、multisortのexecCallbackFuncを呼び出す.
 *
 * Argument:
 *   (コールバック関数と同じ引数）
 *
 * Return:
 *   処理に成功した場合:0、失敗した場合:1
 *
 *************************************************/
int pssExecutor::execCallbackFunc(const char *seqid1, const char *seqid2, TYPE_INDEX index_of_seq1, TYPE_INDEX index_of_seq2, char* aln1, char* aln2, double dist, int aln_size)
{
    omp_set_lock(myLock_ptr);
    int ret = grfptr(seqid1, seqid2, index_of_seq1, index_of_seq2, aln1, aln2, dist, aln_size);
    omp_unset_lock(myLock_ptr);
    return ret;
}


/**************************************************
 *
 * ブロック組合せ情報設定関数.
 *
 * ブロック組合せ情報を設定する.
 *
 * Argument:
 *   *rp_bc : I : ブロック組合せ情報
 *
 * Return:
 *
 *************************************************/
void pssExecutor::setBC(BlockCombination *rp_bc)
{
    p_bc = rp_bc;
}
