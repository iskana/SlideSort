#include<limits.h>
#include "parallelslidesort.h"

/**************************************************
 *
 * ブロック情報計算・取得関数.
 *
 * ブロック情報(ブロックあたりの配列数,データセット内配列開始位置)の設定
 * に必要なブロック毎の配列数を計算しブロック情報をVectorに格納する.
 *
 * Argument:
 *   proc_num       : I : マルチプロセス数
 *   seq_head       : I : データセットごとの開始配列位置
 *                        [0]:1st DataSet,[1]:2nd DataSet
 *   seq_num        : I : データセットごとの配列数
 *                        [0]:1st DataSet,[1]:2nd DataSet
 *   fbi            : O : 1st DataSetのブロック情報格納用Vector
 *   sbi            : O : 2nd DataSetのブロック情報格納用Vector
 *   msg            : O : エラーメッセージ(エラー時のみ設定)
 *
 * Return:
 *   エラーが無い場合:true, エラーがある場合:false
 *
 *************************************************/
int blockutil::getBlockInfo
(
    unsigned int proc_num,
    TYPE_INDEX *seq_head,
    TYPE_INDEX *seq_num,
    vector<BlockInformation> &fbi,
    vector<BlockInformation> &sbi,
    char *msg
)
{

    TYPE_INDEX  block_num[2];

    if (proc_num > LIMIT_PROCS) {
        sprintf(msg, "ERROR: limit of proccess exceeded:%d > limit %d", (int)proc_num, LIMIT_PROCS);
        return false;
    }

    // 各データセットのブロック数を計算する。
	TYPE_INDEX real_proc_num = calcBlockNum(proc_num, seq_num, block_num);

    TYPE_INDEX max_proc_num;
    if (isMapMode) {
        max_proc_num = seq_num[0] * seq_num[1];
    } else {
        max_proc_num = seq_num[0];
    } 

    if (real_proc_num > max_proc_num) {
        strcpy(msg, "ERROR: Too many proccess.");
        return false;
    }

    // Comparing two datasets の場合
    if (isMapMode) {
        // 1st Dataset のブロック情報を Vector に格納する。
        if (!pushBlockInfo(fbi, seq_head[0], seq_num[0], block_num[0])) {
            // ブロック分割数が配列数よりも多い場合はエラーメッセージを表示して false を返す。
            strcpy(msg, "ERROR: Too many proccess.");
			return false;
        }

        // 2nd Dataset のブロック情報を Vector に格納する。
        if (!pushBlockInfo(sbi, seq_head[1], seq_num[1], block_num[1])) {
            // ブロック分割数が配列数よりも多い場合はエラーメッセージを表示して false を返す。
            strcpy(msg, "ERROR: Too many proccess.");
            return false;
        }
    }
    // Finding pairs from some part of an input sequences の場合
    else {
        // 1st Dataset のブロック情報を Vector に格納する。
        if (!pushBlockInfo(fbi, seq_head[0], seq_num[0], block_num[0])) {
            // ブロック分割数が配列数よりも多い場合はエラーメッセージを表示して false を返す。
            strcpy(msg, "ERROR: Too many proccess.");
            return false;
        }
    }

    return true;
}

/**************************************************
 *
 * ブロック情報設定関数.
 *
 * ブロック情報(ブロックあたりの配列数,データセット内配列開始位置)を設定する.
 *
 * Argument:
 *   BlockInfoVector: O : ブロック情報格納用Vector
 *   seq_head       : I : データセット内ブロック開始位置
 *   seq_num        : I : 配列数
 *   block_num      : I : ブロック数
 *
 * Return:
 *   正常:1
 *   ブロックあたりの概算配列数が確保できない(=0)場合:0
 *
 *************************************************/
int blockutil::pushBlockInfo
(
    vector<BlockInformation> &BlockInfoVector,
    TYPE_INDEX seq_head,
    TYPE_INDEX seq_num,
    TYPE_INDEX block_num
)
{
    TYPE_INDEX  seqnum_quotient;
    TYPE_INDEX  seqnum_remainder;
    TYPE_INDEX  offset_of_seq = 0;

    // ブロックあたりの概算配列数を求める。
    seqnum_quotient  = int (seq_num / block_num);

    // ブロックあたりの概算配列数が確保できない場合は 0 を返す。
    if (seqnum_quotient < 1) return 0;

    // ブロックあたりの概算配列数を求めた余りの配列数を求める。
    seqnum_remainder = seq_num % block_num;

    // 各ブロックに、ブロックあたりの配列数を格納する。
    for(unsigned int i=0;i<block_num;i++) {
        BlockInformation bi;

        bi.size  = seqnum_quotient;
        bi.head  = seq_head;

        BlockInfoVector.push_back(bi);
    }

    vector<BlockInformation>::iterator it_bi = BlockInfoVector.begin();

    // ブロック番号の若いブロックから、余りの配列分を均等に割り当てる。
    for(unsigned int i=0;i<seqnum_remainder;i++) {
        (*it_bi).size ++;
        (*it_bi).head += seqnum_quotient * i + offset_of_seq;
        offset_of_seq ++;
        it_bi ++;
    }

    // 余りの配列が割り当てられなかったブロックの配列開始位置を格納する。
    for(TYPE_INDEX i=seqnum_remainder;i<block_num;i++) {
        (*it_bi).head += seqnum_quotient * i + offset_of_seq;
        it_bi ++;
    }

    return 1;
}

/**************************************************
 *
 * ブロック数計算関数
 *
 * マルチプロセス数と配列数からブロック数を計算する.
 *
 * Argument:
 *   proc_num       : I : マルチプロセス数
 *   seq_num        : I : データセットごとの配列数
 *                        [0]:1st DataSet,[1]:2nd DataSet
 *   block_num      : O : ブロック数
 *                        [0]:1st DataSet,[1]:2nd DataSet
 *
 * Return:
 *   マルチ処理数
 *
 *************************************************/
int blockutil::calcBlockNum
(
    unsigned int proc_num,
    TYPE_INDEX *seq_num,
    TYPE_INDEX *block_num
)
{
    double     r_ratio_seqnum;
    TYPE_INDEX ratio_seqnum;
    double     r_block_num[2];
    int idx;
    unsigned int real_proc_num = proc_num;
    unsigned int virt_proc_num = proc_num;

    // ブロック数を初期化する。
    block_num[0]   = 0;
    block_num[1]   = 0;
    r_block_num[0] = 0;
    r_block_num[1] = 0;

    // ブロック数を求める。
    if (isMapMode) {
        if (seq_num[0] > seq_num[1]) {
            r_ratio_seqnum = seq_num[0] / seq_num[1];
            idx = 0;
        }
        else {
            r_ratio_seqnum = seq_num[1] / seq_num[0];
            idx = 1;
        }

        // 配列数比の端数を切り上げる。
        ratio_seqnum = int(ceil(r_ratio_seqnum));

        // 配列数の比がプロセス数以上の場合は、配列数が少ない方のブロック数=1,
        // 配列数が大きい方のブロック数=プロセス数とする。
        if (ratio_seqnum >= proc_num) {
            block_num[!idx] = 1;
            block_num[idx]  = proc_num;
        }
        else {
            vector<unsigned int> ds;
            calcDivisor(proc_num, &ds);
            calcBlockNumForMap(proc_num, &ds, seq_num, block_num);
        }
    }
    // isMapMode=FALSEのとき
    else {
        calcBlockForNormal(proc_num, &real_proc_num, &virt_proc_num, &block_num[0]);
    }

    return real_proc_num;
}

/**************************************************
 *
 * １つのデータセットの場合（オプションなし、-p）の場合でのブロック計算関数
 *
 * 指定したマルチ処理数からブロック数を算出する.
 *
 * Argument:
 *   proc_num       : I : マルチ処理数
 *   real_proc_num  : O : 計算時に実行されるマルチ処理数
 *   virt_proc_num  : O : 複数ブロック組合せの同時処理を実行しないとした場合のマルチ処理数
 *   block_num   : O : ブロック数
 *
 * Return:
 *   常に0
 *
 *************************************************/
int blockutil::calcBlockForNormal(unsigned int proc_num, unsigned int *real_proc_num, unsigned int *virt_proc_num, TYPE_INDEX *block_num)
{

    unsigned int tmp_real_proc_num = 0;
    unsigned int tmp_virt_proc_num = 0;
    unsigned int tmp_block_num = 1;

    while(1) {
        tmp_virt_proc_num = tmp_block_num * (tmp_block_num + 1)/2;

        if (isMergeBlock) {
            tmp_real_proc_num = tmp_virt_proc_num - (int)(tmp_block_num/2) * 2;
        } else {
            tmp_real_proc_num = tmp_virt_proc_num ;
        }
        //cerr << tmp_real_proc_num << ":" << tmp_virt_proc_num << endl;

        if (proc_num <= tmp_real_proc_num) {
            break;
        }

        tmp_block_num++;

        if (tmp_real_proc_num >= LIMIT_PROCS) {
            break;
        }
    }

    *real_proc_num = tmp_real_proc_num;
    *virt_proc_num = tmp_virt_proc_num;
    *block_num =tmp_block_num;

    return 0;
}
/**************************************************
 *
 * 約数計算関数
 *
 * 指定した数の約数を計算する.
 *
 * Argument:
 *   num       : I : 数値
 *   ds        : O : 約数
 *
 * Return:
 *   常に0
 *
 *************************************************/
int blockutil::calcDivisor(unsigned int num, vector<unsigned int> * ds)
{
    for (unsigned int i = 1; i <= num; i++) {
        if (num % i == 0)
        {
            ds->push_back(i);
        }
    }
    return 0;
}

/**************************************************
 *
 * ２つのデータセットの場合（オプションなし、-p）の場合でのブロック計算関数
 *
 * 指定したマルチ処理数からブロック数を算出する.
 *
 * Argument:
 *   proc_num       : I : マルチ処理数
 *   ds             : I : ブロック数の候補値リスト
 *   seq_num        : I : データセットごとの配列数
 *                        [0]:1st DataSet,[1]:2nd DataSet
 *   block_num      : O : ブロック数
 *                        [0]:1st DataSet,[1]:2nd DataSet
 *
 * Return:
 *   常に0
 *
 *************************************************/
int blockutil::calcBlockNumForMap(unsigned int proc_num, vector<unsigned int> *ds, TYPE_INDEX *seq_num, TYPE_INDEX *block_num)
{
    unsigned int block_num1 = 0;
    unsigned int block_num2 = 0;

    TYPE_INDEX num_of_seq_in_block1 = 0;
    TYPE_INDEX num_of_seq_in_block2 = 0;

    unsigned int diff = INT_MAX;
    unsigned int tmp_diff = INT_MAX;

    vector<unsigned int>::iterator it = ds->begin();
    while( it != ds->end() )
    {
        //cout << *it << endl;

        block_num1 = *it;
        block_num2 = proc_num / block_num1;

        num_of_seq_in_block1 = seq_num[0]/block_num1;
        num_of_seq_in_block2 = seq_num[1]/block_num2;

        tmp_diff = labs((int)num_of_seq_in_block1 - (int)num_of_seq_in_block2);

//        cout << diff << ":" << tmp_diff << ":" << num_of_seq_in_block1 << ":" << num_of_seq_in_block2 << ":" << block_num1 << ":" << block_num2 << endl;
        if (tmp_diff < diff) {
            block_num[0] = block_num1;
            block_num[1] = block_num2;
            diff = tmp_diff;
        }
/*
        if (block_num1 > block_num2)
        {
            break;
        }
        cout << block_num1 << ":" << block_num2 << endl;
*/
        ++it;
    }

    return 0;
}


/**************************************************
 *
 * multisortオプション設定関数
 *
 * multisortに渡すオプション情報を設定する.
 *
 * Argument:
 *   fbi              : I : 1st DataSet のブロック情報格納Vector
 *   sbi              : I : 2nd DataSet のブロック情報格納Vector
 *   BlockCombiVector : O : 分割ブロックの組合せ別オプション情報
 *
 * Return:
 *   常に0
 *
 *************************************************/
int blockutil::setBlockCombination
(
    vector<BlockInformation> &fbi,
    vector<BlockInformation> &sbi,
    vector<BlockCombination> &BlockCombiVector
)
{
    // Comparing two datasets で指定された場合
    if (isMapMode) {
        vector<BlockInformation>::iterator it_fbi = fbi.begin();

        while (it_fbi != fbi.end()) {
            TYPE_INDEX fst_head  = (*it_fbi).head;
            TYPE_INDEX fst_size  = (*it_fbi).size;
            TYPE_BYTE_SIZE fst_byte_head  = (*it_fbi).byte_head;
            TYPE_BYTE_SIZE fst_byte_size  = (*it_fbi).byte_size;

            vector<BlockInformation>::iterator it_sbi = sbi.begin();

            while (it_sbi != sbi.end()) {
                BlockCombination bc;

                bc.snd_head      = (*it_sbi).head;
                bc.snd_size      = (*it_sbi).size;
                bc.snd_byte_head = (*it_sbi).byte_head;
                bc.snd_byte_size = (*it_sbi).byte_size;

                bc.fst_head      = fst_head;
                bc.fst_size      = fst_size;
                bc.fst_byte_head = fst_byte_head;
                bc.fst_byte_size = fst_byte_size;

                bc.isMapMode     = true;   // co.isMapMode;
                bc.isPartialMode = false; // !co.isPartialMode

                // multisort に渡すブロック情報オプション値を設定する。
                BlockCombiVector.push_back(bc);

                // 2nd DataSet のブロック情報のイテーションをインクリメントして
                // Vector のポインタを次のブロック情報に移動する。
                it_sbi ++;
            }

            // 1st DataSet のブロック情報のイテーションをインクリメントして
            // Vector のポインタを次のブロック情報に移動する。
            it_fbi ++;
        }
    }

    // Finding pairs from some part of an input sequences で指定された場合は総当たり
    else {
        // fst_block と snd_block のブロックの組合せに対応したインデックス情報を作成する。
        vector < vector <TYPE_INDEX> > bc_info;
        bc_info.resize( (TYPE_INDEX)fbi.size() );
        for(TYPE_INDEX i=0 ; i < (TYPE_INDEX)bc_info.size(); i++ ) bc_info[i].resize((TYPE_INDEX)fbi.size());

        TYPE_INDEX  bc_index  = 0;

        for(TYPE_INDEX i=0 ; i<(TYPE_INDEX)fbi.size() ; ++i) {

            // 最初のブロック情報を設定する。
            TYPE_INDEX fst_head = fbi[i].head;
            TYPE_INDEX fst_size = fbi[i].size;
            TYPE_BYTE_SIZE fst_byte_head = fbi[i].byte_head;
            TYPE_BYTE_SIZE fst_byte_size = fbi[i].byte_size;

            TYPE_INDEX snd_head = 0;
            TYPE_INDEX snd_size = 0;
            TYPE_BYTE_SIZE snd_byte_head = 0;
            TYPE_BYTE_SIZE snd_byte_size = 0;

            // 最初のブロック同士の組合せは、2nd DataSet情報を指定しない。
            BlockCombination bc;

            bc.fst_head      = fst_head;
            bc.fst_size      = fst_size;
            bc.fst_byte_head      = fst_byte_head;
            bc.fst_byte_size      = fst_byte_size;
            bc.snd_head      = snd_head;
            bc.snd_size      = snd_size;
            bc.snd_byte_head      = snd_byte_head;
            bc.snd_byte_size      = snd_byte_size;

            bc.isMapMode     = false;  // !co.isMapMode;
            bc.isPartialMode = true;  // co.isPartialMode;

            // multisort に渡すブロック情報オプション値を設定する。
            BlockCombiVector.push_back(bc);

            // ブロックの組合せに対応したインデックス情報にインデックス値を設定する.
            bc_info[i][i] = bc_index;

            bc_index ++;

            for(TYPE_INDEX j=i+1 ; j < (TYPE_INDEX)fbi.size() ; ++j) {
                bc.snd_head      = fbi[j].head;
                bc.snd_size      = fbi[j].size;
                bc.snd_byte_head      = fbi[j].byte_head;
                bc.snd_byte_size      = fbi[j].byte_size;

                bc.isMapMode     = true;  // co.isMapMode;
                bc.isPartialMode = false; // !co.isMapMode;

                // multisort に渡すブロック情報オプション値を設定する。
                BlockCombiVector.push_back(bc);

                // ブロックの組合せに対応したインデックス情報にインデックス値を設定する.
                bc_info[i][j]    = bc_index;

                bc_index ++;
            }
        }
//cerr << "Before BC" << endl;
//        debugprint(fbi, sbi, BlockCombiVector);
//cerr << "----" << endl;
        // 隣接するブロック同士のブロック情報を集約する。
        if (isMergeBlock) {
            mergeBlockCombination((TYPE_INDEX)fbi.size(), bc_info, BlockCombiVector);
        }
    }

    // 常に 0 を返す。
    return 0;
}

/**************************************************
 *
 * 隣接ブロック間領域集約関数
 *
 * 隣接するブロック間の領域を集約し、１つのブロックとして再設定する.
 *
 * b=1stDataSetのブロック数
 *
 * X:対象としない組合せ, +:集約対象としない(隣接しない)組合せ
 * *:集約対象となる(隣接する)組合せ, -:集約後削除される組合せ
 *
 *  b=2,3    i       i
 *          01      01
 *        j0** -> j0*-
 *         1X*     1X-
 *
 *  b=4,5    i         i
 *          0123      0123
 *        j0**++ -> j0*-++
 *         1X*++     1X-++
 *         2XX**     2XX*-
 *         3XXX*     3XXX-
 *
 * Argument:
 *   block_num        : I   : ブロック数
 *   bc_info          : I   : 隣接ブロック情報のインデックス値格納二次元配列
 *   BlockCombiVector : I/O : 分割ブロック組合せ情報
 *
 * Return:
 *   常に0
 *
 *************************************************/
int blockutil::mergeBlockCombination
(
    TYPE_INDEX                     block_num,
    vector < vector <TYPE_INDEX> > &bc_info,
    vector<BlockCombination>       &BlockCombiVector
)
{
    // 基準ブロック数
    TYPE_INDEX   std_block_num = block_num;

    // 奇数ならマイナス１で偶数にする。
    if (std_block_num%2) std_block_num --;

    // 1ブロックにまとめる領域の個数
    TYPE_INDEX  mrg_area_num = std_block_num / 2;

    // 隣接するブロックの配列数
    TYPE_INDEX  nbr_fst_size;
    TYPE_BYTE_SIZE  nbr_fst_byte_size;

    for (TYPE_INDEX k=mrg_area_num-1 ; k>=0 ; k--) {
        for(TYPE_INDEX i=1 ; i>=0 ; i--) {
            for (TYPE_INDEX j=1 ; j>=0 ; j-- ) {
                if (i>j) continue;
                // 基準ブロックに該当する場合は、隣接する2つのブロックを集約する。
                if ((i==0) && (j==0)) {
                    TYPE_INDEX mrg_block_idx = bc_info[k*2][k*2];

                    // 基準ブロックのfst_sizeに,隣接するブロックのfst_sizeを加算する。
                    BlockCombiVector[mrg_block_idx].fst_size += nbr_fst_size;
                    BlockCombiVector[mrg_block_idx].fst_byte_size += nbr_fst_byte_size;
                }
                else {
                    if ((i==1) && (j==1)) {
                        // 隣接するブロックのfst_sizeを取り出す。
                        TYPE_INDEX  nbr_block_idx = bc_info[k*2+1][k*2+1];
                        nbr_fst_size  = BlockCombiVector[nbr_block_idx].fst_size;
                        nbr_fst_byte_size  = BlockCombiVector[nbr_block_idx].fst_byte_size;
                    }
                    // 隣接するブロックのVector要素を削除する。
                    BlockCombiVector.erase(BlockCombiVector.begin() + bc_info[k*2+i][k*2+j]);
                }
            }
        }
    }

    // 常に0を返す。
    return 0;
}

/**************************************************
 *
 * 配列数カウント関数（Fasta形式）
 *
 * Fasta形式のファイルに対し、配列すう　をカウントする.
 *
 * Argument:
 *   *filename : I : ファイルパス
 *
 * Return:
 *   配列数
 *
 *************************************************/
int blockutil::countSeqForFasta(
    string *filename
)
{
    FILE *infp;
    char line[512];
    TYPE_INDEX seq_count = 0;
    infp = fopen(filename->c_str(), "r");
    while(fgets(line, sizeof(line), infp) ) {
        if(strncmp(line, ">", 1) != 0) {
            continue;
        }
        seq_count++;
    }

    fclose(infp);
    return seq_count;
}

/**************************************************
 *
 * 配列数カウント関数（Fastq形式）
 *
 * Fastq形式のファイルに対し、配列すう　をカウントする.
 *
 * Argument:
 *   *filename : I : ファイルパス
 *
 * Return:
 *   配列数
 *
 *************************************************/
int blockutil::countSeqForFastq(
    string *filename
)
{
    TYPE_BYTE_SIZE byte_pos = 0;
    TYPE_BYTE_SIZE byte_size = 0;

    FILE *infp;
    char line[512];
    TYPE_INDEX seq_count = 0;
    TYPE_INDEX row_count = 0;
    infp = fopen(filename->c_str(), "r");

    while(fgets(line, sizeof(line), infp) ) {

        if (strrchr(line, '\n')) {
            row_count++;
        } else {
            continue;
        }

        if(row_count % 4 != 0) {
            continue;
        }

        seq_count++;
    }

    fclose(infp);
    return seq_count;
}

/**************************************************
 *
 * ブロック情報配列数->バイト数変換関数（Fasta形式）
 *
 * Fasta形式のファイルに対し、ブロック情報に格納されている配列位置、配列数から
 * ファイル内のバイト位置、バイトサイズを返す.
 *
 * Argument:
 *   fbi       : I : 1st DataSet のブロック情報格納Vector
 *   sbi       : I : 2nd DataSet のブロック情報格納Vector
 *   *filename : I : ファイルパス
 *
 * Return:
 *   配列数
 *
 *************************************************/
int blockutil::convertBytesForFasta(
    vector<BlockInformation> &fbi,
    vector<BlockInformation> &sbi,
    string *filename
)
{
    TYPE_BYTE_SIZE byte_pos = 0;
    TYPE_BYTE_SIZE offset_byte_pos = 0;
    TYPE_BYTE_SIZE byte_size = 0;
    FILE *infp;
    char line[512];
    TYPE_INDEX seq_count = 0;
    TYPE_INDEX offset_seq_count = -1;
    TYPE_INDEX bc_index = 0;
    infp = fopen(filename->c_str(), "r");

    bool isFst = true;
    bool isSnd = sbi.size() != 0;

    while(fgets(line, sizeof(line), infp) ) {

        int len = strlen(line);

        if(strncmp(line, ">", 1) != 0) {
            byte_pos += len;
            continue;
        }

        // 最初のブロックのオフセットまで
        if (isFst && fbi[0].head == seq_count) {
            fbi[bc_index].byte_head = byte_pos;
            if (bc_index < fbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 中間ブロック
        else if (isFst && fbi[bc_index].head == seq_count) {
            fbi[bc_index].byte_head = byte_pos;
            fbi[bc_index-1].byte_size = byte_pos - fbi[bc_index-1].byte_head;
            if (bc_index < fbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 最後のブロックの終端
        else if (isFst && bc_index == fbi.size()-1 && fbi[bc_index].head + fbi[bc_index].size == seq_count) {
            fbi[bc_index].byte_size = byte_pos - fbi[bc_index].byte_head;
            if (!isSnd) {
                seq_count++;
                break;
            }
            isFst = false;
            bc_index = 0;
            byte_pos += len;
        }

        // 最初のブロックのオフセットまで
        else if (isSnd && sbi[0].head == seq_count) {
            sbi[bc_index].byte_head = byte_pos;
            if (bc_index < sbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 中間ブロック
        else if (isSnd && sbi[bc_index].head == seq_count) {
            sbi[bc_index].byte_head = byte_pos;
            sbi[bc_index-1].byte_size = byte_pos - sbi[bc_index-1].byte_head;
            if (bc_index < sbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 最後のブロックの終端
        else if (isSnd && bc_index == sbi.size()-1 && sbi[bc_index].head + sbi[bc_index].size == seq_count) {
            sbi[bc_index].byte_size = byte_pos - sbi[bc_index].byte_head;
            seq_count++;
            break;
        } else {
            byte_pos += len;
        }

        seq_count++;

    }

    // 最後のブロックの終端
    if (isFst && bc_index == fbi.size()-1 && fbi[bc_index].head + fbi[bc_index].size == seq_count) {
        fbi[bc_index].byte_size = byte_pos - fbi[bc_index].byte_head;
    }

    if (isSnd && bc_index == sbi.size()-1 && sbi[bc_index].head + sbi[bc_index].size == seq_count) {
        sbi[bc_index].byte_size = byte_pos - sbi[bc_index].byte_head;
    }

    fclose(infp);
    return seq_count;

}

/**************************************************
 *
 * ブロック情報配列数->バイト数変換関数（Fastq形式）
 *
 * Fastq形式のファイルに対し、ブロック情報に格納されている配列位置、配列数から
 * ファイル内のバイト位置、バイトサイズを返す.
 *
 * Argument:
 *   fbi       : I : 1st DataSet のブロック情報格納Vector
 *   sbi       : I : 2nd DataSet のブロック情報格納Vector
 *   *filename : I : ファイルパス
 *
 * Return:
 *   配列数
 *
 *************************************************/
int blockutil::convertBytesForFastq(
    vector<BlockInformation> &fbi,
    vector<BlockInformation> &sbi,
    string *filename
)
{
    TYPE_BYTE_SIZE byte_pos = 0;
    TYPE_BYTE_SIZE offset_byte_pos = 0;
    TYPE_BYTE_SIZE byte_size = 0;
    FILE *infp;
    char line[512];
    TYPE_INDEX seq_count = 0;
    TYPE_INDEX row_count = 0;
    TYPE_INDEX offset_seq_count = -1;
    TYPE_INDEX bc_index = 0;
    infp = fopen(filename->c_str(), "r");

    bool isFst = true;
    bool isSnd = sbi.size() != 0;

    while(fgets(line, sizeof(line), infp) ) {

        int len = strlen(line);

        if(row_count % 4 != 0) {
        //if(strncmp(line, ">", 1) != 0) {
            byte_pos += len;
            if (strrchr(line, '\n')) {
                row_count++;
            }
            continue;
        }

        if (strrchr(line, '\n')) {
            row_count++;
        }

        // 最初のブロックのオフセットまで
        if (isFst && fbi[0].head == seq_count) {
            fbi[bc_index].byte_head = byte_pos;
            if (bc_index < fbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 中間ブロック
        else if (isFst && fbi[bc_index].head == seq_count) {
            fbi[bc_index].byte_head = byte_pos;
            fbi[bc_index-1].byte_size = byte_pos - fbi[bc_index-1].byte_head;
            if (bc_index < fbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 最後のブロックの終端
        else if (isFst && bc_index == fbi.size()-1 && fbi[bc_index].head + fbi[bc_index].size == seq_count) {
            fbi[bc_index].byte_size = byte_pos - fbi[bc_index].byte_head;
            if (!isSnd) {
                seq_count++;
                break;
            }
            isFst = false;
            bc_index = 0;
            byte_pos += len;
        }

        // 最初のブロックのオフセットまで
        else if (isSnd && sbi[0].head == seq_count) {
            sbi[bc_index].byte_head = byte_pos;
            if (bc_index < sbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 中間ブロック
        else if (isSnd && sbi[bc_index].head == seq_count) {
            sbi[bc_index].byte_head = byte_pos;
            sbi[bc_index-1].byte_size = byte_pos - sbi[bc_index-1].byte_head;
            if (bc_index < sbi.size()-1) {
                bc_index++;
            }
            byte_pos += len;
        }

        // 最後のブロックの終端
        else if (isSnd && bc_index == sbi.size()-1 && sbi[bc_index].head + sbi[bc_index].size == seq_count) {
            sbi[bc_index].byte_size = byte_pos - sbi[bc_index].byte_head;
            seq_count++;
            break;
        } else {
            byte_pos += len;
        }

        seq_count++;

    }

    // 最後のブロックの終端
    if (isFst && bc_index == fbi.size()-1 && fbi[bc_index].head + fbi[bc_index].size == seq_count) {
        fbi[bc_index].byte_size = byte_pos - fbi[bc_index].byte_head;
    }

    if (isSnd && bc_index == sbi.size()-1 && sbi[bc_index].head + sbi[bc_index].size == seq_count) {
        sbi[bc_index].byte_size = byte_pos - sbi[bc_index].byte_head;
    }

    fclose(infp);
    return seq_count;

}

/**************************************************
 *
 * デバッグ：ブロック情報出力関数
 *
 * ブロック情報を出力する.
 *
 * Argument:
 *   bi       : I : ブロック情報格納Vector
 *
 * Return:
 *   
 *
 *************************************************/
void blockutil::debugPrintBlockInformation
(
    vector<BlockInformation> &bi
)
{
    int blk_no = 0;

    vector<BlockInformation>::iterator it_bi = bi.begin();

    // 1st DataSetのBlock分割データを表示する。
    cerr << "BlockInfoVector:"  << endl;
    cerr << "BlockNo\t" << "head\t" << "size\t" << "byte_head\t" << "byte_size" << endl;

    while (it_bi != bi.end()) {
        cerr << blk_no << "\t" << (*it_bi).head << "\t" << (*it_bi).size << "\t" << (*it_bi).byte_head << "\t" << (*it_bi).byte_size << endl;

        it_bi++;
        blk_no ++;
    }

    cerr << "---" << endl;
}
