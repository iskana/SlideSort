#include <time.h>
#include "parallelslidesort.h"

/**************************************************
 *
 * parallelslidesort class
 *
 * multisort並列実行クラス.
 *
 *  Created on: 2012/12/05
 *
 *************************************************/

/**************************************************
 *
 * コールバック関数設定関数.
 *
 * コールバック関数を設定する.
 *
 * Argument:
 *   fptr:コールバック関数
 *
 * Return:
 *   処理に成功した場合:0、失敗した場合:1
 *
 *************************************************/

int parallelslidesort::setResGetFuncPtr(GETRESFUNC fptr) {
    return multisort::setResGetFuncPtr(fptr);
}

/**************************************************
 *
 * コマンドライン引数取得関数
 *
 * multisort の getParam をオーバーロードして関数を実行する.
 *
 * Argument:
 *   argc : I : コマンドライン引数の数
 *   argv : I : コマンドライン引数
 *
 *
 *************************************************/
int parallelslidesort::getParam(int argc, char **argv)
{
    char  msg[512];

    memset(msg, 0, sizeof(msg));

    // parallelslidesort用のオプション値を取得
    if (!setPSSOptions(argc, argv, msg)) return dieWithMsg(msg);

    cmlOptions inco;
    inco.setCmlOptions(argc, argv);
    inco.isSortOrgSeq=!inco.isMapMode;
    int ret = multisort::getParam(inco);
    return ret;
}

/**************************************************
 *
 * 処理実行関数.
 *
 * pssExecutorを並列実行する.
 *
 * Argument:
 *
 * Return:
 *   処理に成功した場合:0、失敗した場合:1
 *
 *************************************************/
int parallelslidesort::exec()
{
    time_t start, end;
    start = time(NULL);

    // ブロック情報の Vector を定義
    vector<BlockInformation> fbi;
    vector<BlockInformation> sbi;
    vector<BlockCombination> bc;

    char    msg[512];

    memset(msg, 0, sizeof(msg));

    // -mt または -mp オプションが設定されている場合は、parallelslidesort を実行する。
    if (isMTOptSet||isMPOptSet) {
        // オプション値をもとにブロック分割を行い、ブロック情報を取得する。
        if (!getBlockInfo(fbi, sbi, msg)) return dieWithMsg(msg);

        // 分割したブロック情報の組合せ毎に、multisort に渡す引数を設定する。
        // fbi : Input  : 1st DataSetのブロック情報
        // sbi : Input  : 2nd DataSetのブロック情報
        // bc  : Output : 組合せ別ブロック情報
        blockutil bu;
        bu.isMapMode = co.isMapMode;
        bu.setBlockCombination(fbi, sbi, bc);

#ifdef DEBUG_PSS
cerr << "mt:" << mt << endl;
cerr << "mr:" << mr << endl;
cerr << "mp:" << mp << endl;
cerr << "bc.size():" << (int)bc.size() << endl;
        debugprint(fbi, sbi, bc);
#endif
//printHeader(&sq);

		omp_set_num_threads(mt); // マルチスレッド数（オプションで指定）を設定する
		omp_init_lock(&myLock); // ロックの生成

		if (!co.isMapMode) {
			initMapRanges();
			//先頭よりn番目から、s個の配列を一つの領域として設定する。
			//setMapRanges(n,s);

			for (int i = 0; i < (int)bc.size(); ++i) {
				setMapRanges(bc[i].fst_head, bc[i].fst_size, SET_REGION_START+i);
			}
		}

		int num_of_st_blc = co.distance+1;

#pragma omp parallel for
		for (int cnt = 0; cnt < (int)bc.size()*num_of_st_blc; ++cnt) // Vector bc の要素数分ループ
		{

			int i = cnt/num_of_st_blc;
			int st_blc = cnt%num_of_st_blc;

			//            printf("Thread ID:%d, index:%d\n", omp_get_thread_num(), i);
			//            fprintf(stderr, "parallelsort::exec ; Thread ID:%d, index:%d\n", omp_get_thread_num(), i);
			pssExecutor psse;

			psse.setBC(&bc[i]);
			psse.setResGetFuncPtr(grfptr);
			psse.getCmlOptions(co);

			psse.attachSeq(sq);
			psse.attachBox(bx);
			psse.attachCharTable(ct);
			psse.attachBitstream(bs, block_mask, lower_block_mask);
			psse.free_vals_automatic=false;

			// for parallelizing by sort
			psse.co.isDevSort = true;
			psse.co.startBlc = st_blc;

			// -o オプションを無効にする
			//            psse.co.outputfile    = false;
			psse.co.isMapMode = bc[i].isMapMode;
			if(bc[i].isMapMode){
				//先頭よりn1番目からs1個の配列と、先頭よりn2番目からs2個の配列を一つのデータセットとみなし、それに対して計算をする。
				//setSortRange(n,s)
				psse.setSortRange(bc[i].fst_head,bc[i].fst_size,bc[i].snd_head,bc[i].snd_size);
			}else{
				//先頭よりn番目から、s個の配列に対して計算をする。
				//setSortRange(n,s)
				psse.setSortRange(bc[i].fst_head,bc[i].fst_size);
			}

			psse.setMyLock(&myLock);
			psse.exec();

			psse.detachBitstream();
			psse.detachCharTable();
			psse.detachBox();
			psse.detachSeq();
			// スレッド毎の類似配列ペア数を類似配列ペアの合計値に加算
			num_of_similar_pairs += psse.num_of_similar_pairs;
		}
		cerr<<"--\n"<<num_of_similar_pairs<<" pairs were found in total.\n";
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

/**************************************************
 *
 * シングルスレッドでの解析処理実行関数.
 *
 * Argument:
 *
 * Return:
 *   multisort::exec()の戻り値
 *
 *************************************************/
int parallelslidesort::ss_main()
{
    free_vals_automatic=false;
    int ret = multisort::exec();
    return (ret);
}

/**************************************************
 *
 * デバッグ：ブロック情報出力関数.
 *
 * 引数で指定したブロック情報の内容を出力する.
 *
 * Argument:
 *   fbi  : I : 1st DataSetのブロック情報格納用Vector
 *   sbi  : I : 2nd DataSetのブロック情報格納用Vector
 *   bc   : I : ブロック組合せい情報格納用Vector
 *
 * Return:
 *
 *************************************************/
void  parallelslidesort::debugprint
(
    vector<BlockInformation> &fbi,
    vector<BlockInformation> &sbi,
    vector<BlockCombination> &bc
)
{

    int blk_no = 0;

    vector<BlockInformation>::iterator it_fbi = fbi.begin();

    // 1st DataSetのBlock分割データを表示する。
    cerr << "FstBlockInfoVector:"  << endl;
    cerr << "FstBlockNo\t" << "fst_head\t" << "fst_size" << endl;

    while (it_fbi != fbi.end()) {
        cerr << blk_no << "\t" << (*it_fbi).head << "\t" << (*it_fbi).size << endl;

        it_fbi++;
        blk_no ++;
    }

    cerr << "---" << endl;

    vector<BlockInformation>::iterator it_sbi = sbi.begin();

    blk_no = 0;

    // 2nd DataSetのBlock分割データを表示する。
    cerr << "SndBlockInfoVector:"  << endl;
    cerr << "SndBlockNo\t" << "snd_head\t" << "snd_size" << endl;

    while (it_sbi != sbi.end()) {
        cerr << blk_no << "\t" << (*it_sbi).head << "\t" << (*it_sbi).size  << endl;

        it_sbi++;
        blk_no ++;
    }

    cerr << "---" << endl;

    vector<BlockCombination>::iterator it_bc = bc.begin();

    int combi_no = 0;

    // ブロックの組合せデータを表示する。
    cerr << "BlockCombination:"  << endl;
    cerr << "Combi No\t"
         << "isMapMode\t"
         << "isPartialMode\t"
         << "fst_head\t"
         << "fst_size\t"
         << "snd_head\t"
         << "snd_size" << endl;

    while (it_bc != bc.end()) {
        cerr << combi_no << "\t"
             << (*it_bc).isMapMode << "\t"
             << (*it_bc).isPartialMode << "\t"
             << (*it_bc).fst_head << "\t"
             << (*it_bc).fst_size << "\t"
             << (*it_bc).snd_head << "\t"
             << (*it_bc).snd_size << endl;

        combi_no ++;
        it_bc ++;
    }
    cerr << "---" << endl;
}

/**************************************************
 *
 * setPSSOptions
 *
 * コマンドライン引数からのオプション取得関数
 *
 * コマンドライン引数から、-mt,-mr,-mp の各オプション値
 * を取得し、クラス変数に格納する.
 *
 * Argument:
 *   argc           : I : 引数の数
 *   argv           : I : コマンドライン引数
 *   msg            : O : エラー発生時のメッセージ
 *                        (エラー発生時のみ設定)
 *
 * Return:
 *   エラーが無い場合:true, エラーがある場合:false
 *
 *************************************************/
int parallelslidesort::setPSSOptions
(
    int argc,
    char **argv,
    char *msg
)
{
    char *key = "\0";
    TYPE_LABEL val;
    bool status = false;  // 初期状態は、FALSE(=引数名待ち)

    // -mt, -mr, -mp の各クラス変数を初期化する。
    isMTOptSet = false;
    isMROptSet = false;
    isMPOptSet = false;

    // -mt, -mr, -mp の各引数の値を 0 にする。
    mt = 0;
    mr = 0;
    mp = 0;

    // 引数を順番に調べ、-mt, -mr, -mp の各引数名と値の対応データを取得する。
    // i=0 はプログラム名
    int missing=0;
    string tmp;
    for(int cnt=0;cnt<argc;cnt++){
        int next = cnt+1;
        tmp.clear();
        if(argv[cnt][0]=='-'){
            switch(argv[cnt][1]){
                case 'm':
                    tmp = argv[cnt];
                    if(tmp=="-mt"){
                        if(next>=argc || argv[next][0]=='-'){
                            missing=cnt;
                            goto err_msg;
                        }
                        isMTOptSet = true;
                        mt = atoi(argv[next]);
                        strncpy(argv[cnt], "-|", 3);
                        argv[next][0] = 0;
                    }else if(tmp=="-mp"){
                        if(next>=argc || argv[next][0]=='-'){
                            missing=cnt;
                            goto err_msg;
                        }
                        isMPOptSet = true;
                        mp = atoi(argv[next]);
                        strncpy(argv[cnt], "-|", 3);
                        argv[next][0] = 0;
                    }else if(tmp=="-mr"){
                        if(next>=argc || argv[next][0]=='-'){
                            missing=cnt;
                            goto err_msg;
                        }
                        isMROptSet = true;
                        mr = atoi(argv[next]);
                        strncpy(argv[cnt], "-|", 3);
                        argv[next][0] = 0;
                    }
                    break;
            }
        }
    }

    // key に対応した値が設定されていない場合はエラーを表示する。
    if (strlen(key)) {
        sprintf(msg, "ERROR: There is a missing value for -%s", key);
        // key 変数のメモリを開放する。
        free(key);
    }

    // コマンドライン引数の整合性をチェックする。
    // -mt , -mr および -mp が指定されている場合はエラーとする。
    if (isMTOptSet && isMROptSet && isMPOptSet) {
        sprintf(msg, "ERROR: -%s, -%s and -%s option must not be set at the same time", OPT_MT, OPT_MP, OPT_MR);
    }
    else if (isMPOptSet || isMROptSet) {
        if (!isMTOptSet) sprintf(msg, "ERROR: -%s option must be set", OPT_MT);
    }

    return (strlen(msg)==0);

err_msg:
    cerr<<"ERROR: There is a missing value for "<<argv[missing]<<"\n";
    exit(1);
}

/**************************************************
 *
 * オプション指定チェック関数
 *
 **************************************************/
bool parallelslidesort::chkOpts(char *opt)
{
    bool  status;

    if ((strcmp(opt, OPT_MP)==0) ||
        (strcmp(opt, OPT_MT)==0) ||
        (strcmp(opt, OPT_MR)==0)) {
        status = true;
    }
    else {
        status = false;
    }

    return(status);
}
 
/**************************************************
 *
 * オプション取得関数 
 *
 * -mt, -mr, -mp の各オプションの値を返す。
 *
 * Arguments.
 *   char *opt : 引数名の文字列
 *
 * Return:
 *   OPT_MP, OPT_MT, OPT_MR で定義された文字列:オプションの値
 *   それ以外:0
 *
 **************************************************/
TYPE_INDEX  parallelslidesort::getPSSOptions
(
    char *opt
)
{
    TYPE_INDEX  val;

    val = 0;

    if (strcmp(opt, OPT_MP)==0)      val = mp;
    else if (strcmp(opt, OPT_MT)==0) val = mt;
    else if (strcmp(opt, OPT_MR)==0) val = mr;

    return(val);
}

/**************************************************
 *
 * メッセージ表示＆プログラム終了関数
 *
 * 引数で渡されたメッセージ文字列を表示してプログラムを終了する。
 *
 * Argument:
 *   *msg  : I : メッセージ
 *
 * Return:
 *   -1
 *
 *************************************************/
int parallelslidesort::dieWithMsg(char *msg) {
    cerr << msg << endl;
//    return -1;
    exit(1);
}

/**************************************************
 *
 * 数値判定関数
 *
 * 入力値が数値(入力文字列がすべて数字)かどうかチェックし、
 * 数字でない文字が含まれる場合は false を返す。
 * (入力値は、0以上の値であることが前提)
 *
 * Argument:
 *   str            : I : 入力文字列
 *
 * Return:
 *   入力値が数値の場合:true, 数字でない文字が含まれる場合:false
 *
 *************************************************/
bool parallelslidesort::isNumber(char *str)
{
    unsigned int j;
    bool status = true;  // 初期状態は TRUE (=数字である)

//cout << "isNumber in " << endl;

    for(j=0;j<strlen(str);j++) {
        if (isdigit(str[j])==false) {
            status = false;
            break;
        }
    }

    return(status);
}

/**************************************************
 *
 * 配列数・配列開始位置取得関数
 *
 * ブロック分割に用いる配列数と配列開始位置を取得する.
 *
 * Argument:
 *  TYPE_INDEX seq_head[2] : O : ファイル内配列開始位置
 *  TYPE_INDEX seq_num[2]  : O : DataSet別配列数
 *
 * Return:
 *   常に0
 *
 *************************************************/
int parallelslidesort::getSeqSize( TYPE_INDEX *seq_head, TYPE_INDEX *seq_num)
{
//cerr << "getSeqSize in" << endl;

    // -m オプション(Multi Mode)が指定されている場合
    if (co.isMapMode) {
        // 各データセットの配列数を格納する。
        seq_num[0] = parallelslidesort::getFstSizeForMap();
        seq_num[1] = sq.num_of_seq - seq_num[0];

        // 各データセットの配列開始位置を格納する。
        seq_head[0] = 0;
        seq_head[1] = seq_num[0];

#ifdef DEBUG_PSS
cerr << "MapMode!" << endl;
#endif
    }
    // -p オプション(Partial Mode)が指定されている場合または
    // -p, -m のいずれのオプションも指定されていない場合は
    // 配列数をメモリに取り込んだ配列数として指定する。
    else {
        seq_num[0]  = sq.num_of_seq;
        seq_head[0] = 0;
#ifdef DEBUG_PSS
        if (co.isPartialMode) cerr << "PartialMode!" << endl;
        else                  cerr << "OtherMode!" << endl;
#endif
    }

    return(0);
}

/**************************************************
 *
 * Map Mode時でのfst_size取得関数.
 *
 * seq.org_seq_mapIDのフラグを参照して、fst_sizeに相当する
 * 数を算出する.seqを生成する段階で、Nを含む配列等が
 * が削除され、コマンドオプションで与えた、fst_sizeから
 * 減少する可能性があるため、この関数を用いてseqでのfst_sizeを
 * 算出する
 *
 * Argument:
 *
 * Return:
 *   fst_size
 *
 *************************************************/
int parallelslidesort::getFstSizeForMap()
{
    int fst_size = 0;
    for (int i = 0; i < sq.num_of_seq; i++) {
        if (sq.org_seq_mapID[i] == FST_DATASET) {
            fst_size++;
        } else {
            break;
        }
    }
    return fst_size;
}

/**************************************************
 *
 * ブロック情報計算・取得関数.
 *
 * ブロック情報(ブロックあたりの配列数,データセット内配列開始位置)の設定
 * に必要なブロック毎の配列数を計算しブロック情報をVectorに格納する.
 *
 * Argument:
 *   fbi            : O : 1st DataSetのブロック情報格納用Vector
 *   sbi            : O : 2nd DataSetのブロック情報格納用Vector
 *   msg            : O : エラーメッセージ(エラー時のみ設定)
 *
 * Return:
 *   エラーが無い場合:true, エラーがある場合:false
 *
 *************************************************/
int parallelslidesort::getBlockInfo
(
    vector<BlockInformation> &fbi,
    vector<BlockInformation> &sbi,
    char *msg
)
{
    TYPE_INDEX  proc_num;
    TYPE_INDEX  seq_num[2];
    TYPE_INDEX  seq_head[2];
    TYPE_INDEX  block_num[2];

//cerr << "getBlockInfo in " << endl;

    // プロセス数を求める。
    // スレッド数が 0 の場合(-mt オプションが未指定の場合を含む)は 1 とする。
    // 倍率が 0 の場合(-mr オプションが未指定の場合を含む)は 1 とする。
    if (mt<1) mt = 1;
    if (mr<1) mr = 1;
    proc_num = (mp>0) ? mp : mt * mr;

	// for parallelization of sort
	proc_num /= (co.distance+1);

    if (proc_num > LIMIT_PROCS) {
        sprintf(msg, "ERROR: limit of proccess exceeded:%d > limit %d", (int)proc_num, LIMIT_PROCS);
        return false;
    }

    // 各データセットの指定配列数を格納する。
    getSeqSize(seq_head, seq_num);

    if (seq_num[0]<1) {
        strcpy(msg, "ERROR: no valid input sequence.\n       use -u option if sequence including unknown characters.");
        return false;
    }
    else if (co.isMapMode && (seq_num[1]<1)) {
        strcpy(msg, "ERROR: no valid input sequence in the second data set.\n       use -u option if sequence including unknown characters.");
        return false;
    }

    blockutil bu;
    bu.isMapMode = co.isMapMode;
    return bu.getBlockInfo(proc_num, seq_head, seq_num, fbi, sbi, msg);
}

