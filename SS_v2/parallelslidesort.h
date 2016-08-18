#include "mscls.h"
#include <omp.h>

class BlockCombination
{
public:
    TYPE_INDEX  fst_head;
    TYPE_INDEX  fst_size;
    TYPE_INDEX  snd_head;
    TYPE_INDEX  snd_size;
    bool isMapMode;
    bool isPartialMode;

    // for Hadoop
    TYPE_INDEX  fst_byte_head;
    TYPE_INDEX  fst_byte_size;
    TYPE_INDEX  snd_byte_head;
    TYPE_INDEX  snd_byte_size;

    BlockCombination() {}
    ~BlockCombination() {}

    // BlockCombination のデフォルトコンストラクタ
    BlockCombination(const BlockCombination& rbc) {
        fst_head = rbc.fst_head;
        fst_size = rbc.fst_size;
        snd_head = rbc.snd_head;
        snd_size = rbc.snd_size;
        isMapMode = rbc.isMapMode;
        isPartialMode = rbc.isPartialMode;
        fst_byte_head = rbc.fst_byte_head;
        fst_byte_size = rbc.fst_byte_size;
        snd_byte_head = rbc.snd_byte_head;
        snd_byte_size = rbc.snd_byte_size;
    }

    // BlockCombination のコピーコンストラクタ
    BlockCombination& Copy(const BlockCombination& rbc) {
        if(this == &rbc) {
            return *this;
        }
        fst_head = rbc.fst_head;
        fst_size = rbc.fst_size;
        snd_head = rbc.snd_head;
        snd_size = rbc.snd_size;
        isMapMode = rbc.isMapMode;
        isPartialMode = rbc.isPartialMode;

        fst_byte_head = rbc.fst_byte_head;
        fst_byte_size = rbc.fst_byte_size;
        snd_byte_head = rbc.snd_byte_head;
        snd_byte_size = rbc.snd_byte_size;
        return *this;
    }

    BlockCombination& operator =(const BlockCombination& rbc) {
		return Copy(rbc);
    }
};

class BlockInformation
{
public:
    TYPE_INDEX  head;
    TYPE_INDEX  size;

    // for Hadoop
    unsigned long long int byte_head;
    unsigned long long int byte_size;

    BlockInformation() {}
    ~BlockInformation() {}

    // BlockInformation のデフォルトコンストラクタ
    BlockInformation(const BlockInformation& rbi) {
        head = rbi.head;
        size = rbi.size;
        byte_head = rbi.byte_head;
        byte_size = rbi.byte_size;
    }

    // BlockInformation のコピーコンストラクタ
    BlockInformation& Copy(const BlockInformation& rbi) {
        if(this == &rbi) {
            return *this;
        }
        head = rbi.head;
        size = rbi.size;
        byte_head = rbi.byte_head;
        byte_size = rbi.byte_size;
        return *this;
    }

    BlockInformation& operator =(const BlockInformation& rbi) {
        return Copy(rbi);
    }
};

class blockutil
{
    public:
        bool isMapMode;
        bool isMergeBlock;
        int getBlockInfo(unsigned int proc_num, TYPE_INDEX *seq_head, TYPE_INDEX *seq_num, vector<BlockInformation> &fbi, vector<BlockInformation> &sbi, char *msg);
        int pushBlockInfo(vector<BlockInformation> &bi, TYPE_INDEX seq_head, TYPE_INDEX seq_num, TYPE_INDEX block_num);
        int calcBlockNum(unsigned int proc_num, TYPE_INDEX *seq_num, TYPE_INDEX *block_num);
        int setBlockCombination(vector<BlockInformation> &fbi, vector<BlockInformation> &sbi, vector<BlockCombination> &bc);
        int calcDivisor(unsigned int num, vector<unsigned int> * ds);
        int calcBlockNumForMap(unsigned int proc_num, vector<unsigned int> *ds, TYPE_INDEX *seq_num, TYPE_INDEX *block_num);
        int calcBlockForNormal(unsigned int proc_num, unsigned int *real_proc_num, unsigned int *virt_proc_num, TYPE_INDEX *block_num);
        int mergeBlockCombination(TYPE_INDEX block_num, vector < vector <TYPE_INDEX> > &bc_info, vector<BlockCombination> &bc);
        int countSeqForFasta(string *filename);
        int countSeqForFastq(string *filename);
        int convertBytesForFasta(vector<BlockInformation> &fbi, vector<BlockInformation> &sbi, string *filename);
        int convertBytesForFastq(vector<BlockInformation> &fbi, vector<BlockInformation> &sbi, string *filename);
        void debugPrintBlockInformation(vector<BlockInformation> &bi);
    blockutil() {
        isMergeBlock = true;
    }

};

class blockOptions
{
    private:
        void debugPrintOptions();
        void debugPrintGetSeqSize(TYPE_INDEX *seq_head, TYPE_INDEX *seq_num);
    public:
        TYPE_INDEX num_of_seq;
        TYPE_INDEX proc_num;
        TYPE_INDEX fst_head;
        TYPE_INDEX snd_head;
        TYPE_INDEX fst_size;
        TYPE_INDEX snd_size;
        long long int fst_bz;
        long long int snd_bz;
        int inputFileType;
        bool isMapMode;
        bool isPartialMode;
        string inputFileName;
        int setOptions(int argc, char **argv);
        int getSeqSize(TYPE_INDEX num_of_seq, TYPE_INDEX *seq_head, TYPE_INDEX *seq_num);

        blockOptions(){
            isMapMode = false;
            inputFileType=FORMAT_FASTA;
            fst_bz = 0;
            snd_bz = 0;
        }
};

class blockCountOptions
{
    private:
        void debugPrintOptions();
    public:
        TYPE_INDEX proc_num;
        long long int fst_bz;
        long long int snd_bz;
        int inputFileType;
        bool isMapMode;
        bool isPartialMode;
        string inputFileName;
        int setOptions(int argc, char **argv);

        blockCountOptions(){
            inputFileType=FORMAT_FASTA;
            isMapMode = false;
            fst_bz = 0;
            snd_bz = 0;
        }
};


class parallelslidesort : public multisort
{
    private:
        omp_lock_t myLock;
        double gettimeofday_sec();
        int dieWithMsg(char *msg);
        void debugprint(vector<BlockInformation> &fbi, vector<BlockInformation> &sbi, vector<BlockCombination> &bc);
    public:
        TYPE_INDEX mt;    // -mt オプションの値
        TYPE_INDEX mr;    // -mr オプションの値
        TYPE_INDEX mp;    // -mp オプションの値
        bool       isMTOptSet;  // -mt オプション設定:1, なし:false
        bool       isMROptSet;  // -mr オプション設定:1, なし:false
        bool       isMPOptSet;  // -mp オプション設定:1, なし:false
        BlockCombination *p_bc;
        int setResGetFuncPtr(GETRESFUNC fptr);
        int getParam(int argc, char **argv);
        int getParam(cmlOptions co);
        int exec();
		int exec2(); 		// test parallel
        int ss_main();
        int setPSSOptions(int argc, char **argv, char *msg);
        TYPE_INDEX getPSSOptions(char *argv);
        bool chkOpts(char *opt);
        bool isNumber(char *str);
        int getSeqSize(TYPE_INDEX *seq_head, TYPE_INDEX *seq_num);
        int getFstSizeForMap();
        int getBlockInfo(vector<BlockInformation> &fbi, vector<BlockInformation> &sbi, char *msg);
        int printHeader(seq *sq);
};

class pssExecutor : public multisort
{
    private:
        omp_lock_t *myLock_ptr;
        BlockCombination *p_bc;
    public:
        int execCallbackFunc(const char *seqid1, const char *seqid2, TYPE_INDEX index_of_seq1, TYPE_INDEX index_of_seq2, char* aln1, char* aln2, double dist, int aln_size);
        void setMyLock(omp_lock_t *myLock);
        void setBC(BlockCombination *p_bc);
};

// マルチ処理数の上限値
#define LIMIT_PROCS 1225 // ブロック数50の場合のマルチ処理数
// parallelslidesort のオプション名
#define OPT_MT       "mt"
#define OPT_MR       "mr"
#define OPT_MP       "mp"

#define TYPE_BYTE_SIZE unsigned long long

//#define DEBUG_PSS
