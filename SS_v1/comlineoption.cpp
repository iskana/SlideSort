/***********************************************/
//     Copyright (c) 2010-2013, Kana Shimizu
//         All rights reserved.
//              comlineoption.cpp
/***********************************************/

#include"mscls.h"

int cmlOptions::setCmlOptions(int &argc, char **&argv)
{
	// help msg

	string HELP_MSG = "\
slidesort -d <distance>  [options]\
\n\n\
Typical usage:\n\
./slidesort -d <distance>  -i <input_fasta_file> -o <output_file>\n\
\n\
help: \n\
  -h show this message.\n \
\n\
Important options: \n\
  -d  distance threshold \n\
  -i  input filename  (default=input.txt)\n\
  -j  input file format; F:fasta  Q:fastq \n\
  -o  output filename  (default=stdout)\n\
  -t  distance type;  E: edit-distance  H: hamming-distance (default=E)\n\
  -v  search with both original seq and reverse complement seq.\n\
  -mt number of logical processors used for calculation.\n\
\n\
Advanced options:\n\
  -a  output alignment\n\
  -c  type of input string.\n\
      DNA: DNA seq, PROTEIN: protein seq, INT: integer seq (default=DNA)\n\
  -g  gap extention cost (default=1, must be positive value. it is better to use larger value to avoid slow-down of the search.)\n\
  -G  gap open cost (default=0, must be positive value)\n\
  -k  size of sorting key\n\
  -l  do not output results. (for library version)\n\
  -m  use cross search mode. (find pairs between two datasets. first set: from fst_head to fst_head+fst_size. second set snd_head to snd_head+snd_size)\n\
  -p  use partial mode. (find pairs from fst_head to fst_head+fst_size) \n\
  -u  do not exclude sequences with unknown character. ex) n, N, Z, etc... \n\
  -V  output a same pair twice if dist(A,B)<=d and dist(A,B')<=d. B' is reverse complement of B. (use with -v)\n\
  -I  input file from STDIN. If possible, add approximate filesize by Megabyte. ex) -I 1000M (default=100M)\n\
\n\
SlideSort developed by Kana Shimizu\n\
Contact: slidesort@m.aist.go.jp \n";

	//decomposing elemetns.
	string tmp;
	key_size = 0;
	distance = -1;
	gap_limit = -1;
	gap_cost = 1.0;
	gap_open = 0;

	outputaln = false;
	outputfile=true;
	exclude_unknown_character=true;
	isMapMode=false;
	isPartialMode=false;
	isRevComp=false;
	isOutputBothStrand=false;

	inputFileName = "input.txt";
	outputFileName = "-";

	charType=DNA;
	distanceType=EDIT_DISTANCE;

	int FLAGS_fst_head=0;
	int FLAGS_fst_size=0;
	int FLAGS_snd_head=0;
	int FLAGS_snd_size=0;

	int missing=0;

	for(int cnt=0;cnt<argc;cnt++){
		int next = cnt+1;
		tmp.clear();
		if(argv[cnt][0]=='-'){
			switch(argv[cnt][1]){
				case 'a':
					outputaln=true;
					break;
				case 'c':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					tmp = argv[next];
					if(tmp=="DNA"){
						charType=DNA;
					}else if(tmp=="PROTEIN"){
						charType=PROTEIN;
					}else if(tmp=="INT"){
						charType=INPUT_INT;
					}else{
						cerr<<"ERROR: unknown character type ("<<tmp<<")\n";
						exit(1);
					}
					break;
				case 'j':
					if(argv[cnt+1][0]=='F'){
						inputFileType= FORMAT_FASTA;
					}else if(argv[cnt+1][0]=='Q'){
						inputFileType = FORMAT_FASTQ;	
					}else{
						cerr<<"ERROR: unknown inputfile type ("<<argv[cnt+1]<<")\n";
						exit(1);
					}
					break;					
				case 'd':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					distance = atoi(argv[next]);
					break;
				case 'g':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					gap_cost = atof(argv[next]);
					gap_limit = 0;
					break;
				case 'G':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					gap_open = atof(argv[next]);
					gap_limit = 0;
					break;
				case 'h':
					cerr<<HELP_MSG;
					exit(1);
					break;
				case 'i':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					tmp = argv[next]; 
					if(tmp.size()>0)
						inputFileName=tmp;
					break;
				case 'k':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					key_size=atoi(argv[cnt+1]);
					break;
				case 'l':
					outputfile=false;
					break;
				case 'o':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					tmp = argv[next]; 
					if(tmp.size()>0)
						outputFileName=tmp;
					break;
				case 't':
					if(next>=argc || argv[next][0]=='-'){
						missing=cnt;
						goto err_msg;
					}
					if(argv[cnt+1][0]=='H'){
						distanceType=HAMMING_DISTANCE;
					}else if(argv[cnt+1][0]=='E'){
						distanceType=EDIT_DISTANCE;	
					}else{
						cerr<<"ERROR: unknown distance type ("<<argv[cnt+1]<<")\n";
						exit(1);
					}
					break;
				case 'u':
					exclude_unknown_character=false;
					break;
				case 'm':
					isMapMode=true;
					break;
				case 'p':
					isPartialMode=true;
					break;
				case 'f':
					tmp = argv[cnt];
					if(tmp=="-fst_head"){
						if(next>=argc || argv[next][0]=='-'){
							missing=cnt;
							goto err_msg;
						}
						FLAGS_fst_head = atoi(argv[next]);
					}else if(tmp=="-fst_size"){
						if(next>=argc || argv[next][0]=='-'){
							missing=cnt;
							goto err_msg;
						}
						FLAGS_fst_size = atoi(argv[next]);	
					}
					break;
				case 's':
					tmp = argv[cnt];
					if(tmp=="-snd_head"){
						if(next>=argc || argv[next][0]=='-'){
							missing=cnt;
							goto err_msg;
						}
						FLAGS_snd_head = atoi(argv[next]);
					}else if(tmp=="-snd_size"){
						if(next>=argc || argv[next][0]=='-'){
							missing=cnt;
							goto err_msg;
						}
						FLAGS_snd_size = atoi(argv[next]);					
					}
					break;
				case 'v':
					isRevComp=true;
					break;
				case 'V':
					isOutputBothStrand=true;
					break;
				case 'I':
					isInputFromStdin=true;
					if(next>=argc || argv[next][0]=='-'){
					}else{
						string tmp = argv[next];
						tmp_file_size = atoi(tmp.substr(0,tmp.find("M")).c_str());
						tmp_file_size *= 1000000;
						if(tmp_file_size<=0){
							tmp_file_size = TMP_FILE_SIZE;
						}
					}
					break;
			}
		}
	}

	if(distance==-1){
		cerr<<"ERROR: set distance\n";		
		exit(1);	
	}

	if(gap_limit==0){
		/**
		if(gap_cost<1.0){
			cerr<<"ERROR: Gap extention cost must be larger than mismatch cost (1.0)\n";
			exit(1);
		}
		**/
		if(gap_cost<=0){
			cerr<<"ERROR: Gap extention cost must be larger than 0.\n";
			exit(1);
		}
		if(gap_open<0){
			cerr<<"ERROR: Gap open cost must be positive value\n";
			exit(1);
		}
#ifdef INTERNAL_GAP_OPEN
		//gap openは、リード内のみ追加。オーバーラップに関しては追加しないので。
		if(distance==(distance/gap_cost)*gap_cost){
			gap_limit = distance/gap_cost;
		}else{
			gap_limit = distance/gap_cost;
		}
#else
//	define in other part
#endif
		cerr<<gap_limit<<"  "<<gap_cost<<"\n";
	}else{
		gap_limit=distance;
	}

	if(isMapMode && isPartialMode){
		cerr<<"partial mode and map mode is exclusive.\n";
		exit(1);
	}

	if(isMapMode){
		if(FLAGS_fst_head<0 || FLAGS_snd_head<0 || FLAGS_fst_size<0 || FLAGS_snd_size<0){
			cerr<<"ERROR: values should not be minus.\n";
			exit(1);
		}
		if(FLAGS_fst_head + FLAGS_fst_size > FLAGS_snd_head){
			cerr<<"ERROR: selected data are overlapped.  (fst_head + fst_size <= snd_head) \n";
			exit(1);
		}else{
			fst_head = FLAGS_fst_head;
			snd_head = FLAGS_snd_head;
			fst_size = FLAGS_fst_size;
			snd_size = FLAGS_snd_size;
		}
	}

	if(isPartialMode){
		if(FLAGS_fst_head<0 || FLAGS_fst_size<0 ){
			cerr<<"ERROR: values should not be minus.\n";
			exit(1);
		}else{
			fst_head = FLAGS_fst_head;
			fst_size = FLAGS_fst_size;
		}
	}
	if(charType!=DNA){
		isRevComp = false;
	}

	return(0);
err_msg:
	cerr<<"ERROR: There is a missing value for "<<argv[missing]<<"\n";
	exit(1);
}

#if 0

#include"gflags.h"
#include"mscls.h"

DEFINE_bool(usage, false, "show usage of mscls");
DEFINE_string(c,"DNA","character type, DNA, PROTEIN or INT");
DEFINE_string(i,"input.txt","input filename");
DEFINE_string(o,"output.txt","output filename");
DEFINE_bool(l,false,"do not output result (for lib version).");
DEFINE_string(t,"H","H:hamming distance or E:edit distance");
DEFINE_bool(u, false, "include sequence with unknown character. ex:n, N, Z, etc...");
DEFINE_int32(d,-1,"distance of similar pairs");
DEFINE_int32(g,-1,"maximum gap length");
DEFINE_int32(k,0,"size of sorting key");
DEFINE_bool(m,false,"use mapping mode");
DEFINE_bool(p,false,"use partial mode");
DEFINE_int32(fst_head,0,"set data point of the first dataset (starting from 0)");
DEFINE_int32(snd_head,0,"set data point of the second dataset (starting from 0)");
DEFINE_int32(fst_size,0,"set data size of the first dataset (automatically select all data before snd_head if 0 is sat)");
DEFINE_int32(snd_size,0,"set data size of the second dataset (automatically select end of the file if 0 is sat)");
DEFINE_bool(a,false,"output alignment");

void cmlOptions::setCmlOptions(int &argc, char **&argv)
{
	google::SetUsageMessage(usage);
	google::ParseCommandLineFlags(&argc, &argv, true);	

	if(FLAGS_usage){
		cerr<<usage;
	}

	if(FLAGS_i.length()>0){
		inputFileName=FLAGS_i;
	}

	if(FLAGS_o.length()>0){
		outputFileName=FLAGS_o;
	}	

	if(FLAGS_t=="H"){
		distanceType=HAMMING_DISTANCE;
	}else if(FLAGS_t=="E"){
		distanceType=EDIT_DISTANCE;	
	}else{
		cerr<<"ERROR: unknown distance type ("<<FLAGS_t<<")\n";
		exit(1);
	}

	if(FLAGS_l){
		outputfile=false;
	}else{
		outputfile=true;
	}

	if(FLAGS_d==-1){
		cerr<<"ERROR: set distance\n";		
		exit(1);
	}
	distance = FLAGS_d;

	if(FLAGS_g>=0){
		if(FLAGS_g>FLAGS_d){
			cerr<<"ERROR: Gap length should be smaller than distance d\n";
			exit(1);			
		}else{
			gap_limit=FLAGS_g;
		}
	}else{
		gap_limit=FLAGS_d;
	}

	if(FLAGS_u){
		exclude_unknown_character=false;
	}else{
		exclude_unknown_character=true;
	}

	if(FLAGS_c=="DNA"){
		charType=DNA;
	}else if(FLAGS_c=="PROTEIN"){
		charType=PROTEIN;
	}else if(FLAGS_c=="INT"){
		charType=INPUT_INT;
	}else{
		cerr<<"ERROR: unknown character type ("<<FLAGS_c<<")\n";
		exit(1);
	}

	if(FLAGS_m==true && FLAGS_p){
		cerr<<"paritla mode and map mode is exclusive.\n";
		exit(1);
	}

	isMapMode=false;
	if(FLAGS_m){
		isMapMode=true;
		if(FLAGS_fst_head<0 || FLAGS_snd_head<0 || FLAGS_fst_size<0 || FLAGS_snd_size<0){
			cerr<<"ERROR: values should not be minus.\n";
			exit(1);
		}
		if(FLAGS_fst_head + FLAGS_fst_size > FLAGS_snd_head){
			cerr<<"ERROR: selected data are overlapped.  (fst_head + fst_size <= snd_head) \n";
			exit(1);
		}else{
			fst_head = FLAGS_fst_head;
			snd_head = FLAGS_snd_head;
			fst_size = FLAGS_fst_size;
			snd_size = FLAGS_snd_size;
		}
	}
	isPartialMode=false;
	if(FLAGS_p){
		isPartialMode=true;
		if(FLAGS_fst_head<0 || FLAGS_fst_size<0 ){
			cerr<<"ERROR: values should not be minus.\n";
			exit(1);
		}else{
			fst_head = FLAGS_fst_head;
			fst_size = FLAGS_fst_size;
		}
	}

	key_size=FLAGS_k;

	outputaln=false;
	if(FLAGS_a){
		outputaln=true;
	}

	//	if(allowedgap==0){
	//		distanceType=HAMMING_DISTANCE;
	//	}

}

#endif


