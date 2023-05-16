#include "olmst.h"

//#define DEBUG

#define FORMAT_SIF 0
#define FORMAT_ORG 1
#define top_N 30

Node *forest;
TYPE_SIGNED_INDEX nodesize;
int dist_thr;
TYPE_SIGNED_INDEX debug_cnt=0;
TYPE_SIGNED_INDEX max_cls;
TYPE_SIGNED_INDEX max_mst_in_fst;
int isPairOut;
int isTreeOut;
int isAlnOut;
int isDegreeOut;
int isVisualOut;
int isSTDOUT;
int isSingleOut;
ofstream pFile;
ofstream dFile;
ofstream vFile;
ofstream tFile;
char *fname=NULL;
char *pairFname=NULL;
char *degreeFname=NULL;
char *visualFname=NULL;
TYPE_SIGNED_INDEX member, cls, single;
TYPE_SIGNED_INDEX largest_tree;
char *resource_path=NULL;

#ifdef SINGLE_CORE
multisort ml;
#else
parallelslidesort ml;
#endif

TYPE_SIGNED_INDEX *invindex;

bool reverse_complement;

/**
  -a  output alignment\n\
  -f  output pair\n\
  -F  output pair file name\n\
**/

int MSTree::OutPair(const char* head1, const char* head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size){

	if(isPairOut){
		pFile<<forest[seqid1].id<<" , "<<forest[seqid2].id<<" , "<<dist<<"\n";
	}
	if(isAlnOut){
		pFile<<aln1<<"\n";
		pFile<<aln2<<"\n";
	}
	return(0);
}

int MSTree::mkMST_main(int argc, char** argv){

isPairOut=false;
isAlnOut=false;
isTreeOut=true;
isDegreeOut=false;
isVisualOut = false;
isSTDOUT = false;
isSingleOut = false;
int fmtType = FORMAT_ORG;
int pairs=0;

	        string HELP_MSG = "\
ssmst -d <distance>  [options]\
\n\n\
Typical usage:\n\
./ssmst -d <distance>  -i <input_fasta_file> -o <output_file>\n\
\n\
help: \n\
  -h show this message.\n \
\n\
Important options: \n\
  -d  distance threshold \n\
  -t  distance type  E: edit-distance  H: hamming-distance (default=E)\n\
  -i  input filename  (default=input.txt)\n\
  -j  input file format; F:fasta  Q:fastq \n\
  -o  output filename  (default=tree.txt)\n\
  -v  search with both original seq and reverse complement seq.\n\
  -mt number of logical processors used for calculation.\n\
\n\
Advanced options:\n\
  -T  type of tree format. (S:SIF, O:SlideSort original)\n\
  -c  type of input string.\n\
      DNA: DNA seq, PROTEIN: protein seq, INT: integer seq (default=DNA)\n\
  -u  do not exclude sequence with unknown character. ex) n, N, Z, etc... \n\
  -g  gap extention cost (default=1, same value as mismatch cost)\n\
  -G  gap open cost (default=0, must be positive value)\n\
  -k  size of sorting key\n\
  -p  use partial mode. (find pairs from fst_head to fst_head+fst_size) \n\
  -m  use cross search mode. (find pairs between two datasets. first set: from fst_head to fst_head+fst_size. second set snd_head to snd_head+snd_size)\n\
  -V  output a same pair twice if dist(A,B)<=d and dist(A,B')<=d. B' is reverse complement of B. (use with -v)\n\
  -I  input from STDIN. If possible, add approximate filesize by Megabyte. ex) -I 1000M (default=100M)\n\
  -Z  output to STDOUT.\n\
  -U  output singleton ID\n\
\n\
SSMST developed by Kana Shimizu\n\
Contact: shimizu-kana@aist.go.jp \n";

	char **argv_t =(char**)malloc(sizeof(char*)*(argc+1));
	string tmp;

	for(int cnt=0;cnt<argc;cnt++){
		argv_t[cnt] = argv[cnt];
		int next = cnt+1;
		tmp.clear();
		if(argv[cnt][0]=='-'){
			switch(argv[cnt][1]){
				case 'U':
					argv_t[cnt][0] = ' ';
					argv_t[cnt][1] = ' ';
					isSingleOut = true;
				break; 
				case 'Z':
					argv_t[cnt][0] = ' ';
					argv_t[cnt][1] = ' ';
					isSTDOUT = true;
					break;
				case 'h':
					cerr<<HELP_MSG;
					argv_t[cnt][1]='l';
					exit(1);
				case 'o':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					fname = argv[cnt+1];
					break;
				case 'F':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					pairFname = argv[cnt+1];
					break;
				case 'R':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					degreeFname = argv[cnt+1];
					break;
				case 'D':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					isDegreeOut=true;
					break;
				case 'r':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					isPairOut=true;
					break;
				case 'L':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					isVisualOut=true;
					break;
				case 'O':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					visualFname = argv[cnt+1];
					break;
				case 'S':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					resource_path = argv[cnt+1];
					break;
				case 'X':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					isTreeOut = false;
					break;
				case 'a':
					isAlnOut = true;
					break;
				case 'T':
					argv_t[cnt][0]=' ';
					argv_t[cnt][1]=' ';
					if(argv[cnt+1][0]=='S'){
						fmtType = FORMAT_SIF;
					}else if(argv[cnt+1][0]=='O'){
						fmtType = FORMAT_ORG;
					}else{
						cerr<<"ERROR: unknow file format of '"<<argv[cnt+1][0]<<"' \n";
						exit(1);
					}
					break;
			}
		}
	}
	resource_path="app://";
		argv_t[argc] = (char *)malloc(sizeof(char)*2);
		argv_t[argc][0] = '-';		
		argv_t[argc][1] = 'l';		
// bug:sprintf(buffer overflow)		sprintf(argv_t[argc],"-l");
		argc++;

	debug_cnt=0;
	GETRESFUNC infptr;
	infptr = &MSTree::mkForest;

	/**
	if(isTreeOut){
	//	GETRESFUNC infptr = &MSTree::VOID;
	}else if(isPairOut||isAlnOut){
		infptr = &MSTree::OutPair;
	}else{
		infptr = &MSTree::OutPair;
		cerr<<"WARRNING: No Output is set.\n";
	}
	**/

	if(isPairOut || isAlnOut){
		pFile.open(pairFname);
		if(!pFile){
			cerr<<"ERROR: output file for alignment is not found\nPlease use -F filename\n";
			exit(1);
		}
	}

	if(isDegreeOut){
		dFile.open(degreeFname);
		if(!dFile){
			cerr<<"ERROR: output file for degree is not found\nPlease use -R filename\n";
			exit(1);
		}
	}
	ml.free_vals_automatic=false;

	ml.setResGetFuncPtr(infptr);
	ml.getParam(argc, argv_t);	

	int num_of_seq = ml.sq.num_of_valid_seq;
	if(ml.co.isRevComp==true){
		reverse_complement = true;
		num_of_seq /=2;
	}else{
		reverse_complement=false;
	}
	//forest
	nodesize = num_of_seq + 1;
	//tree
	//nodesize = ml.sq.seqName.size();

	forest = (Node *)malloc(sizeof(Node)*(nodesize));
	dist_thr = ml.co.distance;
	for(TYPE_SIGNED_INDEX i=0;i<num_of_seq;i++){

		if(ml.co.isRevComp==true){
			forest[i].id = (char*)malloc(sizeof(char)*ml.sq.seqName[i*2].size()+1);
			sprintf(forest[i].id,"%s",ml.sq.seqName[i*2].c_str());
		}else{
			forest[i].id = (char*)malloc(sizeof(char)*ml.sq.seqName[i].size()+1);
			sprintf(forest[i].id,"%s",ml.sq.seqName[i].c_str());
		}

		forest[i].parent=i;
		forest[i].edge = -10;
		forest[i].rank=0;
		forest[i].member=1;
		forest[i].num_of_sod=0;
		forest[i].sod_ptr=0;
		forest[i].degree=0;
	}

	ml.exec();
	pairs = ml.num_of_similar_pairs;

	// for walrus output
	linkChildrenAsForest();
	//	outputLargeTreeForWalrus("walrus_test.txt");

	cerr<<"LINK CHILDREN ..";
	//	linkChildren();
	cerr<<"done.\nOUTPUT..";
	//	outputRoot();
	cerr<<"done\n";
	if(fname==NULL && !isSTDOUT){
		cerr<<"<WARNING> tree info will be written in tree.txt\n";
		cerr<<"use -o option for setting a filename.\n";
		fname=(char*)malloc(sizeof(char)*9);
		sprintf(fname,"tree.txt\0");
	}
	if(isTreeOut){
		if(fmtType==FORMAT_SIF){
			outputLargeTreeForSIF(fname);
		}else if(fmtType==FORMAT_ORG){
			outputLargeTreeForText(fname);
		}
	}else{
		outputLargeTreeForText(fname);
	}
	if(isPairOut || isAlnOut){
		pFile.close();
	}
	if(isDegreeOut){
		for(TYPE_SIGNED_INDEX i=0;i<num_of_seq;i++){
			dFile<<forest[i].id<<": "<<forest[i].degree<<"\n";
		}
		dFile.close();
	}

	if(isVisualOut){
		outputVisualResult(visualFname);
	}
	return(0);
}

int MSTree::outputVisualResult(char *ofname){

	TYPE_SIGNED_INDEX *toptop_Nchar;
	vFile.open(ofname);
	string tfilename(ofname);
	tfilename += ".txt";
	tFile.open(tfilename.c_str());
	
//	vFile<<"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n";
	vFile<<"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\">\n";
	vFile<<"<html>";
	vFile<<"<head>";
	vFile<<"<link rel=\"stylesheet\" href=\""<<resource_path<<"css/msg.css\" type=\"text/css\">";
	vFile<<"</head>";

	vFile<<"<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">\n";
	vFile<<"<script type=\"text/javascript\" src=\""<<resource_path<<"/import/jquery-1.4.2.js\"></script>\n";
	vFile<<"<script type=\"text/javascript\" src=\""<<resource_path<<"/import/highcharts.js\"></script>\n";

//	/** for earthquake 

	vFile<<"<script type=\"text/javascript\" src=\"http://ajax.googleapis.com/ajax/libs/jquery/1.4.2/jquery.min.js\"></script>\n";
	vFile<<"<script type=\"text/javascript\" src=\"http://www.cbrc.jp/~shimizu/SS/import/highcharts.js\"></script>\n";
// **/


	//	vFile<<"<script type=\"text/javascript\" src=\""<<resource_path<<"/highcharts/themes/gray.js\"></script>\n";
//	vFile<<"<script type=\"text/javascript\" src=\""<<resource_path<<"/highcharts/modules/exporting.js\"></script>\n";
//	vFile<<"<script src=\"app:///script/prop.js\"></script>\n";
//	vFile<<"<script type=\"text/python\" src=\"app:///script/gui.py\"></script>\n";
	vFile<<"<script src=\"app:///script/save.js\"></script>\n";

	vFile<<"<script type=\"text/javascript\">\n";

	TYPE_SIGNED_INDEX hist[10];
	TYPE_SIGNED_INDEX t_hist[10];
	TYPE_SIGNED_INDEX toptop_N[top_N];
	vector<TYPE_SIGNED_INDEX> tc;
	memset(hist,0x00,sizeof(TYPE_SIGNED_INDEX)*10);
	for(TYPE_SIGNED_INDEX cnt=0;cnt<nodesize-1;cnt++){
		if(forest[cnt].degree==0){
			hist[0]++;
		}else if(forest[cnt].degree==1){
			hist[1]++;
		}else{
			TYPE_SIGNED_INDEX h = (TYPE_SIGNED_INDEX)log10((double)forest[cnt].degree);
			if(h+2>9){
				hist[9]++;
			}else{
				hist[h+2]++;
			}
		}
	}

vFile<<"			var chart;\n";
vFile<<"			$(document).ready(function() {\n";
vFile<<"				chart = new Highcharts.Chart({\n";
vFile<<"					chart: {\n";
vFile<<"						renderTo: 'degree'\n";
vFile<<"					},\n";
vFile<<"					title: {\n";
vFile<<"						text: 'Degree of sequences.'\n";
vFile<<"					},\n";
vFile<<"					plotArea: {\n";
vFile<<"						shadow: null,\n";
vFile<<"						borderWidth: null,\n";
vFile<<"						backgroundColor: null\n";
vFile<<"					},\n";
vFile<<"					tooltip: {\n";
vFile<<"						formatter: function() {\n";
vFile<<"							return '<b>'+ this.point.name +'</b>: '+ this.y +' seqs';\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"					plotOptions: {\n";
vFile<<"						pie: {\n";
vFile<<"							allowPointSelect: true,\n";
vFile<<"							cursor: 'pointer',\n";
vFile<<"							dataLabels: {\n";
vFile<<"								enabled: true,\n";
vFile<<"								color: '#000000',\n";
vFile<<"								connectorColor: '#000000',\n";
vFile<<"								formatter: function() {\n";
vFile<<"									return '<b>'+ this.point.name +'</b>: '+ Math.round(this.y*10000/"<<(nodesize-1)<<")/100 +' %';\n";
vFile<<"								}\n";
vFile<<"							}\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"				    series: [{\n";
vFile<<"						type: 'pie',\n";
vFile<<"						name: 'Degree distribution',\n";
vFile<<"						data: [\n";
for(int cnt=9;cnt>2;cnt--){
	if(hist[cnt]>0){
		vFile.precision(2);
		vFile<<"['"<<(int)pow(10,cnt-2)<<" - "<<(int)pow(10,cnt-1)-1<<"',"<<hist[cnt]<<"],\n";
	}
}
if(hist[2]>0){
	vFile.precision(2);
	vFile<<"['"<<2<<" - "<<9<<"',"<<hist[2]<<"],\n";
}
if(hist[1]>0)
	vFile<<"['"<<"1"<<"',"<<hist[1]<<"],\n";
vFile<<"['"<<"0"<<"',"<<hist[0]<<"]\n";
vFile<<"						]\n";
vFile<<"					}]\n";
vFile<<"				});\n";
vFile<<"			});\n";


memset(t_hist,0x00,sizeof(TYPE_SIGNED_INDEX)*10);
for(int i=0;i<top_N;i++){
	toptop_N[i]=-1;
}
bool isTreeZero=true;
for(TYPE_SIGNED_INDEX cnt=0;cnt<nodesize-1;cnt++){
	if(forest[cnt].parent == nodesize-1 && forest[cnt].member>1){
		isTreeZero=false;
		if(toptop_N[top_N-1]==-1 || forest[cnt].member>forest[ toptop_N[top_N-1] ].member){
			toptop_N[top_N-1] = cnt;
			for(int i=top_N-1;i>0;i--){
				if(toptop_N[i-1]==-1 || forest[ toptop_N[i] ].member > forest[ toptop_N[i-1] ].member){
					int tmp = toptop_N[i];
					toptop_N[i] = toptop_N[i-1];
					toptop_N[i-1] = tmp; 
				}else{
					break;
				}
			}
		}
		TYPE_SIGNED_INDEX h = (TYPE_SIGNED_INDEX)log10((double)forest[cnt].member);
		if(h>9){
			t_hist[9]++;
		}else{
			t_hist[h]++;
		}
	}
}

if(isTreeZero){
	goto end_script_output;
}

vFile<<"			var chart2;\n";
vFile<<"			$(document).ready(function() {\n";
vFile<<"				chart2 = new Highcharts.Chart({\n";
vFile<<"					chart: {\n";
vFile<<"						renderTo: 'tree'\n";
vFile<<"					},\n";
vFile<<"					title: {\n";
vFile<<"						text: 'Distribution of tree size.'\n";
vFile<<"					},\n";
vFile<<"					plotArea: {\n";
vFile<<"						shadow: null,\n";
vFile<<"						borderWidth: null,\n";
vFile<<"						backgroundColor: null\n";
vFile<<"					},\n";
vFile<<"					tooltip: {\n";
vFile<<"						formatter: function() {\n";
vFile<<"							return '<b>'+ this.point.name +'</b>: '+ this.y +' trees';\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"					plotOptions: {\n";
vFile<<"						pie: {\n";
vFile<<"							allowPointSelect: true,\n";
vFile<<"							cursor: 'pointer',\n";
vFile<<"							dataLabels: {\n";
vFile<<"								enabled: true,\n";
vFile<<"								color: '#000000',\n";
vFile<<"								connectorColor: '#000000',\n";
vFile<<"								formatter: function() {\n";
vFile<<"									return '<b>'+ this.point.name +'</b>: '+ Math.round(this.y*10000/"<<cls<<")/100 +' %';\n";
vFile<<"								}\n";
vFile<<"							}\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"				    series: [{\n";
vFile<<"						type: 'pie',\n";
vFile<<"						name: 'tree size distribution',\n";
vFile<<"						data: [\n";
for(int cnt=9;cnt>0;cnt--){
	if(t_hist[cnt]>0){
		vFile<<"['"<<(int)pow(10,cnt)<<" - "<<(int)pow(10,cnt+1)-1<<"',"<<t_hist[cnt]<<"],\n";
	}
}
if(t_hist[0]>0){
	vFile<<"['"<<(int)pow(10,0)+1<<" - "<<(int)pow(10,1)-1<<"',"<<t_hist[0]<<"],\n";
}
vFile<<"						]\n";
vFile<<"					}]\n";
vFile<<"				});\n";
vFile<<"			});\n";


invindex = (TYPE_SIGNED_INDEX *)malloc(sizeof(TYPE_SIGNED_INDEX)*(nodesize-1));
for(TYPE_SIGNED_INDEX i=0;i<ml.sq.num_of_seq;i++){
	for(TYPE_SIGNED_INDEX j=0;j<ml.sq.toOrgSeqIndex[i]->size;j++){
		if(ml.co.isRevComp){
			invindex[ ml.sq.toOrgSeqIndex[i]->val[j]/2 ] = i;
		}else{
			invindex[ ml.sq.toOrgSeqIndex[i]->val[j] ] = i;		
		}
	}
}
toptop_Nchar = (TYPE_SIGNED_INDEX *)malloc(sizeof(TYPE_SIGNED_INDEX)*top_N*ml.ct.num_of_character);
memset(toptop_Nchar,0x00,sizeof(TYPE_SIGNED_INDEX)*top_N*ml.ct.num_of_character);
TYPE_SIGNED_INDEX toptop_N_z[top_N];
for(TYPE_SIGNED_INDEX cnt=0;cnt<top_N;cnt++){
	toptop_N_z[cnt]=0;
	if(toptop_N[cnt]!=-1){
		tc.push_back(toptop_N[cnt]);
		while(tc.size()>0){
			TYPE_SIGNED_INDEX next = tc[tc.size()-1];
			tc.pop_back();
			for(TYPE_SIGNED_INDEX i=ml.sq.head[ invindex[ next ] ]+ml.bx.box_set_center;i<ml.sq.head[ invindex[ next ] ]+ml.sq.max_seq_length+ml.bx.box_set_center;i++){
				if(ml.sq.nr_seq[i]<ml.ct.num_of_character){
					toptop_Nchar[ cnt*ml.ct.num_of_character + ml.sq.nr_seq[i] ]++;
					toptop_N_z[cnt]++;
				}
			}
			for(TYPE_SIGNED_INDEX i=0;i<forest[next].num_of_sod;i++){
				tc.push_back(forest[next].sod[i]);
			}
		}
	}
}

vFile<<"			var chart3;\n";
vFile<<"			$(document).ready(function() {\n";
vFile<<"				chart3 = new Highcharts.Chart({\n";
vFile<<"					chart: {\n";
vFile<<"						renderTo: 'composition',\n";
vFile<<"						defaultSeriesType: 'bar'\n";
vFile<<"					},\n";
vFile<<"					title: {\n";
vFile<<"						text: 'Composition  of top ";
if(top_N<cls){
vFile<<top_N;
}else{
vFile<<cls;
}
vFile<<" largest trees.'\n";
vFile<<"					},\n";
vFile<<"					xAxis: {\n";
vFile<<"						categories: [\n";
for(int cnt=0;cnt<top_N-1;cnt++){
	if(toptop_N[cnt]!=-1)
		vFile<<"'"<<"id:"<<toptop_N[cnt]<<" size:"<<forest[toptop_N[cnt]].member<<"',";
}
	if(toptop_N[top_N-1]!=-1)
		vFile<<"'"<<"id:"<<toptop_N[top_N-1]<<" size:"<<forest[toptop_N[top_N-1]].member<<"'";
vFile<<"]\n";

vFile<<"					},\n";
vFile<<"					yAxis: {\n";
vFile<<"						max: 100,\n";
vFile<<"						title: {\n";
vFile<<"							text: 'Composition (%)'\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"					legend: {\n";
vFile<<"						backgroundColor: '#FFFFFF',\n";
vFile<<"						reversed: true\n";
vFile<<"					},\n";
vFile<<"					tooltip: {\n";
vFile<<"						formatter: function() {\n";
vFile<<"							return ''+\n";
vFile<<"								 this.series.name +': '+ this.y +'%';\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"					plotOptions: {\n";
vFile<<"						series: {\n";
vFile<<"							stacking: 'normal'\n";
vFile<<"						}\n";
vFile<<"					},\n";
vFile<<"				        series: [\n";

for(int cnt=ml.ct.num_of_character-1;cnt>=0;cnt--){
	vFile<<"{name: '"<<ml.ct.toChar[cnt]<<"',\ndata:[ ";
	for(int i=0;i<top_N;i++){
		if(toptop_N[i]!=-1){
			vFile.precision(3);
			vFile<<(double)toptop_Nchar[i*ml.ct.num_of_character + cnt]*100/toptop_N_z[i]<<",";
		}
	}
	vFile<<"]},";
}
vFile<<"					]\n";
vFile<<"				});\n";					
vFile<<"			});\n";

end_script_output:;

	vFile<<"window.document.onkeydown = Function('if(event.keyCode == 0008) return false;');";
	vFile<<"</script>\n";

	vFile<<"<body class=\"wgbg\">";

	vFile<<"<div class=\"wlabel\">";
	vFile<<"InputFile : "<<ml.co.inputFileName<<"<br>\n";
	tFile<<"InputFile : "<<ml.co.inputFileName<<"\n";

	vFile<<"Distance : "<<ml.co.distance<<"<br>\n";
	tFile<<"Distance : "<<ml.co.distance<<"\n";

	if(ml.co.distanceType==HAMMING_DISTANCE){
		vFile<<"Distance type : Hamming<br>\n";
		tFile<<"Distance type : Hamming\n";
	}else{
		vFile<<"Distance type : Edit<br>\n";	
		tFile<<"Distance type : Edit\n";
	}

	vFile<<"Gap : extension="<<ml.co.gap_cost<<" open="<<ml.co.gap_open<<"<br>\n";
	tFile<<"Gap : extension="<<ml.co.gap_cost<<" open="<<ml.co.gap_open<<"\n";


	vFile<<"Num of Input Seqs : "<<ml.sq.num_of_seq_of_all_input<<"<br>\n";	
	tFile<<"Num of Input Seqs : "<<ml.sq.num_of_seq_of_all_input<<"\n";
	if(ml.co.exclude_unknown_character && ml.co.isRevComp==false){
		vFile<<"Num of Input Seqs with unknown_character : "<<ml.sq.num_of_seq_of_all_input-nodesize+1<<"<br>\n";	
		tFile<<"Num of Input Seqs with unknown_character : "<<ml.sq.num_of_seq_of_all_input-nodesize+1<<"\n";
	}

	if(ml.ct.isUseWildCard){
		vFile<<"Max sequence length : "<<ml.sq.max_seq_length<<"<br>\n";
		vFile<<"Min sequence length : "<<ml.sq.min_seq_length<<"<br>\n";
		tFile<<"Max sequence length : "<<ml.sq.max_seq_length<<"\n";
		tFile<<"Min sequence length : "<<ml.sq.min_seq_length<<"\n";
	}else{
		vFile<<"Sequence length : "<<ml.sq.max_seq_length<<"<br>\n";	
		tFile<<"Sequence length : "<<ml.sq.max_seq_length<<"\n";	
	}

	if(fname!=NULL){
		vFile<<"Trees File : "<<fname<<"<br>\n";
		tFile<<"Trees File : "<<fname<<"\n";
	}
	if(pairFname!=NULL){
		vFile<<"Pairs File : "<<pairFname<<"<br>\n";
		tFile<<"Pairs File : "<<pairFname<<"\n";
	}
	if(degreeFname!=NULL){
		vFile<<"Degrees File : "<<degreeFname<<"<br>\n";
		tFile<<"Degrees File : "<<degreeFname<<"\n";
	}
	vFile<<"<br>\n";
	vFile<<"Total number of Similar Pairs : "<<ml.num_of_similar_pairs<<"<br>\n";
	vFile<<"Total number of Trees : "<<cls<<"<br>\n";
	vFile<<"Total number of Seqs included in Minimum Spanning Trees : "<<member<<"<br>\n";
	vFile<<"Total number of Singletons : "<<single<<"<br>\n";
	vFile<<"Size of the largest tree : "<<largest_tree<<"<br>\n";

	tFile<<"Total number of Similar Pairs : "<<ml.num_of_similar_pairs<<"\n";
	tFile<<"Total number of Trees : "<<cls<<"\n";
	tFile<<"Total number of Seqs included in Minimum Spanning Trees : "<<member<<"\n";
	tFile<<"Total number of Singletons : "<<single<<"\n";
	tFile<<"Size of the largest tree : "<<largest_tree<<"\n";

	if(cls==0){
		vFile<<"Average tree size : "<<0<<"<br>\n";
		tFile<<"Average tree size : "<<0<<"\n";
	}else{
		vFile<<"Average tree size : "<<(double)(member)/cls<<"<br>\n";		
		tFile<<"Average tree size : "<<(double)(member)/cls<<"\n";		
	}

//	vFile<<"<div class="button" onclick=window.print()>print</div>";

	vFile<<"<div class=\"button\" align=\"center\" style=\"TEXT-ALIGN:center;WIDTH:120px;FLOAT:left;MARGIN-LEFT:320px;CURSOR:pointer\"";
	vFile<<"onMouseOut=\"this.style.color='#363636';\" onMouseOver=\"this.style.color='#FF1493';\"";
	vFile<<"onClick=\"saveSummaryText();\">";
	vFile<<"<div style=\"MARGIN-TOP:3px\">Save as Text</div></div><BR>";

	vFile<<"<div class=\"button\" align=\"center\" style=\"TEXT-ALIGN:center;WIDTH:120px;FLOAT:left;MARGIN-LEFT:320px;MARGIN-TOP:5px;CURSOR:pointer\"";
	vFile<<"onMouseOut=\"this.style.color='#363636';\" onMouseOver=\"this.style.color='#FF1493';\"";
	vFile<<"onClick=\"window.print()\">";
	vFile<<"<div style=\"MARGIN-TOP:3px\">Print</div></div>";

	vFile<<"<br><br><br><br>";

	vFile<<"<div id=\"degree\" style=\"width: 450px; height: 300px; \"></div>\n";
	vFile<<"<br><br>";
	vFile<<"<div id=\"tree\" style=\"width: 450px; height: 300px; \"></div>\n";
	vFile<<"<br><br>";
	vFile<<"<div id=\"composition\" style=\"width: 450px; height:";
	if(top_N<cls){
		vFile.precision(2);
		vFile<<200+20*top_N;
	}else{
		vFile.precision(2);
		vFile<<200+20*cls;
	}
	vFile<<"px; \"></div>\n";
	vFile<<"</div>";
//	if(tfilename.find_last_of("/")>=0){
//		vFile<<"<div id=\"job_id\" >"<<tfilename.substr(tfilename.find_last_of("/")+1)<<"</div>\n";
//	}
	vFile<<"<div id=\"job_id\" style=\"visibility:hidden\">"<<tfilename<<"</div>\n";
	vFile<<"</body>";
	vFile<<"</html>";
	vFile.close();

	tFile<<"\ndegree distribution\n";	
	tFile<<"0"<<","<<hist[0]<<"\n";
	if(hist[1]>0)
		tFile<<"1"<<","<<hist[1]<<"\n";
	if(hist[2]>0)
		tFile<<"2-9"<<","<<hist[2]<<"\n";
	for(int cnt=3;cnt<10;cnt++){
		if(hist[cnt]>0){
			tFile.precision(2);
			tFile<<(int)pow(10,cnt-2)<<" - "<<(int)pow(10,cnt-1)-1<<","<<hist[cnt]<<"\n";
		}
	}

	tFile<<"\ntree size distribution\n";
	if(t_hist[0]>0){
		tFile<<2<<" - "<<9<<","<<t_hist[0]<<"\n";
	}
	for(int cnt=1;cnt<10;cnt++){
		if(t_hist[cnt]>0){
			tFile<<(int)pow(10,cnt)<<" - "<<(int)pow(10,cnt+1)-1<<","<<t_hist[cnt]<<"\n";
		}
	}
	tFile<<"\nComposition of large tree (%)\n";
	tFile<<"tree_id,tree_size,";
	for(int cnt=0;cnt<ml.ct.num_of_character;cnt++){
		tFile<<ml.ct.toChar[cnt]<<",";
	}
	tFile<<"\n";
	for(int i=0;i<top_N;i++){
		if(toptop_N[i]!=-1){
			tFile<<toptop_N[i]<<","<<forest[toptop_N[i]].member<<",";
		}
		for(int cnt=0;cnt<ml.ct.num_of_character;cnt++){
			if(toptop_N[i]!=-1){
				tFile.precision(3);
				tFile<<(double)toptop_Nchar[i*ml.ct.num_of_character + cnt]*100/toptop_N_z[i]<<",";
			}
		}
		tFile<<"\n";
	}

	tFile.close();
	cerr<<"outputdone.";
	return(0);
}


int MSTree::linkChildrenAsForest(){

	//for root
	forest[nodesize-1].parent=nodesize-1;
	forest[nodesize-1].member=1;
	forest[nodesize-1].num_of_sod=0;
	forest[nodesize-1].sod_ptr=0;

	TYPE_SIGNED_INDEX n0m=0;
	TYPE_SIGNED_INDEX min_cls=MIN_CLS; // 最小クラスタサイズ
	for(TYPE_SIGNED_INDEX i=0;i<nodesize -1;i++){
		TYPE_SIGNED_INDEX p = forest[i].parent;
		if(p==i && forest[i].member>min_cls){
			n0m++;
		}
	}

//	forest[0].sod = (int *)malloc(sizeof(int)*nodesize-1);
	forest[nodesize-1].sod = (TYPE_SIGNED_INDEX *)malloc(sizeof(TYPE_SIGNED_INDEX)*(forest[nodesize-1].member+n0m));
	for(TYPE_SIGNED_INDEX i=0;i<nodesize-1;i++){
		if(forest[i].num_of_sod!=0){
			forest[i].sod = (TYPE_SIGNED_INDEX *)malloc(sizeof(TYPE_SIGNED_INDEX)*forest[i].num_of_sod);
		}
	}

	max_mst_in_fst=0; //max cluster のgid
	TYPE_SIGNED_INDEX mm = min_cls;
	for(TYPE_SIGNED_INDEX i=0;i<nodesize-1;i++){
		TYPE_SIGNED_INDEX p = forest[i].parent;
		if(p==i){
			if(forest[i].member>min_cls){
				if(forest[i].member>mm){
					max_mst_in_fst = forest[nodesize-1].sod_ptr;
					mm = forest[i].member;
				}
				//cerr<<forest[nodesize-1].sod_ptr<<"  "<<forest[i].member<<"\n";
				forest[i].parent=nodesize-1;
				forest[i].edge = dist_thr+1;
				forest[nodesize-1].member += forest[i].member;
				forest[nodesize-1].num_of_sod++;
				forest[nodesize-1].sod[ forest[nodesize-1].sod_ptr++ ] = i;
			}
		}else{
			forest[p].sod[ forest[p].sod_ptr++ ]=i;	
		}
	}
	forest[nodesize-1].parent=nodesize-1;
	max_cls=nodesize-1;

	return(0);
}

int MSTree::linkChildren(){
	for(TYPE_SIGNED_INDEX i=0;i<nodesize;i++){
		if(forest[i].num_of_sod!=0){
			forest[i].sod = (TYPE_SIGNED_INDEX *)malloc(sizeof(TYPE_SIGNED_INDEX)*forest[i].num_of_sod);
		}
	}
	max_cls = 0;
	int mm = 0;
	for(int i=0;i<nodesize;i++){
		int p = forest[i].parent;
		if(p==i){
			if(forest[i].member>mm){
				max_cls = i;
				mm = forest[i].member;
			}
		}else{
			forest[p].sod[ forest[p].sod_ptr++ ]=i;	
		}
	}
	return(0);
}

int MSTree::outputLargeTreeForSIF(char *ofname){
	vector<int> tc;
	int min_cls=MIN_CLS;
	member=0; cls=0; single=0;
	largest_tree=0;


	FILE *outFile;
	if (isSTDOUT) outFile = stdout;
	else outFile = fopen(ofname, "w");

	for(int i=0;i<nodesize;i++){
		int p = forest[i].parent;
		if(p==i){
			if(forest[i].member==1){
				single++;
				if (isTreeOut&&isSingleOut)
					fprintf(outFile, "%s\n", forest[i].id);
			}else if(forest[i].member>min_cls){
				if(i == max_cls){
					member = forest[i].member-1; // for linkchildasforest
					for(int cnt=0;cnt<forest[max_cls].num_of_sod;cnt++){
						tc.push_back(forest[max_cls].sod[cnt]);
						if( largest_tree < forest[ forest[max_cls].sod[cnt] ].member){
							largest_tree = forest[ forest[max_cls].sod[cnt] ].member;
						}
						cls++;
					}
					int rn = tc.size();
					while(tc.size()>0){
						int pos = tc[tc.size()-1];
						if(forest[pos].num_of_sod==0){
							tc.pop_back();
						}else{
							if(tc.size()-1<rn){
								rn--;
							}
							tc.pop_back();
							for(int cnt=0;cnt<forest[pos].num_of_sod;cnt++){
								tc.push_back(forest[pos].sod[cnt]);
								fprintf(outFile,"%s\t-\t%s\t%0.1f\n",forest[pos].id, forest[ forest[pos].sod[cnt] ].id, forest[ forest[pos].sod[cnt] ].edge);
							}
						}
					}

				}
			}
		}
	}
	if(member == -1){
		member = 0;
	}
	cerr<<"singleton="<<single<<" num of trees="<<cls<<" mem in MSF="<<member<<"\n";
//	cerr<<"size of largest tree="<<largest_tree<<" "<<"average tree size="<<(member-1)/cls<<"\n";
	fclose(outFile);

	return(0);
}

int MSTree::outputLargeTreeForText(char *ofname){
	vector<int> tc;
	 member=0; cls=0; single=0;
	int min_cls=MIN_CLS;
	largest_tree=0;
	FILE *outFile;
	if (isTreeOut) {		
		if (isSTDOUT) outFile = stdout;
		else outFile = fopen(ofname, "w");
	}
	for(int i=0;i<nodesize;i++){
		int p = forest[i].parent;
		if(p==i){
			if(forest[i].member==1){
				single++;
				if (isTreeOut&&isSingleOut)
					fprintf(outFile, "singleton: %s\n", forest[i].id);
			}else if(forest[i].member>min_cls){
				if(i == max_cls){
					member = forest[i].member-1; // for linkchildasforest
					for(int cnt=0;cnt<forest[max_cls].num_of_sod;cnt++){
						tc.push_back(forest[max_cls].sod[cnt]);
						if( largest_tree < forest[ forest[max_cls].sod[cnt] ].member){
							largest_tree = forest[ forest[max_cls].sod[cnt] ].member;
						}
						cls++;
					}
					int rn = tc.size();
					while(tc.size()>0){
						int pos = tc[tc.size()-1];
						if(forest[pos].num_of_sod==0){
							tc.pop_back();
						}else{
							if(tc.size()-1<rn){
								if(isTreeOut){
									fprintf(outFile,"id:%d ",tc[tc.size()-1]);
									fprintf(outFile,"--\n");
								}
								rn--;
							}
							if(isTreeOut)
								fprintf(outFile,"%s: ",forest[pos].id);
							tc.pop_back();
							for(int cnt=0;cnt<forest[pos].num_of_sod;cnt++){
								tc.push_back(forest[pos].sod[cnt]);
								if(isTreeOut)
									fprintf(outFile,"%s(%0.1f),",forest[ forest[pos].sod[cnt] ].id,forest[ forest[pos].sod[cnt] ].edge);
							}
							if(isTreeOut)
								fprintf(outFile,"\n");						
						}
					}

				}
			}
		}
	}
	if(member == -1){
		member = 0;
	}
	cerr<<"singleton="<<single<<" num of trees="<<cls<<" mem in MSF="<<member<<"\n";
//	cerr<<"size of largest tree="<<largest_tree<<" "<<"average tree size="<<(member-1)/cls<<"\n";
	if(isTreeOut)
		fclose(outFile);
	return(0);
}

int MSTree::outputLargeTreeForWalrus(char *ofname){
	member=0; cls=0; single=0;
	vector<int> tc;
	vector<int> td;
	vector<int> ed;
	vector<int> e_color;
	string walrus_ct;
	char tmp;

	FILE *outFile;
	if (isSTDOUT) outFile = stdout;
	else outFile = fopen(ofname, "w");

	int ect=0;
	int r=255,g=255,b=0;

	/**
	//	for MST
	for(int cnt=0;cnt<dist_thr+1;cnt++){
		g = cnt*255/(dist_thr+1);
		b = g;
		r = 255-g;
		ect = r*256*256 + g*256 + b;
		e_color.push_back(ect);
	}
	e_color.push_back(0);
	// for MST
**/

// for MSF
	e_color.push_back(0);
	cerr<<"MAX CLS ID"<<max_mst_in_fst<<"\n";
	for(int cnt=0;cnt<forest[nodesize-1].num_of_sod;cnt++){
		g = cnt*255/forest[nodesize-1].num_of_sod;
		b = g;
		r = 255-g;
		ect = r*256*256 + g*256 + b;
		e_color.push_back(ect);
	}
	r=255;
	g=255;
	b=255;
	ect = r*256*256 + g*256 + b;
	e_color[forest[nodesize-1].num_of_sod - max_mst_in_fst] = ect;

// for MSF

	int min_cls=MIN_CLS;
	for(int i=0;i<nodesize;i++){
		int p = forest[i].parent;
		if(p==i){
			if(forest[i].member==1){
				single++;
			}else if(forest[i].member>min_cls){
				cls++;
				member += forest[i].member;
				if(i == max_cls){
					fprintf(outFile,"Graph \n {\n \"MST_BY_SLIDESORT\" ;\n;\n  ");
					fprintf(outFile,"%d; %d; ",forest[i].member, forest[i].member-1);
					fprintf(outFile,"0; 0;\n\n@links=[\n");
					int t=0;
					int nn=0;
					int gid=0;
					while(1){
						for(int j=0;j<forest[p].num_of_sod;j++){
							nn++;
							fprintf(outFile,"{%d; %d;}",t,nn);
							tc.push_back(forest[p].sod[j]);
							td.push_back(nn);
							//for MST
//							ed.push_back(forest[forest[p].sod[j]].edge);
							
							//for MSF
							//cerr<<gid<<" ";
							ed.push_back(gid);
							//for MSF

							if(ed.size()!=forest[i].member-1){
								fprintf(outFile,",\n");
							}else{
								fprintf(outFile,"\n");
							}
						}
						//for MSF
						if(t==0){
							gid++;
						}
						//for MSF
						if(tc.size()==0){
							break;
						}
						p = tc[ tc.size()-1 ];
						tc.pop_back();
						t = td[ td.size()-1 ];
						td.pop_back();
						if(t<forest[nodesize-1].num_of_sod){
							gid++;
						}
					}
					fprintf(outFile,"];\n\n");
					fprintf(outFile," @paths=;\n  @enumerations=; \n \
						@attributeDefinitions=[ \n \
					{\n \
					@name=$root; \n \
					@type=bool;\n \
					@default=|| false ||;\n \
					@nodeValues=[ { @id=0; @value=T; } ];\n \
					@linkValues=;\n \
					@pathValues=;\n \
					},\n \
					{\n \
					@name=$tree_link;\n \
					@type=bool;\n \
					@default=|| false ||;\n \
					@nodeValues=;\n \
					@linkValues=[	\n ");
					for(TYPE_SIGNED_INDEX cnt=0;cnt<forest[max_cls].member-2;cnt++){
						fprintf(outFile,"{%d; T;},\n",cnt);
					}
					fprintf(outFile,"{%d; T;}\n ];\n @pathValues=; \n},\n\n",forest[max_cls].member-2);
					fprintf(outFile,"{ \n \
						@name=$tree_link_color; \n \
						@type=int; \n \
						@default=|| 0 ||; \n \
						@nodeValues=; \n \
						@linkValues=[ \n ");
						for(TYPE_SIGNED_INDEX cnt=0;cnt<forest[max_cls].member-2;cnt++){
							fprintf(outFile,"{%d; %d;},\n",cnt,e_color[ ed[cnt] ]);
						}
						fprintf(outFile,"{%d; %d;}\n ];\n @pathValues=; \n}\n ]; \n",forest[max_cls].member-2, e_color[ ed[forest[max_cls].member-2] ]);

						fprintf(outFile," \
							@qualifiers=[ \n \
						{ \n \
						@type=$spanning_tree; \n \
						@name=$sample_spanning_tree; \n \
						@description=; \n \
						@attributes=[ \n \
						{ @attribute=0; @alias=$root; }, \n \
						{ @attribute=1; @alias=$tree_link; } \n \
						]; \n \
						} \n \
							]; \n \
							@filters=; \n \
							@selectors=; \n \
							@display=[ \n \
						{ \n \
						\"test2\"; \n \
						[ \n \
						{2;\"color\";F;T;F;} \n \
						]; \n \
						} \n \
							]; \n \
							@presentations=; \n \
							@presentationMenus=; \n \
							@displayMenus=; \n \
							@selectorMenus=; \n \
							@filterMenus=; \n \
							@attributeMenus=; \n \
						}\n ");
					}
				}
			}
		}
		cerr<<"singleton="<<single<<" others="<<cls<<" mem of the largest cls="<<forest[max_cls].member<<" mem in MSF="<<member<<"\n";
		return(0);
	}

	int MSTree::outputRoot(){
		member=0; cls=0; single=0; 
		vector<TYPE_SIGNED_INDEX> tc;
		vector<TYPE_SIGNED_INDEX> td;
		for(TYPE_SIGNED_INDEX i=0;i<nodesize;i++){
			TYPE_SIGNED_INDEX p = forest[i].parent;
			if(p==i){
				if(forest[i].member==1){
					single++;
				}else{
					if(i == max_cls){
						cerr<<"-----\n";
						cerr<<p<<" = root\n";
						TYPE_SIGNED_INDEX t=0;
						t=1;
						while(1){
							for(TYPE_SIGNED_INDEX j=0;j<forest[p].num_of_sod;j++){
								//						cerr<<forest[p].sod[j]<<","<<forest[ forest[p].sod[j] ].edge<<"=>"<<p<<"    "<<t<<"\n";
								cout<<p<<"\t"<<forest[p].sod[j]<<"\t"<<forest[forest[p].sod[j]].edge<<"\n";
								tc.push_back(forest[p].sod[j]);
								td.push_back(t+1);
							}
							if(tc.size()==0){
								break;
							}
							p = tc[ tc.size()-1 ];
							tc.pop_back();
							t = td[ td.size()-1 ];
							td.pop_back();
						}
						cerr<<"\t\t"<<i<<","<<forest[i].member<<"\n";
					}
				}
				member += forest[i].member;
				cls++;
			}else{				
			}	
		}
		cerr<<"cls = "<<cls<<"  single="<<single<<"\n";
		//		cerr<<cls<<"  "<<member<<"\n";
		//	cerr<<debug_cnt<<"  ind pairs.\n";
		return(0);
	}

	int MSTree::VOID(const char* head1, const char* head2, TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size){
		return(0);
	}

	int MSTree::mkForest(const char* head1, const char* head2,TYPE_INDEX seqid1, TYPE_INDEX seqid2, char* aln1, char* aln2, double dist,int aln_size){
		if(isPairOut){
			pFile<<forest[seqid1].id<<" , "<<forest[seqid2].id<<" , "<<dist<<"\n";
		}
		if(isAlnOut){
			pFile<<aln1<<"\n";
			pFile<<aln2<<"\n";
		}

		if(reverse_complement){
			seqid1 /= 2;
			seqid2 /= 2;
		}

		forest[seqid1].degree++;
		forest[seqid2].degree++;


#ifdef DEBUG
		cerr<<"--\n";
		cerr<<forest[seqid1].id<<" - "<<forest[seqid2].id<<"\n";
#endif
		int rank[2]={0,0};
		int rdist[2]={0,0};	
		TYPE_SIGNED_INDEX next=0, prev=0, cur=0;
		double prev_e=0, next_e=0;
		TYPE_SIGNED_INDEX mr;
		TYPE_SIGNED_INDEX far, near;
		TYPE_SIGNED_INDEX ft, nt;
		//	int A, B1, B2, L1=0, L2=0, R=0;
		TYPE_SIGNED_INDEX APL=0, B2=0, R=0, orgR=0;
		TYPE_SIGNED_INDEX ct,nct;
		double max;
		TYPE_SIGNED_INDEX max_n;
		TYPE_SIGNED_INDEX grp[2]={0,0};

		for(int i=0;i<2;i++){
			if(i==0){
				next = seqid1;
			}else{
				next = seqid2;
			}
			while(forest[next].parent!=next){
				next = forest[next].parent;
				rdist[i]++;
			}
			rank[i]=forest[next].rank; 
			grp[i]=next;
		}
		orgR = next;

		if( grp[0]==grp[1] ){
#ifdef DEBUG
			cerr<<"grp="<<orgR<<"\n";
#endif
			if(dist==dist_thr){
				return(0);
			}
#ifdef DEBUG
			cerr<<"grp=grp\n";
#endif
			//processing circuit
			max=dist;
			max_n=(int)NULL;
			if(forest[seqid1].parent == forest[seqid2].parent){
				far = seqid1;
				near = seqid2;
				ft=far; nt=near;
#ifdef DEBUG
				cerr<<"p=p\n";
#endif
			}else{
#ifdef DEBUG
				cerr<<"p != p\n";
#endif
				if(rdist[0]<rdist[1]){
					//				L1 = rdist[0];
					far = seqid2;
					near = seqid1;
				}else{
					//				L1 = rdist[1];
					far = seqid1;
					near = seqid2;
				}
				ft=far; nt=near;
				for(TYPE_SIGNED_INDEX i=0;i<abs(rdist[0]-rdist[1]);i++){
					if(max<forest[ft].edge){
						max = forest[ft].edge;
						ct = far;
						max_n = ft;
					}
					ft = forest[ft].parent;			
				}
			}
			while(ft!=nt){
				//			L1--;
				if(max<forest[ft].edge){
					max_n = ft;
					ct = far;
					max = forest[ft].edge;
				}
				if(max<forest[nt].edge){
					max_n = nt;
					ct =near;
					max = forest[nt].edge;
				}
				ft = forest[ft].parent;
				nt = forest[nt].parent;
			}
#ifdef DEBUG
			if(max_n==NULL)
				cerr<<"max null\n";
#endif
			B2=0;
			if(max_n!=(int)NULL){
				TYPE_SIGNED_INDEX mr = ft;
				next = ct;
				while(next!=max_n){
					next = forest[next].parent;	
					B2++;
				}
				if(ct==seqid1){
					nct=seqid2;
					APL = rdist[1];
				}else{
					nct=seqid1;
					APL = rdist[0];
				}

				next=ct;
				prev=nct;
				prev_e = dist;
				forest[ct].num_of_sod++;
				forest[nct].num_of_sod++;

				while(1){
					cur = next;
					next_e = forest[cur].edge;
					next = forest[cur].parent;
					forest[cur].parent = prev;
					forest[cur].edge = prev_e;
					prev = cur;
					prev_e = next_e;
#ifdef DEBUG
					cerr<<"fst loop="<<cur<<"   "<<forest[ cur ].parent<<"\n";
#endif
					if(cur==max_n){
						forest[cur].num_of_sod--;
						forest[next].num_of_sod--;
						break;
					}
				}

				if(B2>0){
					next = nct;
					for(TYPE_SIGNED_INDEX i=0;i<APL-B2;i++){
						next = forest[next].parent;	
					}
					forest[next].rank = forest[orgR].rank + B2 + 1;
					forest[next].member = forest[orgR].member;
					forest[next].num_of_sod++;
					prev = next;
					prev_e = -10;
					while(1){
						cur = next;
						next_e = forest[cur].edge;
						next = forest[cur].parent;
						forest[cur].parent = prev;
						forest[cur].edge = prev_e;
						prev = cur;
						prev_e = next_e;
#ifdef DEBUG
						cerr<<"snd loop="<<cur<<"   "<<forest[ cur ].parent<<"\n";
#endif
						if(cur == next){
							forest[cur].num_of_sod--;
							break;
						}
					}
				}else{
					forest[orgR].rank++;
#ifdef DEBUG
					cerr<<"root="<<orgR<<"\n";
#endif
				}
			}
		}else{
#ifdef DEBUG
			cerr<<"grp != grp\n";
			cerr<<"grp="<<grp[0]<<"  "<<grp[1]<<"\n";
#endif

			TYPE_SIGNED_INDEX md = (rdist[0]+rank[0]+rdist[1]+rank[1])/2;

			if(rdist[0]+rank[0]<rdist[1]+rank[1]){
				far = seqid1;
				near = seqid2;
				md -= (rdist[0]+rank[0]);
			}else{
				far = seqid2;
				near = seqid1;
				md -= (rdist[1]+rank[1]);
			}

			mr = near;
			for(TYPE_SIGNED_INDEX i=0;i<md;i++){
				if(forest[mr].parent==mr){
					break;
				}else{
					mr = forest[mr].parent;
				}
			}
			forest[mr].member = forest[grp[0]].member + forest[grp[1]].member;
			forest[mr].rank = md+1;
			forest[near].num_of_sod++;

			//rootまでGo
			forest[mr].num_of_sod++;
			forest[far].num_of_sod++;
			for(int i=0;i<2;i++){
				if(i==0){
					next = mr;
					prev = mr;
					prev_e = -10;
				}else{
					next = far;
					prev = near;
					prev_e = dist;
				}
				while(1){
					cur = next;
					next = forest[cur].parent;
					next_e = forest[cur].edge;
					forest[cur].parent = prev;
					forest[cur].edge = prev_e;
					prev = cur;
					prev_e = next_e;
					if(cur == next){
						rank[i] = forest[cur].rank;
						forest[cur].num_of_sod--;
						break;
					}
				}
			}
		}

#ifdef DEBUG
		cerr<<"id id pr rk grp\n";
		for(TYPE_SIGNED_INDEX i=0;i<nodesize;i++){
			cerr<<i<<"  ";
			cerr<<forest[i].id<<"  ";
			cerr<<forest[i].parent<<"  ";
			cerr<<forest[i].rank<<"\n";
		}
		cerr<<"--\n";
#endif	
		return(0);
	}
