#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Hmm.h"
//#include "nrutil.h"

#define	MAX_WORD_LEN 20
#define MAX_TOTAL_WORD_NUM 300000
#define MAX_WORD_NUM 50000
#define MAX_TAG_NUM 44
#define NOT_FOUND -1
#define END_TAG MAX_TAG_NUM+1

typedef  struct  Word_Tag_Pair {
	int wordid;
	int	 tag;
}  Word_Tag_Pair ;

Word_Tag_Pair wordstream[MAX_TOTAL_WORD_NUM];

typedef struct HeadWordItem {
//	int						wordid;
	char		   word_string[MAX_WORD_LEN];//the Word String
} HeadWordItem;

typedef struct Lexicon{
	int word_number_of_lexicon;
	HeadWordItem lexicon[MAX_WORD_NUM];
} Lexicon;

char * tags[]={"#",
	"$",
	"''",
	"(",
	")",
	",",
	".",
	":",
	"CC",
	"CD",
	"DT",
	"EX",
	"FW",
	"IN",
	"JJ",
	"JJR",
	"JJS",
	"MD",
	"NN",
	"NNP",
	"NNPS",
	"NNS",
	"PDT",
	"POS",
	"PRP",
	"PRP$",
	"RB",
	"RBR",
	"RBS",
	"RP",
	"SYM",
	"TO",
	"UH",
	"VB",
	"VBD",
	"VBG",
	"VBN",
	"VBP",
	"VBZ",
	"WDT",
	"WP",
	"WP$",
	"WRB",
	"``",
	};
int bsearch_tag(char ** A,int low,int high , char* key);
int bsearch_wordid(HeadWordItem *A,int low,int high , char* key);

//#include <stdio.h> 
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//#include "hmm.h"
//#include "nrutil.h"
//#include "trainer.h"

//#include "Hmm.h"

int main (int argc, char **argv)
{

	FILE * fp;
	FILE * flexicon;
	Lexicon *lexicon; 

	char s_line[256];
	char *p;
	int word_num;
	int nlexicon_word_num;
	int tag_no;
	char temp[256];
	int i,j,k;
	int tag_freq[MAX_TAG_NUM+1];
	int tag_tag_freq[MAX_TAG_NUM+1][MAX_TAG_NUM+1];
	int ** tag_word_freq;
	int tagsum;

	//HMM  *phmm;
	Hmm  phmm;

	int   T;
	int  *O;
	int	 *q;
	double **delta;
	int **psi;
	double proba;
	
	char *test_string;
/****************************************************************/

//Initialization

//	if (argc != 3) {
//		printf("Usage error \n");
//		printf("Usage: training <lexicon file> <training text>  \n");
//		exit (1);
//	}

//	flexicon = fopen(argv[1], "r");
	flexicon = fopen("lexicon.txt", "r");
	if (flexicon == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[1]);
		exit (1);
	}

	/*
//	fp = fopen(argv[2], "r");
	fp = fopen("train01.txt", "r");
	if (fp == NULL) {
		fprintf(stderr, "Error: File %s not found\n", argv[1]);
		exit (1);
	}
	*/
	for (i=0;i<MAX_TAG_NUM+1;i++){
		tag_freq[i]=0;
	};
	for (i=0;i<MAX_TAG_NUM+1;i++)
		for (j=0;j<MAX_TAG_NUM+1;j++){
			tag_tag_freq[i][j]=0;
		}
	//tag_word_freq=(int **)imatrix(0,MAX_TAG_NUM-1,0,MAX_WORD_NUM-1);
	tag_word_freq=(int **)iMatrix(0,MAX_TAG_NUM-1,0,MAX_WORD_NUM-1);

	for (i=0;i<MAX_TAG_NUM;i++)
		for (j=0;j<MAX_WORD_NUM;j++){
			tag_word_freq[i][j]=0;
		};

/****************************************************************/
/**************************************************************/
// Step 1
//load lexicon into memory

	lexicon=(Lexicon *)malloc(sizeof(Lexicon));
	nlexicon_word_num=0;
	while( !feof(flexicon) )
	{
		fgets(s_line,256,flexicon);
		p=strchr((char *)s_line,'\n');
		if( p )
		{
			*p='\0';
		}
		strcpy((char *)lexicon->lexicon[nlexicon_word_num].word_string,s_line);
		nlexicon_word_num++;
	}
	nlexicon_word_num--;
	lexicon->word_number_of_lexicon=nlexicon_word_num;
	fclose(flexicon);
/****************************************************************/	
/****************************************************************
// Step 2
//Load the training text into memory

	word_num=0;
	while(!feof(fp)){
		fgets(s_line,256,fp);
		p=strchr((char *)s_line,'\n');
		if(p){
			*p='\0';
		};
		p=s_line;
		if (strlen(p)!=0){
			i=0;
			while (*p!=' '){
				temp[i++]=*p;
				p++;
			}
			temp[i]=0;
			wordstream[word_num].wordid=bsearch_wordid(lexicon->lexicon,0,lexicon->word_number_of_lexicon-1,temp);
			while (*p==' ') p++;
			i=0;
			while (*p!=' '){
				temp[i++]=*p;
				p++;
			}
			temp[i]=0;
			tag_no=bsearch_tag(tags,0,MAX_TAG_NUM-1,temp);
			wordstream[word_num].tag=tag_no;
			word_num++;
		}
	}
	wordstream[word_num].tag=END_TAG;
	word_num--;
	fclose(fp);
****************************************************************/
/****************************************************************
// Step 3
//statistics M,N,A,B,pi

//Step 3.1 The construction of A[i][j]
	for (i=0;i<word_num;i++){
		tag_freq[wordstream[i].tag]++;
		tag_tag_freq[wordstream[i].tag][wordstream[i+1].tag]++;
		tag_word_freq[wordstream[i].tag][wordstream[i].wordid]++;
	}
	tagsum=0;
	for (i=0;i<MAX_TAG_NUM;i++){
		tagsum+=tag_freq[i];
	}

	//phmm=(HMM  *)malloc(sizeof(HMM));
	//phmm=(Hmm  *)malloc(sizeof(Hmm));
	phmm.M=lexicon->word_number_of_lexicon;
	phmm.N=MAX_TAG_NUM;
	//phmm->A=(double **) dmatrix(1, phmm->N, 1, phmm->N);
	phmm.A=(double **) dMatrix(0, phmm.N, 0, phmm.N);
	for (i = 0; i < phmm.N; i++) { 
		for (j = 0; j < phmm.N; j++) {
			phmm.A[i][j]=(double)tag_tag_freq[i][j]/tag_freq[i]; 
		}
	}
	
	//phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
	phmm.B = (double **) dMatrix(0, phmm.N, 0, phmm.M);
	for (j = 0; j < phmm.N; j++) { 
		for (k = 0; k < phmm.M; k++) {
			phmm.B[j][k]=(double)tag_word_freq[j][k]/tag_freq[j];
		}
	}
	//phmm->pi = (double *) dvector(1, phmm->N);
	phmm.pi.resize(MAX_TAG_NUM ) ;
	for (i =0 ; i < phmm.N; i++) 
		phmm.pi[i]=(double)tag_freq[i]/tagsum;
****************************************************************/
//to write hmm to file
	//phmm.WriteHmm( "HmmData.txt" ) ;
	phmm.ReadHmm( "HmmData.txt" ) ;
	
/**********************************************************************/
	T=3;
	test_string="I love you";

	//O = ivector(1,T);
	O = (int*)malloc( T*sizeof(int) ) ;
	p=test_string;
	i=0;
	while (i<T){
		j=0;
		while (*p!=' '){
			temp[j++] = *p;
			p++;
		}
		temp[j]=0;
		O[i]=bsearch_wordid(lexicon->lexicon,0,lexicon->word_number_of_lexicon-1,temp);
		i++;
		while ((*p==' ')&&i<=T) p++;
	};
	//q = ivector(1,T);
	//q = (int*)malloc( T*sizeof(int) ) ;

	//delta = dmatrix(1, T, 1, phmm->N);
	//delta = dMatrix(0, T, 0, phmm.N);
	//psi = imatrix(1, T, 1, phmm->N);
	//psi = iMatrix(0, T, 0, phmm.N);
	//Viterbi(phmm, T, O, delta, psi, q, &proba); 
	
	proba = phmm.Viterbi( T, O, q ) ;
/*Output*/
    for (i=0; i < T; i++) 
        fprintf(stdout,"%s ", tags[q[i]]);
	printf( "\n" ) ;
	printf( "the probability from viterbi is : %e\n", proba ) ;

	double probaForward ;
	probaForward = phmm.Forward( T, O ) ;
	printf("the probability from forward algorithm is : %e\n", probaForward ) ;

	double probaForwardNormal ;
	probaForwardNormal = phmm.ForwardNormalized( T, O ) ;
	printf("the probability from forward Normalized algorithm is : %e\n", probaForwardNormal ) ;

	double probaBackward ;
	probaBackward = phmm.Backward( T, O ) ;
	printf("the probability from backward algorithm is : %e\n", probaBackward ) ;

	double probaBackwardNormal ;
	probaBackwardNormal = phmm.BackwardNormalized( T, O ) ;
	printf("the probability from backward Normalized algorithm is : %e\n", probaBackwardNormal ) ;

	printf("\n");
/*
A DT I
further JJ I
decline NN I
in IN O
prices NNS I
will MD O
lead VB O
to TO O
mine VB I
production NN I
cuts NNS I
in IN O
the DT I
U.S. NNP I
*/
 
	return 0;
}
int bsearch_tag(char ** A,int low,int high , char* key)
{
	int mid;
	char * midvalue;
	if (low >high)
		return NOT_FOUND;
	else {
		mid=(low+high)/2;
		midvalue=A[mid];
		if (strcmp((char *)key,(char *)midvalue)==0){
			return mid;
		}
		else if (strcmp((char *)key,(char *)midvalue)<0)
			return bsearch_tag(A,low,mid-1,key);
		else 
			return bsearch_tag(A,mid+1,high,key);
	}
}
int bsearch_wordid(HeadWordItem * A,int low,int high , char* key)
{
	int mid;
	char * midvalue;
	if (low >high)
		return NOT_FOUND;
	else {
		mid=(low+high)/2;
		midvalue=A[mid].word_string;
		if (strcmp((char *)key,(char *)midvalue)==0){
			return mid;
		}
		else if (strcmp((char *)key,(char *)midvalue)<0)
			return bsearch_wordid(A,low,mid-1,key);
		else 
			return bsearch_wordid(A,mid+1,high,key);
	}
}