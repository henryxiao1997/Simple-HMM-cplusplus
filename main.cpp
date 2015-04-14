#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <sstream>

#include "Hmm.h"

using namespace std ;

//some constants
const int MaxTagNum = 44 ;

struct WordTagPair
{
	int iWordID ;
	int iTagID ;
} ;

string Tags[] = { "#",
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
} ;

int main()
{
	cout<<"Hello world for HmmOne!"<<endl ;
	cout<<"Test for viterbi algorithm!"<<endl ;

	Hmm theHmm ;/////////////////////////////
	//contain the whole word list
	vector<string> WordList ;
	//initialize it from file of Lexicon.txt
	ifstream in ;
	in.open( "Lexicon.txt" ) ;
	if( !in.is_open() )
	{
		cerr<<"Can not open the file of Lexicon.txt to load the word list"<<endl ;
		exit( EXIT_FAILURE ) ;
	}
	string sTemp ;
	while( in.good() )
	{
		in>>sTemp ;
		WordList.push_back( sTemp ) ;
	}
	in.close() ;
/*
	//load the training file to memory struct
	vector<WordTagPair> WordTagPairArray ;
	in.clear() ;
	in.open( "train01.txt" ) ;
	//in.clear() ;
	if( !in.is_open() )
	{
		cerr<<"Can not open the file of train01.txt to load the word-tag pair list"<<endl ;
		exit( EXIT_FAILURE ) ;
	}
	string word ;
	string tag ;
	WordTagPair tempWordTagPair ;
	int WordTagPairNum = 0 ;
	while( in.good() )
	{
		getline( in, sTemp,'\n' ) ;
		istringstream inStream( sTemp ) ;
		inStream>>word>>tag ;
		//cout<<"the word is : "<<word<<endl ;
		//cout<<"the tag is : "<<tag<<endl ;
		vector<string>::iterator Iter = lower_bound( WordList.begin(), WordList.end(), word ) ;
		tempWordTagPair.iWordID = Iter - WordList.begin() ;
		int iTag = -1 ;
		string* pBegin = Tags ;
		string* pEnd = Tags + MaxTagNum ;
		string* pResult = NULL ;
		pResult = lower_bound( pBegin, pEnd, tag ) ;
		if( pResult )
			iTag = pResult - pBegin ;
		//////////////////////////// the binary search using c programme///////////////////////////
		int iLow = 0 ; 
		int iHigh = MaxTagNum ;
		while( iLow<=iHigh )
		{
			int iMid = ( iLow + iHigh ) / 2 ;
			if( Tags[iMid] == tag )
			{
				iTag = iMid ;
				break ;
			}
			else if( Tags[iMid] < tag )
				iLow = iMid + 1 ;
			else
				iHigh = iMid - 1 ;
		}
		///////////////////////////////////////////////////////////////////////////////////////////
		tempWordTagPair.iTagID = iTag ;
		WordTagPairArray.push_back( tempWordTagPair ) ;
		WordTagPairNum++ ;
		//cout<<"the number is : "<<WordTagPairNum <<endl ;
	}
	in.close() ;

	//initialize some variables
	int N = MaxTagNum ;
	int M = WordList.size() ;
	int iSeed = 1 ;
	Hmm theHmm( N, M, iSeed ) ;

	//statistic the frequency
	int* TagFreq = new int[N+1] ;
	int** TagTagFreq = iMatrix( 0, N, 0, N ) ;
	int** TagWordFreq ;
	TagWordFreq = iMatrix( 0, N, 0 , M ) ;
	//firstly initialize them
	for( int i=0 ; i<N+1 ; i++ )
		for( int j=0 ; j<N+1 ; j++ )
			TagTagFreq[i][j] = 0 ;
	for( i=0 ; i<N+1 ; i++ )
		TagFreq[i] = 0 ;
	for( i=0 ; i<N+1 ; i++ )
		for( int j=0 ; j<M+1 ; j++ )
			TagWordFreq[i][j] = 0 ;
	//secondly statistic
	double tagFreqSum = 0.0 ;
	for( i=0 ; i<(int)WordTagPairArray.size()-1 ; i++ )
	{
		TagFreq[WordTagPairArray[i].iTagID]++;
		TagTagFreq[WordTagPairArray[i].iTagID][WordTagPairArray[i+1].iTagID]++;
		TagWordFreq[WordTagPairArray[i].iTagID][WordTagPairArray[i].iWordID]++;
	}
	for( i=0 ; i<N ; i++ )
		tagFreqSum += TagFreq[i] ;

	//set the value of the matrix
	for (i = 0; i < theHmm.N; i++) 
	{ 
		for ( int j = 0; j < theHmm.N; j++) 
		{
			theHmm.A[i][j]=(double)TagTagFreq[i][j]/TagFreq[i]; 
		}
	}
	for ( int j = 0; j < theHmm.N ; j++) 
	{ 
		for ( int k = 0; k < theHmm.M ; k++) 
		{
			theHmm.B[j][k]=(double)TagWordFreq[j][k]/TagFreq[j];
		}
	}
	for ( i = 0; i < theHmm.N; i++) 
		theHmm.pi[i]=(double)TagFreq[i]/tagFreqSum;
*/
	//string filename( "HmmData.txt" ) ;
	theHmm.ReadHmm( "HmmData.txt" ) ;
	//theHmm.ReadHmm( filename) ;
	//test viterbi
	int T=3;
	char test_string[] ="I love you";
	vector<int> O ;
	vector<int> S ;
	char* p=test_string;
int	i=1;/////////////////
	char temp[256] ;
	while (i<=T){
		int j=0;/////////////////
		while ( (*p!=' ')&&(*p)){
			temp[j++] = *p;
			p++;
		}
		temp[j]=0;
		O.push_back( (find( WordList.begin(), WordList.end(), temp )-WordList.begin() ) ) ;
		i++;
		while ((*p==' ')&&i<=T) p++;
	};

	////////test generate sequence//////////////
	int t = 15 ;
	int seed = 108 ;
	vector<int> theSymbol ;
	vector<int> theState ;
	theHmm.GenerateSequence( seed, t, theSymbol, theState ) ;
	for( i=0 ; i<t ; i++ )
	{
		cout<<"the "<<(i+1)<<" word is : "<<WordList[theSymbol[i]]<<" according to state "<<theState[i]<<endl ;
	}
	////////////////end/////////////////////////

	/////////test baum-welch here//////////////
	double probInit = 0.0 ;
	double probFinal = 0.0 ;
	theHmm.BaumWelch( T, O, probInit, probFinal ) ;
	cout<<"the initial probability is : "<<probInit<<endl ;
	cout<<"the final probability is : "<<probFinal<<endl ;
	///////////////end ////////////////////////
	
	double dProb ;
	dProb = theHmm.Viterbi( T, O, S ) ;
	vector<int>::iterator SBegin = S.begin() ;
	vector<int>::iterator SEnd = S.end() ;
	while( SBegin != SEnd )
	{
		cout<<Tags[(*SBegin)]<<" " ;
		SBegin++ ;
	}
	cout<<endl ;

	double probForward = theHmm.Forward( T, O ) ;
	cout<<"the probability of forward algorithm is : "<<probForward<<endl ;
	double probBackward = theHmm.Backward( T, O ) ;
	cout<<"the probability of backward algorithm is : "<<probBackward<<endl ;

	//string filename( "HmmData.txt" ) ;
	//theHmm.WriteHmm( filename ) ;

	return 0 ;
}