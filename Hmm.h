#pragma once

// Author : Jinghui Xiao   xiaojinghui@insun.hit.edu.cn, xceman@hotmail.com
// Organization : the natual language processing center of computer department of HIT
// Data : from 5 September 2002 to 14 september 2002
// File : Hmm.h and Hmm.cpp
// Develop evirnment : windows xp, vs.net and vc++6.0
// Purpose : 1. to learn hmm
//			 2. to develop a general hmm class to be shared on the net and updated by others
//			 3. to provide a hmm source code compiling on windows system, 
//              and provide it to the beginner of learning hmm who can not find these source on the net, just like me
//           4. to learn and practise stl
// Note1 : any one can use this class and update it, but please write down why and when and by who do this, 
//        and what exactly has been done. that is convenint to others who read this code
// Note2 : the code of Tapas Kanungo writen by c is the reference, but there are some mistakes in it
//         one is the backward procedure and the normalized backward procedure is even not finished yet! 
// Note3 : there maybe mistake in the normalized backward procedure
// Note4 : if cause any damage by using this, i'm not responsibility at all!

#include <vector>

using namespace std ;

//the default false value
const double FALSEVALUE = -1.0 ;
//in baum-welch algorithm, the stop-condition is set using ratio( that is (probFinal-probInit)/probInit )
const double RATIOLIMIT = 0.001 ;
//in baum-welch algorithm, when the circle of do...while exceed certain number, break it and return the result
const int ROUNDLIMIT = 5 ;
//in baum-welch algorithm, in case the denominator is zero, so set it to a very very small double number
const double DENOMINATORLIMIT = 10e-70 ;

class Hmm
{
	//some variables in the class
//private:
public:
	//******** these are variables in Hmm itself **********************//
	//the number of states : Q={1,2,...,N}
	int N ;
	//the number of observation symbols : V={1,2,...,M}
	int M ;
	//the probability transition matrix : A[1..N][1..N] 
	//A[i][j] is the probability from state i to j
	double** A ;
	//the probablity emition matrix : B[1..N][1..M]
	//B[j][k] is the probabiltiy of the observation symbol k in the state j
	double** B ;
	//the initial state distribution
	//pi[i] is the initial state probability
	vector<double> pi ;
private:
	//********** these are variables used(created) in Hmm ****************//
	//the T value of the last time
	int iPreT ;
	//the alpha variable used in forward procedure and Baum-Welch procedure
	//alpha[t][i] is the probability of the observation sequence O1,O2....Ot and reach state i given modul u
	double** alpha ;
	//the beta variable used in forward procedure and Baum-Welch procedure
	//beta[t][i] is the probability of the observation sequence Ot,...OT given modul u and state i
	double** beta ;
	//the delta variable used in viterbi algorithm
	//delat[t][j] = max P(X1...Xt,O1...Ot,Xt=j | u )
	double** delta ;
	//the gamma variable used in baum-welch algorithm
	//gamma[t][i] = P( Xt=i | O,u )
	double** gamma ;

/********* move these out of the class**********************
	//some service programmes needed by others
private:
	//some allocating memory function to 2-dimonsion array
	void nrerror( string errStr ) ;
	int** iMatrix( int irLow, int irHigh, int icLow, int icHigh ) ;
	float** fMatrix( int irLow, int irHigh, int icLow, int icHigh ) ;
	double** dMatrix( int irLow, int irHigh, int icLow, int icHigh ) ;
	//the functions to free memory
	void iFreeMatrix( int** iMatrix, int irLow, int irHigh, int icLow, int icHigh ) ;
	void fFreeMatrix( float** fMatrix, int irLow, int irHigh, int icLow, int icHigh ) ;
	void dFreeMatrix( double** dMatrix, int irLow, int irHigh, int icLow, int icHigh ) ;
*********************************************************/

private:
	//the initialize function
	//random set the  matrixes to the value between 0 and 1
	void InitHmm( int iSeed ) ;

	//some public input-output functions
public:
	//output Hmm to the file,espcially for c++
	void WriteHmm( string strFileName ) ;
	//output Hmm to the file,espcially for c
	void WriteHmm( char* sFileName ) ;
	//read a Hmm from the file named "strFileName", for c++
	void ReadHmm( string strFileName ) ;
	//read a Hmm from the file named "strFileName", for c
	void ReadHmm( char* sFileName ) ;
	//read the observation sequence to the vector sequence, espcially for c++
	void ReadSequence( ifstream& in, vector<int>& sequence ) ;
	//read the observation sequence to the array sequence, espcially for c
	void ReadSequence( FILE* pFile, int*& sequence, int& iNum ) ;

	//generate sequence function, it's not much use in fact, but funny
public:
	//generate the sequence according Hmm for c++
	//parameters : iSeed : the seed to generate random number
	//             O : the returned observation sequence
	//             S : the returned state sequence
	void GenerateSequence( int iSeed, int T, vector<int>& O, vector<int>& S ) ;
	//generate the sequence according Hmm for c
	//parameters : iSeed : the seed to generate random number
	//             O : the returned observation sequence
	//             S : the returned state sequence
	void GenerateSequence( int iSeed, int T, int* O, int* S ) ;

	//the functions serve GenerateSequence( int iSeed, vector<int>& O, vector<int>& S )
private:
	//get the initial state using iSeed
	int GetInitialState( int iSeed ) ;
	//generate an output symbol according to the given state
	int GetSymbol( int iState ) ;
	//get the next state according to the transittion probability and given state
	int GetNextState( int iState ) ;

	//the core functions of Hmm
	//such as Forward and Backward functions and viterbi function and Baum-Welch function
public:
	//the forward function for c++ interface
	double Forward( int T, vector<int>& O ) ;
	//the forward function for c interface
	double Forward( int T, int* O ) ;
	//the forward function with scale ( been normalized ) and return the log value of the probability
	//the reason using a log value, i think, is to smooth the original value which is much small
	//this function is designed especially for c++ interface
	double ForwardNormalized( int T, vector<int>& O) ;
	//the forward function with scale ( been normalized ) and return the log value of the probability
	//the reason using a log value, i think, is to smooth the original value which is much small
	//this function is designed especially for c interface
	double ForwardNormalized( int T, int* O) ;

	//the backward function for c++ interface
	double Backward( int T, vector<int>& O ) ;
	//the backward function for c interface
	double Backward( int T, int* O ) ;
	//the backward function with scale to normalize and return the log value of the probability
	//the reason using a log value, is the same with forwardnormalized procedure
	//it's very similar to the forwardnormalized procedure
	//this function is designed especially for c++ interface
	double BackwardNormalized( int T, vector<int>& O ) ;
	//the backward function with scale to normalize and return the log value of the probability
	//the reason using a log value, is the same with forwardnormalized procedure
	//it's very similar to the forwardnormalized procedure
	//this function is designed especially for c interface
	double BackwardNormalized( int T, int* O ) ;

	//the viterbi algorithm for c++ interface
	//parameter: T : the number of the observation sequence
	//			 O : the observation sequence
	//			 S : the most likely state sequence as return value
	//return value : the probability of the state sequence
	double Viterbi( int T, vector<int>& O, vector<int>& S ) ;
	//the viterbi algorithm for c interface
	//parameter: T : the number of the observation sequence
	//			 O : the observation sequence
	//			 S : the most likely state sequence as return value
	//return value : the probability of the state sequence
	double Viterbi( int T, int* O, int*& S ) ;

	//the baum-welch algorithm, for c++ interface
	//parameter: T: the number of the observation sequence
	//           O: the observation sequence
	//           probInit : the probability of the observation sequence calculated by the origin model
	//           probFinal : the probability of the observation sequence calculated by the adjusted model
	//note: this algorithm is often not used and it' also the most complex one in hmm
	void BaumWelch( int T, vector<int>& O, double& probInit, double& probFinal ) ;
	//the baum-welch algorithm, for c interface
	//parameter: T: the number of the observation sequence
	//           O: the observation sequence
	//           probInit : the probability of the observation sequence calculated by the origin model
	//           probFinal : the probability of the observation sequence calculated by the adjusted model
	//note: this algorithm is often not used and it' also the most complex one in hmm
	void BaumWelch( int T, int* O, double& probInit, double& probFinal ) ;

	//the function called by baum-welch algorithm
private:
	//allocate gamma memory and calculate it
	double** ComputeGamma( int T, int N, double probTemp ) ;

	//allocate p memory and calculate it, for c++ interface
	//p[t][i][j] = p( Xt=i,Xt+1=j | O,u ), is the probability of transitting from state i to state j
	double*** ComputeP( double*** p,int T, int N, vector<int>& O, double prob ) ;
	//allocate p memory and calculate it, for c interface
	//p[t][i][j] = p( Xt=i,Xt+1=j | O,u ), is the probability of transitting from state i to state j
	double*** ComputeP( double*** p,int T, int N, int* O, double prob ) ;

	//recalculate the probability of hmm parameters( pi, A and B ), for c++ interface
	void RecomputeParameter( double***& p, int T, vector<int>& O, double probInit, double& probFinal ) ;
	//recalculate the probability of hmm parameters( pi, A and B ), for c++ interface
	void RecomputeParameter( double***& p, int T, int* O, double probInit, double& probFinal ) ;

	//free the memory of 3-diminson matrix
	void PFreeMatrix( double*** p, int T ) ;

public:
	Hmm(void);
	~Hmm(void);

public:
	//There are two methods to initialize a Hmm:
	//One way is to read a Hmm from a file which store it
	//  ( of course, you can state a Hmm object and explicitly call the function of Hmm::ReadHmm(char* sFileName) )
	//The other is to set the size of state set and output alphabet set and initialize it with random value
	Hmm( string strFileName ) ; // for c++ interface
	Hmm( char* sFileName ) ; // for c interface
	Hmm( int NHmm, int MHmm, int iSeed ) ;

	//and there is the third one : to initialize Hmm with every variable in it
	Hmm( int NHmm, int MHmm, double** AHmm, double** BHmm, vector<double>& piHmm ) ;
};

//some allocating memory function to 2-dimonsion array
void nrerror( string errStr ) ;
int** iMatrix( int irLow, int irHigh, int icLow, int icHigh ) ;
float** fMatrix( int irLow, int irHigh, int icLow, int icHigh ) ;
double** dMatrix( int irLow, int irHigh, int icLow, int icHigh ) ;
//the functions to free memory
void iFreeMatrix( int** iMatrix, int irLow, int irHigh, int icLow, int icHigh ) ;
void fFreeMatrix( float** fMatrix, int irLow, int irHigh, int icLow, int icHigh ) ;
void dFreeMatrix( double** dMatrix, int irLow, int irHigh, int icLow, int icHigh ) ;
