#include "hmm.h"
//#using <mscorlib.dll>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>

using namespace std ;

/********************************the default construction and destruction functions****************************/
Hmm::Hmm(void)
{
	cout<<"This is the construction function of class Hmm "<<endl ;

	cout<<"initialize some variables in Hmm "<<endl ;
	//initialize the variables used in Hmm
	M = 0  ;
	N = 0 ;
	iPreT = 0 ;

	A = NULL ;
	B = NULL ;

	alpha = NULL ;
	beta = NULL ;
	delta = NULL ;
	gamma = NULL ;
}

Hmm::~Hmm(void)
{
	cout<<"This is the destruction function of class Hmm "<<endl ;

	//free all the memory of Hmm
	cout<<"Now free all the memory of the matrix A and B "<<endl ;
	if( A )
		dFreeMatrix( A, 0, N-1, 0, N-1 ) ;
	if( B )
		dFreeMatrix( B, 0, N-1, 0, M-1 ) ;

	cout<<"Now free all the memory of the matrix alpha and beta, delta and gamma"<<endl ;
	if( alpha )
		dFreeMatrix( alpha, 0, iPreT-1, 0, N-1 ) ;
	if( beta )
		dFreeMatrix( beta, 0, iPreT-1, 0, N-1 ) ;
	if( delta )
		dFreeMatrix( delta, 0, iPreT-1, 0, N-1 ) ;
	if( gamma )
		dFreeMatrix( gamma, 0, iPreT-1, 0, N-1 ) ;
}
/********************************the default construction and destruction functions end************************/

/***********************the adjusted construction and destruction functions****************************/
Hmm::Hmm( string strFileName )
{
	cout<<"This is the construction function of class Hmm "<<endl ;

	cout<<"initialize some variables in Hmm "<<endl ;
	//initialize the variables used in Hmm
	alpha = NULL ;
	beta = NULL ;
	delta = NULL ;
	gamma = NULL ;

	iPreT = 0 ;

	ReadHmm( strFileName.c_str() ) ;
}

Hmm::Hmm( char* sFileName)
{
	cout<<"This is the construction function of class Hmm "<<endl ;

	cout<<"initialize some variables in Hmm "<<endl ;
	//initialize the variables used in Hmm
	alpha = NULL ;
	beta = NULL ;
	delta = NULL ;
	gamma = NULL ;

	iPreT = 0 ;

	ReadHmm( sFileName ) ;
}

Hmm::Hmm( int NHmm, int MHmm, int iSeed )
{
	cout<<"This is the construction function of class Hmm "<<endl ;

	cout<<"allocate the memory needed for matrix A and B"<<endl ;
	//allocate the memory for the matrix A
	A = dMatrix( 0, NHmm-1, 0, NHmm-1 ) ;
	//allocate the memory for the matrix B
	B = dMatrix( 0, NHmm-1, 0, MHmm-1 ) ;
	//allocate the memory for the initial probability vector
	pi.resize( NHmm ) ;

	cout<<"initialize some variables in Hmm "<<endl ;
	N = NHmm ;
	M = MHmm ;
	iPreT = 0 ;

	//initialize the variables used in Hmm
	alpha = NULL ;
	beta = NULL ;
	delta = NULL ;
	gamma = NULL ;
	//initialize the Hmm with random value
	InitHmm( iSeed ) ;
}


Hmm::Hmm( int NHmm, int MHmm, double** AHmm, double** BHmm, vector<double>& piHmm )
{
	cout<<"This is the construction function of class Hmm "<<endl ;

	cout<<"initialize some variables in Hmm "<<endl ;
	N = NHmm ;
	M = MHmm ;
	iPreT = 0 ;

	A = AHmm ;
	B = BHmm ;
	pi = piHmm ;

	//initialize the variables used in Hmm
	alpha = NULL ;
	beta = NULL ;
	delta = NULL ;
	gamma = NULL ;
}

/***********************the adjusted construction and destruction functions****************************/

/***********************  Globle some service programmes needed by others******************************/
//the function output some messages and exit the programe when there is error allocating memory
//void Hmm::nrerror( string errStr )
void nrerror( string errStr )
{
	cerr<<"Numerical Recipes run-time error..."<<endl ;
	cerr<<errStr<<endl ;
	cerr<<"...now exit the programme"<<endl ;

	exit( EXIT_FAILURE ) ;
}

//some allocating memory function to 2-dimonsion array
//int** Hmm::iMatrix( int irLow, int irHigh, int icLow, int icHigh ) 
int** iMatrix( int irLow, int irHigh, int icLow, int icHigh ) 
{
	int** p ;
	p = new int* [irHigh-irLow+1] ;
	if( p == NULL )
		nrerror("allocate failure in fMatrix") ;
	p -= irLow ;

	for( int i=irLow ; i<irHigh+1 ; i++ )
	{
		p[i] = new int[icHigh-icLow+1] ;
		if( p[i] == NULL )
			nrerror("allocate failure in fMatrix") ;
		p[i] -= icLow ;
	}
	return p ;
}

//float** Hmm::fMatrix( int irLow, int irHigh, int icLow, int icHigh ) 
float** fMatrix( int irLow, int irHigh, int icLow, int icHigh ) 
{
	float** p ;
	p = new float* [irHigh-irLow+1] ;
	if( p == NULL )
		nrerror("allocate failure in fMatrix") ;
	p -= irLow ;

	for( int i=irLow ; i<irHigh+1 ; i++ )
	{
		p[i] = new float[icHigh-icLow+1] ;
		if( p[i] == NULL )
			nrerror("allocate failure in fMatrix") ;
		p[i] -= icLow ;
	}
	return p ;
}

//double** Hmm::dMatrix( int irLow, int irHigh, int icLow, int icHigh ) 
double** dMatrix( int irLow, int irHigh, int icLow, int icHigh ) 
{
	double** p ;
	p = new double* [irHigh-irLow+1] ;
	if( p == NULL )
		nrerror("allocate failure in fMatrix") ;
	p -= irLow ;

	for( int i=irLow ; i<irHigh+1 ; i++ )
	{
		p[i] = new double[icHigh-icLow+1] ;
		if( p[i] == NULL )
			nrerror("allocate failure in fMatrix") ;
		p[i] -= icLow ;
	}
	return p ;
}

//void Hmm::iFreeMatrix( int** iMatrix, int irLow, int irHigh, int icLow, int icHigh )
void iFreeMatrix( int** iMatrix, int irLow, int irHigh, int icLow, int icHigh )
{
	for( int i=irLow ; i<=irHigh ; i++ )
		delete[] (iMatrix[i]+icLow) ;

	delete[] (iMatrix+irLow) ;
}

//void Hmm::fFreeMatrix( float** fMatrix, int irLow, int irHigh, int icLow, int icHigh )
void fFreeMatrix( float** fMatrix, int irLow, int irHigh, int icLow, int icHigh )
{
	for( int i=irLow ; i<=irHigh ; i++ )
		delete[] (fMatrix[i]+icLow) ;

	delete[] (fMatrix+irLow) ;
}

//void Hmm::dFreeMatrix( double** dMatrix, int irLow, int irHigh, int icLow, int icHigh )
void dFreeMatrix( double** dMatrix, int irLow, int irHigh, int icLow, int icHigh )
{
	for( int i=irLow ; i<=irHigh ; i++ )
		delete[] (dMatrix[i]+icLow) ;

	delete[] (dMatrix+irLow) ;
}
/************************ Globle some service programmes needed by others end**************************/

/********************************some initialize function**********************************************/
void Hmm::InitHmm( int iSeed ) 
{
	cout<<"Initialize the Hmm with random values"<<endl ;

	//set the seed 
	srand( iSeed ) ;

	//for normalize
	double sum = 0.0 ;

	//set the random value to pi
	for( int i=0 ; i<N ; i++ )
	{
		pi[i] = rand()/RAND_MAX ;
		sum += pi[i] ;
	}
	//normalize
	for( i=0 ; i<N ; i++ )
		pi[i] /= sum ;

	//set the random value to the matrix A
	for( i=0 ; i<N ; i++ )
	{
		sum = 0.0 ;
		for( int j=0 ; j<N ; j++ )
		{
			A[i][j] = rand()/RAND_MAX ;
			sum += A[i][j] ;
		}
		//normalize
		for( j=0 ; j<N ; j++ )
			A[i][j] /= sum ;
	}

	//set the random value to the matrix B
	for( int j=0 ; j<N ; j++ )
	{
		sum = 0.0 ;
		for( int k=0 ; k<M ; k++ )
		{
			B[j][k] = rand()/RAND_MAX ;
			sum += B[j][k] ;
		}
		//normalize
		for( k=0 ; k<M ; k++ )
			B[j][k] /= sum ;
	}

}

/********************************some initialize function end******************************************/

/********************************the input-output function*********************************************/
void Hmm::ReadHmm( string strFileName )
{
	cout<<"Now read Hmm from file "<<strFileName<<endl ;

	ifstream in ; 
	in.open( strFileName.c_str() ) ;
	if( !in.is_open() )
	{
		cerr<<"Can not open the file of Hmm to read!"<<endl ;
		exit( EXIT_FAILURE ) ;
	}

	//read the matrix A
	cout<<"Now read the matrix A "<<endl ;
	in>>N ;
	//allocate memory first
	A = dMatrix( 0 , N, 0, N ) ;
	//initialize
	for( int i=0 ; i<N ; i++ )
		for( int j=0 ; j<N ; j++ )
			A[i][j] = 0.0 ;
	//read data from file
	for( i=0 ; i<N ; i++ )
		for( int j=0 ; j<N ; j++ )
			in>>A[i][j] ;

	//read the matrix B
	cout<<"Now read the matrix B "<<endl ;
	in>>M ;
	//allocate memory first
	B = dMatrix( 0, N, 0, M ) ;
	//initialize
	for( int j=0 ; j<N ; j++ )
		for( int k=0 ; k<M ; k++ )
			B[j][k] = 0.0 ;
	//read data from data
	for( j=0 ; j<N ; j++ )
		for( int k=0 ; k<M ; k++ )
			in>>B[j][k] ;

	//read  the initial probability vector pi
	cout<<"Read vector pi"<<endl ;
	double dTemp = 0.0 ;
	for( i=0 ; i<N ;i++ )
	{
		in>>dTemp ;
		pi.push_back( dTemp ) ;
	}
		
	in.close() ;
}

/*
void Hmm::ReadHmm( char* sFileName )
{
	cout<<"Now read Hmm from file "<<sFileName<<endl ;

	FILE* pFile ;
	pFile = fopen( sFileName, "r" ) ;
	if( pFile == NULL )
	{
		cerr<<"Can not open the file of Hmm to read"<<endl ;
		exit( EXIT_FAILURE ) ;
	}

	//read the matrix A
	cout<<"Now read the matrix A "<<endl ;
	fscanf( pFile, "%d\n", &N ) ;
	//allocate memory first
	A = dMatrix( 0 , N, 0, N ) ;
	//initialize
	for( int i=0 ; i<N ; i++ )
		for( int j=0 ; j<N ; j++ )
			A[i][j] = 0.0 ;
	//read data from file
	for( i=0 ; i<N ; i++ )
	{
		for( int j=0 ; j<N ; j++ )
			fscanf( pFile, "%f ", &A[i][j] ) ;

		fscanf( pFile, "\n" ) ;
	}

	//read the matrix B
	cout<<"Now read the matrix B "<<endl ;
	fscanf( pFile, "%d\n", &M ) ;
	//allocate memory first
	B = dMatrix( 0, N, 0, M ) ;
	//initialize
	for( int j=0 ; j<N ; j++ )
		for( int k=0 ; k<M ; k++ )
			B[j][k] = 0.0 ;
	//read data from data
	for( j=0 ; j<N ; j++ )
	{
		for( int k=0 ; k<M ; k++ )
			fscanf( pFile, "%f ", &B[j][k] ) ;

		fscanf( pFile, "\n" ) ;
	}

	//read  the initial probability vector pi
	cout<<"Read vector pi"<<endl ;
	double dTemp = 0.0 ;
	for( i=0 ; i<N ;i++ )
	{
		fscanf( pFile, "%f ", &dTemp ) ;
		pi.push_back( dTemp ) ;
	}
		
	fclose( pFile ) ;
}
*/

void Hmm::ReadHmm( char* sFileName ) 
{
	cout<<"Now read Hmm from file "<<sFileName<<endl ;

	ifstream in ; 
	in.open( sFileName ) ;
	if( !in.is_open() )
	{
		cerr<<"Can not open the file of Hmm to read!"<<endl ;
		exit( EXIT_FAILURE ) ;
	}

	//read the matrix A
	cout<<"Now read the matrix A "<<endl ;
	in>>N ;
	//allocate memory first
	A = dMatrix( 0 , N, 0, N ) ;
	//initialize
	for( int i=0 ; i<N ; i++ )
		for( int j=0 ; j<N ; j++ )
			A[i][j] = 0.0 ;
	//read data from file
	for( i=0 ; i<N ; i++ )
		for( int j=0 ; j<N ; j++ )
			in>>A[i][j] ;

	//read the matrix B
	cout<<"Now read the matrix B "<<endl ;
	in>>M ;
	//allocate memory first
	B = dMatrix( 0, N, 0, M ) ;
	//initialize
	for( int j=0 ; j<N ; j++ )
		for( int k=0 ; k<M ; k++ )
			B[j][k] = 0.0 ;
	//read data from data
	for( j=0 ; j<N ; j++ )
		for( int k=0 ; k<M ; k++ )
			in>>B[j][k] ;

	//read  the initial probability vector pi
	cout<<"Read vector pi"<<endl ;
	double dTemp = 0.0 ;
	for( i=0 ; i<N ;i++ )
	{
		in>>dTemp ;
		pi.push_back( dTemp ) ;
	}
		
	in.close() ;
}

void Hmm::WriteHmm( string strFileName )
{
	cout<<"Now write Hmm to file "<<strFileName<<endl ;

	ofstream out ;
	out.open( strFileName.c_str() ) ;
	if( !out.is_open() )
	{
		cerr<<"Can not open the file of Hmm to write"<<endl ;
		exit( EXIT_FAILURE ) ;
	}
	//write the matrix A
	cout<<"Now output the matrix A "<<endl ;
	out<<N<<endl ;
	for( int i=0 ; i<N ; i++ )
	{
		for( int j=0 ; j<N ; j++ )
			out<<A[i][j]<<" " ;

		out<<endl ;
	}
	//write the matrix B
	cout<<"Now output the matrix B "<<endl ;
	out<<M<<endl ;
	for( int j=0 ; j<N ; j++ )
	{
		for( int k=0 ; k<M ; k++ )
			out<<B[j][k]<<" " ;

		out<<endl ;
	}
	//write the initial probability vector pi
	cout<<"Now output vector pi"<<endl ;
	for( i=0 ; i<N ; i++ )
	{
		out<<pi[i]<<" " ;
	}
	out<<endl ;

	out.close() ;
}

/*
void Hmm::WriteHmm( char* sFileName )
{
	cout<<"Now write Hmm to file "<<sFileName<<endl ;

	FILE* pFile ;
	pFile = fopen( sFileName, "w" ) ;
	if( pFile == NULL )
	{
		cerr<<"Can not open the file of Hmm to write"<<endl ;
		exit( EXIT_FAILURE ) ;
	}
	//write the matrix A
	cout<<"Now output the matrix A "<<endl ;
	fprintf( pFile, "%d\n", N ) ;
	for( int i=0 ; i<N ; i++ )
	{
		for( int j=0 ; j<N ; j++ )
			fprintf( pFile, "%f ", A[i][j] ) ;

		fprintf( pFile, "\n" ) ;
	}
	cout<<"Now output the matrix B "<<endl ;
	fprintf( pFile, "%d\n", M ) ;
	for( int j=0 ; j<N ; j++ )
	{
		for( int k=0 ; k<M ; k++ )
			fprintf( pFile, "%f ", B[j][k] ) ;

		fprintf( pFile, "\n" ) ;
	}
	cout<<"Now output the vector pi"<<endl ;
	for( i=0 ; i<N ; i++ )
		fprintf( pFile, "%f ", pi[i] ) ;
	fprintf( pFile, "\n" ) ;

	fclose( pFile ) ;
}
*/

void Hmm::WriteHmm( char* sFileName ) 
{
	cout<<"Now write Hmm to file "<<sFileName<<endl ;

	ofstream out ;
	out.open( sFileName ) ;
	if( !out.is_open() )
	{
		cerr<<"Can not open the file of Hmm to write"<<endl ;
		exit( EXIT_FAILURE ) ;
	}
	//write the matrix A
	cout<<"Now output the matrix A "<<endl ;
	out<<N<<endl ;
	for( int i=0 ; i<N ; i++ )
	{
		for( int j=0 ; j<N ; j++ )
			out<<A[i][j]<<" " ;

		out<<endl ;
	}
	//write the matrix B
	cout<<"Now output the matrix B "<<endl ;
	out<<M<<endl ;
	for( int j=0 ; j<N ; j++ )
	{
		for( int k=0 ; k<M ; k++ )
			out<<B[j][k]<<" " ;

		out<<endl ;
	}
	//write the initial probability vector pi
	cout<<"Now output vector pi"<<endl ;
	for( i=0 ; i<N ; i++ )
	{
		out<<pi[i]<<" " ;
	}
	out<<endl ;

	out.close() ;
}

void Hmm::ReadSequence( ifstream& in, vector<int>& sequence )
{
	cout<<"Now read the observation sequence from file"<<endl ;

	int iTemp ;
	if( !in.is_open() )
	{
		cerr<<"Wrong file stream"<<endl ;
		exit( EXIT_FAILURE ) ;
	}

	while( !in.eof() )
	{
		in>>iTemp ;
		sequence.push_back( iTemp ) ;
	}

}

void Hmm::ReadSequence( FILE* pFile, int*& sequence, int& iNum )
{
	if( pFile == NULL )
	{
		cerr<<"Wrong file stream"<<endl ;
		exit( EXIT_FAILURE ) ;
	}

	iNum = 0 ;
	int iSize = 0 ;
	sequence = (int*)malloc( sizeof(int) ) ;//just allocate memory for one item
	while( !feof(pFile) )
	{
		if( iNum > iSize )
		{
			iSize = 2*iNum ;
			sequence = (int*)realloc( sequence,iSize*sizeof(int) ) ;
		}
		if( sequence )
		{
			fscanf( pFile, "%d", &sequence[iNum] ) ;
			iNum++ ;
		}
		else
		{
			cerr<<"Error in reallcate memory"<<endl ;
		}
	}

}
/********************************the input-output function end ******************************************/

/********************************the core functions of Hmm **********************************************/
double Hmm::Forward( int T, vector<int>& O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( alpha != NULL )
	{
		//alpha is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of alpha
		dFreeMatrix( alpha, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	alpha = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int j=0 ; j<N ; j++ )
			alpha[t][j] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//***************this is the forward precedure below****************************************/

	double sum = 0.0 ;
	double dProbability = 0.0 ;

	//step one: initialize
	for( int i=0 ; i<N ; i++ )
	{
		alpha[0][i] = pi[i] * B[i][O[0]] ;
	}

	//step two ; induction
	for( t=1 ; t<T ; t++ )
	{
		for( int j=0 ; j<N ; j++ )
		{
			sum = 0.0 ;
			for( int i=0 ; i<N ; i++ )
			{
				sum += alpha[t-1][i] * A[i][j] ;
			}
			alpha[t][j] = sum * B[j][O[t]] ;
		}
	}

	//step three : termination
	for( i=0 ; i<N ; i++ )
	{
		dProbability += alpha[T-1][i] ;
	}

	return dProbability ;
}

double Hmm::Forward( int T, int* O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( alpha != NULL )
	{
		//alpha is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of alpha
		dFreeMatrix( alpha, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	alpha = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int j=0 ; j<N ; j++ )
			alpha[t][j] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//***************this is the forward precedure below****************************************/

	double sum = 0.0 ;
	double dProbability = 0.0 ;

	//step one: initialize
	for( int i=0 ; i<N ; i++ )
	{
		alpha[0][i] = pi[i] * B[i][O[0]] ;
	}

	//step two ; induction
	for( t=1 ; t<T ; t++ )
	{
		for( int j=0 ; j<N ; j++ )
		{
			sum = 0.0 ;
			for( int i=0 ; i<N ; i++ )
			{
				sum += alpha[t-1][i] * A[i][j] ;
			}
			alpha[t][j] = sum * B[j][O[t]] ;
		}
	}

	//step three : termination
	for( i=0 ; i<N ; i++ )
	{
		dProbability += alpha[T-1][i] ;
	}

	return dProbability ;
}

double Hmm::ForwardNormalized( int T, vector<int>& O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( alpha != NULL )
	{
		//alpha is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of alpha
		dFreeMatrix( alpha, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	alpha = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int j=0 ; j<N ; j++ )
			alpha[t][j] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//***************this is the forward normalized precedure below*****************************/

	double sum = 0.0 ;
	double dLogProbability = 0.0 ;
	//the scale variable to normalize the matrix alpha and return the log value of the sum of the matrix value
	vector<double> Scale ;
	Scale.resize( T ) ;

	//step one: initialize
	for( int i=0 ; i<N ; i++ )
	{
		alpha[0][i] = pi[i] * B[i][O[0]] ;
		Scale[0] += alpha[0][i] ;
	}
	//normalize
	for( i=0 ; i<N ; i++ )
		alpha[0][i] /= Scale[0] ;

	//step two ; induction
	for( t=1 ; t<T ; t++ )
	{
		Scale[t] = 0.0 ;
		for( int j=0 ; j<N ; j++ )
		{
			sum = 0.0 ;
			for( int i=0 ; i<N ; i++ )
			{
				sum += alpha[t-1][i] * A[i][j] ;
			}
			alpha[t][j] = sum * B[j][O[t]] ;
			Scale[t] += alpha[t][j] ;
		}
		//normalize
		for( j=0 ;j<N ; j++ )
			alpha[t][j] /= Scale[t] ;
	}

	//step three : termination
	/***************************************************************************************************/
	//the reference code return the log value of all the product of the probabiltiy of matrix alpha
	//that is, the sum of the vector of scale[t],why dose do so, i still can not understand.
	//i regard it's wrong. i think it should directly return the value of of log( Scale[T-1] )
	//note: the log function is to smooth the harss orignal value which is much small
	/***************************************************************************************************/
	dLogProbability = log( Scale[T-1] ) ;

	return dLogProbability ;
}

double Hmm::ForwardNormalized( int T, int* O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( alpha != NULL )
	{
		//alpha is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of alpha
		dFreeMatrix( alpha, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	alpha = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int j=0 ; j<N ; j++ )
			alpha[t][j] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//***************this is the forward normalized precedure below*****************************/

	double sum = 0.0 ;
	double dLogProbability = 0.0 ;
	//the scale variable to normalize the matrix alpha and return the log value of the sum of the matrix value
	vector<double> Scale ;
	Scale.resize( T ) ;

	//step one: initialize
	for( int i=0 ; i<N ; i++ )
	{
		alpha[0][i] = pi[i] * B[i][O[0]] ;
		Scale[0] += alpha[0][i] ;
	}
	//normalize
	for( i=0 ; i<N ; i++ )
		alpha[0][i] /= Scale[0] ;

	//step two : induction
	for( t=1 ; t<T ; t++ )
	{
		Scale[t] = 0.0 ;
		for( int j=0 ; j<N ; j++ )
		{
			sum = 0.0 ;
			for( int i=0 ; i<N ; i++ )
			{
				sum += alpha[t-1][i] * A[i][j] ;
			}
			alpha[t][j] = sum * B[j][O[t]] ;
			Scale[t] += alpha[t][j] ;
		}
		//normalize
		for( j=0 ;j<N ; j++ )
			alpha[t][j] /= Scale[t] ;
	}

	//step three : termination
	/***************************************************************************************************/
	//the reference code return the log value of all the product of the probabiltiy of matrix alpha
	//that is, the sum of the vector of scale[t],why dose do so, i still can not understand.
	//i regard it's wrong. i think it should directly return the value of of log( Scale[T-1] )
	//note: the log function is to smooth the harss orignal value which is much small
	/***************************************************************************************************/
	dLogProbability = log( Scale[T-1] ) ;

	return dLogProbability ;
}

double Hmm::Backward( int T, vector<int>& O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( beta != NULL )
	{
		//beta is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of beta
		dFreeMatrix( beta, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	beta = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int i=0 ; i<N ; i++ )
			beta[t][i] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//*********************** this is the backward procedure below***************************//

	double sum = 0.0 ;
	//step one : initialize
	for( int i=0 ; i<N ; i++ )
		beta[T-1][i] = 1 ;

	//step two : induction
	for( t=T-2 ; t>=0 ; t-- )
	{
		for( int i=0 ; i<N ; i++ )
		{
			sum = 0.0 ;
			for( int j=0 ; j<N ; j++ )
			{
				sum += A[i][j] * B[j][O[t+1]] * beta[t+1][j] ;
			}
			beta[t][i] = sum ;
		}
	}

	//step three : termination
	double dProbability = 0.0 ;
	for( i=0 ; i<N ; i++ )
		dProbability += pi[i] * B[i][O[0]] * beta[0][i] ;

	return dProbability ;
}

double Hmm::Backward( int T, int* O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( beta != NULL )
	{
		//beta is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of beta
		dFreeMatrix( beta, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	beta = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int i=0 ; i<N ; i++ )
			beta[t][i] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//*********************** this is the backward procedure below***************************//

	double sum = 0.0 ;
	//step one : initialize
	for( int i=0 ; i<N ; i++ )
		beta[T-1][i] = 1 ;

	//step two : induction
	for( t=T-2 ; t>=0 ; t-- )
	{
		for( int i=0 ; i<N ; i++ )
		{
			sum = 0.0 ;
			for( int j=0 ; j<N ; j++ )
			{
				sum += A[i][j] * B[j][O[t+1]] * beta[t+1][j] ;
			}
			beta[t][i] = sum ;
		}
	}

	//step three : termination
	double dProbability = 0.0 ;
	for( i=0 ; i<N ; i++ )
		dProbability += pi[i] * B[i][O[0]] * beta[0][i] ; 

	return dProbability ;
}

double Hmm::BackwardNormalized( int T, vector<int>& O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( beta != NULL )
	{
		//beta is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of beta
		dFreeMatrix( beta, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	beta = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int i=0 ; i<N ; i++ )
			beta[t][i] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//*********************** this is the backward procedure below***************************//

	double sum = 0.0 ;
	//the scale variable to normalize the matrix beta and return the log value of the sum of the matrix value
	vector<double> Scale ;
	Scale.resize( T ) ;

	//step one : initialize
	for( int i=0 ; i<N ; i++ )
		beta[T-1][i] = 1 ;
	//normalize
	for( i=0 ; i<N ; i++ )
		beta[T-1][i] /= N ;

	//step two : induction
	for( t=T-2 ; t>=0 ; t-- )
	{
		Scale[t] = 0.0 ;
		for( int i=0 ; i<N ; i++ )
		{
			sum = 0.0 ;
			for( int j=0 ; j<N ; j++ )
			{
				sum += A[i][j] * B[j][O[t+1]] * beta[t+1][j] ;
			}
			beta[t][i] = sum ;
			Scale[t] += beta[t][i] ;
		}
		//normalize
		for( i=0 ; i<N ; i++ )
			beta[t][i] /= Scale[t] ;
	}

	//step three : termination
	double dLogProbability = 0.0 ;
	double temp = 0.0 ;
	for( i=0 ; i<N ; i++ )
		temp += pi[0] * B[i][O[0]] * beta[0][i] ;
	dLogProbability = log( temp ) ;

	return dLogProbability ;
}

double Hmm::BackwardNormalized( int T, int* O )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( beta != NULL )
	{
		//beta is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of beta
		dFreeMatrix( beta, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
	}
	//allocate memory
	beta = dMatrix( 0, T-1, 0, N-1 ) ;
	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int i=0 ; i<N ; i++ )
			beta[t][i] = 0.0 ;
	//set the iPreT with the current value
	iPreT = T ;

	//*********************** this is the backward procedure below***************************//

	double sum = 0.0 ;
	//the scale variable to normalize the matrix beta and return the log value of the sum of the matrix value
	vector<double> Scale ;
	Scale.resize( T ) ;

	//step one : initialize
	for( int i=0 ; i<N ; i++ )
		beta[T-1][i] = 1 ;
	//normalize
	for( i=0 ; i<N ; i++ )
		beta[T-1][i] /= N ;

	//step two : induction
	for( t=T-2 ; t>=0 ; t-- )
	{
		Scale[t] = 0.0 ;
		for( int i=0 ; i<N ; i++ )
		{
			sum = 0.0 ;
			for( int j=0 ; j<N ; j++ )
			{
				sum += A[i][j] * B[j][O[t+1]] * beta[t+1][j] ;
			}
			beta[t][i] = sum ;
			Scale[t] += beta[t][i] ;
		}
		//normalize
		for( i=0 ; i<N ; i++ )
			beta[t][i] /= Scale[t] ;
	}

	//step three : termination
	double dLogProbability = 0.0 ;
	double temp = 0.0 ;
	for( i=0 ; i<N ; i++ )
		temp += pi[0] * B[i][O[0]] * beta[0][i] ;
	//dLogProbability = log( temp ) ;
	dLogProbability = log( Scale[0] ) ;

	return dLogProbability ;
}

double Hmm::Viterbi( int T, vector<int>& O, vector<int>& S )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( delta != NULL )
	{
		//delta is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of beta
		dFreeMatrix( delta, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
		//set the iPreT with the current value
		iPreT = T ;
	}
	//allocate memory
	delta = dMatrix( 0, T-1, 0, N-1 ) ;

	//allocate the memory for the state sequence S
	S.resize( T ) ;

	//***************** this is the viterbi algorithm********************************/

	//the matrix for back trace
	int** Psi ;
	//allocate memory for it
	Psi = iMatrix( 0, T-1, 0, N-1 ) ;

	//step one : initialization
	for( int i=0 ; i<N ; i++ )
	{
		delta[0][i] = pi[i] * B[i][O[0]] ;
		Psi[0][i] = 0 ;
	}

	//step two : induction
	double dMax = 0.0 ;
	double dValue = 0.0 ;
	for( int t=1 ; t<T ; t++ )
	{
		for( int j=0 ; j<N ; j++ )
		{
			dMax = 0.0 ;
			dValue = 0.0 ;
			for( int i=0 ; i<N ; i++ )
			{
				dValue = delta[t-1][i] * A[i][j] ;
				if( dMax < dValue )
				{
					dMax = dValue ;
					Psi[t][j] = i ;
				}
			}
			delta[t][j] = dMax * B[j][O[t]] ;
		}
	}

	//step three : termination
	dMax = 0.0 ;
	int iIndex = 0 ;
	for( i=0 ; i<N ; i++ )
	{
		if( dMax < delta[T-1][i] )
		{
			dMax = delta[T-1][i] ;
			iIndex = i ;
		}
	}
	S[T-1] = iIndex ;
	for( t=T-2 ; t>=0 ; t-- )
		S[t] = Psi[t+1][S[t+1]] ;

	//free the memory of Psi
	iFreeMatrix( Psi, 0, T-1, 0, N-1 ) ;

	return dMax ;

}

double Hmm::Viterbi( int T, int* O, int*& S )
{
	//to check the observation sequence
	if( T == 0 )
		return FALSEVALUE ;
	if( delta != NULL )
	{
		//delta is not empty indicate that this is a new observation sequence
		//then reallocate the memory is needed
		//fisrt free the memory of beta
		dFreeMatrix( delta, 0, iPreT-1, 0, N-1 ) ;/*the fourth parament is no use here since the third one is always 0*/
		//set the iPreT with the current value
		iPreT = T ;
	}
	//allocate memory
	delta = dMatrix( 0, T-1, 0, N-1 ) ;

	//allocate the memory for the state sequence S
	S = (int*)malloc( T*sizeof(int) ) ;

	//***************** this is the viterbi algorithm********************************/

	//the matrix for back trace
	int** Psi ;
	//allocate memory for it
	Psi = iMatrix( 0, T-1, 0, N-1 ) ;

	//step one : initialization
	for( int i=0 ; i<N ; i++ )
	{
		delta[0][i] = pi[i] * B[i][O[0]] ;
		Psi[0][i] = 0 ;
	}

	//step two : induction
	double dMax = 0.0 ;
	double dValue = 0.0 ;
	for( int t=1 ; t<T ; t++ )
	{
		for( int j=0 ; j<N ; j++ )
		{
			dMax = 0.0 ;
			dValue = 0.0 ;
			for( int i=0 ; i<N ; i++ )
			{
				dValue = delta[t-1][i] * A[i][j] ;
				if( dMax < dValue )
				{
					dMax = dValue ;
					Psi[t][j] = i ;
				}
			}
			delta[t][j] = dMax * B[j][O[t]] ;
		}
	}

	//step three : termination
	dMax = 0.0 ;
	int iIndex = 0 ;
	for( i=0 ; i<N ; i++ )
	{
		if( dMax < delta[T-1][i] )
		{
			dMax = delta[T-1][i] ;
			iIndex = i ;
		}
	}
	S[T-1] = iIndex ;
	for( t=T-2 ; t>=0 ; t-- )
		S[t] = Psi[t+1][S[t+1]] ;

	//free the memory of Psi
	iFreeMatrix( Psi, 0, T-1, 0, N-1 ) ;

	return dMax ;

}

void Hmm::BaumWelch( int T, vector<int>& O, double& probInit, double& probFinal )
{
	//step one : calculate alpha and initial probability using forward procedure
	probInit = Forward( T, O ) ;

	//step two : calculate beta using backward procedure
	Backward( T, O ) ;//first it's equal to probInit

	//step three : allocate memory for gamma and calculate it
	gamma = ComputeGamma( T, N, probInit ) ;

	//step four : allocate memory for p[t][i][j] and calculate it
	//p[t][i][j] = p( Xt=i,Xt+1=j | O,u ), is the probability of transitting from state i to state j
	double*** p = NULL ;
	p = ComputeP( p, T, N, O, probInit ) ;

	//step five : recalculate all the parameters of hmm( pi, A and B, and also alpha,beta,gamma )
	//untill some consition is abtained
	RecomputeParameter( p, T, O, probInit, probFinal ) ;

	//step six : free all the memory that allocated before
	PFreeMatrix( p, T ) ;

}

void Hmm::BaumWelch( int T, int* O, double& probInit, double& probFinal )
{
	//step one : calculate alpha and initial probability using forward procedure
	probInit = Forward( T, O ) ;

	//step two : calculate beta using backward procedure
	Backward( T, O ) ;//first it's equal to probInit

	//step three : allocate memory for gamma and calculate it
	gamma = ComputeGamma( T, N, probInit ) ;

	//step four : allocate memory for p[t][i][j] and calculate it
	//p[t][i][j] = p( Xt=i,Xt+1=j | O,u ), is the probability of transitting from state i to state j
	double*** p = NULL ;
	p = ComputeP( p, T, N, O, probInit ) ;

	//step five : recalculate all the parameters of hmm( pi, A and B, and also alpha,beta,gamma )
	//untill some consition is abtained
	RecomputeParameter( p, T, O, probInit, probFinal ) ;

	//step six : free all the memory that allocated before
	PFreeMatrix( p, T ) ;
}

/********************************the core functions of Hmm end*******************************************/

/******************************** the function called by baumwelch ***************************************/
double** Hmm::ComputeGamma( int T, int N, double prob ) 
{
	//firstly allocate memory for gamma
	if( gamma )
	{
		//if gamma is not null, then we should free the memory that gamma ocuppied
		dFreeMatrix( gamma, 0, iPreT-1, 0, N-1 ) ;
	}
	gamma = dMatrix( 0, T-1, 0, N-1 ) ;

	//initialize
	for( int t=0 ; t<T ; t++ )
		for( int i=0 ; i<N ; i++ )
			gamma[t][i] = 0.0 ;

	//calculate gamma using alpha and beta
	double dDenominator = DENOMINATORLIMIT ;
	for( t=0 ; t<T ; t++ )
	{
		for( int i=0 ; i<N ; i++ )
		{
			gamma[t][i] = alpha[t][i] * beta[t][i] ;
			gamma[t][i] /= prob ;
		}
	}

	return gamma ;
}

double*** Hmm::ComputeP( double*** p ,int T, int N, vector<int>& O, double prob )
{
	//firstly allocate memory for p
	if( p )
	{
		//if p is not null, then we should free the memory that p ocuppied
		for( int t=0 ; t<T ; t++ )
			dFreeMatrix( p[t], 0, N-1, 0, N-1 ) ;
		delete[] p ;
	}
	p = new double** [T] ;
	for( int t=0 ; t<T ; t++ )
		p[t] = dMatrix( 0, N-1, 0, N-1 ) ;

	//initialize the 3-dimision matrix
	for( t=0 ; t<T ;t++ )
		for( int i=0 ; i<N ; i++ )
			for( int j=0 ; j<N ; j++ )
				p[t][i][j] = 0.0 ;

	//calculate it
	for( t=0 ; t<T-1 ; t++ ) 
	//attension : the time variable is from 0 to T-2, together t-1 time, because there are only t-1 transittions!
	{
		for( int i=0 ; i<N ; i++ )
		{
			for( int j=0 ; j<N ; j++ )
			{
				p[t][i][j] = alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j] ;
				p[t][i][j] /= prob ;
			}
		}
	}

	return p ;
}

double*** Hmm::ComputeP( double*** p, int T, int N, int* O, double prob )
{
	//firstly allocate memory for p
	if( p )
	{
		//if p is not null, then we should free the memory that p ocuppied
		for( int t=0 ; t<T ; t++ )
			dFreeMatrix( p[t], 0, N-1, 0, N-1 ) ;
		delete[] p ;
	}
	p = new double** [T] ;
	for( int t=0 ; t<T ; t++ )
		p[t] = dMatrix( 0, N-1, 0, N-1 ) ;

	//initialize the 3-dimision matrix
	for( t=0 ; t<T ;t++ )
		for( int i=0 ; i<N ; i++ )
			for( int j=0 ; j<N ; j++ )
				p[t][i][j] = 0.0 ;

	//calculate it
	for( t=0 ; t<T-1 ; t++ ) 
	//attension : the time variable is from 0 to T-2, together t-1 time, because there are only t-1 transittions!
	{
		for( int i=0 ; i<N ; i++ )
		{
			for( int j=0 ; j<N ; j++ )
			{
				p[t][i][j] = alpha[t][i] * A[i][j] * B[j][O[t+1]] * beta[t+1][j] ;
				p[t][i][j] /= prob ;
			}
		}
	}

	return p ;
}

void Hmm::RecomputeParameter( double***& p, int T, vector<int>& O, double probInit, double& probFinal )
{
	//record the round number of do..while circle and jump out of it when certain number exceed
	int iRound = 0 ;

	double dRatio = 0.0 ;
	double probTemp = 0.0 ;
	double probOld = probInit ;
	do{
		//recalculate initial probability
		pi.clear() ;
		for( int i=0 ; i<N ; i++ )
			pi.push_back( gamma[0][i] ) ;

		//recalculate the transition matrix A probability
		//double dDenominatorA = 0.0 ;
		double dDenominatorA = DENOMINATORLIMIT ;
		double dNumeratorA = 0.0 ;
		for( i=0 ; i<N ; i++ )
		{
			//first calculate the denominator of A[i][j], that's the sum of gamma[t][i]
			//(excluded the time T-1 because there is no transition from the time T-1)
			dDenominatorA = DENOMINATORLIMIT ;
			for( int t=0 ; t<T-1 ; t++ )
				dDenominatorA += gamma[t][i] ;
			for( int j=0 ; j<N ; j++ )
			{
				//second calculate the numerator of A[i][j], that's the sum of p[t][i][j]
				//(excluded the time T-1 because there is no transition from the time T-1)
				dNumeratorA = 0.0 ;
				for( t=0 ; t<T-1 ; t++ )
					dNumeratorA += p[t][i][j] ;

				//recalculate the probability of A[i][j]
				A[i][j] = dNumeratorA / dDenominatorA ;
			}
		}

		//recalculate the emit probability matrix of B
		double dDenominatorB = DENOMINATORLIMIT ;
		double dNumeratorB = 0.0 ;
		for( int j=0 ; j<N ; j++ )
		{
			dDenominatorB = DENOMINATORLIMIT ;
			for( int t=0 ; t<T ; t++ )
				dDenominatorB += gamma[t][j] ;
			for( int k=0 ; k<M ; k++ )
			{
				dNumeratorB = 0.0 ;
				for( int t=0 ; t<T ;t++ )
				{
					if( O[t] == k )
						dNumeratorB += gamma[t][j] ;
				}

				//recalculate the probability of B[j][k]
				B[j][k] = dNumeratorB / dDenominatorB ;
			}
		}

		//calculate alpha and the probability once again, repeat the "step one"
		probTemp = Forward( T, O ) ;
		//calculate beta , repeat the "step two"
		Backward( T, O ) ;
		//allocate memory for gamma and calculate it, repeat the "step three"
		gamma = ComputeGamma( T, N, probTemp ) ;
		//recalculate p[t][i][j]
		p = ComputeP( p, T, N, O, probTemp ) ;

		dRatio = ( probTemp - probOld ) / probOld ;
		probOld = probTemp ;

		//when the circle number of do..while exceed ROUNDLIMIT times, break it
		//thus, we may not get stuck in a local maximum
		iRound++ ;
		if( iRound > ROUNDLIMIT )
			break ;

	}while( dRatio>RATIOLIMIT ) ;

	probFinal = probTemp ;

}

void Hmm::RecomputeParameter( double***& p, int T, int* O, double probInit, double& probFinal )
{
	//record the round number of do..while circle and jump out of it when certain number exceed
	int iRound = 0 ;

	double dRatio = 0.0 ;
	double probTemp = 0.0 ;
	double probOld = probInit ;
	do{
		//recalculate initial probability
		pi.clear() ;
		for( int i=0 ; i<N ; i++ )
			pi.push_back( gamma[0][i] ) ;

		//recalculate the transition matrix A probability
		//double dDenominatorA = 0.0 ;
		double dDenominatorA = DENOMINATORLIMIT ;
		double dNumeratorA = 0.0 ;
		for( i=0 ; i<N ; i++ )
		{
			//first calculate the denominator of A[i][j], that's the sum of gamma[t][i]
			//(excluded the time T-1 because there is no transition from the time T-1)
			dDenominatorA = DENOMINATORLIMIT ;
			for( int t=0 ; t<T-1 ; t++ )
				dDenominatorA += gamma[t][i] ;
			for( int j=0 ; j<N ; j++ )
			{
				//second calculate the numerator of A[i][j], that's the sum of p[t][i][j]
				//(excluded the time T-1 because there is no transition from the time T-1)
				dNumeratorA = 0.0 ;
				for( t=0 ; t<T-1 ; t++ )
					dNumeratorA += p[t][i][j] ;

				//recalculate the probability of A[i][j]
				A[i][j] = dNumeratorA / dDenominatorA ;
			}
		}

		//recalculate the emit probability matrix of B
		double dDenominatorB = DENOMINATORLIMIT ;
		double dNumeratorB = 0.0 ;
		for( int j=0 ; j<N ; j++ )
		{
			dDenominatorB = DENOMINATORLIMIT ;
			for( int t=0 ; t<T ; t++ )
				dDenominatorB += gamma[t][j] ;
			for( int k=0 ; k<M ; k++ )
			{
				dNumeratorB = 0.0 ;
				for( int t=0 ; t<T ;t++ )
				{
					if( O[t] == k )
						dNumeratorB += gamma[t][j] ;
				}

				//recalculate the probability of B[j][k]
				B[j][k] = dNumeratorB / dDenominatorB ;
			}
		}

		//calculate alpha and the probability once again, repeat the "step one"
		probTemp = Forward( T, O ) ;
		//calculate beta , repeat the "step two"
		Backward( T, O ) ;
		//allocate memory for gamma and calculate it, repeat the "step three"
		gamma = ComputeGamma( T, N, probTemp ) ;
		//recalculate p[t][i][j]
		p = ComputeP( p, T, N, O, probTemp ) ;

		dRatio = ( probTemp - probOld ) / probOld ;
		probOld = probTemp ;

		//when the circle number of do..while exceed ROUNDLIMIT times, break it
		//thus, we may not get stuck in a local maximum
		iRound++ ;
		if( iRound > ROUNDLIMIT )
			break ;

	}while( dRatio>RATIOLIMIT ) ;

	probFinal = probTemp ;
}

void Hmm::PFreeMatrix( double*** p, int T )
{
	for( int t=0 ; t<T ; t++ )
		dFreeMatrix( p[t], 0, N-1, 0, N-1 ) ;

	delete[] p ;
}
/****************** the function called by baumwelch end **********************************/

/************************the generating sequence function*******************************/
void Hmm::GenerateSequence( int iSeed, int T, vector<int>& O, vector<int>& S )
{
	//allocate memory for vector O and S
	O.resize( T ) ;
	S.resize( T ) ;

	//get the initial state
	S[0] = GetInitialState( iSeed ) ;
	O[0] = GetSymbol( S[0] ) ;

	//get the rest state and symbol
	for( int t=1 ; t<T ; t++ )
	{
		S[t] = GetNextState( S[t-1] ) ;
		O[t] = GetSymbol( S[t] ) ;
	}
}

void Hmm::GenerateSequence( int iSeed, int T, int* O, int*S )
{
	//get the initial state
	S[0] = GetInitialState( iSeed ) ;
	O[0] = GetSymbol( S[0] ) ;

	//get the rest state and symbol
	for( int t=1 ; t<T ; t++ )
	{
		S[t] = GetNextState( S[t-1] ) ;
		O[t] = GetSymbol( S[t] ) ;
	}
}

int Hmm::GetInitialState( int iSeed ) 
{
	double dVal = 0.0 ;
	double dAccum = 0.0 ;
	int iState = -1 ;

	//set the seed
	srand( iSeed ) ;
	//generate a value between 0 and 1
	dVal = rand() / RAND_MAX ;
	
	//get a state according the state probability
	for( int i=0 ; i<N ; i++ )
	{
		if( dVal < (pi[i]+dAccum) )
		{
			iState = i ;
			break ;
		}
		else
			dAccum += pi[i] ;
	}

	return iState ;
}

int Hmm::GetSymbol( int iState )
{
	double dVal = 0.0 ;
	double dAccum = 0.0 ;

	int iSymbol = -1 ;

	int iSeed = 0 ;
	//get the current time as the seed
	long ltime = 0 ;
	ltime = time( NULL ) ;
	iSeed = (unsigned)ltime/2 ;

	//set the seed
	srand( iSeed ) ;
	//generate a value between 0 and 1
	dVal = rand() / RAND_MAX ;

	//get an output symbol according the emition probability
	for( int k=0 ; k<M ; k++ )
	{
		if( dVal < (B[iState][k]+dAccum) )
		{
			iSymbol = k ;
			break ;
		}
		else
			dAccum += B[iState][k] ;
	}

	return iSymbol ;
}

int Hmm::GetNextState( int iState )
{
	double dVal = 0.0 ;
	double dAccum = 0.0 ;

	int iStateNext = -1 ;

	int iSeed = 0 ;
	//get the current time as the seed
	long ltime = 0 ;
	ltime = time( NULL ) ;
	iSeed = (unsigned)ltime/2 ;
	//set the seed
	srand( iSeed ) ;
	//generate a value between 0 and 1
	dVal = rand() / RAND_MAX ;

	//get the next state according the state transittion probability
	for( int i=0 ; i<M ; i++ )
	{
		if( dVal < (A[iState][i]+dAccum) )
		{
			iStateNext = i ;
			break ;
		}
		else
			dAccum += A[iState][i] ;
	}

	return iStateNext ;
}

/**********************the generating sequence function end **************************/