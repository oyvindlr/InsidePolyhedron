

/** \file insidepoly_mexfunction.cpp
*	\brief Matlab mex interface for FindInsideOfPolyhedron
*	\author Oyvind L Rortveit
*	\date 2020
*/


#include "pch.h"

#ifdef MATLAB_MEX_FILE
#include "FindInsideOfPolyhedron.h"
#include "mex.h"


template <class T>

static void transpose( T* data, size_t m, size_t n )
{
	T* newData = new T[m * n];
	for ( int i = 0; i < m; i++ )
		for ( int j = 0; j < n; j++ )
			newData[i * n + j] = data[j * m + i];
	for ( int i = 0; i < m * n; i++ )
		data[i] = newData[i];
	delete[] newData;
}

/** In Matlab, all array-indices are doubles, and indexing starts at 1. 
This function converts the doubles to ints and subtracts 1 to 
get 0-based indexing */
static void convertToIndex( const double* from, int* to, size_t n )
{
	for ( int i = 0; i < n; i++ )
		to[i] = ((int)from[i]) - 1;
}

/** Entry point for Matlab function InsidePolyhedron. See InsidePolyhedron.m for documentation */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs,  const mxArray* prhs[] )
{
	double* vertices = nullptr;
	int* faceIndices = nullptr;
	const double* faces = nullptr;
	size_t nFaces;

	if (nrhs != 5)
		mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument", "5 input arguments required");
	if (mxGetN(prhs[0]) != 3)
		mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument2", "Incorrect dimensions of first argument (should be n by 3)");
	if (mxGetN(prhs[1]) != 3)
		mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument", "Incorrect dimensions of second argument (should be n by 3)");
	const double* vertices_t = mxGetDoubles( prhs[0] );
	size_t nVertices = mxGetM( prhs[0] );
	const double* faces_d = mxGetDoubles( prhs[1] );
	nFaces = mxGetNumberOfElements( prhs[1] ) / 3;
	faceIndices = new int[nFaces * 3];
	convertToIndex( faces_d, faceIndices, nFaces * 3 );
	transpose( faceIndices, nFaces, 3 );
	vertices = new double[nVertices * 3];
	memcpy( vertices, vertices_t, nVertices * 3 * sizeof( double ) );
	transpose( vertices, nVertices, 3 );
	
	size_t nx = mxGetNumberOfElements( prhs[2] );
	const double* x = mxGetDoubles( prhs[2] );
	size_t ny = mxGetNumberOfElements( prhs[3] );
	const double* y = mxGetDoubles( prhs[3] );
	size_t nz = mxGetNumberOfElements( prhs[4] );
	const double* z = mxGetDoubles( prhs[4] );
	
	mwSize outdims[3] = { ny, nx, nz };
	plhs[0] = mxCreateLogicalArray( 3, outdims );
	bool* inside = mxGetLogicals( plhs[0] );

	insidePolyhedron( inside, ( double( * )[3] )vertices, ( int( * )[3] )faceIndices, nFaces, x, nx, y, ny, z, nz );
	delete[] vertices;
	delete[] faceIndices;
}
#endif
