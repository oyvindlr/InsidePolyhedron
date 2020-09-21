#include "pch.h"
#include "mex.h"
#include "FindInsideOfPolyhedron.h"

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
	int argNo = 0;
	double* vertices = nullptr;
	int* faceIndices = nullptr;
	const double* faces = nullptr;
	size_t nFaces;
	bool areFacesGiven = false;
	if (nrhs != 4 && nrhs != 5)
		mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument", "4 or 5 input arguments required");
	if ( mxGetNumberOfDimensions( prhs[0] ) == 3 )
	{
		if (nrhs != 4)
			mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument", "4 input arguments required when first argument is face matrix");
		if (mxGetN(prhs[0]) != 9)
			mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument", "Incorrect dimensions of first argument (should be n by 3 by 3 or n by 3)");
		faces = mxGetDoubles( prhs[argNo] );
		nFaces = mxGetNumberOfElements( prhs[argNo] ) / 9;
		argNo++;
		areFacesGiven = true;
	}
	else
	{
		if (nrhs != 5)
			mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument1", "5 input arguments required when first argument is vertex list");
		if (mxGetN(prhs[0]) != 3)
			mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument2", "Incorrect dimensions of first argument (should be n by 3 by 3 or n by 3)");
		if (mxGetN(prhs[1]) != 3)
			mexErrMsgIdAndTxt("InsidePolyhedron:IncorrectArgument", "Incorrect dimensions of second argument (faceIndices) (should be n by 3)");
		const double* vertices_t = mxGetDoubles( prhs[argNo] );
		size_t nVertices = mxGetM( prhs[argNo] );
		argNo++;
		const double* faces_d = mxGetDoubles( prhs[argNo] );
		nFaces = mxGetNumberOfElements( prhs[argNo] ) / 3;
		argNo++;
		faceIndices = new int[nFaces * 3];
		convertToIndex( faces_d, faceIndices, nFaces * 3 );
		transpose( faceIndices, nFaces, 3 );
		vertices = new double[nVertices * 3];
		memcpy( vertices, vertices_t, nVertices * 3 * sizeof( double ) );
		transpose( vertices, nVertices, 3 );
		areFacesGiven = false;
	}
	
	size_t nx = mxGetNumberOfElements( prhs[argNo] );
	const double* x = mxGetDoubles( prhs[argNo++] );
	size_t ny = mxGetNumberOfElements( prhs[argNo] );
	const double* y = mxGetDoubles( prhs[argNo++] );
	size_t nz = mxGetNumberOfElements( prhs[argNo] );
	const double* z = mxGetDoubles( prhs[argNo++] );
	
	mwSize outdims[3] = { ny, nx, nz };
	plhs[0] = mxCreateLogicalArray( 3, outdims );
	bool* inside = mxGetLogicals( plhs[0] );

	if ( areFacesGiven )
		insidePolyhedron( inside, (double (*)[3][3]) faces, nFaces, x, nx, y, ny, z, nz );
	else
	{
		insidePolyhedron( inside, ( double( * )[3] )vertices, ( int( * )[3] )faceIndices, nFaces, x, nx, y, ny, z, nz );
		delete[] vertices;
		delete[] faceIndices;
	}
}