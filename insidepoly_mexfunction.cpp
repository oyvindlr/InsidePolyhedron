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

static void convertToInt( const double* from, int* to, int n )
{
	for ( int i = 0; i < n; i++ )
		to[i] = (int)from[i];
}

void mexFunction( int nlhs, mxArray* plhs[], int nrhs,  const mxArray* prhs[] )
{
	double* vertices_t = mxGetDoubles( prhs[0] );
	size_t nVertices = mxGetM( prhs[0] );

	double* faces_d = mxGetDoubles( prhs[1] );
	size_t nFaces = mxGetNumberOfElements( prhs[1] ) / 3;
	int *faces = new int[nFaces * 3];
	convertToInt( faces_d, faces, nFaces * 3 );
	transpose( faces, nFaces, 3 );
	for ( int i = 0; i < nFaces * 3; i++ )
		faces[i]--;//Convert from matlab 1-based index to C 0-based.
	double * vertices = new double[nVertices * 3];
	memcpy( vertices, vertices_t, nVertices * 3 * sizeof( double ) );
	transpose( vertices, nVertices, 3 );
	
	double* x = mxGetDoubles( prhs[2] );
	double* y = mxGetDoubles( prhs[3] );
	double* z = mxGetDoubles( prhs[4] );
	int nx = mxGetNumberOfElements( prhs[2] );
	int ny = mxGetNumberOfElements( prhs[3] );
	int nz = mxGetNumberOfElements( prhs[4] );
	
	mwSize outdims[3] = { ny, nx, nz };
	plhs[0] = mxCreateLogicalArray( 3, outdims );
	bool* inside = mxGetLogicals( plhs[0] );


	insidePolyhedron( inside, (double(*)[3])vertices, (int(*)[3])faces, nFaces, x, nx, y, ny, z, nz );

	delete[] faces;
	delete[] vertices;
}