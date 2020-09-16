#include "pch.h"
#include "FindInsideOfPolyhedron.h"
#include <cmath>
#include <algorithm>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define printf mexPrintf
#endif

static int dim0, dim1, dim2;

static void findExtremeCoords( const double faces[][3][3], double minCoords[][3], double maxCoords[][3], size_t nFaces )
{
	for ( int i = 0; i < nFaces; i++ )
	{
		for ( int dim = 0; dim < 3; dim++ )
		{
			minCoords[i][dim] = faces[i][0][dim];
			if ( faces[i][1][dim] < minCoords[i][dim] )
				minCoords[i][dim] = faces[i][1][dim];
			if ( faces[i][2][dim] < minCoords[i][dim] )
				minCoords[i][dim] = faces[i][2][dim];
			maxCoords[i][dim] = faces[i][0][dim];
			if ( faces[i][1][dim] > maxCoords[i][dim] )
				maxCoords[i][dim] = faces[i][1][dim];
			if ( faces[i][2][dim] > maxCoords[i][dim] )
				maxCoords[i][dim] = faces[i][2][dim];
		}
	}
}

static int findFacesInDim( int *facesIndex, const double minCoords[][3], const double maxCoords[][3], double value, int dim, size_t nFaces )
{
	int count = 0;
	for ( int i = 0; i < nFaces; i++ )
	{
		if ( minCoords[i][dim] < value && maxCoords[i][dim] > value )
		{
			facesIndex[count] = i;
			count++;
		}
	}
	return count;
}

static void selectFaces( double facesZ[][3][3], const double faces[][3][3], const int facesIndex[], int nFacesZ )
{
	for ( int i = 0; i < nFacesZ; i++ )
		for ( int j = 0; j < 3; j++ )
			for ( int k = 0; k < 3; k++ )
				facesZ[i][j][k] = faces[facesIndex[i]][j][k];
}

static void selectCoords( double minCoordsZ[][3], double maxCoordsZ[][3], const double minCoords[][3], const double maxCoords[][3], const int facesIndex[], size_t nFacesZ )
{
	for ( int i = 0; i < nFacesZ; i++ )
	{
		for ( int j = 0; j < 3; j++ )
		{
			minCoordsZ[i][j] = minCoords[facesIndex[i]][j];
			maxCoordsZ[i][j] = maxCoords[facesIndex[i]][j];
		}
	}
}

static void solve2by2otherway( double A[2][2], double b[2] )
{
	double fac = A[0][0] / A[1][0];
	double a12 = A[0][1] - A[1][1] * fac;

	if ( abs( a12 ) < 1e-10 || abs( A[1][0] ) < 1e-10 )
		printf( "Near singular\n" );
	double b1 = ( b[0] - b[1] * fac ) / a12;
	b[0] = ( b[1] - A[1][1] * b1 ) / A[1][0];
	b[1] = b1;
}


//Solve a 2 by 2 linear equation by Gaussian elimination with partial pivoting
static void solve2by2( double A[2][2], double b[2] )
{
	if ( abs( A[0][0] ) < abs( A[1][0] ) )//Partial pivoting
	{
		solve2by2otherway( A, b );
		return;
	}
	double fac = A[1][0] / A[0][0];
	double a22 = A[1][1] - A[0][1] * fac;
	if ( abs( a22 ) < 1e-10 || abs( A[0][0] ) < 1e-10 )
		printf( "Near singular\n" );
	b[1] = ( b[1] - b[0] * fac ) / a22;
	b[0] = ( b[0] - A[0][1] * b[1] ) / A[0][0];
}

static int getCrossings( double crossings[], const double faces[][3][3], int nFaces, double dim1Value, double dim2Value )
{
	double b[2];
	double A[2][2];
	int nCrossings = 0;

	for ( int i = 0; i < nFaces; i++ )
	{
		b[0] = dim1Value - faces[i][0][dim1];
		b[1] = dim2Value - faces[i][0][dim2];
		A[0][0] = faces[i][1][dim1] - faces[i][0][dim1];
		A[1][0] = faces[i][1][dim2] - faces[i][0][dim2];
		A[0][1] = faces[i][2][dim1] - faces[i][0][dim1];
		A[1][1] = faces[i][2][dim2] - faces[i][0][dim2];
		solve2by2( A, b );
		if ( b[0] > 0 && b[1] > 0 && ( ( b[0] + b[1] ) < 1 ) )
		{
			crossings[nCrossings] = faces[i][0][dim0] + b[0] * ( faces[i][1][dim0] - faces[i][0][dim0]) + b[1] * ( faces[i][2][dim0] - faces[i][0][dim0] );
			nCrossings++;
		}
	}
	if ( nCrossings > 0 )
		std::sort( crossings, crossings + nCrossings );
	return nCrossings;
}


static inline bool isOdd( int n )
{
	return ( n % 2 ) == 1;
}

static void buildFaceMatrix( double faces[][3][3], const double vertices[][3], const int faceIndices[][3], size_t nFaces )
{
	for ( int i = 0; i < nFaces; i++ )
		for ( int j = 0; j < 3; j++ )
			for ( int k = 0; k < 3; k++ )
				faces[i][j][k] = vertices[faceIndices[i][j]][k];
}

static void selectDimensionsForFastestProcessing( size_t *nx, size_t *ny, size_t *nz, size_t dimStep[3])
{
	size_t sizes[] = { *nx, *ny, *nz };
	int dims[3] = { 0, 1, 2 };
	size_t dimMult[] = { *ny * *nx, 1, *ny };

	for ( int i = 0; i < 2; i++ )
	{
		for ( int j = 0; j < 2 - i; j++ )
		{
			if ( sizes[j] > sizes[j + 1] )
			{
				size_t temp = sizes[j];
				sizes[j] = sizes[j + 1];
				sizes[j + 1] = temp;
				int tempi = dims[j];
				dims[j] = dims[j + 1];
				dims[j + 1] = tempi;
			}
		}
	}
	dim2 = dims[0];
	dim1 = dims[1];
	dim0 = dims[2];
	dimStep[dim0] = dimMult[0];
	dimStep[dim1] = dimMult[1];
	dimStep[dim2] = dimMult[2];
	*nz = sizes[0];
	*ny = sizes[1];
	*nx = sizes[2];
}


void insidePolyhedron( bool inside[], const double vertices [][3], const int faceIndices[][3], size_t nFaces, const double x[], size_t nx, const double y[], size_t ny, const double z[], size_t nz )
{
	double( *faces )[3][3] = new double[nFaces][3][3];
	double( *maxCoords )[3] = new double[nFaces][3];
	double( *minCoords )[3] = new double[nFaces][3];
	int* facesIndex = new int[nFaces];
	double( *facesD2 )[3][3] = new double[nFaces][3][3];
	double( *facesD1 )[3][3] = new double[nFaces][3][3];
	double( *minCoordsD2 )[3] = new double[nFaces][3];
	double( *maxCoordsD2 )[3] = new double[nFaces][3];
	double* crossings = new double[nFaces];
	size_t dimSteps[3];

	selectDimensionsForFastestProcessing( &nx, &ny, &nz, dimSteps );

	size_t nxny = nx * ny;

	buildFaceMatrix( faces, vertices, faceIndices, nFaces );

	findExtremeCoords( faces, minCoords, maxCoords, nFaces );

	const double* xyz[] = { x, y, z };

	for ( int i = 0; i < nz; i++ )
	{
		int nFacesD2 = findFacesInDim( facesIndex, minCoords, maxCoords, xyz[dim2][i], dim2, nFaces );
		if ( nFacesD2 == 0 )
			continue;
		selectFaces( facesD2, faces, facesIndex, nFacesD2 );
		selectCoords( minCoordsD2, maxCoordsD2, minCoords, maxCoords, facesIndex, nFacesD2 );
		for ( int j = 0; j < ny; j++ )
		{
			int nFacesD1 = findFacesInDim( facesIndex, minCoordsD2, maxCoordsD2, xyz[dim1][j], dim1, nFacesD2 );
			if ( nFacesD1 == 0 )
				continue;
			selectFaces( facesD1, facesD2, facesIndex, nFacesD1 );
			int nCrossings = getCrossings( crossings, facesD1, nFacesD1, xyz[dim1][j], xyz[dim2][i] );
			if ( nCrossings == 0 )
				continue;
			if ( isOdd( nCrossings ) )
			{
				printf( "Odd number of crossings found. Is the polyhedron not closed?\n" );
				return;
			}
			bool isInside = false;
			int crossingsPassed = 0;
			for ( int k = 0; k < nx; k++ )
			{
				while ( ( crossingsPassed < nCrossings ) && ( crossings[crossingsPassed] < xyz[dim0][k] ) )
				{
					crossingsPassed++;
					isInside = !isInside;
				}
				inside[i * dimSteps[0] + j * dimSteps[1]+ k * dimSteps[2]] = isInside;
			}
		}
	}
	delete[] maxCoords;
	delete[] minCoords;
	delete[] facesIndex;
	delete[] facesD2;
	delete[] facesD1;
	delete[] minCoordsD2;
	delete[] maxCoordsD2;
	delete[] crossings;
	delete[] faces;
}

