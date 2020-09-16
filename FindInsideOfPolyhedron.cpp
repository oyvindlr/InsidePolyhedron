#include "pch.h"
#include "FindInsideOfPolyhedron.h"
#include <cmath>
#include <algorithm>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define printf mexPrintf
#endif



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

static int getCrossings( double crossings[], const double faces[][3][3], int nFaces, double y, double z )
{
	double v1[3];
	double v2[3];
	double p[3];
	double b[2];
	double A[2][2];
	int nCrossings = 0;

	for ( int i = 0; i < nFaces; i++ )
	{
		for ( int j = 0; j < 3; j++ )
		{
			p[j] = faces[i][0][j];
			v1[j] = faces[i][1][j] - p[j];
			v2[j] = faces[i][2][j] - p[j];
		}
		b[0] = y - p[1];
		b[1] = z - p[2];
		A[0][0] = v1[1];
		A[1][0] = v1[2];
		A[0][1] = v2[1];
		A[1][1] = v2[2];
		solve2by2( A, b );
		if ( b[0] > 0 && b[1] > 0 && ( ( b[0] + b[1] ) < 1 ) )
		{
			crossings[nCrossings] = p[0] + b[0] * v1[0] + b[1] * v2[0];
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


void insidePolyhedron( bool inside[], const double vertices [][3], const int faceIndices[][3], size_t nFaces, const double x[], size_t nx, const double y[], size_t ny, const double z[], size_t nz )
{
	double( *faces )[3][3] = new double[nFaces][3][3];
	double( *maxCoords )[3] = new double[nFaces][3];
	double( *minCoords )[3] = new double[nFaces][3];
	int* facesIndex = new int[nFaces];
	double( *facesZ )[3][3] = new double[nFaces][3][3];
	double( *facesY )[3][3] = new double[nFaces][3][3];
	double( *minCoordsZ )[3] = new double[nFaces][3];
	double( *maxCoordsZ )[3] = new double[nFaces][3];
	double* crossings = new double[nFaces];

	int nxny = nx * ny;

	buildFaceMatrix( faces, vertices, faceIndices, nFaces );

	findExtremeCoords( faces, minCoords, maxCoords, nFaces );

	for ( int i = 0; i < nz; i++ )
	{
		int nFacesZ = findFacesInDim( facesIndex, minCoords, maxCoords, z[i], 2, nFaces );
		if ( nFacesZ == 0 )
			continue;
		selectFaces( facesZ, faces, facesIndex, nFacesZ );
		selectCoords( minCoordsZ, maxCoordsZ, minCoords, maxCoords, facesIndex, nFacesZ );
		for ( int j = 0; j < ny; j++ )
		{
			int nFacesY = findFacesInDim( facesIndex, minCoordsZ, maxCoordsZ, y[j], 1, nFacesZ );
			if ( nFacesY == 0 )
				continue;
			selectFaces( facesY, facesZ, facesIndex, nFacesY );
			int nCrossings = getCrossings( crossings, facesY, nFacesY, y[j], z[i] );
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
				while ( ( crossingsPassed < nCrossings ) && ( crossings[crossingsPassed] < x[k] ) )
				{
					crossingsPassed++;
					isInside = !isInside;
				}
				inside[i * nxny + k * ny + j] = isInside;
			}
		}
	}
	delete[] maxCoords;
	delete[] minCoords;
	delete[] facesIndex;
	delete[] facesZ;
	delete[] facesY;
	delete[] minCoordsZ;
	delete[] maxCoordsZ;
	delete[] crossings;
	delete[] faces;
}

