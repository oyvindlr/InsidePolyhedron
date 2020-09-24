/** \file FindInsideOfPolyhedron.cpp
*	\brief Functions to determine whether points on a grid lie inside or outside a surface in 3D-space defined by triangular faces.
*	\author Oyvind L Rortveit
*	\date 2020
*/

#include "pch.h"
#include "FindInsideOfPolyhedron.h"
#include <cmath>
#include <algorithm>
#include "DynamicArray.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define warning(msg) mexWarnMsgIdAndTxt("InsidePolyhedron:LogicError", msg)
#else
#define warning(msg) printf(msg)
#endif

//Remember to change this if you ever rewrite to use single precision 
#define EPS 1e-14

using namespace std;

typedef vector<array<double, 3>> nBy3Array;


static const int NWARNINGS = 2;

static bool hasWarned[NWARNINGS];

static void resetWarnings()
{
	for (int i = 0; i < NWARNINGS; i++)
		hasWarned[i] = false;
}

static void warnOnce(const char *msg, int warningNo)
{
	if (!hasWarned[warningNo])
	{
		warning(msg);
		hasWarned[warningNo] = true;
	}
}


static void findBoundingBox(const nBy3By3Array& faces, nBy3Array &minCoords, nBy3Array &maxCoords)
{
	minCoords.resize(faces.size());
	maxCoords.resize(faces.size());
	for (int i = 0; i < faces.size(); i++)
	{
		for (int dim = 0; dim < 3; dim++)
		{
			minCoords[i][dim] = faces[i][0][dim];
			if (faces[i][1][dim] < minCoords[i][dim])
				minCoords[i][dim] = faces[i][1][dim];
			if (faces[i][2][dim] < minCoords[i][dim])
				minCoords[i][dim] = faces[i][2][dim];
			maxCoords[i][dim] = faces[i][0][dim];
			if (faces[i][1][dim] > maxCoords[i][dim])
				maxCoords[i][dim] = faces[i][1][dim];
			if (faces[i][2][dim] > maxCoords[i][dim])
				maxCoords[i][dim] = faces[i][2][dim];
		}
	}
}

static void findFacesInDim(DynamicArray<int> &facesIndex, const nBy3Array &minCoords, const nBy3Array &maxCoords, double value, int dim)
{
	facesIndex.clear();
	size_t nFaces = minCoords.size();
	for (int i = 0; i < nFaces; i++)
		if (minCoords[i][dim] <= value && maxCoords[i][dim] >= value)
			facesIndex.insertLast_unsafe(i); //Actually safe because the array has a capacity of nFaces elements
}

static void findFacesInDim(DynamicArray<int>& facesIndex, const nBy3Array &minCoords, const nBy3Array &maxCoords, double value, int dim, const DynamicArray<int>& facesIndexOld)
{
	facesIndex.clear();
	size_t nFaces = facesIndexOld.size();
	for (int i = 0; i < nFaces; i++)
	{
		int k = facesIndexOld[i];
		if (minCoords[k][dim] <= value && maxCoords[k][dim] >= value)
			facesIndex.insertLast_unsafe(k); //Actually safe because the array has a capacity of nFaces elements
	}
}


/** @brief  Invert a 2 by 2 matrix by Gaussian elimination with partial pivoting
 * @param[in] A Matrix to be inverted
 * @param[out] Ainv Inverted matrix
*/
static void invert2by2(const double A[2][2], double Ainv[2][2])
{

	if (abs(A[0][0]) > abs(A[1][0]))//Partial pivoting
	{
		double fac = A[1][0] / A[0][0];
		double a22i = 1.0/ (A[1][1] - A[0][1] * fac);
		Ainv[1][0] = -fac * a22i;
		Ainv[0][0] = (1.0 - A[0][1] * Ainv[1][0]) / A[0][0];
		Ainv[1][1] = a22i;
		Ainv[0][1] = - A[0][1] * a22i / A[0][0];
	}
	else
	{
		double fac = A[0][0] / A[1][0];
		double a12i = 1.0 / (A[0][1] - A[1][1] * fac);
		Ainv[0][0] = -A[1][1] * a12i / A[1][0];
		Ainv[1][0] = a12i;
		Ainv[1][1] = -fac * a12i;
		Ainv[0][1] = (1.0 - A[1][1] * Ainv[1][1]) / A[1][0];
	}
}


static void multiply2by2(const double A[2][2], double b[2])
{
	double b0 = A[0][0]*b[0] + A[0][1]*b[1];
	b[1] = A[1][0] * b[0] + A[1][1] * b[1];
	b[0] = b0;
}


/** Get all crossings of triangular faces by a line in a specific direction */
static void getCrossings(vector<double> &crossings, const nBy3By3Array& faces, const double coords[2], const int dimOrder[3], 
	const DynamicArray<int> &facesIndex, double faceMatrices[][2][2])
{
	double b[2];
	double A[2][2];
	double Ainv[2][2];
	int dim2 = dimOrder[2];
	size_t nFaces = facesIndex.size();
	crossings.clear();

	for (int k = 0; k < nFaces; k++)
	{
		int i = facesIndex[k];
		for (int dimNo = 0; dimNo < 2; dimNo++)
		{
			int dim = dimOrder[dimNo];
			b[dimNo] = coords[dimNo] - faces[i][0][dim];
		}

		multiply2by2(faceMatrices[i], b);

		if (b[0] > 0 && b[1] > 0 && ((b[0] + b[1]) < 1))
		{
			//No pun intended
			double crossing = faces[i][0][dim2] + b[0] * (faces[i][1][dim2] - faces[i][0][dim2]) + b[1] * (faces[i][2][dim2] - faces[i][0][dim2]);
			crossings.push_back(crossing);
		}
	}
	sort(crossings.begin(), crossings.end());
}

static void computeFaceMatrices(const nBy3By3Array &faces, double matrices[][2][2], int dimOrder[3])
{
	double A[2][2];
	for (int i = 0; i < faces.size(); i++)
	{
		for (int dimNo = 0; dimNo < 2; dimNo++)
		{
			int dim = dimOrder[dimNo];
			A[dimNo][0] = faces[i][1][dim] - faces[i][0][dim];
			A[dimNo][1] = faces[i][2][dim] - faces[i][0][dim];
		}
		invert2by2(A, matrices[i]);
	}
}

static inline bool isOdd(size_t n)
{
	return (n % 2) == 1;
}

static void buildFaceArray(nBy3By3Array& faces, const double vertices[][3], const int faceIndices[][3], size_t nFaces)
{
	faces.resize(nFaces);
	for (int i = 0; i < nFaces; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				faces[i][j][k] = vertices[faceIndices[i][j]][k];
}


//Selects the order in which the dimensions are processed
static void selectDimensionsForFastestProcessing(int dimOrder[3], size_t dimSize[3], size_t dimStep[3])
{
	size_t nx = dimSize[0];
	size_t ny = dimSize[1];

	//Bubble sort
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2 - i; j++)
		{
			if (dimSize[j] > dimSize[j + 1])
			{
				size_t temp = dimSize[j];
				dimSize[j] = dimSize[j + 1];
				dimSize[j + 1] = temp;
				int tempi = dimOrder[j];
				dimOrder[j] = dimOrder[j + 1];
				dimOrder[j + 1] = tempi;
			}
		}
	}

	dimStep[dimOrder[0]] = ny;
	dimStep[dimOrder[1]] = 1;
	dimStep[dimOrder[2]] = ny * nx;
}


/**
\brief Check whether a set of points on a 3D-grid is inside or outside a surface defined by a polyhedron.
This function uses ray-tracing to determine whether or not a point is inside the surface. Since the points
to be checked are aligned on a grid, we can reuse information for each point to perform the calculation 
significantly faster than if we were to check each point individually.
\param[out] inside Boolean array of output values, must be large enough to contain nx*ny*nz values. The result corresponding to the coordinate (x[i], y[j], z[k]) is found
in inside[(j * nx * ny) + (i * ny) +  k]. The reason for this configuration is to align with Matlab's meshgrid(x, y, z) function.
\param vertices[in] Array of vertices in the polyhedron. Each vertex consists of 3 coordinates, x, y and z, therefore this is an n x 3 array.
\param faceIndices[in] Definition of the triangular faces of the surface. Each row of this matrix consists of three indices into the vertex-list, which together define
a triangular face.
\param nFaces Number of faces in the surface (size of faceIndices)
\param x X-coordinate values on the grid to be checked. 
\param nx Number of X-coordinates
\param y Y-coordinate values on the grid to be checked.
\param ny Number of Y-coordinates
\param z Z-coordinate values on the grid to be checked.
\param nz Number of Z-coordinates
*/
void insidePolyhedron(bool inside[], const double vertices[][3], const int faceIndices[][3], size_t nFaces, const double x[], size_t nx, const double y[], size_t ny, const double z[], size_t nz)
{
	nBy3By3Array faces;
	buildFaceArray(faces, vertices, faceIndices, nFaces);
	insidePolyhedron(inside, faces, x, nx, y, ny, z, nz);
}

/**
\brief Check whether a set of points on a 3D-grid is inside or outside a surface defined by a polyhedron.
This function uses ray-tracing to determine whether or not a point is inside the surface. Since the points
to be checked are aligned on a grid, we can reuse information for each point to perform the calculation
significantly faster than if we were to check each point individually.
This function is exactly the same as the other insidePolyhedron function, except that the surface is defined in a single list of triangular faces instead of a separate
list of vertices and faces.
\param[out] inside Boolean array of output values, must be large enough to contain nx*ny*nz values. The result corresponding to the coordinate (x[i], y[j], z[k]) is found
in inside[(j * nx * ny) + (i * ny) +  k]. The reason for this configuration is to align with Matlab's meshgrid(x, y, z) function.
\param vertices[in] Array of vertices in the polyhedron. Each vertex consists of 3 coordinates, x, y and z, therefore this is an n x 3 array.
\param faces Array (n by 3 by 3) of triangular faces, such that faces[i][j][k] represents the k-coordinate (where x=0, y=1, z= 2) of the j'th vertex of the i'th face
of the polyhedron.
\param x X-coordinate values on the grid to be checked.
\param nx Number of X-coordinates
\param y Y-coordinate values on the grid to be checked.
\param ny Number of Y-coordinates
\param z Z-coordinate values on the grid to be checked.
\param nz Number of Z-coordinates
*/
void insidePolyhedron(bool inside[], const nBy3By3Array& faces, const double x[], size_t nx, const double y[], size_t ny, const double z[], size_t nz)
{
	nBy3Array minCoords;
	nBy3Array maxCoords;	
	vector<double> crossings;
	size_t dimSteps[3];

	
	double(*faceMatrices)[2][2] = new double[faces.size()][2][2];

	DynamicArray<int> facesIndex{faces.size()};
	DynamicArray<int> facesIndex2{faces.size()};

	resetWarnings();


	size_t dimSize[3] = {nx, ny, nz};
	int dimOrder[3] = {0, 1, 2};
	
	selectDimensionsForFastestProcessing(dimOrder, dimSize, dimSteps);
	const int dim0 = dimOrder[0], dim1 = dimOrder[1], dim2 = dimOrder[2];

	computeFaceMatrices(faces, faceMatrices, dimOrder);
	

	findBoundingBox(faces, minCoords, maxCoords);

	const double *gridCoords[] = {x, y, z};

	for (int i = 0; i < dimSize[0]; i++)
	{
		findFacesInDim(facesIndex, minCoords, maxCoords, gridCoords[dim0][i], dim0);
		if (facesIndex.size() == 0)
			continue;
		for (int j = 0; j < dimSize[1]; j++)
		{
			const double coords[2] = {gridCoords[dim0][i], gridCoords[dim1][j]};
			findFacesInDim(facesIndex2, minCoords, maxCoords, coords[1], dim1, facesIndex);
			if (facesIndex2.size() == 0)
				continue;
			getCrossings(crossings, faces, coords, dimOrder, facesIndex2, faceMatrices);
			size_t nCrossings = crossings.size();
			if (nCrossings == 0)
				continue;
			if (isOdd(nCrossings))
				warnOnce("Odd number of crossings found. The polyhedron may not be closed, or one of the triangular faces may lie in the exact direction of the traced ray.", 0);
			bool isInside = false;
			int crossingsPassed = 0;
			for (int k = 0; k < dimSize[2]; k++)
			{
				while ((crossingsPassed < nCrossings) && (crossings[crossingsPassed] < gridCoords[dim2][k]))
				{
					crossingsPassed++;
					isInside = !isInside;
				}
				inside[i * dimSteps[0] + j * dimSteps[1] + k * dimSteps[2]] = isInside;
			}
		}
	}
	delete[] faceMatrices;
}


