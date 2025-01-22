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
#include <cstdio>
#define warning(msg) printf("%s\n",msg)
#endif

//Remember to change this if you ever rewrite to use single precision 
#define EPS 1e-14

using namespace std;

typedef vector<array<double, 3>> nBy3Array;


static const int NWARNINGS = 2;

static bool hasWarned[NWARNINGS];

static void resetWarnings()
{
	for (size_t i = 0; i < NWARNINGS; i++)
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


static void findExtremeCoords(const nBy3By3Array& faces, nBy3Array &minCoords, nBy3Array &maxCoords)
{
	minCoords.resize(faces.size());
	maxCoords.resize(faces.size());
	for (size_t i = 0; i < faces.size(); i++)
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
	for (size_t i = 0; i < nFaces; i++)
		if (minCoords[i][dim] <= value && maxCoords[i][dim] >= value)
			facesIndex.insertLast_unsafe(i); //Actually safe because the array has a capacity of nFaces elements
}


static void selectFaces(nBy3By3Array &selectedFaces, const nBy3By3Array& faces, const DynamicArray<int> &facesIndex)
{
	selectedFaces.clear();
	for (size_t i = 0; i < facesIndex.size(); i++)
		selectedFaces.push_back(faces[facesIndex[i]]);
}

static void selectCoords(nBy3Array& minCoordsDim, nBy3Array& maxCoordsDim, const nBy3Array& minCoords, 
	const nBy3Array& maxCoords, const DynamicArray<int>& facesIndex)
{
	minCoordsDim.clear();
	maxCoordsDim.clear();
	for (size_t i = 0; i < facesIndex.size(); i++)
	{
		minCoordsDim.push_back(minCoords[facesIndex[i]]);
		maxCoordsDim.push_back(maxCoords[facesIndex[i]]);
	}
}


/******
* Removed this part for now. It used to check if the matrix is singular, which means that the plane defined
* by the triangle is parallel to the ray. However, this can mean either that the ray doesn't cross the triangle
* (which is fine) or that the triangular face lies along the ray, in which case we cannot determine whether or
* not there is a crossing. Since the two situations are difficult to discern, I have chosen for now to ignore the 
* problem. In most cases, if there is a real problem, we will at some point end up with an odd number of crossings
* for a single ray, and therefore get a warning.
******/
#if 0
static const char *singularWarning = "The plane defined by one of the triangle faces is along the line used in ray tracing. "
	"Try adding random noise to your vertex coordinates to avoid this problem.";


/** Solve a 2 by 2 linear equation A*x = b by Gaussian elimination with partial pivoting.
*   b is overwritten by the result x.
*  \returns False if the matrix is (numerically) singular, true otherwise.
*/
static bool solve2by2(const double A[2][2], double b[2])
{
	if (abs(A[0][0]) > abs(A[1][0]))//Partial pivoting
	{
		double fac = A[1][0] / A[0][0];
		double a22 = A[1][1] - A[0][1] * fac;		
		b[1] = (b[1] - b[0] * fac) / a22;
		b[0] = (b[0] - A[0][1] * b[1]) / A[0][0];
		return abs(a22) > EPS && abs(A[0][0]) > EPS;
	}
	else
	{
		double fac = A[0][0] / A[1][0];
		double a12 = A[0][1] - A[1][1] * fac;
		double b1 = (b[0] - b[1] * fac) / a12;
		b[0] = (b[1] - A[1][1] * b1) / A[1][0];
		b[1] = b1;
		return abs(a12) > EPS && abs(A[1][0]) > EPS;
	}
}
#endif

/** Solve a 2 by 2 linear equation A*x = b by Gaussian elimination with partial pivoting.
*   b is overwritten by the result x.
*/
static void solve2by2(const double A[2][2], double b[2])
{
	if (abs(A[0][0]) > abs(A[1][0]))//Partial pivoting
	{
		double fac = A[1][0] / A[0][0];
		double a22 = A[1][1] - A[0][1] * fac;
		b[1] = (b[1] - b[0] * fac) / a22;
		b[0] = (b[0] - A[0][1] * b[1]) / A[0][0];
	}
	else
	{
		double fac = A[0][0] / A[1][0];
		double a12 = A[0][1] - A[1][1] * fac;
		double b1 = (b[0] - b[1] * fac) / a12;
		b[0] = (b[1] - A[1][1] * b1) / A[1][0];
		b[1] = b1;
	}
}

/** Get all crossings of triangular faces by a line in a specific direction */
static void getCrossings(vector<double> &crossings, const nBy3By3Array& faces, const double coords[2], const int dimOrder[3])
{
	double b[2];
	double A[2][2];
	int dim2 = dimOrder[2];
	size_t nFaces = faces.size();
	crossings.clear();

	for (size_t i = 0; i < nFaces; i++)
	{
		for (int dimNo = 0; dimNo < 2; dimNo++)
		{
			int dim = dimOrder[dimNo];
			b[dimNo] = coords[dimNo] - faces[i][0][dim];
			A[dimNo][0] = faces[i][1][dim] - faces[i][0][dim];
			A[dimNo][1] = faces[i][2][dim] - faces[i][0][dim];
		}
		solve2by2(A, b);
		if (b[0] > 0 && b[1] > 0 && ((b[0] + b[1]) < 1))
		{
			//No pun intended
			double crossing = faces[i][0][dim2] + b[0] * (faces[i][1][dim2] - faces[i][0][dim2]) + b[1] * (faces[i][2][dim2] - faces[i][0][dim2]);
			crossings.push_back(crossing);
		}
	}
	sort(crossings.begin(), crossings.end());
}


static inline bool isOdd(size_t n)
{
	return (n % 2) == 1;
}

static void buildFaceMatrix(nBy3By3Array& faces, const double vertices[][3], const int faceIndices[][3], size_t nFaces)
{
	faces.resize(nFaces);
	for (size_t i = 0; i < nFaces; i++)
		for (size_t j = 0; j < 3; j++)
			for (size_t k = 0; k < 3; k++)
				faces[i][j][k] = vertices[faceIndices[i][j]][k];
}


//Selects the order in which the dimensions are processed
static void selectDimensionsForFastestProcessing(int dimOrder[3], size_t dimSize[3], size_t dimStep[3])
{
	size_t nx = dimSize[0];
	size_t ny = dimSize[1];

	//Bubble sort
	for (size_t i = 0; i < 2; i++)
	{
		for (size_t j = 0; j < 2 - i; j++)
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

	size_t dimSteps_orig[3] = {ny, 1, ny*nx};
	dimStep[0] = dimSteps_orig[dimOrder[0]];
	dimStep[1] = dimSteps_orig[dimOrder[1]];
	dimStep[2] = dimSteps_orig[dimOrder[2]];
}


/**
\brief Check whether a set of points on a 3D-grid is inside or outside a surface defined by a polyhedron.
This function uses ray-tracing to determine whether or not a posize_t is inside the surface. Since the points
to be checked are aligned on a grid, we can reuse information for each point to perform the calculation 
significantly faster than if we were to check each posize_t individually.
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
	buildFaceMatrix(faces, vertices, faceIndices, nFaces);
	insidePolyhedron(inside, faces, x, nx, y, ny, z, nz);
}

/**
\brief Check whether a set of points on a 3D-grid is inside or outside a surface defined by a polyhedron.
This function uses ray-tracing to determine whether or not a posize_t is inside the surface. Since the points
to be checked are aligned on a grid, we can reuse information for each point to perform the calculation
significantly faster than if we were to check each posize_t individually.
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
	nBy3Array minCoordsD2;
	nBy3Array maxCoordsD2;
	nBy3By3Array facesD2;
	nBy3By3Array facesD1;
	size_t dimSteps[3];

	resetWarnings();

	DynamicArray<int> facesIndex {faces.size()};

	size_t dimSize[3] = {nx, ny, nz};
	int dimOrder[3] = {0, 1, 2};
	
	selectDimensionsForFastestProcessing(dimOrder, dimSize, dimSteps);
	const int dim0 = dimOrder[0], dim1 = dimOrder[1], dim2 = dimOrder[2];

	findExtremeCoords(faces, minCoords, maxCoords);

	const double *gridCoords[] = {x, y, z};

	for (size_t i = 0; i < dimSize[0]; i++)
	{
		findFacesInDim(facesIndex, minCoords, maxCoords, gridCoords[dim0][i], dim0);
		if (facesIndex.size() == 0)
			continue;
		selectFaces(facesD2, faces, facesIndex);
		selectCoords(minCoordsD2, maxCoordsD2, minCoords, maxCoords, facesIndex);
		for (size_t j = 0; j < dimSize[1]; j++)
		{
			const double coords[2] = {gridCoords[dim0][i], gridCoords[dim1][j]};
			findFacesInDim(facesIndex, minCoordsD2, maxCoordsD2, coords[1], dim1);
			if (facesIndex.size() == 0)
				continue;
			selectFaces(facesD1, facesD2, facesIndex);
			getCrossings(crossings, facesD1, coords, dimOrder);
			size_t nCrossings = crossings.size();
			if (nCrossings == 0)
				continue;
			if (isOdd(nCrossings))
				warnOnce("Odd number of crossings found. The polyhedron may not be closed, or one of the triangular faces may lie in the exact direction of the traced ray.", 0);
			bool isInside = false;
			size_t crossingsPassed = 0;
			for (size_t k = 0; k < dimSize[2]; k++)
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
}


