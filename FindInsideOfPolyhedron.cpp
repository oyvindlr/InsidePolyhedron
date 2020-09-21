#include "pch.h"
#include "FindInsideOfPolyhedron.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define printf mexPrintf
#define warning(msg) mexWarnMsgIdAndTxt("InsidePolyhedron:LogicError", msg)
#else
#define warning(msg) printf(msg)
#endif

using namespace std;

typedef vector<array<double, 3>> nBy3Array;
typedef vector<array<array<double, 3>, 3>> nBy3By3Array;


static int dim0, dim1, dim2;

static const char *singularWarning = "The plane defined by one of the triangle faces is along the line used in ray tracing. Try adding random noise to your vertex coordinates to avoid this problem.";


static void findExtremeCoords(const nBy3By3Array &faces, nBy3Array &minCoords, nBy3Array &maxCoords, size_t nFaces)
{
	minCoords.resize(nFaces);
	maxCoords.resize(nFaces);
	for (int i = 0; i < nFaces; i++)
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

static int findFacesInDim(int *facesIndex, const nBy3Array &minCoords, const nBy3Array &maxCoords, double value, int dim, size_t nFaces)
{
	int count = 0;
	for (int i = 0; i < nFaces; i++)
	{
		if (minCoords[i][dim] < value && maxCoords[i][dim] > value)
		{
			facesIndex[count] = i;
			count++;
		}
	}
	return count;
}

static void selectFaces(nBy3By3Array& selectedFaces, const nBy3By3Array& faces, const int facesIndex[], int nFacesZ)
{
	selectedFaces.clear();
	for (int i = 0; i < nFacesZ; i++)
		selectedFaces.push_back(faces[facesIndex[i]]);
				
}

static void selectCoords(nBy3Array& minCoordsDim, nBy3Array& maxCoordsDim, const nBy3Array& minCoords, 
	const nBy3Array& maxCoords, const int facesIndex[], size_t nFacesZ)
{
	minCoordsDim.clear();
	maxCoordsDim.clear();
	for (int i = 0; i < nFacesZ; i++)
	{
		minCoordsDim.push_back(minCoords[facesIndex[i]]);
		maxCoordsDim.push_back(maxCoords[facesIndex[i]]);
	}
}

static void solve2by2otherway(double A[2][2], double b[2])
{
	double fac = A[0][0] / A[1][0];
	double a12 = A[0][1] - A[1][1] * fac;

	if (abs(a12) < 1e-10 || abs(A[1][0]) < 1e-10)
		warning(singularWarning);
	double b1 = (b[0] - b[1] * fac) / a12;
	b[0] = (b[1] - A[1][1] * b1) / A[1][0];
	b[1] = b1;
}


//Solve a 2 by 2 linear equation by Gaussian elimination with partial pivoting
static void solve2by2(double A[2][2], double b[2])
{
	if (abs(A[0][0]) < abs(A[1][0]))//Partial pivoting
	{
		solve2by2otherway(A, b);
		return;
	}
	double fac = A[1][0] / A[0][0];
	double a22 = A[1][1] - A[0][1] * fac;
	if (abs(a22) < 1e-10 || abs(A[0][0]) < 1e-10)
		warning(singularWarning);
	b[1] = (b[1] - b[0] * fac) / a22;
	b[0] = (b[0] - A[0][1] * b[1]) / A[0][0];
}

/** Get all crossings of triangular faces by a line in a specific direction */
static void getCrossings(vector<double> &crossings, const nBy3By3Array& faces, int nFaces, double dim1Value, double dim2Value)
{
	double b[2];
	double A[2][2];
	crossings.clear();

	for (int i = 0; i < nFaces; i++)
	{
		b[0] = dim1Value - faces[i][0][dim1];
		b[1] = dim2Value - faces[i][0][dim2];
		A[0][0] = faces[i][1][dim1] - faces[i][0][dim1];
		A[1][0] = faces[i][1][dim2] - faces[i][0][dim2];
		A[0][1] = faces[i][2][dim1] - faces[i][0][dim1];
		A[1][1] = faces[i][2][dim2] - faces[i][0][dim2];
		solve2by2(A, b);
		if (b[0] > 0 && b[1] > 0 && ((b[0] + b[1]) < 1))
		{
			//No pun intended
			double crossing = faces[i][0][dim0] + b[0] * (faces[i][1][dim0] - faces[i][0][dim0]) + b[1] * (faces[i][2][dim0] - faces[i][0][dim0]);
			crossings.push_back(crossing);
		}
	}
	sort(crossings.begin(), crossings.end());
}


static inline bool isOdd(int n)
{
	return (n % 2) == 1;
}

static void buildFaceMatrix(double faces[][3][3], const double vertices[][3], const int faceIndices[][3], size_t nFaces)
{
	for (int i = 0; i < nFaces; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				faces[i][j][k] = vertices[faceIndices[i][j]][k];
}

static void selectDimensionsForFastestProcessing(size_t *nx, size_t *ny, size_t *nz, size_t dimStep[3])
{
	size_t sizes[] = {*nx, *ny, *nz};
	int dims[3] = {0, 1, 2};
	size_t dimMult[] = {*ny * *nx, 1, *ny};

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2 - i; j++)
		{
			if (sizes[j] > sizes[j + 1])
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


/**
\brief Check whether a set of points on a 3D-grid is inside or outside a surface defined by a polyhedron.
This function uses ray-tracing to determine whether or not a point is inside the surface. Since the points
to be checked are aligned on a grid, we can reuse information for each point to perform the calculation 
significantly faster than if we were to check each point individually.
\param inside[out] Boolean array of output values, must be large enough to contain nx*ny*nz values. The result corresponding to the coordinate (x[i], y[j], z[k]) is found
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
	double (*faces)[3][3] = new double[nFaces][3][3];
	buildFaceMatrix(faces, vertices, faceIndices, nFaces);
	insidePolyhedron(inside, faces, nFaces, x, nx, y, ny, z, nz);
	delete[] faces;
}

/**
\brief Check whether a set of points on a 3D-grid is inside or outside a surface defined by a polyhedron.
This function uses ray-tracing to determine whether or not a point is inside the surface. Since the points
to be checked are aligned on a grid, we can reuse information for each point to perform the calculation
significantly faster than if we were to check each point individually.
This function is exactly the same as the other insidePolyhedron function, except that the surface is defined in a single list of triangular faces instead of a separate
list of vertices and faces.
\param inside[out] Boolean array of output values, must be large enough to contain nx*ny*nz values. The result corresponding to the coordinate (x[i], y[j], z[k]) is found
in inside[(j * nx * ny) + (i * ny) +  k]. The reason for this configuration is to align with Matlab's meshgrid(x, y, z) function.
\param vertices[in] Array of vertices in the polyhedron. Each vertex consists of 3 coordinates, x, y and z, therefore this is an n x 3 array.
\param faces[in] Definition of the triangular faces of the surface. faces[i][j][k] represents the k-coordinate (where x=0, y=1, z= 2) of the j'th vertex of the i'th face 
of the polyhedron.
\param nFaces Number of faces in the surface (size of faces)
\param x X-coordinate values on the grid to be checked.
\param nx Number of X-coordinates
\param y Y-coordinate values on the grid to be checked.
\param ny Number of Y-coordinates
\param z Z-coordinate values on the grid to be checked.
\param nz Number of Z-coordinates
*/
void insidePolyhedron(bool inside[], const double faces_[][3][3], size_t nFaces, const double x[], size_t nx, const double y[], size_t ny, const double z[], size_t nz)
{
	nBy3By3Array faces;
	nBy3Array minCoords;
	nBy3Array maxCoords;
	int *facesIndex = new int[nFaces];
	vector<double> crossings;
	nBy3Array minCoordsD2;
	nBy3Array maxCoordsD2;
	nBy3By3Array facesD2;
	nBy3By3Array facesD1;
	size_t dimSteps[3];

	faces.resize(nFaces);
	for (int i = 0; i < nFaces; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				faces[i][j][k] = faces_[i][j][k];

	selectDimensionsForFastestProcessing(&nx, &ny, &nz, dimSteps);

	size_t nxny = nx * ny;


	findExtremeCoords(faces, minCoords, maxCoords, nFaces);

	const double *gridCoords[] = {x, y, z};

	for (int i = 0; i < nz; i++)
	{
		int nFacesD2 = findFacesInDim(facesIndex, minCoords, maxCoords, gridCoords[dim2][i], dim2, nFaces);
		if (nFacesD2 == 0)
			continue;
		selectFaces(facesD2, faces, facesIndex, nFacesD2);
		selectCoords(minCoordsD2, maxCoordsD2, minCoords, maxCoords, facesIndex, nFacesD2);
		for (int j = 0; j < ny; j++)
		{
			int nFacesD1 = findFacesInDim(facesIndex, minCoordsD2, maxCoordsD2, gridCoords[dim1][j], dim1, nFacesD2);
			if (nFacesD1 == 0)
				continue;
			selectFaces(facesD1, facesD2, facesIndex, nFacesD1);
			getCrossings(crossings, facesD1, nFacesD1, gridCoords[dim1][j], gridCoords[dim2][i]);
			size_t nCrossings = crossings.size();
			if (nCrossings == 0)
				continue;
			if (isOdd(nCrossings))
				warning("Odd number of crossings found. The polyhedron may not be closed, or one of the triangular faces may lie in the exact direction of the traced ray.");
			bool isInside = false;
			int crossingsPassed = 0;
			for (int k = 0; k < nx; k++)
			{
				while ((crossingsPassed < nCrossings) && (crossings[crossingsPassed] < gridCoords[dim0][k]))
				{
					crossingsPassed++;
					isInside = !isInside;
				}
				inside[i * dimSteps[0] + j * dimSteps[1] + k * dimSteps[2]] = isInside;
			}
		}
	}
	delete[] facesIndex;
}


