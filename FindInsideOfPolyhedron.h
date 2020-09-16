#pragma once

void insidePolyhedron( bool inside[], const double vertices[][3], const int faceIndices[][3], size_t nFaces, const double x[], size_t nx, const double y[], size_t ny, const double z[], size_t nz );
