#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include "FindInsideOfPolyhedron.h"

namespace py = pybind11;

py::array_t<bool> wrap_insidePolyhedron(
    py::array_t<double> vertices,
    py::array_t<int> faceIndices,
    py::array_t<double> x,
    py::array_t<double> y,
    py::array_t<double> z) {
    // Validate and get input buffers
    auto vertices_info = vertices.request();
    auto faceIndices_info = faceIndices.request();
    auto x_info = x.request();
    auto y_info = y.request();
    auto z_info = z.request();

    size_t nVertices = vertices_info.shape[0];
    size_t nFaces = faceIndices_info.shape[0];
    size_t nx = x_info.shape[0];
    size_t ny = y_info.shape[0];
    size_t nz = z_info.shape[0];

    if (vertices_info.ndim != 2 || vertices_info.shape[1] != 3) {
        throw std::runtime_error("vertices must be an n x 3 array");
    }
    if (faceIndices_info.ndim != 2 || faceIndices_info.shape[1] != 3) {
        throw std::runtime_error("faceIndices must be an n x 3 array");
    }
    if (x_info.ndim != 1 || y_info.ndim != 1 || z_info.ndim != 1) {
        throw std::runtime_error("x, y, z must be 1-dimensional arrays");
    }

    // Prepare output array
    size_t total_points = nx * ny * nz;
    py::array_t<bool> inside(total_points);
    auto inside_info = inside.request();

    // Call the original function
    insidePolyhedron(static_cast<bool*>(inside_info.ptr),
                     reinterpret_cast<const double(*)[3]>(vertices_info.ptr),
                     reinterpret_cast<const int(*)[3]>(faceIndices_info.ptr),
                     nFaces,
                     static_cast<const double*>(x_info.ptr), nx,
                     static_cast<const double*>(y_info.ptr), ny,
                     static_cast<const double*>(z_info.ptr), nz);

    // Reshape output to (nx, ny, nz) for convenience
    inside.resize({nx, ny, nz});
    return inside;
}

PYBIND11_MODULE(polyhedrontools, m) {
    m.def("inside_polyhedron", &wrap_insidePolyhedron, "Check if points lie inside a polyhedron");
}
