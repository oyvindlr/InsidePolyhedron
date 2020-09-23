%InsidePolyhedron  
%   Check which of a set of 3D-points on a grid are inside 
%   and which are outside of one or more closed surfaces defined by a 
%   polyhedron. Implemented as a mex function for additional speedup.
% 
%   I = InsidePolyhedron(V, F, x, y, z) returns a 
%   3-dimensional boolean array where I(i, j, k) is true if and only 
%   if the point (y(i), x(j), z(k)) is inside the polyhedron surface. The
%   dimensions of I are (length(y), length(x), length(z)), which corresponds
%   to the size of the arrays returned by meshgrid(x, y, z). 
%   
%   I = InsidePolyhedron(V, F, x, y, z, useDithering) adds a tiny bit of
%   noise to the vertex coordinates if useDithering is True. This can 
%   usually solve issues where one or more of the triangular faces are
%   parallel to the traced ray, which in other circumstances would make the
%   ray-tracing algorithm fail. This problem may occur when the vertices 
%   of the surface are themselves aligned on a grid. Points that lie 
%   exactly on (or numerically indistingushable from) the surface will 
%   end up being randomly assigned as inside or outside the surface.
% 
%   V is an n by 3 matrix of vertices on the surface. Each row consists of 
%   the x, y and z coordinates of a vertex.
%   F is an m by 3 matrix of indices into the rows of V. Each row of F d
%   defines a triangular face defined by three vertices.
%   This is the same format as returned by surf = isosurface(...), which 
%   returns a struct with the fields surf.vertices and surf.faces.
%
%   The arguments x, y, and z are vectors of x, y and z coordinates
%   respectively. This function exploits the fact that the points to
%   be checked are on a grid in order to speed up processing significantly.
%   The points checked by InsidePolyhedron are those that are on the grid 
%   defined by [X, Y, Z] = meshgrid(x, y, z). Note, however, that the input
%   arguments are the vectors x, y and z, not the 3-dimensional arrays 
%   X, Y and Z.
%
%   This function relies on ray-tracing. For every point, a ray extending 
%   from this point in any direction will cross an odd number of faces if
%   the point is inside the surface, and an even number of faces if it is
%   outside. This technique can work even if the structure defines multiple
%   surfaces.
%   
%   The key to the speedup of checking points on a grid compared to
%   checking points individually comes from reducing the number of faces
%   we need to check for each point. For every surface [x, y, z0], we need
%   check only those faces that cross z0, and for every line [x, y0, z0]
%   in that plane we need only to check those faces that also cross y0.
%   Finally, we only need to check where each of these faces crosses the
%   line once in order to check each point on the line. 
%   
%   Example:
%
%   dx = 0.1;
%   x = -7:dx:7;y = -7:dx:7; z = -7:dx:7;
%   [X, Y, Z] = meshgrid(x, y, z);
%
%   %Create a logical array where points inside a radius of 5 are assigned 
%   %true.
%   I0 = sqrt(X.^2 + Y.^2 + Z.^2) < 5;
%   %Turn it into a smooth function for use in isosurface
%   I0smooth = smooth3(I0, 'box', 5);
%   %Find the isosurface as a polyhedron
%   surface = isosurface(X, Y, Z, I0smooth, 0.5);
%
%   %Use InsidePolyhedron to find points inside and outside the surface.
%   %Since the grid we check is the same that we used to generate the
%   %surface, we will need to use dithering to get a good result.
%   I1 = InsidePolyhedron(surface.vertices, surface.faces, x, y, z, true);
%
%   %Show a slice of the result
%   imagesc(squeeze(I1(:,50,:)));
%   axis equal;
