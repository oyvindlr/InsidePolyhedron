%Compile
mex -R2018a FindInsideOfPolyhedron.cpp insidepoly_mexfunction.cpp -output InsidePolyhedron
%%
%Run tests

%Resolution
dx = 0.1;
dy = 0.12;
dz = 0.08;
x = -7:dx:7.3;y = -7:dy:7; z = -7:dz:7.0;
%Create a logical array I0 of points inside and outside a structure
%First, create a sphere of radius 5
[X, Y, Z] = meshgrid(x, y, z);
I0 = sqrt(X.^2 + Y.^2 + Z.^2) < 5;
%Create a hollow cylinder inside the sphere
I0(sqrt(X.^2 + Z.^2) < 2.5) = false;
%Remove one half of the sphere
I0(y>0,:,:) = false;



%Find the isosurface as a polyhedron
surface = isosurface(X, Y, Z, I0, 0.5);


xyz = {x, y,z};

%Find all possible permutations of the dimensions
P = perms([1 2 3]);

for i = 1:size(P, 2) %Try all permutations of the three axes
    
    %Pick a permutation
    perm = P(i, :);
    x0 = xyz{perm(1)};
    y0 = xyz{perm(2)};
    z0 = xyz{perm(3)};
    
    %Use InsidePolyhedron to find points inside and outside the surface.
    %Since the grid we check is the same that we used to generate the
    %surface, we will need to use dithering to get a good result.
    I1 = InsidePolyhedron(surface.vertices(:, perm), surface.faces,x0, y0, z0, true);
    
    %Check for  correct dimensions of output
    assert(size(I1, 1) == length(y0));
    assert(size(I1, 2) == length(x0));
    assert(size(I1, 3) == length(z0));
    
    
    %Fix for a weirdness of meshgrid, (and inherited by InsidePolyhedron),
    %that y-dim is first and x-dim is second
    I0a = permute(I0, [2 1 3]);
    I1a = permute(I1, [2 1 3]);
    
    %Check that I1 reproduced I0
    assert(all(I0a == ipermute(I1a, perm), 'all'), 'Logical array not perfectly reproduced');
end
disp('All logical arrays perfectly reproduced');
%%
%Create a figure showing the result


p = patch(surface, 'EdgeColor', 'black', 'FaceColor', [0.7 0.7 0.7]);
p.FaceAlpha = 0.4;
%camlight
camlight(-60,-10)
camlight(-30,-10)
lighting gouraud
view(3)
axis off
hold on
reducepatch(p, 0.01);

x2i = 11:10:length(x)-11;
x2 = x(x2i);
y2i=30:15:61;
y2= y(y2i);
z2i = 35:15:141;
z2 = z(z2i);
I3 = InsidePolyhedron(surface.vertices, surface.faces, x2, y2, z2, true);

hold on
for i = 1:length(y2)
    for j = 1:length(z2)
        plot3(x2, y2(i)*ones(length(x2), 1), z2(j)*ones(length(x2), 1), 'k*-', 'MarkerEdgeColor', 'r');
        ins = x2(squeeze(I3(i, :,j)));
        plot3(ins, y2(i)*ones(length(ins), 1), z2(j)*ones(length(ins), 1), 'g*');
        plot3(ins, y2(i)*ones(length(ins), 1), z2(j)*ones(length(ins), 1), 'go');
    end
end
view(3);
axis equal;
grid on
xticklabels([]);
yticklabels([]);
zticklabels([]);