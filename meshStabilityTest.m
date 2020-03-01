
%--------------------------------------------------------------------------
%choose mesh
if 0 %square mesh into 8 right angled triangles
    [X,Y] = meshgrid([0:0.5:1],[0:0.5:1]);
    X = reshape(X,[],1);Y = reshape(Y,[],1);
    coordinates = [X,Y];
    DT = delaunayTriangulation(coordinates);
    elements = DT.ConnectivityList;
    dirichlet = DT.freeBoundary;
    neumann = [];
end

if 0 %acute randomized square
    coordinates = [	0, 0; 0.4, 0; 0.75,	0; 1, 0; 0.78, 0.22; 0.65, 0.2; 0.35, 0.2; 1, 0.36; 0.65, 0.5;
        0.3, 0.55; 0, 0.4; 0, 0.6; 0.2, 0.8; 0.4, 0.8; 0.75, 0.8; 1, 0.7; 1, 1; 0.6, 1; 0.25, 1; 0, 1];
    elements = [20 19 13; 13 19 14; 14 19 18; 14 18 15; 15 18 17; 15 17 16; 12 20 13; 12 13 10; 10 13 14; 10 14 9; 9 14 15; 
               9 15 16; 9 16 8; 9 8 5; 5 8 4; 5 4 3; 3 6 5; 2 6 3; 6 9 5; 7 9 6; 7 6 2; 7 10 9; 7 2 1; 1 11 7; 7 11 10; 11 12 10]; 
    dirichlet = [1 11; 11 12; 12 20; 20 19; 19 18; 18 17; 17 16; 16 8; 8 4; 4 3; 3 2; 2 1];
    neumann = [];
end

if 0 %copy of mesh from TW496 
    coordinates = [ 0, 0; 1, 0; 0.75, 0.25; 0.25, 0.5; 0, 0.5; 1, 0.5; 1, 1; 0.75, 1; 0, 1];
    elements = [1 3 2; 1 4 3; 1 5 4; 5 9 4; 4 9 8; 4 8 3; 8 7 6; 3 8 6; 3 6 2];
    dirichlet = [1 5; 5 9; 9 8; 8 7; 7 6; 6 2; 2 1];
    neumann = [];
end

if 0 %two very acute triangles on the main diagonal
    coordinates = [ 0, 0; 1, 0; 0.51, 0.49; 1, 1; 0.49, 0.51; 0, 1];
    elements = [1 3 2; 2 3 4; 3 5 4; 1 5 3; 1 6 5; 5 6 4];
    dirichlet = [1 6; 6 4; 4 2; 2 1];
    neumann = [];    
end

if 0 %slim obtuse triangles on all 4 borders
    coordinates = [ 0, 0; 1, 0; 0.5, 0.01; 0.99, 0.5; 1, 1; 0.5, 0.99; 0.5, 0.5; 0.01, 0.5; 0, 1];
    elements = [1 3 2; 2 3 4; 2 4 5; 5 4 6; 4 7 6; 4 3 7; 3 8 7; 3 1 8; 1 9 8; 8 9 6; 9 5 6; 8 6 7];
    dirichlet = [1 9; 9 5; 5 2; 2 1];
    neumann = [];     
end

if 1 %4 pointed star
    coordinates = [ 0,0; 0.5,0; 1,0; 1,0.5; 0.51,0.51; 0.51,0.49; 0.49,0.49; 0.5,0.5; 0.49,0.51; 0,0.5; 0,1; 0.5,1; 1,1];
    elements = [7 1 10; 1 7 2; 7 6 2; 2 6 3; 3 6 4; 4 6 5; 4 5 13; 5 12 13; 5 9 12; 9 11 12; 9 10 11; 9 7 10; 7 9 8; 9 5 8; 6 8 5; 6 7 8];
    dirichlet = [1 10; 10 11; 11 12; 12 13; 13 4; 4 3; 3 2; 2 1];
    neumann = [];
end


%--------------------------------------------------------------------------
%plot initial mesh
figure
triplot(elements, coordinates(:,1), coordinates(:,2),'linewidth',1)
pbaspect([1,1,1])
axis off



%--------------------------------------------------------------------------
%refine the mesh 
%
%I recommend adding break points so the computer isn't overwhelmed.Each 
%iteration increases triangles by factor of 3, so i=9 will usually be 
%around 100k triangles depending on the initial mesh. Increase I at your 
%own peril!
triangleType = 0;
for i = 1:9
    marked = (1:size(elements,1))';
    
    
    [coordinates,elements,dirichlet,neumann, triangleType] = refineCOM4(coordinates,elements,dirichlet,neumann,marked, triangleType);
    %[coordinates,elements,dirichlet,neumann] = refineRGB(coordinates,elements,dirichletelements,neumann,marked);
    figure
    triplot(elements, coordinates(:,1), coordinates(:,2),'linewidth',1)
    pbaspect([1,1,1])
    axis off
end


