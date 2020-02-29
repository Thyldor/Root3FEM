
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

if 1
    coordinates = [ 0, 0; 1, 0; 0.51, 0.49; 1, 1; 0.49, 0.51; 0, 1];
    elements = [1 3 2; 2 3 4; 3 5 4; 1 5 3; 1 6 5; 5 6 4];
    dirichlet = [1 6; 6 4; 4 2; 2 1];
    neumann = [];
    
end

if 0 %square mesh into 8 right angled triangles
    [X,Y] = meshgrid([0:0.5:1],[0:0.5:1]);
    X = reshape(X,[],1);Y = reshape(Y,[],1);
    coordinates = [X,Y];
    DT = delaunayTriangulation(coordinates);
    elements = DT.ConnectivityList;
    dirichlet = DT.freeBoundary;
    neumann = [];
    %triplot(DT)
end


%figure
%triplot(elements, coordinates(:,1), coordinates(:,2),'linewidth',1)
%pbaspect([1,1,1])
%axis off



%triangleType = zeros(size(elements,1),1);
triangleType = 0;

for i = 1:9
    marked = (1:size(elements,1))';
    
    
    [coordinates,elements,dirichlet,neumann, triangleType] = refineCOM4(coordinates,elements,dirichlet,neumann,marked, triangleType);
    %[coordinates,elements,dirichlet,neumann] = refineRGB(coordinates,elements,dirichlet,neumann,marked);
    %figure
    %triplot(elements, coordinates(:,1), coordinates(:,2),'linewidth',1)
    %pbaspect([1,1,1])
    %axis off
end

