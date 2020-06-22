%% write test gmsh file
% write to gmsh file with cell data
x = rand(50,1);
y = rand(50,1);
z= rand(50,1);
dt = delaunayTriangulation(x,y,z);
dataex=1:size(dt.ConnectivityList,1);
meshio.write('example.msh',dt.Points,dt.ConnectivityList,dataex);

%% read gmsh file

% example gmsh file is conversion from .stl file 
% this contains vertcies, lines, triangles and tetrahedra which are all read
% order is always vert/line/tri/tetra
fname=('finger_3D.msh');

P=meshio.read(fname);
meshio.plot(P);

% extract just the tetrahedra
Ptet=P;
Ptet.cells=Ptet.cells(end);
meshio.plot(Ptet);
%% read vtu and write new data
fname=('NNexample.vtu');
fnameoutcell='NNexampleNewDataC.vtu';
fnameoutpoint='NNexampleNewDataP.vtu';

% read mesh with subdomain as cell data (2 = skull, 1 = scalp/brain)
R=meshio.read(fname);
meshio.plot(R);

% create example test data
datacell=1:size(R.cells.tri,1);
datapoint=1:size(R.vtx,1);

% write to vtu - default is binary compressed
meshio.write(fnameoutcell,R.vtx,R.cells.tri,datacell,'Celldata');
meshio.write(fnameoutpoint,R.vtx,R.cells.tri,datapoint,'Pointdata');











