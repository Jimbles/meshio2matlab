%% write test gmsh file

% create small tetra triangulation
rng(1337);
x = rand(20,1);
y = rand(20,1);
z= rand(20,1);
dt = delaunayTriangulation(x,y,z);

% cell data 
celldataex=1:size(dt.ConnectivityList,1);
celldata={celldataex,celldataex+1,celldataex+2};
celldataname={'blah','blaah','blaaah'};

%pointdata
pointdataex=1:size(dt.Points,1);
pointdata={pointdataex,pointdataex+1,pointdataex+2};
pointdataname={'test','teest','tesest'};

%write with no data
meshio.write('example1.msh',dt.Points,dt.ConnectivityList);
A=meshio.read('example1.msh');

assert(isequal(A.vtx,dt.Points));


%write with cell data
meshio.write('examplecell1.msh',dt.Points,dt.ConnectivityList,celldata(1),celldataname(1));
B=meshio.read('examplecell1.msh');

assert(isequal(celldata(1),B.cell_data));

meshio.write('examplecell2.msh',dt.Points,dt.ConnectivityList,celldata,celldataname);
C=meshio.read('examplecell2.msh');

assert(isequal(celldata,C.cell_data));


% write with point data
meshio.write('examplepoint1.msh',dt.Points,dt.ConnectivityList,[],[],pointdata(1),pointdataname(1));
D=meshio.read('examplepoint1.msh');

assert(isequal(pointdata(1),D.point_data));

meshio.write('examplepoint2.msh',dt.Points,dt.ConnectivityList,[],[],pointdata,pointdataname);
E=meshio.read('examplepoint2.msh');

assert(isequal(pointdata,E.point_data));

%write with both point and cell data
meshio.write('examplepointcell.msh',dt.Points,dt.ConnectivityList,celldata,celldataname,pointdata,pointdataname);
F=meshio.read('examplepointcell.msh');

assert(isequal(pointdata,F.point_data) && isequal(celldata,F.cell_data));

%convert file using structwrite
meshio.structwrite('examplecell.vtu',C);
meshio.structwrite('examplepoint.vtu',E);
meshio.structwrite('examplepointcell.vtu',F);

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
meshio.write(fnameoutcell,R.vtx,R.cells.tri,{datacell},{'Celldata'});
meshio.write(fnameoutpoint,R.vtx,R.cells.tri,[],[],{datapoint},{'Pointdata'});











