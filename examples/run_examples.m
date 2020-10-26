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
meshio.write('testout/example1.msh',dt.Points,dt.ConnectivityList);
A=meshio.read('testout/example1.msh');

assert(isequal(A.vtx,dt.Points));


%write with cell data
meshio.write('testout/examplecell1.msh',dt.Points,dt.ConnectivityList,celldata(1),celldataname(1));
B=meshio.read('testout/examplecell1.msh');

assert(isequal(celldata(1),B.cell_data));

meshio.write('testout/examplecell2.msh',dt.Points,dt.ConnectivityList,celldata,celldataname);
C=meshio.read('testout/examplecell2.msh');

assert(isequal(celldata,C.cell_data));


% write with point data
meshio.write('testout/examplepoint1.msh',dt.Points,dt.ConnectivityList,[],[],pointdata(1),pointdataname(1));
D=meshio.read('testout/examplepoint1.msh');

assert(isequal(pointdata(1),D.point_data));

meshio.write('testout/examplepoint2.msh',dt.Points,dt.ConnectivityList,[],[],pointdata,pointdataname);
E=meshio.read('testout/examplepoint2.msh');

assert(isequal(pointdata,E.point_data));

%write with both point and cell data
meshio.write('testout/examplepointcell.msh',dt.Points,dt.ConnectivityList,celldata,celldataname,pointdata,pointdataname);
F=meshio.read('testout/examplepointcell.msh');

assert(isequal(pointdata,F.point_data) && isequal(celldata,F.cell_data));

%convert file using structwrite
meshio.structwrite('testout/examplecell.vtu',C);
meshio.structwrite('testout/examplepoint.vtu',E);
meshio.structwrite('testout/examplepointcell.vtu',F);

%% read gmsh file

% example gmsh file is conversion from .stl file 
% this contains vertcies, lines, triangles and tetrahedra which are all read
% order is always vert/line/tri/tetra
fname=('finger_3D');

P=meshio.read([fname '.msh']);
meshio.plot(P);

% extract just the tetrahedra
Ptet=P;
Ptet.cells=Ptet.cells(end);
meshio.plot(Ptet);

% write file back to vtu N.B. meshio merges cells of the same type
meshio.structwrite([fname '.vtu'],P);
Pvtu=meshio.read([fname '.vtu']);
%Ptvu only has 4 cells, each for all vertex,line,tri,tetra

% meshio cannot write .msh files with different cell types in at the moment
% (meshio v 4.0.16). So write each type to a separate file and merge them
% in gmsh

try
meshio.structwrite('testout/error.msh',P);

catch
    fprintf(2,'\n -----------\n meshio can only write 1 cell type (e.g. tetra) to msh file \n -----------\n');
end


%% read vtu and write new data
fname='NNexample.vtu';
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











