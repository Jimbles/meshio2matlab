classdef meshio
    %MESHIO Matlab2Meshio functions
    %   Use the python meshio library to read and write unstructured mesh
    %   files .msh .vtu .vtk .xdmf tetgen ansys .stl .obj
    %
    %   more formats at:
    %   https://github.com/nschloe/meshio
    %
    % meshio.read        - read mesh file to struct
    % meshio.write       - write matlab mesh to file
    % meshio.plot        - plot all contents of a meshfile
    % meshio.structwrite - write meshio.read output struct to file
    %
    % utilities
    % meshio.np2mat - np array to matlab array
    % meshio.mat2nparray - matlab to np array
    
    
    properties
        
    end
    
    methods(Static)
        
        function arrayout = np2mat(nparrayin)
            %np2mat convert python nparray to matlab
            % matlab doesnt like nparray, so first convert to standard
            % py.array.array
            %
            % small modification of
            % https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double#comment_437274
            %
            
            % Add order='F' to get data in column-major order (as in
            % Fortran 'F' and Matlab)
            data_size = cellfun(@int64,cell(nparrayin.shape));
            if any(numel(data_size) == [0,1]) % 1D or scalars
                % This is a simple operation
                arrayout=double(py.array.array('d', py.numpy.nditer(nparrayin)));
            elseif length(data_size)==2
                % order='F' is used to get data in column-major order (as
                % in Fortran 'F' and Matlab)
                arrayout=reshape(double(py.array.array('d', ...
                    py.numpy.nditer(nparrayin, pyargs('order', 'F')))), ...
                    data_size);
            else
                % For multidimensional arrays more manipulation is
                % required First recover in python order (C contiguous
                % order)
                arrayout=double(py.array.array('d', ...
                    py.numpy.nditer(nparrayin, pyargs('order', 'C'))));
                % Switch the order of the dimensions (as Python views
                % this in the opposite order to Matlab) and reshape to
                % the corresponding C-like array
                arrayout=reshape(arrayout,fliplr(data_size));
                % Now transpose rows and columns of the 2D sub-arrays to
                % arrive at the correct Matlab structuring
                arrayout=permute(arrayout,[length(data_size):-1:1]);
            end
        end
        
        
        function objout = read(filename)
            %meshio.read read mesh file using python meshio
            %   Calls python meshio library and processes output to a
            %   matlab struct
            %   
            % Inputs:
            % filename - needs extension .msh .vtu .vtk etc.
            % 
            % Output:
            % objout.vtx - verticies
            % objout.Cells - structure array for each geometry saved in file
            %       .Cells.tri - Trigangulation connectivity list for this cell
            %       .Cells.type - 'Tetra','Triangle','Line','Vertex'
            % objout.cell_data - data for each element {cell array}
            % objout.cell_data_name - {cell array}
            % objout.point_data - data for each vertex {cell array}
            % objout.point_data_name - {cell array}
            %
            % Usage:
            % %read data 
            % M=meshio.read('example.msh');
            % %plot contents
            % meshio.plot(M);
            
            
            % ------ load meshes ------
            
            fprintf('Meshio reading meshfile : %s\n',filename);
            
            % read mesh into meshio python format
            pymesh=py.meshio.read(filename);
            
            % vtx are saved as 1 array
            vtx =meshio.np2mat(pymesh.points);
            
            objout.vtx=vtx;
            
            % files can have mix of data tetra/tri/line/point
            numCells=size(pymesh.cells,2);
            
            % for each cell data array save the type and data
            % many file types can have mix of tetra/tri/lines etc
            
            for iCell = 1 : numCells
                % find cell type
                C(iCell).type=char(pymesh.cells{iCell}.type);
                
                % convert cell data to matlab array
                tritmp=pymesh.cells{iCell}.data;
                tri =meshio.np2mat(tritmp);
                
                % python uses 0 indexing
                tri= tri+1;
                
                C(iCell).tri=tri;
            end
            
            % ------ load data ------
            
            % cannot convert from dict to struct directly as
            % python allows spaces in dict keys but matlab doesnt in
            % struct fields
            
            % cell data
            cell_data_names_py=py.list(pymesh.cell_data.keys);
            numCelldata=double(py.len(cell_data_names_py));
            
            if numCelldata > 0
                for iCelldata = 1:numCelldata
                    
                    % find name in dict
                    cell_data_name=char(cell_data_names_py{iCelldata});
                    
                    % get value for this key, data is stored in list so get
                    % the first (only?) array in this list
                    cell_data_py=pymesh.cell_data{cell_data_name}{1};
                    cell_data=meshio.np2mat(cell_data_py);
                    
                    objout.cell_data{iCelldata}=cell_data;
                    objout.cell_data_name{iCelldata}=cell_data_name;
                end
            else
                objout.cell_data=cell(1);
                objout.cell_data_name=cell(1);
            end
            
            % point data
            point_data_names_py=py.list(pymesh.point_data.keys);
            numPointdata=double(py.len(point_data_names_py));
            
            if numPointdata > 0
                for iPointdata = 1:numPointdata
                    
                    %find name in dict
                    point_data_name=char(point_data_names_py{iPointdata});
                    
                    %get value for this key, usually stored as an array,
                    %but sometimes stored as list so check for that
                    if length(pymesh.point_data{point_data_name}) > 1
                        point_data_py=pymesh.point_data{point_data_name}{1};
                    else
                        point_data_py=pymesh.point_data{point_data_name};
                    end
                    
                    point_data=meshio.np2mat(point_data_py);
                    
                    objout.point_data{iPointdata}=point_data;
                    objout.point_data_name{iPointdata}=point_data_name;
                end
            else
                objout.point_data=cell(1);
                objout.point_data_name=cell(1);
            end
            
            objout.cells=C;
            objout.pymesh=pymesh;
            
            % ------ print output ------
            numvtx=size(vtx,1);
            
            fprintf('Verticies: %d and Cells: %d\n',numvtx,numCells);
            
            for iCell=1:numCells
                fprintf('Cell %d: %s %d elements\n',iCell,C(iCell).type,size(C(iCell).tri,1));
            end
            
            if (numPointdata > 0)
                fprintf('Point data: %s\n',strjoin(objout.point_data_name))
            end
            
            if (numCelldata > 0)
                fprintf('Cell data: %s\n',strjoin(objout.cell_data_name))
            end
            
            
        end
        
        function fileout = write(filename,points,cells,cell_data,cell_data_name,point_data,point_data_name)
            %meshio.write write mesh to file using meshio library
            %
            % Inputs:
            % filename  -   needs extension .msh .vtu .vtk etc.
            % points    -   vtx coordinates
            % cells     -   connectivity array OR struct of cells.tri cells.type
            %               as obtained from meshio.read
            % [cell_data]         -  element data in {cell array}
            % [cell_data_name]    -  strings in {cell array}
            % [point_data]        -  optional must match size of points or nodes
            % [point_data_name]   -  optional string
            %
            %
            % Usage - simple example
            % x = rand(20,1);
            % y = rand(20,1);
            % z= rand(20,1);
            % dt = delaunayTriangulation(x,y,z);
            % meshio.write('test.msh',dt.Points,dt.ConnectivityList);
            %
            % With cell and point data
            % celldata={1:size(dt.ConnectivityList,1)};
            % celldataname={'celldata'};
            % meshio.write('test.msh',dt.Points,dt.ConnectivityList,celldata,celldataname);
            % 
            % pointdata={1:size(dt.Points,1)};
            % pointdataname={'pointdata'};
            % meshio.write('test.msh',dt.Points,dt.ConnectivityList,[],[],pointdata,pointdataname);

            
            fprintf('Meshio writing to meshfile : %s\n',filename);
            
            % ------ convert points and cells ------
            % convert points into array
            pypoints = meshio.mat2nparray(points);
            
            
            %if not structure nodes.tri and nodes.type then make one
            if ~isstruct(cells)
                
                % need to know triangle or tetra
                if size(cells,2) ==4
                    celltype='tetra';
                else
                    if size(cells,2) ==3
                        celltype='triangle';
                    else
                        error('nodes must be Nx3 or Nx4');
                    end
                end
                % convert cells into struct form
                cellstmp=cells;
                clear cells;
                cells.tri=cellstmp;
                cells.type=celltype;
            end
            
            pycellslist=py.list;
            
            nCells=size(cells,2);

            
            for iC=1:nCells
                
                % convert cell array correct for python 0 indexing
                % meshio needs nodes as int32
                pycells = meshio.mat2nparray(int32(cells(iC).tri -1));
                
                %correct for wrong shape if only 1 vertex
                if size(cells(iC).tri) ==1
                    pycells.shape=py.tuple({int32(1),int32(1)});
                end
                                    
                % put nodes into a list of CellBlock type
                pycellblock=py.tuple({cells(iC).type,pycells});
                
                pycellslist.append(pycellblock);
                
                celltypes{iC}=cells(iC).type;
                cellsize(iC)=size(cells(iC).tri,1);
                
            end
            
            % multiple cells of the same type will be merged into single
            % cell automatically by meshio, so warn user this is happening
            if (size(celltypes,2) > size(unique(celltypes),2))
                fprintf(2,'Caution repeated cell types will be merged by meshio \n');
            end
            
            
            
            
            % ------ convert cell and point data if needed ------
            
            % check if cell data exists
            celldata =0;
            pointdata =0;
            
            if (exist('cell_data','var') == 1  && ~isempty(cell_data))
                celldata=1;
            end
            
            if (exist('point_data','var') == 1  && ~isempty(point_data))
                pointdata=1;
            end
            % ------ cell data ------
            
            if celldata
                fprintf('Writing cell data...');
                
                if (exist('cell_data_name','var') == 0  && isempty(cell_data))
                    error('cell_data needs cell_data_name');
                end
                
                if ~iscell(cell_data) || ~iscell(cell_data_name)
                    error('cell_data and cell_data name must be cell arrays');
                end
                
                if size(cell_data,2) ~= size(cell_data_name,2)
                    error('cell_data and cell_data must be same length');
                end
                % meshio needs data in certain way:
                % celldict{'dataname' : list[nparray]}
                
                numCelldata=size(cell_data,2);
                CelldataDict=py.dict;
                
                for iCelldata=1:numCelldata
                    curData=cell_data{iCelldata};
                    curName=cell_data_name{iCelldata};
                    fprintf(' %s ',curName);
                    
                    if ~any(ismember(size(curData),cellsize))
                        warning('Celldata %d %s does not match number of elements/cells',iCelldata,curName);
                    end
                    
                    % convert data into pyton nparray
                    curData_py   = meshio.mat2nparray(curData);
                    
                    % convert into list
                    curDatalist = py.list({curData_py});
                    % update dict
                    CelldataDict.update(pyargs(curName,curDatalist));
                    
                end
                fprintf('\n');
            end
            
            % ------ point data ------
            
            if pointdata
                fprintf('Writing point data...');
                
                if (exist('point_data_name','var') == 0  && isempty(point_data))
                    error('point_data needs point_data_name');
                end
                
                if ~iscell(point_data) || ~iscell(point_data_name)
                    error('point data and point_data name must be cell arrays');
                end
                
                if size(point_data,2) ~= size(point_data_name,2)
                    error('point_data and point_data must be same length');
                end
                
                
                % meshio needs data in certain way:
                % pointdict{'dataname' : nparray}
                
                numPointdata=size(point_data,2);
                PointdataDict=py.dict;
                
                for iPointdata=1:numPointdata
                    curData=point_data{iPointdata};
                    curName=point_data_name{iPointdata};
                    fprintf(' %s ',curName);
                    
                    if ~any(size(curData) == size(points,1))
                        warning('Pointdata %d %s does not match number of elements/cells',iPointdata,curName);
                    end
                    
                    % convert data into pyton nparray
                    curData_py   = meshio.mat2nparray(curData);
                    
                    % update dict
                    PointdataDict.update(pyargs(curName,curData_py));
                    
                end
                fprintf('\n');
            end
            
            % ------ call meshio ------
            
            % need to create python cell_data=data point_data=data etc.
            % cannot combine pyarg types, so using this kludge for now:
            if (celldata && pointdata)
                dataargs=pyargs('cell_data',CelldataDict,'point_data',PointdataDict);
            else
                if celldata
                    dataargs=pyargs('cell_data',CelldataDict);
                end
                
                if pointdata
                    dataargs=pyargs('point_data',PointdataDict);
                end
            end
            
            
            % write the file out adding data arguments i.e. cell_data=data
            % in python
            if (exist('dataargs','var') == 1)
                py.meshio.write_points_cells(filename,pypoints,pycellslist,dataargs);
            else
                %write with no extra arguments
                py.meshio.write_points_cells(filename,pypoints,pycellslist);
            end
        end
        
        function fileout = structwrite(filename,objin)
            %meshio.structwrite write file using struct output from meshio.read
            % this function is a wrapper for meshio.write to make it more
            % convenient to convert files
            %
            % Inputs:
            % filename  - needs extension .msh .vtu .vtk etc.
            % objin     - struct same as meshio.read output
            %
            % Usage:
            % M=meshio.read('example.msh');
            % meshio.structwrite('example.vtu',M);
            
            Docelldata=0;
            Dopointdata=0;
            
            % check if cell data is present
            if all(isfield(objin,{'cell_data','cell_data_name'}))
                if ~isempty(objin.cell_data{1})
                    Docelldata=1;
                end
            end
            
            % check if point data is present
            if all(isfield(objin,{'point_data','point_data_name'}))
                if ~isempty(objin.point_data{1})
                    Dopointdata=1;
                end
            end
            
            % call meshio.write with appropriate inputs. this is a kludge
            % im sorry!
            if (Docelldata && Dopointdata)
                
                meshio.write(filename,objin.vtx,objin.cells,objin.cell_data,objin.cell_data_name,objin.point_data,objin.point_data_name)
            else
                if Docelldata
                    meshio.write(filename,objin.vtx,objin.cells,objin.cell_data,objin.cell_data_name)
                end
                if Dopointdata
                    meshio.write(filename,objin.vtx,objin.cells,[],[],objin.point_data,objin.point_data_name)
                end
                if ~any([Docelldata Dopointdata])
                    meshio.write(filename,objin.vtx,objin.cells);
                end
            end
        end
        
        
        function h = plot(objin)
            %plot plot all elements within a chosen file
            % plots verticies, lines and triangles. plots surface mesh of
            % tetra meshes
            %
            % Use paraview for viewing data as its much better!
            %
            % Inputs:
            % objin     - struct from meshio.write
            
            numCells = size(objin.cells,2);
            figure
            hold on
            
            colours=lines(10);
            
            for iCell = numCells:-1:1 %go backwards so smaller elements show on top of tetra/tri
                curCell=objin.cells(iCell).tri;
                
                switch size(objin.cells(iCell).tri,2)
                    case 1 %verticies plot as points
                        plot3(objin.vtx(curCell,1),objin.vtx(curCell,2),objin.vtx(curCell,3),'*','color',colours(1,:));
                    case 2 % lines plot as separate lines
                        for iLine = 1:size(curCell,1)
                            plot3([objin.vtx(curCell(:,1),1) objin.vtx(curCell(:,2),1)],[objin.vtx(curCell(:,1),2) objin.vtx(curCell(:,2),2)],[objin.vtx(curCell(:,1),3) objin.vtx(curCell(:,2),3)],'color',colours(2,:),'Linewidth',1.5);
                        end
                    case 3 % triangles plot as surface - orange
                        h= trisurf(curCell, objin.vtx(:,1), objin.vtx(:,2), objin.vtx(:,3));
                        %                         set(h,'EdgeColor',colours(3,:)*0.9,'FaceColor',colours(3,:));
                    case 4 % tetra plot surface only as tetramesh is v slow
                        trep = triangulation(curCell, objin.vtx);
                        [Triangle_Boundary, Nodes_Boundary] = freeBoundary(trep);
                        h= trisurf(Triangle_Boundary, Nodes_Boundary(:,1), Nodes_Boundary(:,2), Nodes_Boundary(:,3));
                        set(h,'EdgeColor',[0.3,0.3,0.3],'FaceColor','w','FaceAlpha',1); % make tetra grey so easier to see other geometry
                end
            end
            daspect([1,1,1]);
            title(sprintf('Geometry in meshio file, with %d cells',numCells))
        end
        
        function arrayout = mat2nparray( matarray )
            %mat2nparray Convert a Matlab array into an nparray
            %   Convert an n-dimensional Matlab array into an equivalent
            %   nparray
            % small modification of
            % https://www.mathworks.com/matlabcentral/answers/157347-convert-python-numpy-array-to-double#comment_437274
            
            
            data_size=size(matarray);
            if length(data_size)==2 && any(data_size == 1)
                % 1-D vectors are trivial
                arrayout=py.numpy.array(matarray);
            elseif length(data_size)==2
                % A transpose operation is required either in Matlab, or
                % in Python due to the difference between row major and
                % column major ordering
                transpose=matarray';
                % Pass the array to Python as a vector, and then reshape
                % to the correct size
                arrayout=py.numpy.reshape(transpose(:)', int32(data_size));
            else
                % For an n-dimensional array, transpose the first two
                % dimensions to sort the storage ordering issue
                transpose=permute(matarray,[length(data_size):-1:1]);
                % Pass it to python, and then reshape to the python style
                % of matrix sizing
                arrayout=py.numpy.reshape(transpose(:)', int32(fliplr(size(transpose))));
            end
        end
    end
end

