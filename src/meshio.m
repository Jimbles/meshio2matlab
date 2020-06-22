classdef meshio
    %MESHIO Matlab2Meshio functions
    %   Use the python meshio library to read and write unstructured mesh
    %   files .msh .vtu .vtk .xdmf tetgen ansys .stl .obj
    %
    %   more formats at:
    %   https://github.com/nschloe/meshio
    %
    % meshio.read - read mesh file
    % meshio.write - write matlab mesh to file
    % meshio.plot - plot all contents of a meshfile
    %
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
            %read read mesh file using meshio
            %   Calls python meshio library and processes output to a
            %   matlab struct
            %
            % objout.vtx - verticies
            % objout.Cells - structure array for each geometry saved in file
            % objout.Cells.tri - Trigangulation connectivity list for this cell
            % objout.Cells.type - 'Tetra','Triangle','Line','Vertex'
         
            
            
            
            fprintf('Meshio reading meshfile : %s\n',filename);
            
            % read mesh into meshio python format
            pymesh=py.meshio.read(filename);
            
            % vtx are saved as 1 array
            vtx =meshio.np2mat(pymesh.points);
            
            objout.vtx=vtx;
            
            % files can have mix of data
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
            
            % get data - cannot convert from dict to struct as python
            % allows spaces but matlab doesnt
            
            % cell data
            cell_data_names_py=py.list(pymesh.cell_data.keys);
            numCelldata=double(py.len(cell_data_names_py));
            
            if numCelldata > 0
                cell_data_name=char(cell_data_names_py{1});
                cell_data_py=pymesh.cell_data{cell_data_name}{1};
                cell_data=meshio.np2mat(cell_data_py);
                
                objout.cell_data=cell_data;
                objout.cell_data_name=cell_data_name;
            else
                cell_data=[];
                cell_data_name='';
            end
            
            % point data
            point_data_names_py=py.list(pymesh.point_data.keys);
            numPointdata=double(py.len(point_data_names_py));
            
            if numPointdata > 0
                point_data_name=char(point_data_names_py{1});
                point_data_py=pymesh.point_data{point_data_name}{1};
                point_data=meshio.np2mat(point_data_py);
                
                objout.point_data=point_data;
                objout.point_data_name=cell_data_name;
            else
                point_data=[];
                point_data_name='';
            end
            
            objout.cells=C;
            objout.pymesh=pymesh;
            
            % print output
            numvtx=size(vtx,1);
            
            fprintf('Verticies: %d \n',numvtx);
            fprintf('Cells: %d\n',numCells);
            
            for iCell=1:numCells
                fprintf('Cell %d: %s %d elements\n',iCell,C(iCell).type,size(C(iCell).tri,1));
            end
            
        end
        
        function fileout = write(filename,points,nodes,data,dataname)
            %write write mesh to file using meshio library
            %
            % Inputs
            % filename - needs extension
            % points - vtx
            % nodes - connectivity array
            % [data] - optional must match size of points or nodes
            % [dataname] - optional string
            
            
            fprintf('Meshio writing to meshfile : %s\n',filename);
            
            % ------ convert points and cells ------
            
            % need to know triangle or tetra
            
            if size(nodes,2) ==4
                celltype='tetra';
            else
                if size(nodes,2) ==3
                    celltype='triangle';
                else
                    error('nodes must be Nx3 or Nx4');
                end
            end
            
            % convert cell array correct for python 0 indexing
            % meshio needs nodes as int32
            pycells = meshio.mat2nparray(int32(nodes -1));
            
            % put nodes into a list of CellBlock type
            pycellblock = py.meshio.CellBlock(celltype,pycells);
            pycellslist = py.list({pycellblock});
            
            % convert points into array
            pypoints = meshio.mat2nparray(points);
            
            % ------ convert cell or point data if needed ------
            
            % default dataname
            if exist('dataname','var') == 0  || isempty(dataname)
                dataname = 'Data';
            end
            
            % write data if it exists
            if (exist('data','var') == 1  && ~isempty(data))
                celldata=0;
                pointdata=0;
                
                % check celldata and nodes match, transpose is ok
                if any(size(nodes,1) == (size(data)))
                    celldata=1;
                else
                    if any(size(points,1) == (size(data)))
                        pointdata=1;
                    end
                end
                
                if ~any([celldata pointdata])
                    error('Data does not match number cells for celldata or points for pointdata');
                end
                
                %convert data into pyton nparray
                npdata   = meshio.mat2nparray(data);
                
                %data is different format for cell or point data
                if celldata
                    % meshio needs data in certain way:
                    % celldict{'dataname' : list[nparray]}
                    
                    datalist = py.list({npdata});
                    datadict = py.dict(pyargs(dataname,datalist));
                    dataargs=pyargs('cell_data',datadict);
                else
                    if pointdata
                        % meshio needs data in certain way:
                        % celldict{'dataname' : nparray}
                        datadict = py.dict(pyargs(dataname,npdata));
                        dataargs=pyargs('point_data',datadict);
                    end
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
        
        function h = plot(objin)
            %plot plot all elements within a chosen file
            % plots verticies, lines and triangles. plots surface mesh of
            % tetra meshes
            % Use paraview for viewing data as its much better!
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

