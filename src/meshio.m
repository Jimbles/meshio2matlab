classdef meshio
    %MESHIO Summary of this class goes here
    %   Detailed explanation goes here
    
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
                nodestmp=pymesh.cells{iCell}.data;
                nodes =meshio.np2mat(nodestmp);
                
                % python uses 0 indexing
                nodes= nodes+1;
                
                C(iCell).nodes=nodes;
            end
            
            objout.numCells=numCells;
            objout.Cells=C;
            objout.pymesh=pymesh;
            
            % print output
            numvtx=size(vtx,1);
            
            fprintf('Verticies: %d \n',numvtx);
            fprintf('Cells: %d\n',numCells);
            
            for iCell=1:numCells
                fprintf('Cell %d: %s %d elements\n',iCell,C(iCell).type,size(C(iCell).nodes,1));
            end
            
        end
        
        function fileout = write(filename,points,nodes,celldata,dataname)
            
            % check celldata and nodes match, transpose is ok
            if ~any(size(nodes,1) == (size(celldata)))
                error('nodes and celldata dont match');
            end
            
            
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
            pycells=meshio.mat2nparray(nodes -1);
            
            % put nodes into a list of CellBlock type
            pycellblock=py.meshio.CellBlock(celltype,pycells);
            pycellslist=py.list({pycellblock});
            
            % convert points into array
            pypoints = meshio.mat2nparray(points);
            
            
            % meshio needs data in certain way:
            % celldict{'dataname' : list[nparray]}
            
            npdata = meshio.mat2nparray(celldata);
            datalist = py.list({npdata});
            datadict=py.dict(pyargs(dataname,datalist));
            
            % write the file out
            py.meshio.write_points_cells(filename,pypoints,pycellslist,pyargs('cell_data',datadict));
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

