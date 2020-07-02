# create file with multiple point and cell data, showing correct file types and meshio format

import meshio
import os
import numpy

path = os.path.dirname(os.path.abspath(__file__))

inpath = os.path.join(path,'example.msh')
outpath = os.path.join(path,'pyexample.msh')
outpathvtu = os.path.join(path,'pyexample.vtu')

M = meshio.read(inpath)

c=[(M.cells[-1].type,M.cells[-1].data)]

celldata=numpy.zeros([c[0][1].shape[0]])
celldict={'ex_cd1': [celldata],'ex_cd2':[celldata+1],'ex_cd3':[celldata+2]}


pointdata=numpy.zeros(M.points.shape[0])

pointdict={'ex_pd1': pointdata,'ex_pd2': pointdata +1,'ex_pd3': pointdata +2}


meshio.write_points_cells(outpath,M.points,M.cells,cell_data=celldict,point_data=pointdict)
P=meshio.read(outpath)
meshio.write(outpathvtu,P)

print("Done")