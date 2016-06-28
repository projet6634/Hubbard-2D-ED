import struct
import numpy as np
fh=open("./h_header_01-0000.txt","r")
fv=open("./h_val_01-0000.txt","r")
fo=open("./h_mm.mtx","w")
ncols = int(fh.readline())
colptr = np.zeros(ncols+1,dtype=np.int32)
for i in range(0,ncols):
    colptr[i] = int(fh.readline())
nnztot = int(fh.readline())
colptr[ncols] = nnztot+1

fo.write("%%MatrixMarket matrix coordinate real general\n")
fo.write("%s %s %d\n"%(ncols, ncols, nnztot))

fh.seek(0)
fh.readline()
for colidx in range(1,ncols+1):
    nnz = colptr[colidx] - colptr[colidx-1]
    print colidx,nnz
    for j in range(0,nnz):
        row = fv.readline().split()
        rowidx = int(row[0])
        val = float(row[1])
        fo.write("%5d %5d %12.6f\n"%(rowidx,colidx,val))
fo.close()
fh.close()
fv.close()

def readslice(inputfilename,nx,ny,timeslice):
    f = open(inputfilename,'rb')
    f.seek(8*timeslice*nx*ny)
    field = np.fromfile(f,dtype='float64',count=nx*ny)
    field = np.reshape(field,(nx,ny))
    f.close()
    return field
