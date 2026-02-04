import cupy as np

#Returns the indices of the 2x(2nx)x(2ny) grid triangles with same centroid as the triangles with indices x assuming xy ordering of the triangles
def eta(nx,ny,x):
    index=1+(4*nx+1)*(x%2)+4*((x//2)%ny)+8*nx*(x//(2*nx))
    return index

#returns the indices of the 2x(nx*2)x(ny*2) ... 2x(nx*2^k)x(ny*2^k) grid triangles with same centroid as the triangles with indices x assuming xy ordering of the triangles
def itereta(k,nx,ny,x):
    result=[]
    result.append(eta(nx,ny,x))
    if k==1:
        return result
    for i in range(1,k):
        result.append(eta((2**i)*nx,(2**i)*ny,result[-1]))
    return result

def delta(nx,ny,x):
    row   = x%(2*ny+1)
    col   = x//(2*ny+1)
    nrowp = 4*ny+1

    index = 2*row+nrowp*(2*col+(1-x%2))
    return index

def iterdelta(k,nx,ny,x):
    result=[]
    result.append(delta(nx,ny,x))
    if k==1:
        return result
    for i in range(1,k):
        result.append(delta((2**i)*nx,(2**i)*ny,result[-1]))
    return result