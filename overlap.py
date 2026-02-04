from backend import np

def eta(nx, ny, x):
    """
    Returns the indices of the 2x(2nx)x(2ny) rectangular grid triangles with same centroid as 
    the triangles with indices x, assuming xy ordering of the triangles.
    
    Args:
        nx (int): Number of triangles in x direction.
        ny (int): Number of triangles in y direction.
        x (int or array): Triangle index/indices.
        
    Returns:
        int or array: Refided triangle indices.
    """
    index=1+(4*nx+1)*(x%2)+4*((x//2)%ny)+8*nx*(x//(2*nx))
    return index


#returns the indices of the 2x(nx*2)x(ny*2) ... 2x(nx*2^k)x(ny*2^k) rectangular grid triangles with same centroid as the triangles with indices x assuming xy ordering of the triangles
def itereta(k,nx,ny,x):
    """
    Returns indices of refined rectangular grid triangles with same centroid (recursive eta).
    
    Args:
        k (int): Number of refinements.
        nx (int): Initial number of triangles in x.
        ny (int): Initial number of triangles in y.
        x (int or array): Triangle indices.
    """
    result=[]
    result.append(eta(nx,ny,x))
    if k==1:
        return result
    for i in range(1,k):
        result.append(eta((2**i)*nx,(2**i)*ny,result[-1]))
    return result

def delta(nx,ny,x):
    """
    Returns indices for equilateral grid refinement.
    
    Args:
        nx (int): Number of triangles in x.
        ny (int): Number of triangles in y.
        x (int or array): Triangle indices.
    """
    row   = x%(2*ny+1)
    col   = x//(2*ny+1)
    nrowp = 4*ny+1

    index = 2*row+nrowp*(2*col+(1-x%2))
    return index

def iterdelta(k,nx,ny,x):
    """
    Returns indices of refined equilateral grid triangles with same centroid (recursive delta).
    
    Args:
        k (int): Number of refinements.
        nx (int): Initial number of triangles in x.
        ny (int): Initial number of triangles in y.
        x (int or array): Triangle indices.
    """
    result=[]
    result.append(delta(nx,ny,x))
    if k==1:
        return result
    for i in range(1,k):
        result.append(delta((2**i)*nx,(2**i)*ny,result[-1]))
    return result