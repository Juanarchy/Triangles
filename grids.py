from backend import np

def remove_near_duplicates_2d(array, tolerance=1e-5):
    # Sort the array lexicographically
    ind = np.lexsort(np.array([array[:,1],array[:,0]]))
    sorted_array=array[ind]
    
    # Compute the Euclidean distance between consecutive rows
    diff = np.sqrt(np.sum(np.diff(sorted_array, axis=0)**2, axis=1))
    
    # Create a boolean mask where True indicates a unique row (distance greater than tolerance)
    unique_mask = np.ones(len(sorted_array), dtype=bool)
    unique_mask[1:] = diff > tolerance
    
    # Use the mask to select unique rows
    unique_array = sorted_array[unique_mask]
    
    return unique_array

def trunc(values, decs=0):
    return np.trunc(values*10**decs)/(10**decs)


def rectangular_grid(xmin, xmax, ymin, ymax, dxdy=None, nxny=None):
    """
    Generates a rectangular grid of nodes.
    
    Args:
        xmin (float): Minimum x coordinate.
        xmax (float): Maximum x coordinate.
        ymin (float): Minimum y coordinate.
        ymax (float): Maximum y coordinate.
        dxdy (tuple, optional): Tuple of (dx, dy) spacing.
        nxny (tuple, optional): Tuple of (nx, ny) number of elements.
        
    Returns:
        tuple: (x_mesh, y_mesh) arrays of node coordinates.
    """
    # Domain
    if dxdy is not None:
        nx = int( (xmax-xmin)/dxdy[0]) + 1
        ny = int( (ymax-ymin)/dxdy[1]) + 1
    elif nxny is not None:
        nx=nxny[0]+1
        ny=nxny[1]+1
    else:
        print("Arguments missing! assuming nx=ny=100")
        nx=101
        ny=101

    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(x, y, indexing='xy')
    x_mesh, y_mesh = X.flatten(), Y.flatten()

    return x_mesh, y_mesh

def triangular_grid(xmin,xmax,ymin,ymax,dxdy=None,nxny=None):
    """
    Generates a triangular grid (equilateral-like arrangement).
    
    Args:
        xmin (float): Minimum x coordinate.
        xmax (float): Maximum x coordinate.
        ymin (float): Minimum y coordinate.
        ymax (float): Maximum y coordinate.
        dxdy (tuple, optional): Tuple of (dx, dy) spacing. Note: dxdy takes precedence if both provided.
        nxny (tuple, optional): Tuple of (nx, ny) grid dimensions.
        
    Returns:
        tuple: (x_mesh, y_mesh) arrays of node coordinates.
    """
    if dxdy is not None:
        print("Can't use dXdY option with equilateral mesh! Assuming dL=(dX+dY)/2 instead...\n")
        dL=sum(dxdy)/2
    elif nxny is not None:
        dL=trunc((ymax-ymin)/nxny[1],8)

    dLy = dL
    dLx = dL*np.cos(np.pi/6)
    nx = int( (xmax-xmin)/(2*dLx) )
    ny = int( (ymax-ymin)/dLy )
    # Even points
    nx_even = nx+1
    ny_even = ny+1
    x1 = np.linspace(xmin, xmax, nx_even)
    y1 = np.linspace(ymin, ymax, ny_even)
    dLx = x1[1]-x1[0]
    dLy = y1[1]-y1[0]
    X, Y = np.meshgrid(x1, y1, indexing='xy')
    x1_mesh, y1_mesh = X.flatten(), Y.flatten()
    # Odd points
    nx_odd = nx
    ny_odd = ny
    x2 = np.linspace(xmin+.5*dLx, xmax-.5*dLx, nx_odd)
    y2 = np.linspace(ymin+.5*dLy, ymax-.5*dLy, ny_odd)
    X, Y = np.meshgrid(x2, y2, indexing='xy')
    x2_mesh, y2_mesh = X.flatten(), Y.flatten()
    # Top and bottom
    x_bottom = x1[:-1] + .5*dLx
    y_bottom = y1[0] + 0*x_bottom
    x_top = x1[:-1] + .5*dLx
    y_top = y1[-1] + 0*x_top
    # join the meshes
    x_mesh = np.hstack([x1_mesh,x2_mesh,x_top,x_bottom])
    y_mesh = np.hstack([y1_mesh,y2_mesh,y_top,y_bottom])
    # True mesh
    dLx = x2[0]-x1[0]
    dLy = y2[1]-y2[0]

    print("dlx="+str(dLx))
    print("dly="+str(dLy))

    x_mesh=np.round(trunc(x_mesh,7),6)
    y_mesh=np.round(trunc(y_mesh,6),5)

    mesh=np.array([x_mesh, y_mesh],dtype=float).T

    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))
    mesh=mesh[ind]
    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))
    mesh=mesh[ind]
    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))
    mesh=mesh[ind]
    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))
    mesh=mesh[ind]
    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))
    mesh=mesh[ind]

    x_mesh=mesh[:,0]
    y_mesh=mesh[:,1]

    return x_mesh,y_mesh


def nested_grid(nodes,coords,top_indices,bottom_indices):
    """
    Generates a nested grid by refining specified elements.

    Args:
        nodes (array): 3x(#Triangles) element connectivity array.
        coords (array): 2x(#Vertices) array of node coordinates.
        top_indices (array): Indices of the top border cells.
        bottom_indices (array): Indices of the bottom border cells.
        
    Returns:
        array: New coordinates for the refined grid.
    """
    dy=np.abs(nodes[0][1]-nodes[1][1])

    mask=np.zeros(len(nodes),dtype=bool)

    mask[top_indices]=True
    mask[bottom_indices]=True

    mask=~mask

    nontopbot=nodes[np.where(mask)[0]]

    first=(coords[nontopbot[:,0]]+coords[nontopbot[:,1]])/2
    second=(coords[nontopbot[:,1]]+coords[nontopbot[:,2]])/2
    third=(coords[nontopbot[:,2]]+coords[nontopbot[:,0]])/2

    bot=(coords[nodes[bottom_indices][:,0]]+coords[nodes[bottom_indices][:,1]])/2
    top_evens=(coords[nodes[top_indices[::2]][:,2]]+coords[nodes[top_indices[::2]][:,0]])/2
    top_odds =(coords[nodes[top_indices[1::2]][:,1]]+coords[nodes[top_indices[1::2]][:,2]])/2

    newcoords=np.vstack((coords,first,second,third,bot,top_evens,top_odds))

    newcoords=remove_near_duplicates_2d(newcoords,dy*0.000001)

    return newcoords


def load_grid(path):
    """
    Loads grid coordinates (x, y) from a text file.
    
    Args:
        path (str): Path to the node coordinates file.
        
    Returns:
        tuple: (x_mesh, y_mesh) arrays.
    """

    array=np.loadtxt(path, skiprows=1, dtype=float)
    x_mesh=array[:,0]
    y_mesh=array[:,1]

    return x_mesh,y_mesh