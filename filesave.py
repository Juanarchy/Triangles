import os
from backend import np

def mkdir(folder):
    """
    Creates a directory if it doesn't exist.
    
    Args:
        folder (str): Path to the folder to create.
        
    Returns:
        str: Absolute path to the created folder.
    """
    path = os.path.abspath(folder)
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise
    return path


#saves array data in txt
def save_data(path, filename, data, fmt="%1.8E"):
    """
    Saves numpy/cupy array data to a text file.
    
    Args:
        path (str): Directory path.
        filename (str): Name of the file.
        data (array): Data array to save.
        fmt (str): Format string for saving.
    """
    # Load the file
    file_path = os.path.join(path, filename)
    if data is None:
        return
    if len(data.shape)==1:
        header = "%d 1" %(data.shape)
    else:
        header = "%d %d" %(data.shape)
    np.savetxt(file_path, data, fmt=fmt, header=header)
    return

#writes bathymetry file
def bathymetry(folder_path,B):
    """Saves bathymetry data to Bathymetry.txt"""
    save_data(folder_path,"Bathymetry.txt",B)
    return

#writes config file
def config(folder_path,options):
    """
    Writes the Master Config.swe file.
    
    Args:
        folder_path (str): Output directory.
        options (dict): Simulation options.
    """
    cfl = options['CFL'] if "CFL" in options.keys() else 0.25

    file_path = os.path.join(folder_path, "Config.swe")
    fileswe = open(file_path, "w")
    fileswe.write("Tmax           " + str(options['Tmax']) + "\n")
    fileswe.write("Dry            " + str(options['tol_dry']) + "\n")
    fileswe.write("Tolerance      " + str(options['tol_dry']) + "\n")
    fileswe.write("Gravity        " + str(options['g']) + "\n")
    fileswe.write("Roughness      " + str(options['manning']) + "\n")
    fileswe.write("CFL            " + str(cfl) + "\n")
    fileswe.write("dt_save        " + str(options['dt_save']) + "\n")
    fileswe.write("WaterLevel     WaterLevel.txt\n")
    fileswe.write("Discharge      Discharge.txt\n")
    fileswe.write("Bathymetry     Bathymetry.txt\n")
    fileswe.write("Coordinates    NodeCoords.txt\n")
    fileswe.write("Sides          ElemNeighSides.txt\n")
    fileswe.write("Elements       ElemNodes.txt\n")
    fileswe.write("Neighboors     ElemNeighs.txt\n")
    fileswe.write("GhostCells     GhostCells.txt\n")
    fileswe.close()
    return


#writes initial discharge file
def discharge(folder_path,HUHV):
    """Saves discharge data to Discharge.txt"""
    save_data(folder_path,"Discharge.txt",HUHV)
    return

#writes neighbors file
def elemneighs(folder_path,neighs):
    """Saves element neighbor data to ElemNeighs.txt"""
    save_data(folder_path,"ElemNeighs.txt",neighs,"%d")
    return

#writes connected sides file
def elemneighsides(folder_path,neighsides):
    """Saves neighbor sides data to ElemNeighSides.txt"""
    save_data(folder_path,"ElemNeighSides.txt",neighsides,"%d")
    return

#write triangle nodes file
def elemnodes(folder_path,nodes):
    """Saves element node connectivity to ElemNodes.txt"""
    save_data(folder_path,"ElemNodes.txt",nodes,"%d")
    return

#writes ghost cells file
def ghostcells(folder_path,ghosts):
    """Saves ghost cells data to GhostCells.txt"""
    save_data(folder_path,"GhostCells.txt",ghosts,"%d")
    return

#writes node coordinates file
def nodecoords(folder_path,mesh):
    """Saves node coordinates to NodeCoords.txt"""
    save_data(folder_path,"NodeCoords.txt",mesh)
    return

#writes initial water level file
def waterlevel(folder_path,W):
    """Saves water level data to WaterLevel.txt"""
    save_data(folder_path,"WaterLevel.txt",W)
    return

def refplaces(folder_path,indices):
    """Saves refinement indices to RefPlaces.txt"""
    save_data(folder_path,"RefPlaces.txt",indices,"%d")
    return

def save_all(folder_path,B,options,HUHV,neighs,neighsides,nodes,ghosts,mesh,W):
    """
    Orchestrates saving all mesh and simulation data to the specified folder.
    
    Args:
        folder_path (str): Output directory.
        B (array): Bathymetry.
        options (dict): Simulation options.
        HUHV (array): Discharge.
        neighs (array): Neighbors.
        neighsides (array): Neighbor sides.
        nodes (array): Element nodes.
        ghosts (array): Ghost cells.
        mesh (array): Node coordinates.
        W (array): Water level.
    """
    mkdir(folder_path)
    bathymetry(folder_path,B)
    config(folder_path,options)
    discharge(folder_path,HUHV)
    elemneighs(folder_path,neighs)
    elemneighsides(folder_path,neighsides)
    elemnodes(folder_path,nodes)
    ghostcells(folder_path,ghosts)
    nodecoords(folder_path,mesh)
    waterlevel(folder_path,W)
    return