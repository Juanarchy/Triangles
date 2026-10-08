#This is an example script that shows the usage of triangles within code instead of through CLI.
#It's just calling main.execute() with a dictionary containing all the desired options within elements with keys matching their name

from main import execute

options={
    # "nxny": Grid dimensions [nx, ny] - Number of elements/nodes in X and Y directions
    # "dxdy": Grid spacing [dx, dy] - Alternative to nxny. If provided, it takes precedence in some grids.
    "nxny": [100, 50],
    # "dxdy": [1.0, 1.0],

    # "folder": Name of the output directory where files will be saved
    "folder": "folder_name",

    # "abspath": Path to the location in which directory will be created
    "abspath": ".",

    # "Tmax": Maximum simulation time (in seconds)
    "Tmax": 180,

    # "tol_dry": Tolerance for dry cells (water depth threshold)
    "tol_dry": 0.000001,

    # "g": Acceleration due to gravity (m/s^2)
    "g": 9.81,

    # "manning": Manning's roughness coefficient (n)
    "manning": 0,

    # "CFL": Courant-Friedrichs-Lewy condition number for time step control
    "CFL": 0.25,

    # "dt_save": Time interval for saving simulation results to disk
    "dt_save": 1,

    # "bconds": Boundary conditions for the 4 sides: [West, East, North, South]
    # Options: "soft" (open), "wall" (reflective), "periodic"
    "bconds": ["wall", "wall", "wall", "wall"], 

    # "divisions": Number of times to recursively refine the grid (0 = no refinement)
    "divisions": 0,

    # "triangles": Mesh topology type - "equilateral" or "rectangular"
    "triangles": "equilateral",

    # "test": Name of the test case function in `tests.py` to use for bathymetry/IC setup
    "test": "dambreak_channel2shoebox" 
}

execute(options)