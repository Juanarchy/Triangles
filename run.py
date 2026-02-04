
"""
Master script to generate a mesh for SWEpy (cuSWE).

Usage:
    python run.py --test dambreak_rose --folder test_rose --nx 100 --ny 10
    python run.py --help

This script defines the `options` dictionary with all simulation parameters and then
calls `execute(options)` to generate the necessary files.
You can override parameters using command line arguments.
"""

import argparse
from main import execute

# Default Configuration options
options={
    # "nxny": Grid dimensions [nx, ny] - Number of elements/nodes in X and Y directions
    # "dxdy": Grid spacing [dx, dy] - Alternative to nxny. If provided, it takes precedence in some grids.
    "nxny": [100, 10],
    # "dxdy": [1.0, 1.0],

    # "folder": Name of the output directory where files will be saved
    "folder": "test_rose",

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
    "bconds": ["soft", "soft", "periodic", "periodic"], 

    # "divisions": Number of times to recursively refine the grid (0 = no refinement)
    "divisions": 0,

    # "triangles": Mesh topology type - "equilateral" or "rectangular"
    "triangles": "equilateral",

    # "test": Name of the test case function in `tests.py` to use for bathymetry/IC setup
    "test": "dambreak_rose" 
}

def parse_args():
    parser = argparse.ArgumentParser(description="Triangles Mesh Generator for SWEpy")
    parser.add_argument("--test", type=str, help="Test case name (from cases.py)")
    parser.add_argument("--folder", type=str, help="Output folder name")
    
    # Grid Dimensions
    parser.add_argument("--nx", type=int, help="Number of elements in X")
    parser.add_argument("--ny", type=int, help="Number of elements in Y")
    parser.add_argument("--dx", type=float, help="Grid spacing in X")
    parser.add_argument("--dy", type=float, help="Grid spacing in Y")
    parser.add_argument("--divisions", type=int, help="Number of successive grid refinements. If >0, it will create additional finer grids by midpoint subdivision.")
    parser.add_argument("--triangles", type=str, choices=["equilateral", "rectangular"], help="Mesh topology")

    # Simulation Parameters
    parser.add_argument("--Tmax", type=float, help="Maximum simulation time")
    parser.add_argument("--tol_dry", type=float, help="Tolerance for dry cells")
    parser.add_argument("--g", type=float, help="Gravity acceleration")
    parser.add_argument("--manning", type=float, help="Manning's roughness coefficient")
    parser.add_argument("--CFL", type=float, help="CFL condition number")
    parser.add_argument("--dt_save", type=float, help="Time interval for saving results")

    # Boundary Conditions
    parser.add_argument("--bconds", nargs=4, type=str, metavar=('LEFT', 'RIGHT', 'TOP', 'BOTTOM'), 
                        help="Boundary conditions (e.g., --bconds soft soft periodic periodic)")

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    
    # Override defaults with CLI args
    if args.test:
        options["test"] = args.test
    if args.folder:
        options["folder"] = args.folder
    
    # Grid dimensions handling
    if args.nx:
        options["nxny"][0] = args.nx
    if args.ny:
        options["nxny"][1] = args.ny
        
    # dxdy handling (if provided, it might override nxny in grids.py logic depending on implementation)
    # The default options doesn't have dxdy, but main.py docstring says it's an alternative.
    # main.py passes `options` to `mybench` and `assembly`. 
    # grids.py functions accept `dxdy`. So we should add it to options if present.
    if args.dx or args.dy:
        # Default dx/dy if one is missing but the other is present? 
        # Usually both are needed or assumed equal. Let's just pass what is given.
        dx = args.dx if args.dx else 1.0
        dy = args.dy if args.dy else 1.0
        options["dxdy"] = [dx, dy]

    if args.divisions is not None:
        options["divisions"] = args.divisions
    if args.triangles:
        options["triangles"] = args.triangles

    # Simulation Params
    if args.Tmax: 
        options["Tmax"] = args.Tmax
    if args.tol_dry:
        options["tol_dry"] = args.tol_dry
    if args.g:
        options["g"] = args.g
    if args.manning:
        options["manning"] = args.manning
    if args.CFL:
        options["CFL"] = args.CFL
    if args.dt_save:
        options["dt_save"] = args.dt_save
        
    # Boundary Conditions
    if args.bconds:
        options["bconds"] = args.bconds

    print(f"Generating mesh for test case '{options['test']}' in folder '{options['folder']}'...")
    execute(options)
