# Triangles

**Triangles** is a specialized mesh generation tool designed for the **SWEpy** (Shallow Water Equation Solver) ecosystem. It generates high-quality triangular and rectangular meshes with bathymetry, initial conditions, and boundary definitions.

## Features

- **Flexible Grids**: Generates structured equilateral or rectangular triangular grids.
- **Topography/Bathymetry**: Supports complex bathymetry generation including islands, slopes, and custom functions.
- **Initial Conditions**: define initial water levels, discharges, and velocities.
- **Pipeline Integration**: Automatically generates the `Config.swe` configuration file required by SWEpy.

## Installation

Ensure you have Python 3.8+ installed.

```bash
pip install -r requirements.txt
```

*Note: `cupy` is recommended for performance if you have an NVIDIA GPU, but the tool will fallback to `numpy` if `cupy` is not available.*

## Usage

The main entry point is `run.py`.

### Command Line Interface

You can run the script with default settings or override parameters via command-line arguments:

```bash
# Run with defaults (dambreak_rose test case)
python run.py

# Run a specific test case (from cases.py)
python run.py --test circular_dambreak --folder my_circle_mesh

# Configure simulation parameters
python run.py --test dambreak_rose --nx 200 --ny 20 --Tmax 300 --dt_save 5

# Override Boundary Conditions (Left, Right, Top, Bottom)
python run.py --bconds soft soft periodic periodic
```

### Full Options List
Use `--help` to see all available options:
```bash
python run.py --help
# Options: --nx, --ny, --dx, --dy, --divisions, --triangles, --Tmax, --tol_dry, --g, --manning, --CFL, --dt_save, --bconds
```

## Operations Output

Running the tool will create an output folder (e.g., `test_rose`) containing:
-   `Config.swe`: The master configuration file.
-   `.txt` files: Data for mesh nodes, elements, bathymetry, water levels, etc.

## Defining Custom Tests

You can define new test cases in `cases.py`. A test function receives the `options` dictionary and must return the following 6 arrays:

```python
def my_custom_test(options):
    # 1. Generate Grid
    nxny = options["nxny"]
    #dxdy = options["dxdy"]
    # Use grids.rectangular_grid or grids.triangular_grid
    x, y = grids.rectangular_grid(xmin, xmax, ymin, ymax, nxny=nxny)
    
    # 2. Define Bathymetry (B)
    B = np.zeros_like(x) - depth
    
    # 3. Define Initial Conditions
    # Water Surface Elevation (W)
    W = np.zeros_like(x) + initial_water_level
    
    # Velocities (U, V)
    U = np.zeros_like(x) + initial_velocity_x
    V = np.zeros_like(x) + initial_velocity_y
    
    # Discharge (HU, HV) calculations
    HU = (W - B) * U
    HV = (W - B) * V
    HUHV = np.array([HU, HV]).T
    
    # 4. Create Mesh Array
    mesh = np.array([x, y]).T
    
    return x, y, B, HUHV, W, mesh
```

Once defined, you can run it using the CLI:
```bash
python run.py --test my_custom_test --folder my_test_folder
```

## License

Distributed under the GPL v3 License. See `LICENSE` for more information.
