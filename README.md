# Triangles

**Triangles** is a specialized mesh generation tool designed for the **SWEpy** (Shallow Water Equation Solver) ecosystem. It generates simple structured equilateral or rectangular triangular meshes with bathymetry, initial conditions, and boundary definitions.

## Features

- **Flexible Grids**: Generates structured equilateral or rectangular triangular grids.
- **Topography/Bathymetry**: Supports complex bathymetry generation including islands, slopes, and custom functions.
- **Initial Conditions**: Define initial water levels, discharges, and velocities.
- **Custom Case Support**: Allows for user-defined tests following a straightforward template.
- **Pipeline Integration**: Automatically generates the files required by SWEpy.

## Installation

Ensure you have Python 3.8+ installed.

```bash
pip install .
```

Or for development (editable mode):

```bash
pip install -e .
```

*Note: `cupy` is recommended for performance if you have an NVIDIA GPU, but the tool will fallback to `numpy` if `cupy` is not available.*

## Usage

## Usage

### Command Line Interface (Recommended)

After installation, use the `triangles-mesh` system-wide command:

```bash
# Get help
triangles-mesh --help

# Run a test case
triangles-mesh --test dambreak_rose --folder my_output
```

### Via Python Script

Alternatively, you can run the script directly from the source directory:

```bash
python run.py --test dambreak_rose --folder my_output
```

### Configuration Examples

Both methods support the same configuration flags.

```bash
# Run a specific test case (from cases.py)
triangles-mesh --test circular_dambreak --folder my_circle_mesh

# Configure simulation parameters
triangles-mesh --test dambreak_rose --nx 200 --ny 20 --Tmax 300 --dt_save 5

# Override Boundary Conditions (Left, Right, Top, Bottom)
triangles-mesh --bconds soft soft periodic periodic
```

### Full Options List
```bash
triangles-mesh --help
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
triangles-mesh --test my_custom_test --folder my_test_folder
```

## License

Distributed under the GPL v3 License. See `LICENSE` for more information.

