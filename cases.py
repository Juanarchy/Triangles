from backend import np
import grids
import numpy as op
import scipy as sp

sech = lambda x: np.divide(2.0*np.exp(-1.0*x),1.0 + np.exp(-2.0*x))

def equilateral_plain_grid(options):

    xmax=12141500
    xmin=-xmax
    ymax=6590950
    ymin=-ymax

    nxnyu=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=0*np.exp( -25.*(x-(xmax-xmin)/2)**2 - 50.*(y-(ymax-ymin)/2)**2)

    U=np.zeros_like(x)
    V=np.zeros_like(x)

    W=np.zeros_like(x)

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh


def bryson_example_1(options):

    xmax = 2
    xmin = 0
    ymax = 1
    ymin = 0

    nxnyu=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=0.5*np.exp( -25.*(x-(xmax-xmin)/2)**2 - 50.*(y-(ymax-ymin)/2)**2)

    U=np.zeros_like(x)+0.3
    V=np.zeros_like(x)

    W=np.zeros_like(x)+1

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def bryson_example_1_equi(options):

    xmax = 2
    xmin = 0
    ymax = 1
    ymin = 0

    nxnyu=options["nxny"]

    if "forced_mesh" not in options.keys():
        x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)
    else:
        x=options["forced_mesh"][:,0]
        y=options["forced_mesh"][:,1]

    B=0.5*np.exp( -25.*(x-(xmax-xmin)/2)**2 - 50.*(y-(ymax-ymin)/2)**2)

    U=np.zeros_like(x)+0.3
    V=np.zeros_like(x)

    W=np.zeros_like(x)+1

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def bryson_example_1DX(options):

    xmax = 2
    xmin = 0
    ymax = 1
    ymin = 0

    nxnyu=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=0.5*np.exp( -25.*(x-(xmax-xmin)/2)**2)

    U=np.zeros_like(x)+0.3
    V=np.zeros_like(x)

    W=np.zeros_like(x)+1

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def bryson_example_1DX_equi(options):

    xmax = 2
    xmin = 0
    ymax = 1
    ymin = 0

    nxnyu=options["nxny"]

    if "forced_mesh" not in options.keys():
        x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)
    else:
        x=options["forced_mesh"][:,0]
        y=options["forced_mesh"][:,1]

    B=0.5*np.exp( -25.*(x-(xmax-xmin)/2)**2)

    U=np.zeros_like(x)+0.3
    V=np.zeros_like(x)

    W=np.zeros_like(x)+1

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def bryson_example_1DY(options):

    xmax = 1
    xmin = 0
    ymax = 2
    ymin = 0

    nxnyu=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=0.5*np.exp( -25.*(y-(ymax-ymin)/2)**2)

    U=np.zeros_like(x)
    V=np.zeros_like(x)+0.3

    W=np.zeros_like(x)+1

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def wavesplit_eq(options):
    
    xmax = 500
    xmin = 0
    ymax = 5
    ymin = 0

    H=0.25
    d=1
    x0=(xmax-xmin)/10

    gamma = np.sqrt( 0.75 * H/d**3)

    nxnyu=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=np.zeros_like(x)-d

    U=np.zeros_like(x)
    V=np.zeros_like(x)

    W=H*sech(gamma*(x-x0))**(2)

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def wavesplit_rec(options):
    
    xmax = 500
    xmin = 0
    ymax = 5
    ymin = 0

    H=0.1
    d=1
    x0=(xmax+xmin)/10

    gamma = np.sqrt( 0.75 * H/d**3)

    nxnyu=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=np.zeros_like(x)-d

    U=np.zeros_like(x)
    V=np.zeros_like(x)

    W=H*sech(gamma*(x-x0))**(2)

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh


def maule_reconfig(options):
    
    base_path = options.get("data_path", ".")
    
    # Check if files exist
    node_coords_file = os.path.join(base_path, "NodeCoords.txt")
    if not os.path.exists(node_coords_file):
         raise FileNotFoundError(f"Could not find required data file: {node_coords_file}. Please set 'data_path' in options.")

    x,y=grids.load_grid(node_coords_file)

    B=np.loadtxt(os.path.join(base_path, "Bathymetry.txt"),dtype=float,skiprows=1)

    HUHV=np.loadtxt(os.path.join(base_path, "Discharge.txt"),dtype=float,skiprows=1)

    W=np.loadtxt(os.path.join(base_path, "WaterLevel.txt"),dtype=float,skiprows=1)

    mesh=np.array([x, y]).T

    if "forced_mesh" in options.keys():
        x2=options["forced_mesh"][:,0]
        y2=options["forced_mesh"][:,1]

        splB = sp.interpolate.SmoothBivariateSpline(np.asnumpy(x),np.asnumpy(y),np.asnumpy(B))
        splW = sp.interpolate.SmoothBivariateSpline(np.asnumpy(x),np.asnumpy(y),np.asnumpy(W))
        splHU= sp.interpolate.SmoothBivariateSpline(np.asnumpy(x),np.asnumpy(y),np.asnumpy(HUHV[:,0]))
        splHV= sp.interpolate.SmoothBivariateSpline(np.asnumpy(x),np.asnumpy(y),np.asnumpy(HUHV[:,1]))

        B=splB.ev(np.asnumpy(x2),np.asnumpy(y2))
        W=splW.ev(np.asnumpy(x2),np.asnumpy(y2))
        HU=splHU.ev(np.asnumpy(x2),np.asnumpy(y2))
        HV=splHV.ev(np.asnumpy(x2),np.asnumpy(y2))

        x=np.asarray(x2)
        y=np.asarray(y2)
        mesh=np.array([x, y]).T
        HUHV=np.array([np.asarray(HU),np.asarray(HV)]).T


    return x,y,B,HUHV,W,mesh

def runup_simple(options):

    # Suggested Model Parameters  
    Xo   = 19.85   # [m]
    X1   = 34.53   # [m]
    d    = 1.000   # [m]
    H    = 0.100   # [m]
    g    = 9.800   # [m/s2]
  
    # Discretization parameters
    xmin = -10.0   # [m]
    xmax =  70.0   # [m]
    ymin =  -5.0   # [m]
    ymax =   5.0   # [m]

    # Related variables
    gamma = np.sqrt(0.75 * H/d**3) # [1/m]

    nxnyu=options["nxny"]
    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B = np.where(x <= Xo, -d/Xo*x, -d)

    W = (H/d)*sech(gamma*(x - X1))**2;

    U = -np.sqrt(g/d)*W*(1.00 - 0.25*W/d)
    V =  np.zeros_like(x)

    HU=np.maximum((W-B),0)*U
    HV=np.maximum((W-B),0)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def runup_simple_2(options):

    # Suggested Model Parameters  
    Xo   = 19.85   # [m]
    X1   = 34.53   # [m]
    d    = 1.000   # [m]
    H    = 0.100   # [m]
    g    = 9.800   # [m/s2]
  
    # Discretization parameters
    xmin = -10.0   # [m]
    xmax =  70.0   # [m]
    ymin =  -5.0   # [m]
    ymax =   5.0   # [m]

    # Related variables
    gamma = np.sqrt(0.75 * H/d**3) # [1/m]

    nxnyu=options["nxny"]
    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B = np.where(x <= Xo, -d, -d)

    W = (H/d)*sech(gamma*(x - X1))**2;

    U = -np.sqrt(g/d)*W*(1.00 - 0.25*W/d)
    V =  np.zeros_like(x)

    HU=np.maximum((W-B),0)*U
    HV=np.maximum((W-B),0)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def dambreak(options):

    xmin=0
    xmax=50
    ymin=0
    ymax=5

    nxny=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B = np.zeros_like(x)

    W = np.where(x<(xmin+xmax)/2, 1,0)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh

def dambreak_rose(options):
    xmin=0
    xmax=4000
    ymin=0
    ymax=100

    w0 = 30.5

    xl=1695
    xr=2310

    nxny=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B = np.zeros_like(x)

    W = np.where((xl<x)&(x<xr), w0,1)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh


def circular_dambreak(options):

    xmin=0
    xmax=100
    ymin=0
    ymax=100
    r=25
    h=(xmax+xmin)/2
    k=(ymax+ymin)/2

    nxny=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B = np.zeros_like(x)

    W = np.where((x-h)**2+(y-k)**2<r, 1,0.5)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh

def circular_dambreak_rough(options):

    xmin=0
    xmax=100
    ymin=0
    ymax=100
    r=25
    h=(xmax+xmin)/2
    k=(ymax+ymin)/2

    nxny=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B=np.random.rand(x.size)/10

    W = np.where((x-h)**2+(y-k)**2<r, 1,0.5)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh

def circular_dambreak_parabolic(options):

    xmin=-200
    xmax=200
    ymin=-200
    ymax=200
    wmin=1
    wmax=1.1
    bmin=-2
    r=20
    h=(xmax+xmin)/2
    k=(ymax+ymin)/2

    scale = ((xmax-h)**2+(ymax-k)**2)/(0.5*wmin-bmin)

    nxny=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B=((x-h)**2+(y-k)**2)/scale+bmin

    W = np.where((x-h)**2+(y-k)**2<r, wmax,wmin)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh

def circular_dambreak_eq(options):

    xmin=-200
    xmax=200
    ymin=-200
    ymax=200
    wmin=1
    wmax=1.1
    bmin=-2
    r=20
    h=(xmax+xmin)/2
    k=(ymax+ymin)/2

    nxny=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B=np.zeros_like(x)+bmin

    W = np.where((x-h)**2+(y-k)**2<=r, wmax,wmin)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh

def circular_dambreak_rec(options):

    xmin=-200
    xmax=200
    ymin=-200
    ymax=200
    r=20
    h=(xmax+xmin)/2
    k=(ymax+ymin)/2

    nxny=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,None,nxny)

    B = np.zeros_like(x)-2

    W = np.where((x-h)**2+(y-k)**2<r, 1.1,1)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x,y]).T

    return x,y,B,HUHV,W,mesh

def conical_island(options):

    # Suggested Model Parameters
    g      = options["g"] # [m/s2]
    water_depth       = 0.320 # [m]
    tank_width        = 25 # [m]
    tank_length       = 30 # [m]
    island_height     = 0.625 # [m]
    island_toe_diam   = 7.200 # [m]
    island_crest_diam = 2.200 # [m]
    epsilon = 0.093

    # Water level Parameters
    H  =   epsilon*water_depth # [m]
    X1 = -13  # [m]

    # Discretization parameters
    xmin = -25 # [m]
    xmax =  16 # [m]
    ymin = -tank_length/2.0 # [m]
    ymax =  tank_length/2.0 # [m]
    nxnyu= options["nxny"]

    # Related variables
    gamma = np.sqrt(3*H/(4*water_depth**3)) # [1/m]

    ################################################################################
    # Mesh
    ################################################################################

    if options["triangles"].lower() in "rectangular":
        from grids import rectangular_grid as gridfunc
    else:
        from grids import triangular_grid as gridfunc

    x, y = gridfunc(xmin,xmax,ymin,ymax,nxny=nxnyu)
    mesh = np.array([x, y]).T

    ################################################################################
    # Bathymetry: (B)
    ################################################################################
    R  = np.sqrt(x**2 + y**2)
    Z = 0.25*(island_toe_diam/2.0 - R)*(island_crest_diam/2.0 < R)*(R < island_toe_diam/2.0)
    Z = Z + island_height*(R <= island_crest_diam/2.0)
    Z = Z - water_depth
    B  = Z

    ################################################################################
    # Water Level Surface: (W)
    ################################################################################
    Z = (H)*sech(gamma*(x - X1))**2
    W = Z

    ################################################################################
    # Flux Velocities: (U0,V0)
    ################################################################################
    U = np.sqrt(g/water_depth)*Z
    V = np.zeros_like(x)

    HU=np.maximum((W-B),0)*U
    HV=np.maximum((W-B),0)*V

    HUHV=np.array([HU,HV]).T

    return x,y,B,HUHV,W,mesh

def wet_dry_test(options):

    xmax = 5
    xmin = 0
    ymax = 1
    ymin = 0
    x0 = 2.41
    x1 = 2.59
    y0 = 0.4
    y1 = 0.6
    apex = 10
    d = 3
    xmid=(x0+x1)/2
    ymid=(y0+y1)/2
    slopex=apex/(x1-x0)
    slopey=apex/(y1-y0)

    nxnyu=options["nxny"]

    if options["triangles"].lower() in "rectangular":
        from grids import rectangular_grid as gridfunc
    else:
        from grids import triangular_grid as gridfunc

    x, y = gridfunc(xmin,xmax,ymin,ymax,nxny=nxnyu)
    mesh = np.array([x, y]).T

    line1=(y1-y0)/(x1-x0)*(x-x0)+y0
    line2=(y0-y1)/(x1-x0)*(x-x1)+y0

    B = np.where((x>=x0) & (x<=x1) & (y<=np.minimum(line1,line2)) & (y>=y0),  slopey*(y-y0),0)
    B = np.where((x>=x0) & (x<=x1) & (y>=np.maximum(line1,line2)) & (y<=y1), -slopey*(y-y1),B)
    B = np.where((x>=x0) & (x<=xmid) & (y>=line1) & (y<=line2),  slopex*(x-x0),B)
    B = np.where((x<=x1) & (x>=xmid) & (y<=line1) & (y>=line2), -slopex*(x-x1),B)

    U=np.zeros_like(x)
    V=np.zeros_like(x)

    W=np.zeros_like(x)+d

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh

def wellbalancedtest(options):
    
    xmax = 1
    xmin = 0
    ymax = 1
    ymin = 0

    H=0
    d=2
    x0=(xmax-xmin)/2

    gamma = np.sqrt( 0.75 * H/d**3)

    nxnyu=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=np.random.rand(x.size)/10-d
    #B=np.zeros_like(x)-d

    U=np.zeros_like(x)
    V=np.zeros_like(x)

    W=H*sech(gamma*(x-x0))**(2)

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh


def redo_mesh(options):

    base_path = options.get("data_path", ".")
    
    # Try to find files in relative path if absolute fails, or just assume relative structure
    # For now, we will assume the user provides a 'data_path' in options if they use this test.
    # If not provided, we can default to a 'data' folder or similar.
    
    node_coords_file = os.path.join(base_path, "NodeCoordsWGhosts.txt")
    
    if not os.path.exists(node_coords_file):
        print(f"Warning: Data file {node_coords_file} not found. Using random/generated data or skipping.")
        # Fallback to a simple grid if files are missing, or raise error
        raise FileNotFoundError(f"Could not find required data file: {node_coords_file}. Please set 'data_path' in options.")

    x0,y0=grids.load_grid(node_coords_file)

    x0=x0[0:119149]
    y0=y0[0:119149]

    mesh0=np.array([x0, y0]).T
    
    bath_file = os.path.join(base_path, "BathymetryInterp.txt")
    disch_file = os.path.join(base_path, "DischargeInterp.txt")
    wl_file = os.path.join(base_path, "WaterLevelInterp.txt")

    B0=np.loadtxt(bath_file,dtype=float,skiprows=1)[0:119149]
    HUHV0=np.loadtxt(disch_file,dtype=float,skiprows=1)[0:119149,:]
    W0=np.loadtxt(wl_file,dtype=float,skiprows=1)[0:119149]

    xmin=np.min(x0)
    xmax=np.max(x0)
    ymin=np.min(y0)
    ymax=np.max(y0)

    nxnyu=options["nxny"]

    x,y=grids.rectangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    mesh=np.array([x, y]).T

    B =np.asarray(sp.interpolate.griddata(mesh0.get(),B0.get(),mesh.get(),'nearest',B0.get().min(),True))
    HU=np.asarray(sp.interpolate.griddata(mesh0.get(),HUHV0[:,0].get(),mesh.get(),'nearest',HUHV0[:,0].get().min(),True))
    HV=np.asarray(sp.interpolate.griddata(mesh0.get(),HUHV0[:,1].get(),mesh.get(),'nearest',HUHV0[:,1].get().min(),True))
    W =np.asarray(sp.interpolate.griddata(mesh0.get(),W0.get(),mesh.get(),'nearest',W0.get().mean(),True))

    HUHV=np.array([HU,HV]).T

    return x,y,B,HUHV,W,mesh


def wavesplit_bump(options):
    
    xmax = 150
    xmin = 0
    ymax = 5
    ymin = 0

    H=0.25
    d=1
    x0=(xmax-xmin)/10

    xL = 50
    xR = 51

    gamma = np.sqrt( 0.75 * H/d**3)

    nxnyu=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    B=np.zeros_like(x)-(d - np.where((xL<x)&(x<xR),0.5*d,0))

    U=np.zeros_like(x)
    V=np.zeros_like(x)

    W=H*sech(gamma*(x-x0))**(2)

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh