from backend import np
import grids
import numpy as op
import scipy as sp

rng = np.random.default_rng()

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

def dambreak_channel2shoebox(options):
    #WORK IN PROGRESS
    """
    Generates mesh of a channel which discharges into a box from the middle of the left side. ny dictates the number of vertical divisions for the channel, e.g., channel_width=4 and ny=4 means dy=1 for the whole mesh.

                                                            shoebox_length
                                                <----------------------------------->
                                                +-----------------------------------+
                                                |                                   | ^
                                                |                                   | |      
                          channel_length        |                                   | |
                   <--------------------------->|                                   | |
                   +----------------------------+                                   | |
                 ^ |  ~       ~       ~       ~ {                                   | |
                 | |      ~       ~       ~     {                                   | |
   channel_width | |  ~     ~  W0=h  ~      ~   {               W0=0                | | shoebox_width
                 | |    ~       ~       ~       {                                   | | 
                 | |  ~      ~      ~      ~    {                                   | | 
                 V |      ~     ~        ~    ~ {                                   | |
                   +----------------------------+                                   | |
                                                |                                   | |
    y ^                                         |                                   | |
      |                                         |                                   | |
      +-->                                      |                                   | v
          x                                     +-----------------------------------+
    """

    #Geometry parameters
    #To enforece correct gluing of channel and shoebox, shoebox_width should be a k-multiple of channel_width such that ny*k is an integer of the same parity as ny.
    #I.e., if the channel is Lc units wide and is segmented into ny parts vertically, the shoebox's width must be Ls = Lc+Lc*n/ny, with n even. 
    channel_length = 5
    channel_width = 2
    shoebox_length = 10
    shoebox_width = 10

    #Check and correction
    #We solve the problem "Find minimum eps>0 for which Ls+eps is a k-multiple of Lc such that ny*k is an integer of the same parity as ny".
    #By looking at the form of Ls from last comment, we can build a strictly increasing sequence of candidate solutions eps(n)=Lc-Ls+Lc*2n/ny. 
    #Then, solving eps(x0)==0 gives (possibly) non integer x0 such that eps(x)>0 iff x>=x0.
    ny=options['nxny'][1]
    n0=np.ceil((shoebox_width/channel_width-1)*ny/2) #Since eps(n) is strictly increasing,  .
    eps = (1+2*n0/ny)*channel_width-shoebox_width
    shoebox_width+=eps
    if eps>0:
        print('Shoebox elongated ', eps, ' units to ensure matching triangles with channel.')

    #Initial conditions
    h = 1

    #Coordinate limits
    xcmin = 0
    xcmax = channel_length
    ycmin = -channel_width/2.0
    ycmax = channel_width/2.0

    xc,yc = grids.triangular_grid(xcmin,xcmax,ycmin,ycmax,None,options['nxny'])

    meshc = np.array([xc,yc]).T

    xcmax=xc.max()

    xsmin = xcmax
    xsmax = xsmin + shoebox_length
    ysmin = -shoebox_width/2.0
    ysmax = shoebox_width/2.0

    nys = shoebox_width/(channel_width/ny)

    xs,ys=grids.triangular_grid(xsmin,xsmax,ysmin,ysmax,None,[10,nys])

    meshs = np.array([xs,ys]).T

    glue_idxc = np.where(meshc[:,0]==xcmax)
    glue_idxs = np.where((meshs[:,0]==xcmax)&(meshs[:,1]>=ycmin)&(meshs[:,1]<=ycmax))

    if len(glue_idxc)!=len(glue_idxs):
        print("[ERROR] Gluing indices don't match!")

    mesh = np.zeros((len(xc)+len(xs)-len(glue_idxc),2))

    candidate = np.vstack((meshc[0:-len(glue_idxc)],meshs))

    meshg = candidate

    x = meshg[:,0]
    y = meshg[:,1]

    B = np.zeros_like(x)

    W = np.where(x<(xcmin+xcmax)/2, h,0)

    HU = np.zeros_like(x)
    HV = HU

    HUHV=np.array([HU,HV]).T

    mesh = {"global_mesh":meshg,"submeshes":(meshc,meshs),"glue_idx":((-1,(glue_idxc),-1,-1),((glue_idxs),-1,-1,-1)),"gluing":((-1,1,-1,-1),(0,-1,-1,-1))}

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

def kurganov_sine(options):
    """
    python run.py --test kurganov_sine --abspath "../SWEpy-tests" --folder kurganov_sine --nx 100 --ny 16 --triangles equilateral --divisions 4 --Tmax 0.5 --tol_dry 1e-6 --g 9.81 --manning 0 --CFL 0.5 --dt_save 0.5 --bconds soft soft periodic periodic
    """
    xmin=-1
    xmax=11
    ymin=-2
    ymax=2

    if "forced_mesh" not in options.keys():
        x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,None,options["nxny"])
    else:
        x=options["forced_mesh"][:,0]
        y=options["forced_mesh"][:,1]

    mesh=np.array([x,y]).T

    B=np.zeros_like(x)

    u=2*np.sin(np.pi*x/5+np.pi/4)
    
    W=(u+10)**2/(4*options["g"])

    HV=np.zeros_like(x)
    HU=(W-B)*u

    HUHV=np.array([HU,HV]).T


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

def wet_dry_showoff(options):

    """
    Builds a noisy rectangular mesh consisting of a flat bottom and a protrusion(/depression? haven't tested if works the other
    way around) in the form of a truncated ellipsoidal cone (flat cap in the form of an ellipse and a straight slope connecting 
    it with a bigger ellipse on the bottom). Water is 

                                   xmin                   x0                             xmax
                                     |                    |                                |
                                     V                    V                                V
                             ymax--->+-----------------------------------------------------+ 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~ | 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\ ~  ~  ~  ~  ~  ~  ~  ~ .~ -~--|     
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\  ~  ~  ~  ~  ~  ~  ~'´¨wet/dry|
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\~  B=-d  ~  ~  ~  ~/ ~ ~slope~~| 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\  ~  ~  ~  ~  ~  ~' ~  ~ ~     | 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\~  ~  ~  ~  ~  ~ '~  ~~   . ---| 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\ ~  ~  ~  ~  ~  ~   ~   /´  <--+--B=h_cap
                                     |¨¨¨¨¨ W=H ¨¨¨¨¨¨¨¨¨¨^\  ~ W=0 ~  ~  ~ |~  ~  |dry cap| 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\~  ~  ~  ~  ~  ~| ~  ~  \      |
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\ ~  ~  ~  ~  ~  ~. ~  ~   ` ---|            
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\  ~  ~  ~  ~  ~  ~. ~  ~~      |             
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\~  ~wet~ ~  ~  ~  ~  ~  ~  ~~~~| 
                                     |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\ ~ bottom ~  ~  ~  ~ .~  ~  ~  | 
    y ^                              |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\  ~  ~  ~  ~  ~  ~  ~ ¨~' ~--~-| 
      |                              |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\~  ~  ~  ~  ~  ~  ~  ~  ~  ~  ~| 
      +-->                           |¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨^\ ~  ~  ~  ~  ~  ~  ~  ~  ~  ~  | 
          x                  ymin--->+-----------------------------------------------------+ 
                                     l_____________________ĵ---> subsequent bore direction                      
                                                T
                                     bore-producing reserve
    """

    #Bounding box limits
    xmax = 100
    xmin = -300
    ymax = 100
    ymin = -100

    #Depth of non-ellipse bottom
    d = 10

    #Ellipse parameters
    a = 1       #semi axis
    b = 4.6     #semi axis
    h = 100     #horizontal center
    k = 0       #vertical center
    r = 100     #pseudo-radius

    #Flats parameters
    h_cap = 1 #Above-water height
    s = 1.5 #Scale between smaller (top) and bigger (foot) ellipses

    #Noise half-range
    w_hedge = 1.5     #Additive noise (Length range in meters around 0)
    w_ledge = 1.5     
    w_slope = 0.5
    w_cap = 0.125
    w_bottom = 1

    #Wave parameters (eastward bore as wet/wet dambreak)
    H = 2
    g = options["g"]
    x0 = -50

    #Grid construction
    nxnyu=options["nxny"]

    x,y=grids.triangular_grid(xmin,xmax,ymin,ymax,nxny=nxnyu)

    #Fuzzy ellipse borders
    r_top_fuzz = r+rng.uniform(-w_hedge,w_hedge,x.size)
    r_bottom_fuzz = s*r+rng.uniform(-w_ledge,w_ledge,x.size)

    #Indices of points in different regions
    points_in_top = a*(x-h)**2+b*(y-k)**2<=(r_top_fuzz)**2
    points_in_slope = (a*(x-h)**2+b*(y-k)**2<=r_bottom_fuzz**2)&(a*(x-h)**2+b*(y-k)**2>=r_top_fuzz**2)
    points_in_bottom = (a*(x-h)**2+b*(y-k)**2>r_bottom_fuzz**2)

    t = (r_bottom_fuzz-np.sqrt((a*(x-h)**2+b*(y-k)**2)))/(r_bottom_fuzz-r_top_fuzz) #linear parameter ranging from 0 at bottom of foot and 1 at cap

    #Fuzzy bathymetry at different regions
    cap = np.where(points_in_top, h_cap+rng.uniform(-w_cap,w_cap,x.size), 0)
    slope = np.where(points_in_slope,(h_cap+d)*(t)-d+(rng.uniform(-w_slope,w_slope,x.size)*(1-t*(1-w_cap/w_slope))),0)
    bottom = np.where(points_in_bottom, -d, 0)+rng.uniform(-w_bottom,w_bottom,x.size)

    #Bathymetry composition
    B=np.zeros_like(x) + bottom + cap + slope

    #Water variables composition
    W = np.where(x<x0,H,0)

    U = np.zeros_like(W)
    V = np.zeros_like(U)

    HU=(W-B)*U
    HV=(W-B)*V

    HUHV=np.array([HU,HV]).T

    mesh=np.array([x, y]).T

    return x,y,B,HUHV,W,mesh