import os
import cupy as np

import filesave as fs
import triangles as tri
import rhetoric as rh
import tests as tst
from grids import nested_grid

#TO-DO: options->argument parser shit
options={"nxny":[1,1],
         "folder":"equi_test",
         "Tmax":0.07,
         "tol_dry":0.000001,
         "g":1,
         "manning":0,
         "CFL":0.25,
         "dt_save":0.01,
         "bconds":["wall","wall","wall","wall"], #left,right,top,bottom (west, east, north, south)
         "divisions":3,
         "triangles":"equilateral",
         "nestings":3,
         "test":"circular_dambreak_parabolic"}

def execute(options):
    #TO-DO: tests->test selection via argument parser shit
    mybench = getattr(tst,options["test"])

    print("Creating mesh for test "+options["test"]+"...")

    if options["triangles"].lower() in "rectangular":
        from triangles import rectriangles as assembly
    else:
        from triangles import eqtriangles as assembly

    #get path to save files
    folder_path = os.path.join(os.path.abspath("."), options["folder"])

    #get node coords, bathymetry, velocities, and water level
    x,y,B,HUHV,W,oldmesh=mybench(options)

    nodes, neighs, neighsides, ghosts, newmesh, ni, bottom_indices, top_indices = assembly(oldmesh,options["bconds"])

    newB=np.mean(B[nodes[ghosts[:,1]]],axis=1)
    newW=np.mean(W[nodes[ghosts[:,1]]],axis=1)
    newHUHV=np.mean(HUHV[nodes[ghosts[:,1]]],axis=1)

    B=np.hstack((B,newB))
    W=np.hstack((W,newW))
    HUHV=np.vstack((HUHV,newHUHV))

    fs.save_all(folder_path,B,options,HUHV,neighs,neighsides,nodes,ghosts,newmesh,W)

    if options["divisions"]>0:
        these_indices=np.array(range(ni))
        if options["triangles"].lower() in "rectangular":
            ref_indices=rh.itereta(options["divisions"],options["nxny"][0],options["nxny"][1],these_indices)[-1]
        elif options["triangles"].lower() in "equilateral":
            ref_indices=rh.iterdelta(options["divisions"],options["nxny"][0],options["nxny"][1],these_indices)[-1]
            delete_indices=np.hstack((top_indices,bottom_indices))
            ref_indices=np.delete(ref_indices,delete_indices,axis=0)
            these_indices=np.delete(these_indices,delete_indices,axis=0)
        indices=np.array([these_indices,ref_indices]).T
        fs.refplaces(folder_path,indices)
        print("Saved indices of elements in reference grid to compare! Reference grid is "+str(4**(options["divisions"]))+" times larger.")

    print("Finished grid with "+str(len(nodes))+" elements, and "+str(len(newmesh))+" nodes!")

    if options["nestings"]>0:
        for i in range(1,options["divisions"]+1):
            options["divisions"]-=1
            options["forced_mesh"]=nested_grid(nodes[:ghosts[0][0]],oldmesh,top_indices,bottom_indices)
            folder_path = os.path.join(os.path.abspath("."), options["folder"]+"_sub"+str(i))

            x,y,B,HUHV,W,oldmesh=mybench(options)

            nodes, neighs, neighsides, ghosts, newmesh, ni, bottom_indices, top_indices = assembly(oldmesh,options["bconds"])

            newB=np.mean(B[nodes[ghosts[:,1]]],axis=1)
            newW=np.mean(W[nodes[ghosts[:,1]]],axis=1)
            newHUHV=np.mean(HUHV[nodes[ghosts[:,1]]],axis=1)

            B=np.hstack((B,newB))
            W=np.hstack((W,newW))
            HUHV=np.vstack((HUHV,newHUHV))

            fs.save_all(folder_path,B,options,HUHV,neighs,neighsides,nodes,ghosts,newmesh,W)

            if options["divisions"]>0:
                these_indices=np.array(range(ni))
                if options["triangles"].lower() in "rectangular":
                    ref_indices=rh.itereta(options["divisions"],options["nxny"][0]*2**i,options["nxny"][1]*2**i,these_indices)[-1]
                elif options["triangles"].lower() in "equilateral":
                    ref_indices=rh.iterdelta(options["divisions"],options["nxny"][0]*2**i,options["nxny"][1]*2**i,these_indices)[-1]
                    delete_indices=np.hstack((top_indices,bottom_indices))
                    ref_indices=np.delete(ref_indices,delete_indices,axis=0)
                    these_indices=np.delete(these_indices,delete_indices,axis=0)
                indices=np.array([these_indices,ref_indices]).T
                fs.refplaces(folder_path,indices)
                print("Saved indices of elements in reference grid to compare! Reference grid is "+str(4**(options["divisions"]))+" times larger.")
            
            print("Finished grid with "+str(len(nodes))+" elements, and "+str(len(newmesh))+" nodes!")

    print("PROGRAM ENDED CORRECTLY")

    return



#for method in ["RK3","FE","CFE"]:
#    nxny=35
#    for n in range(5):
#        options["folder"]="brysons_equilateral/brysonx_ne"+str(nxny)+method
#        options["nxny"]=[nxny,nxny]
#        options["divisions"]=4-n
#        execute(options)
#        nxny*=2