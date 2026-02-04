from backend import np



def bcond_develop(left_indices,right_indices,top_indices,bottom_indices,mesh,nodes,nelems,bconds,rect=True):
    """
    Develops boundary conditions and ghost cells for the mesh.

    Args:
        left_indices (array): Indices of elements on the left boundary.
        right_indices (array): Indices of elements on the right boundary.
        top_indices (array): Indices of elements on the top boundary.
        bottom_indices (array): Indices of elements on the bottom boundary.
        mesh (array): Array of node coordinates (N, 2).
        nodes (array): Array of element connectivity (nelems, 3).
        nelems (int): Number of elements in the mesh.
        bconds (list): List of boundary condition strings [left, right, top, bottom].
        rect (bool, optional): Flag indicating if the mesh is rectangular (True) or equilateral (False). Defaults to True.

    Returns:
        tuple: (ghosts, ghost_nodes, ghost_neighs, ghost_neighsides, left_neighs, right_neighs, top_neighs, bot_neighs, left_neighsides, right_neighsides, top_neighsides, bot_neighsides, newnodes)
            - ghosts: Array of ghost elements.
            - ghost_nodes: Connectivity for ghost elements.
            - ghost_neighs: Neighbors for ghost elements.
            - ghost_neighsides: Neighbor sides for ghost elements.
            - left_neighs, right_neighs, top_neighs, bot_neighs: Updated boundary neighbors.
            - left_neighsides, right_neighsides, top_neighsides, bot_neighsides: Updated boundary neighbor sides.
            - newnodes: Coordinates of newly created ghost nodes.
    """
    
    nleft=len(left_indices)
    nright=len(right_indices)
    ntop=len(top_indices)
    nbot=len(bottom_indices)
    
    nghosts_left=len(left_indices) if bconds[0] in ["soft","wall"] else 0
    nghosts_right=len(right_indices) if bconds[1] in ["soft","wall"] else 0
    nghosts_top=len(top_indices) if bconds[2] in ["soft","wall"] else 0
    nghosts_bot=len(bottom_indices) if bconds[3] in ["soft","wall"] else 0

    ghosts_left=np.zeros((nghosts_left,3),dtype=int)
    ghosts_right=np.zeros((nghosts_right,3),dtype=int)
    ghosts_top=np.zeros((nghosts_top,3),dtype=int)
    ghosts_bot=np.zeros((nghosts_bot,3),dtype=int)

    nodes_ghosts_left=np.zeros((nghosts_left,3),dtype=int)
    nodes_ghosts_right=np.zeros((nghosts_right,3),dtype=int)
    nodes_ghosts_top=np.zeros((nghosts_top,3),dtype=int)
    nodes_ghosts_bot=np.zeros((nghosts_bot,3),dtype=int)

    newnelems=nelems
    newncoords=mesh.shape[0]

    if bconds[0] in ["soft","wall"]:
        oldnelems=nelems
        newnelems=oldnelems+nghosts_left
        oldncoords=mesh.shape[0]
        newncoords=oldncoords+nghosts_left
        ghosts_left[:,0]=np.array(range(oldnelems,newnelems),dtype=int)
        ghosts_left[:,1]=left_indices
        ghosts_left[:,2]=1 if bconds[0]=="soft" else -1
        ghosts_left_nodes=np.array([nodes[left_indices][:,2],nodes[left_indices][:,0],np.array(range(oldncoords,newncoords))],dtype=int).T
        ghosts_left_sides=np.tile(np.array([2,-1,-1],dtype=int),(nghosts_left,1))
        left_neighsides=np.repeat(np.array([0],dtype=int),nghosts_left)
    else:
        ghosts_left=np.array([],dtype=int).reshape(0,3)
        ghosts_left_nodes=np.array([],dtype=int).reshape(0,3)
        ghosts_left_sides=np.array([],dtype=int).reshape(0,3)

    if bconds[1] in ["soft","wall"]:
        oldnelems=newnelems
        newnelems=oldnelems+nghosts_right
        oldncoords=newncoords
        newncoords=oldncoords+nghosts_right
        ghosts_right[:,0]=np.array(range(oldnelems,newnelems),dtype=int)
        ghosts_right[:,1]=right_indices
        ghosts_right[:,2]=1 if bconds[1]=="soft" else -1
        ghosts_right_nodes=np.array([nodes[right_indices][:,0],nodes[right_indices][:,1],np.array(range(oldncoords,newncoords))],dtype=int).T if rect else np.array([nodes[right_indices][:,1],nodes[right_indices][:,2],np.array(range(oldncoords,newncoords))],dtype=int).T
        ghosts_right_sides=np.tile(np.array([0,-1,-1],dtype=int),(nghosts_right,1)) if rect else np.tile(np.array([1,-1,-1],dtype=int),(nghosts_right,1))
        right_neighsides=np.repeat(np.array([0],dtype=int),nghosts_right)
    else:
        ghosts_right=np.array([],dtype=int).reshape(0,3)
        ghosts_right_nodes=np.array([],dtype=int).reshape(0,3)
        ghosts_right_sides=np.array([],dtype=int).reshape(0,3)
    

    if bconds[2] in ["soft","wall"]:
        oldnelems=newnelems
        newnelems=oldnelems+nghosts_top
        oldncoords=newncoords
        newncoords=oldncoords+nghosts_top
        ghosts_top[:,0]=np.array(range(oldnelems,newnelems),dtype=int)
        ghosts_top[:,1]=top_indices
        ghosts_top[:,2]=1 if bconds[2]=="soft" else -1
        ghosts_top_nodes=np.array([nodes[top_indices][:,1],nodes[top_indices][:,2],np.array(range(oldncoords,newncoords))],dtype=int).T
        ghosts_top_sides=np.tile(np.array([1,-1,-1],dtype=int),(nghosts_top,1))
        if not rect:
            top_indices_evens=top_indices[np.array(range(0,len(top_indices),2))]
            ghosts_top_nodes[np.array(range(0,len(top_indices),2))]=np.array([nodes[top_indices_evens][:,2],nodes[top_indices_evens][:,0],np.array(range(oldncoords,newncoords))[top_indices_evens]],dtype=int).T
            ghosts_top_sides[np.array(range(0,len(top_indices),2))]=np.tile(np.array([2,-1,-1],dtype=int),(len(ghosts_top_sides[np.array(range(0,len(top_indices),2))]),1))

        top_neighsides=np.repeat(np.array([0],dtype=int),nghosts_top)
    else:
        ghosts_top=np.array([],dtype=int).reshape(0,3)
        ghosts_top_nodes=np.array([],dtype=int).reshape(0,3)
        ghosts_top_sides=np.array([],dtype=int).reshape(0,3)


    if bconds[3] in ["soft","wall"]:
        oldnelems=newnelems
        newnelems=oldnelems+nghosts_bot
        oldncoords=newncoords
        newncoords=oldncoords+nghosts_bot
        ghosts_bot[:,0]=np.array(range(oldnelems,newnelems),dtype=int)
        ghosts_bot[:,1]=bottom_indices
        ghosts_bot[:,2]=1 if bconds[3]=="soft" else -1
        ghosts_bot_nodes=np.array([nodes[bottom_indices][:,0],nodes[bottom_indices][:,1],np.array(range(oldncoords,newncoords))],dtype=int).T
        ghosts_bot_sides=np.tile(np.array([0,-1,-1],dtype=int),(nghosts_bot,1))
        bot_neighsides=np.repeat(np.array([0],dtype=int),nghosts_bot)
    else:
        ghosts_bot=np.array([],dtype=int).reshape(0,3)
        ghosts_bot_nodes=np.array([],dtype=int).reshape(0,3)
        ghosts_bot_sides=np.array([],dtype=int).reshape(0,3)



    ghosts=np.vstack((ghosts_left,ghosts_right,ghosts_top,ghosts_bot))
    ghost_nodes=np.vstack((ghosts_left_nodes,ghosts_right_nodes,ghosts_top_nodes,ghosts_bot_nodes))
    ghost_neighs=np.hstack((ghosts[:,1].reshape(len(ghosts),1),np.tile(np.array([-1,-1]),(len(ghosts),1))))
    ghost_neighsides=np.vstack((ghosts_left_sides,ghosts_right_sides,ghosts_top_sides,ghosts_bot_sides))

    left_neighs=ghosts_left[:,0]
    right_neighs=ghosts_right[:,0]
    top_neighs=ghosts_top[:,0]
    bot_neighs=ghosts_bot[:,0]

    x1=mesh[ghost_nodes[:,0]][:,0]
    y1=mesh[ghost_nodes[:,0]][:,1]
    x2=mesh[ghost_nodes[:,1]][:,0]
    y2=mesh[ghost_nodes[:,1]][:,1]

    newnodes=np.array([(-np.sqrt(3)/2)*(y1-y2)+(x1+x2)/2,(-np.sqrt(3)/2)*(x2-x1)+(y1+y2)/2]).T

    if bconds[0:2]==["periodic","periodic"]:
        left_neighs=right_indices
        right_neighs=left_indices
        left_neighsides=np.repeat(np.array([0],dtype=int),nleft)
        right_neighsides=np.repeat(np.array([2],dtype=int),nright)

    if bconds[2:4]==["periodic","periodic"]:
        top_neighs=bottom_indices
        bot_neighs=top_indices
        top_neighsides=np.repeat(np.array([0],dtype=int),ntop)
        bot_neighsides=np.repeat(np.array([1],dtype=int),nbot)

    return ghosts,ghost_nodes,ghost_neighs,ghost_neighsides,left_neighs,right_neighs,top_neighs,bot_neighs,left_neighsides,right_neighsides,top_neighsides,bot_neighsides,newnodes



def rectriangles(mesh,bconds):
    """
    Assembles a mesh using rectangular-style triangulation.
    
    Args:
        mesh (array): Array of node coordinates (N, 2).
        bconds (list): List of boundary condition strings [left, right, top, bottom].
        
    Returns:
        tuple: (nodes, neighs, neighsides, ghosts, newmesh, ni, bottom_indices, top_indices)
    """

    xmax=np.max(mesh[:,0])
    xmin=np.min(mesh[:,0])
    ymax=np.max(mesh[:,1])
    ymin=np.min(mesh[:,1])

    topp_indices=np.where(mesh[:,1]==ymax)[0]
    leftt_indices=np.where(mesh[:,0]==xmin)[0]

    nx=len(topp_indices)-1
    ny=len(leftt_indices)-1
    ni=2*nx*ny

    nodes=np.zeros((ni,3),dtype=int)

    evens=np.array(range(0,ni,2))
    odds=np.array(range(1,ni,2))

    botleft=evens/2+np.floor(evens/2/nx)
    botright=botleft+1
    topleft=np.floor(evens/2/nx)+(nx+1)+evens/2

    nodes[evens]=np.array([botleft,botright,topleft],dtype=int).T

    botright=(odds-1)/2+np.floor((odds-1)/2/nx)+1
    topleft=(odds-1)/2+nx+1+np.floor((odds-1)/2/nx)
    topright=topleft+1

    nodes[odds]=np.array([botright,topright,topleft],dtype=int).T

    neighs=np.zeros_like(nodes,dtype=int)

    left=evens-1
    right=evens+1
    down=evens-2*nx+1

    neighs[evens]=np.array([down,right,left],dtype=int).T

    left=odds-1
    right=odds+1
    top=odds+2*nx-1

    neighs[odds]=np.array([right,top,left],dtype=int).T

    neighsides=np.zeros_like(neighs,dtype=int)
    neighsides[evens]=np.tile(np.array([1,2,0],dtype=int),(len(evens),1))
    neighsides[odds]=np.tile(np.array([2,0,1],dtype=int),(len(odds),1))

    left_indices=np.array(range(0,2*nx*ny,2*nx),dtype=int)
    right_indices=left_indices+2*nx-1
    bottom_indices=np.array(range(0,2*nx,2),dtype=int)
    top_indices=bottom_indices+2*nx*(ny-1)+1

    ghosts,ghost_nodes,ghost_neighs,ghost_neighsides,left_neighs,right_neighs,top_neighs,bot_neighs,left_neighsides,right_neighsides,top_neighsides,bot_neighsides,newnodes=bcond_develop(left_indices,right_indices,top_indices,bottom_indices,mesh,nodes,ni,bconds)

    neighs[left_indices]=np.hstack((neighs[left_indices][:,0:2],left_neighs.reshape(len(left_neighs),1)))
    neighs[right_indices]=np.hstack((right_neighs.reshape(len(right_neighs),1),neighs[right_indices][:,1:3]))
    neighs[top_indices]=np.vstack((neighs[top_indices][:,0],top_neighs,neighs[top_indices][:,2])).T
    neighs[bottom_indices]=np.hstack((bot_neighs[:,None],neighs[bottom_indices][:,1:3]))

    neighsides[left_indices]=np.hstack((neighsides[left_indices][:,0:2],left_neighsides.reshape(len(left_neighsides),1)))
    neighsides[right_indices]=np.hstack((right_neighsides.reshape(len(right_neighsides),1),neighsides[right_indices][:,1:3]))
    neighsides[top_indices]=np.vstack((neighsides[top_indices][:,0],top_neighsides,neighsides[top_indices][:,2])).T
    neighsides[bottom_indices]=np.hstack((bot_neighsides[:,None],neighsides[bottom_indices][:,1:3]))

    nodes=np.vstack((nodes,ghost_nodes))
    neighs=np.vstack((neighs,ghost_neighs))
    neighsides=np.vstack((neighsides,ghost_neighsides))

    newmesh=np.vstack((mesh,newnodes))

    return nodes, neighs, neighsides, ghosts, newmesh, ni, bottom_indices, top_indices

def eqtriangles(mesh,bconds):
    """
    Assembles a mesh using equilateral-style triangulation.
    
    Args:
        mesh (array): Array of node coordinates (N, 2).
        bconds (list): List of boundary condition strings [left, right, top, bottom].
        
    Returns:
        tuple: (nodes, neighs, neighsides, ghosts, newmesh, ni, bottom_indices, top_indices)
    """

    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))

    mesh=mesh[ind]
    ind = np.lexsort(np.array([mesh[:,1],mesh[:,0]]))
    mesh=mesh[ind]

    xmax=np.max(mesh[:,0])
    xmin=np.min(mesh[:,0])
    ymax=np.max(mesh[:,1])
    ymin=np.min(mesh[:,1])

    topp_indices=np.where(mesh[:,1]==ymax)[0]
    leftt_indices=np.where(mesh[:,0]==xmin)[0]

    ny=len(leftt_indices)-1

    ncols=len(topp_indices)-1
    nrows=2*len(leftt_indices)-1
    ni=ncols*nrows

    nodes=np.zeros((ni,3),dtype=int)

    evens=np.array(range(0,ni,2))
    odds=np.array(range(1,ni,2))

    left=evens/2+evens//(2*ny+1)
    botright=[]
    j=ny+1
    limit=[ny+1,ny]
    count=0
    for i in range(len(evens)):
        botright.append(j)
        j+=1
        count+=1
        if count>=limit[0]:
            j+=1
            limit.reverse()
            count=0
    botright=np.array(botright)
    topright=botright+1

    nodes[evens]=np.array([left,botright,topright],dtype=int).T

    right=[]
    j=ny+2
    limit=[ny,ny+1]
    count=0
    for i in range(len(odds)):
        right.append(j)
        j+=1
        count+=1
        if count>=limit[0]:
            j+=1
            limit.reverse()
            count=0
    right=np.array(right)
    botleft=odds//2+odds//(2*ny+1)
    topleft=botleft+1

    nodes[odds]=np.array([botleft,right,topleft],dtype=int).T

    neighs=np.zeros_like(nodes,dtype=int)

    top=evens+1
    bot=evens-1
    right=evens+nrows

    neighs[evens]=np.array([bot,right,top]).T

    bot=odds-1
    top=odds+1
    left=odds-nrows

    neighs[odds]=np.array([bot,top,left]).T

    neighsides=np.zeros_like(neighs,dtype=int)
    neighsides[evens]=np.tile(np.array([1,2,0],dtype=int),(len(evens),1))
    neighsides[odds]=np.tile(np.array([2,0,1],dtype=int),(len(odds),1))

    left_indices=np.array(range(1,nrows,2),dtype=int)
    right_indices=left_indices+nrows*(ncols-1)
    bottom_indices=np.array(range(0,ni,nrows),dtype=int)
    top_indices=bottom_indices+nrows-1

    ghosts,ghost_nodes,ghost_neighs,ghost_neighsides,left_neighs,right_neighs,top_neighs,bot_neighs,left_neighsides,right_neighsides,top_neighsides,bot_neighsides,newnodes=bcond_develop(left_indices,right_indices,top_indices,bottom_indices,mesh,nodes,ni,bconds,False)

    neighs[left_indices]=np.hstack((neighs[left_indices][:,0:2],left_neighs.reshape(len(left_neighs),1)))
    neighs[right_indices]=np.vstack((neighs[right_indices][:,0],right_neighs,neighs[right_indices][:,2])).T

    top_indices_evens=top_indices[np.array(range(0,len(top_indices),2))]
    top_indices_odds=top_indices[np.array(range(1,len(top_indices),2))]

    top_neighs_evens=top_neighs[np.array(range(0,len(top_neighs),2))]
    top_neighs_odds=top_neighs[np.array(range(1,len(top_neighs),2))]

    neighs[top_indices_evens]=np.hstack((neighs[top_indices_evens][:,0:2],top_neighs_evens.reshape(len(top_neighs_evens),1)))
    neighs[top_indices_odds]=np.vstack((neighs[top_indices_odds][:,0],top_neighs_odds,neighs[top_indices_odds][:,2])).T

    neighs[bottom_indices]=np.hstack((bot_neighs[:,None],neighs[bottom_indices][:,1:3]))

    neighsides[left_indices]=np.hstack((neighsides[left_indices][:,0:2],left_neighsides.reshape(len(left_neighsides),1)))
    neighsides[right_indices]=np.vstack((neighsides[right_indices][:,0],right_neighsides,neighsides[right_indices][:,2])).T
    
    top_neighsides_evens=top_neighsides[np.array(range(0,len(top_neighsides),2))]
    top_neighsides_odds =top_neighsides[np.array(range(1,len(top_neighsides),2))]

    neighsides[top_indices_evens]=np.hstack((neighsides[top_indices_evens][:,0:2],top_neighsides_evens.reshape(len(top_neighsides_evens),1)))
    neighsides[top_indices_odds] =np.vstack((neighsides[top_indices_odds][:,0],   top_neighsides_odds,neighsides[top_indices_odds][:,2])).T

    neighsides[bottom_indices]=np.hstack((bot_neighsides[:,None],neighsides[bottom_indices][:,1:3]))

    nodes=np.vstack((nodes,ghost_nodes))
    neighs=np.vstack((neighs,ghost_neighs))
    neighsides=np.vstack((neighsides,ghost_neighsides))

    newmesh=np.vstack((mesh,newnodes))

    return nodes, neighs, neighsides, ghosts, newmesh, ni, bottom_indices, top_indices