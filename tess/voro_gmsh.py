import numpy as np
import gmsh
import subprocess

from .voro_helper import Timer

def create_grain(dim,conn,verts,gid,lc):
    
    grid=(gid+1)*1000

    i=0
    for vert in verts:
        i+=1
        if not (dim!='3d' and np.abs(vert[2]-1.0)>1.0e-10):
            gmsh.model.geo.addPoint(vert[0],vert[1],vert[2],lc,tag=grid+i)
    
    lines=[]
    face_cents=np.empty((0,3))
    face_tags=np.empty([0],dtype=np.int64)
    edge_cent=[]
    edge_tag=[]
    j=0
    k=0
    for facet in conn:
        n=facet.size
        k+=1
        coord_tot=np.zeros(3)
        edges=[]
        for i in range(n):
            coord_tot=coord_tot+np.array(verts[facet[i]])
            line=[int(grid+facet[(i)%n]+1), int(grid+facet[(i+1)%n]+1)]
            if line in lines:
                edges.append(grid+lines.index(line)+1)
            elif line[::-1] in lines:
                edges.append((grid+lines.index(line[::-1])+1)*-1)
            else:
                j+=1
                edge_tag.append(gmsh.model.geo.addLine(line[0], line[1], grid+j))
                edge_cent.append(np.array(0.5*(np.array(verts[facet[i]])+np.array(verts[facet[(i+1)%n]]))))
                # gmsh.model.geo.mesh.setTransfiniteCurve(grid+j, 5)
                lines.append(line)
                edges.append(grid+j)
        coord_tot=coord_tot/n
        face_cents=np.append(face_cents,np.array([coord_tot]),axis=0)

        gmsh.model.geo.addCurveLoop(edges,grid+k)
        face_id=grid+k if dim=='3d' else gid+1
        face_tags=np.append(face_tags,gmsh.model.geo.addPlaneSurface([grid+k],face_id))
    
    if dim=='3d':
        gmsh.model.geo.addSurfaceLoop(face_tags,gid+1)
        gmsh.model.geo.addVolume([gid+1],gid+1)

    edge_cents=np.array(edge_cent)
    if not dim=='3d':
        edge_cents[:,2]=0.0
    edge_tags=np.array(edge_tag,dtype=np.int64)

    return [face_tags,face_cents,edge_tags,edge_cents]

def create_grains(voro,lc):
    
    dim=3 if voro.dim=='3d' else 2
    conns=voro.conns
    vertss=voro.verts
    gids=voro.gids

    noGrains=len(gids)
    face_cents=np.empty([0,3])
    face_tags=np.empty([0],dtype=np.int64)
    edge_cents=np.empty((0,3))
    edge_tags=np.empty([0],dtype=np.int64)

    for j in range(noGrains):
        grain_faces=create_grain(voro.dim,conns[j], vertss[j], gids[j],lc)
        face_cents=np.vstack((face_cents,grain_faces[1]))
        face_tags=np.hstack((face_tags,grain_faces[0]))
        edge_cents=np.vstack((edge_cents,grain_faces[3]))
        edge_tags=np.hstack((edge_tags,grain_faces[2]))
    
    gmsh.model.geo.removeAllDuplicates()
    gmsh.model.geo.synchronize()

    for j in range(noGrains):
        gmsh.model.addPhysicalGroup(dim, [gids[j]+1])
        gmsh.model.setPhysicalName(dim,gids[j]+1,'grain-'+str(gids[j]+1))

    return [face_tags,face_cents,edge_tags,edge_cents]

def find_hull(voro,faces):

    face=faces[0]
    face_cent=faces[1]
    edge=faces[2]
    edge_cent=faces[3]

    hull=[]
    for hulls in voro.hull:
        side=[]
        side_tag=[]
        if not -1 in hulls:
            for sides in voro.faces[hulls]:
                ind=[]
                for s in sides:
                    if voro.dim=='3d':
                        ind.append(find_match(s,face_cent))
                    else:
                        ind.append(find_match(s,edge_cent))

                side.append(ind)

                if voro.dim=='3d':
                    side_tag.append([face[ind[0]],face[ind[1]]])
                else:
                    side_tag.append([edge[ind[0]],edge[ind[1]]])
                
        hull.append([side,side_tag])

    return hull

def find_match(a,b,tol=1.0-5):
    ind=-1
    for i in range(b.shape[0]):
        if (np.sum(np.abs(a-b[i])))<1.0e-5:
            ind=i
            
    return ind

def make_periodic(voro,hull):

    hull_dist=voro.hull_dist

    dim=2 if voro.dim=='3d' else 1

    for i in range(len(hull)):
        drxn=hull[i]
        left=[s[0] for s in drxn[1]]
        right=[s[1] for s in drxn[1]]
        for j in range(len(left)):
            if left[j]>0:
                translation = [1, 0, 0, hull_dist[i][j][0],
                               0, 1, 0, hull_dist[i][j][1],
                               0, 0, 1, hull_dist[i][j][2],
                               0, 0, 0, 1]
                gmsh.model.mesh.setPeriodic(dim, [left[j]], [right[j]], translation)

def get_allpairs(pair,all_pairs):
    
    p=pair[0:2]

    lst_a=[[el[(el.index(p[0])+1)%2],p[1]] for el in all_pairs if p[0] in el]
    lst_b=[[el[(el.index(p[1])+1)%2],p[0]] for el in all_pairs if p[1] in el]

    chk=p in all_pairs or p[::-1] in all_pairs

    if not chk:
        if len(lst_a)>0:
            for a in lst_a:
                if not (a in all_pairs or a[::-1] in all_pairs):
                    all_pairs.append(a)
                else:
                    chk=False
        if len(lst_b)>0: 
            for b in lst_b:
                if not (b in all_pairs or b[::-1] in all_pairs):
                    all_pairs.append(b)
                else: 
                    chk=False
    
    if not chk: 
        all_pairs.append(p)
        par=pair
    else:
        par=[-1,-1,[0,0,0]]

    return [par,all_pairs]

def get_leftright(pair,left,right):
    
    par=[-1,-1,[0,0,0]]
    if pair[0] in left:
        if not pair[1] in right:
            if not (pair[0] in right or pair[1] in left):
                par[0]=pair[1]
                par[1]=pair[0]
                par[2]=pair[2]*-1.0
                left.append(par[0])
                right.append(par[1])
    else:
        par=pair
        left.append(par[0])
        right.append(par[1])
        
    return [par,left,right]

def get_unique_pairs(pairs):
    
    left=[]
    right=[]
    u_pairs=[]
    all_pairs=[]
    for pair in pairs:
        [par,left,right]=get_leftright(pair,left,right)
        if par[0]>0:
            [par,all_pairs]=get_allpairs(par,all_pairs)
            if par[0]>0:
                u_pairs.append([par[0],par[1],par[2]])
    return u_pairs

def get_periodic(dim,hull):

    pairs=[]
    for i in range(len(hull)):
        drxn=hull[i]
        left=[s[0] for s in drxn[1]]
        for j in range(len(left)):
            [ps,n1,n2,trf]=gmsh.model.mesh.getPeriodicNodes(dim-1,left[j],includeHighOrderNodes=True)
            t=trf[3:12:4]
            for i in range(len(n1)):
                pairs.append([n1[i],n2[i],t])
    
    u_pairs=get_unique_pairs(pairs)

    return u_pairs

def get_insidenode(hull,nodes):

    hl=([h[0] for h in hull])
    hr=([h[1] for h in hull])
    for ni in nodes[0]:
        if not (ni in hl or ni in hr):
            fixed_node=ni
            break

    return fixed_node
 
def gather_mesh(voro,hull):

    nodes=gmsh.model.mesh.getNodes()
    
    dim=3 if voro.dim=='3d' else 2

    elems=[]
    grains=[]
    for e in range(1,voro.n+1):
        elems.append(gmsh.model.mesh.getElements(dim,e))
        grains.append(gmsh.model.getPhysicalName(dim,e))

    pairs=get_periodic(dim,hull)

    inside_node=get_insidenode(pairs,nodes)

    return [nodes,elems,grains,pairs,inside_node]

def set_meshoptions(order,typ,dim):

    options = dict()
    options['Mesh.CharacteristicLengthFromPoints'] = 1
    options['Mesh.Algorithm3D'] = 2
    options['Mesh.CharacteristicLengthExtendFromBoundary'] = 1
    options['Mesh.Smoothing'] = 2
    options['Mesh.ElementOrder'] = order
    options['Mesh.SecondOrderLinear'] = 1
    options['Mesh.OptimizeNetgen'] = 1
    options['Mesh.OptimizeThreshold'] = 0.8
    options['Geometry.Tolerance']=1e-5
    options['General.Verbosity'] = 1
    if typ=='hex':
        options['Mesh.RecombineAll'] = 1
        options['Mesh.Recombine3DAll'] = 1
        options['Mesh.SecondOrderIncomplete']= 1
        options['Mesh.Recombine3DLevel']= 0
        options['Mesh.SubdivisionAlgorithm']=2

    for mo, val in iter(options.items()):
        gmsh.option.setNumber(mo, val)

def create_mesh(voro,mesh_input):
    
    lc=mesh_input.get('mesh_size',0.2)
    order=mesh_input.get('element_order',1)
    typ=mesh_input.get('element_type','tet')
    save_mesh=mesh_input.get('save_mesh',False)
    show_ui=mesh_input.get('show_ui',False)

    t=Timer()
    t.start()
    
    if not gmsh.is_initialized():
        gmsh.initialize()

    gmsh.clear()

    set_meshoptions(order,typ,voro.dim)

    faces=create_grains(voro,lc)
    hull=find_hull(voro,faces)
    make_periodic(voro,hull)

    # print(f"\n Preparing mesh took {t.report():0.4f} seconds")

    gmsh.model.mesh.generate(3 if voro.dim=='3d' else 2)

    [nodes,elems,grains,pairs,inside_node]=gather_mesh(voro,hull)
    
    t.stop('Elapsed time in meshing')
    
    if show_ui: 
        # gmsh.write(voro.fname + ".msh")
        # subprocess.run(["gmsh",voro.fname + ".msh"])
        gmsh.fltk.run()

    if save_mesh: gmsh.write(voro.fname + ".msh")

    gmsh.finalize()
    
    return [pairs,nodes,elems,grains,inside_node,[typ,order]]
