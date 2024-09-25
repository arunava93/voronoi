import numpy as np
import subprocess
import os
from sys import platform

from .voro_helper import Timer,Euler
from .voro_gmsh import create_mesh
from .voro_abaq import abaqus_export

class Voro:

    def __init__(self,voro_input):
    
        self.t=Timer()
        
        if type(voro_input)!=str:
            
            rng=np.random.default_rng()
            
            self.n=voro_input.get('no_grains')
            self.dim=voro_input.get('dim','3d')

            self.sizes=np.array(voro_input.get('size',[1,1,1]))
            self.periodicity=np.array(voro_input.get('periodicity',[1,1,1]))
            
            self.dst=voro_input.get('dst',0)
            self.dst_type=voro_input.get('dst_type','normal')
            
            self.lloyd=voro_input.get('lloyd',50)
            
            cuts=np.array(voro_input.get('cuts',[[]]))
            self.cuts=cuts if cuts.ndim==2 else np.array([cuts])

            self.fname=voro_input.get('voro_name','voro')
            
            self.seed=voro_input.get('seed',rng.integers(100000))

            self.euler=np.empty([0,3])
            self.mesh=[]
    
            self.__main()

        else:

            self.t.start()

            self.__read_voro_vor(voro_input)
            self.__gen_data()
            self.__read_voro_out(f_flag=True)

            self.t.stop('Elapsed time in gathering')

    def __main(self):
        
        self.t.start()

        self.__gen_data()
        self.__run()
        self.__read_voro_out()

        self.t.stop('Elapsed time in tesselation')

    def __gen_data(self):

        self.rng=np.random.default_rng(seed=self.seed)

        if self.lloyd<=1000:
            [points,cnt,ind]=self.__get_points_random()
        else:
            [points,cnt,ind]=self.__get_points_regular()

        self.n=cnt
        self.points=points 
        self.ind=ind.T

        self.v=self.sizes[0]*self.sizes[1]*(self.sizes[2] if self.dim=='3d' else 1)

        w=self.__get_weights()
        self.w=w.T

        if self.euler.size==0:
            self.euler=self.__euler(method='random')

    def __get_dist(self,cnt):

        if self.dst_type=='normal':
            w_dst=self.rng.normal(0.0,0.25,cnt)
        elif self.dst_type=='random':
            w_dst=self.rng.random(cnt)-0.5
        else:
            w_dst=np.zeros(cnt)

        w_dst[w_dst<-1]=-1
        w_dst[w_dst>1]=1

        return self.dst*w_dst

    def __get_weights(self):

        w0=(self.v/self.n*3/4/np.pi)**(1/3) if not self.dim=='2d' else (self.v/self.n/np.pi)**(1/2)
        w=np.array([w0*(1+self.__get_dist(self.n))])

        return w

    def __get_points_random(self):

        points=np.zeros([self.n,3])

        for i in range(self.n):
            points[i]=[self.sizes[0]*(self.rng.random()-1/2),
                       self.sizes[1]*(self.rng.random()-1/2),
                       0 if self.dim=='2d' else self.sizes[2]*(self.rng.random()-1/2)]

        cnt=self.n
        ind=np.array([[i+1 for i in range(cnt)]])

        return [points,cnt,ind]

    def __get_points_regular(self):

        if self.dim=='2d':
            n=round(0.5*self.n**(1/2))*2 
        else:
            if self.lloyd==1001:
                n=round(0.5*self.n**(1/3))*2
            else:
                n=round(self.n**(1/3))

        zn=1 if self.dim=='2d' else n

        points=np.zeros([n**2 if self.dim=='2d' else n**3,3])

        m=0
        for j in range(n):
            for i in range(n):
                for k in range(zn):
                    if self.dim=='2d':
                        points[m]=[self.sizes[0]*(-1/2+i/n+(j%2/n/2 if self.lloyd==1001 else 0)),self.sizes[1]*(-1/2+j/n),0]
                    else:
                        points[m]=[self.sizes[0]*(-1/2+i/n+(k%2/n/2 if self.lloyd==1001 else 0)),self.sizes[1]*(-1/2+j/n+(k%2/n/2 if self.lloyd==1001 else 0)),(-1/2+k/n)*self.sizes[2]]
                    m+=1


        cnt=len(points)
        ind=np.array([[i+1 for i in range(cnt)]])

        return [points,cnt,ind]

    def __run(self):
        
        self.__wrt_vin()
        
        if (platform=="win32"):
            exec_name="tess/voro++/Build/Bin/VoroGen/Windows/VoroGen.exe"
        else:
            exec_name="tess/voro++/Build/Bin/VoroGen/Release/VoroGen"
        
        subprocess.run([exec_name,self.fname+".vin",self.fname+".vol"])
 
        self.__del_vin()

    def __read_voro_out(self,f_flag=False):

        self.__read_voro_vol(f_flag)
        self.__read_voro_hull(f_flag)

    def __read_voro_vol(self,f_flag):
        
        ext='.vor' if f_flag else '.vol'
        
        vor_size=10
        n=vor_size+15
        
        with open(self.fname+ext,"r") as myfile:
            liness=myfile.readlines()
            
            if not f_flag:
                lines=liness
            else:
                lines=liness[n:]
            self.__vol_lines=lines

            conns=[]
            verts=[]
            norms=[]
            cents=[]
            gids=[]
            vols=[]

            for line in lines:

                [ind,vol,cent,conn,norm,vert]=self.__parse_vol(line)
                
                if self.dim=='3d':
                    conns.append(conn)
                    gids.append(ind)
                    cents.append(cent)
                    verts.append(vert)
                    norms.append(norm)
                    vols.append(vol)
                else:
                    conns_dum=[]
                    norms_dum=[]
                    for i in range(len(conn)):
                        if np.sum(abs(norm[i]-np.array([0.,0.,1.])))<1.0e-10:
                            conns_dum.append(conn[i])
                            norms_dum.append(norm[i])
                    
                    if not conns_dum==[]:
                        gids.append(ind)
                        cents.append(cent)
                        verts.append(vert)
                        conns.append(conns_dum)
                        norms.append(norms_dum)
                        if self.dim=='2.5d':
                            vols.append(self.__poly_area(np.array(vert)[conns_dum[0]][:,0:2]))
                        else:
                            vols.append(vol)

        self.cents=cents
        self.conns=conns
        self.verts=verts 
        self.norms=norms
        self.gids=gids
        self.vols=np.array(vols)
        
        if not f_flag and os.path.exists(self.fname+ext):
            os.remove(self.fname+ext)

    def __parse_vol(self,line):
        
        buff=line.rstrip().split(' ')
        
        ind=int(buff[0])

        vol=np.double(buff[1])

        cent=[np.double(buff[2]),np.double(buff[3]),np.double(buff[4])]
        
        n=5

        n_faces=int(buff[n])
        n_verts=int(buff[n+1+2*n_faces])
        
        conn=[]
        for buffy in buff[n+1:n+1+n_faces]:
            conn.append(np.fromstring(buffy[1:-1],dtype=int, sep=","))
        
        norm=[]
        for buffy in buff[n+1+n_faces:n+1+2*n_faces]:
            norm.append(np.round(np.fromstring(buffy[1:-1],dtype=np.double, sep=","),15))

        vert=[]
        for buffy in buff[n+2+2*n_faces:]:
            vert.append(np.round(np.fromstring(buffy[1:-1],dtype=np.double, sep=","),15))
        
        return [ind,vol,cent,conn,norm,vert]

    def __read_voro_hull(self,f_flag):
        
        vor_size=10
        ext='.vor' if f_flag else '.vol.hull'

        with open(self.fname+ext,"r") as myfile:
            liness=myfile.readlines()
            
            if not f_flag:
                lines=liness
            else:
                lines=liness[vor_size:vor_size+15]
            self.__hull_lines=lines

            hulls=[]
            hull_dists=[]
            i=0
            for line in lines:
                if (i==0):
                    buff=line.rstrip()[1:-1].split('] [')
                    faces=np.empty([0,3],dtype=float)
                    for buffy in buff:
                        faces=np.vstack((faces,np.fromstring(buffy,dtype=float,sep=' ')))
                elif (i<8):
                    hull=np.empty([0,2],dtype=np.int32)
                    buff=line.rstrip()[1:-1].split('] [')
                    for buffy in buff:
                        hull=np.vstack((hull,np.fromstring(buffy,dtype=np.int32,sep=' ')))
                    hulls.append(hull)
                else:
                    hull_dist=np.empty([0,3],dtype=float)
                    buff=line.rstrip()[1:-1].split('] [')
                    for buffy in buff:
                        hull_dist=np.vstack((hull_dist,np.fromstring(buffy,dtype=float,sep=' ')))
                    hull_dists.append(hull_dist)
                i+=1;
        
        self.faces=faces
        self.hull=hulls[0:3] if self.dim=='2d' else hulls
        self.hull_dist=hull_dists[0:3] if self.dim=='2d' else hull_dists
        self.hull_stacked=np.vstack(self.faces[np.vstack(self.hull)])
        if self.dim=='2d':
            self.hull_stacked[:,2]=self.hull_stacked[:,2]+1.0

        if not f_flag and os.path.exists(self.fname+ext):
            os.remove(self.fname+ext)

    def __wrt_vin_parameters(self):
        
        line=''
        line+=str(self.lloyd)+' '
        line+=np.array2string(self.sizes)[1:-1]+' '
        line+=np.array2string(self.periodicity)[1:-1]+'\n'
        for cut in self.cuts:
            line+=np.array2string(cut)[1:-1]+' '
        line+='\n'
    
        return line

    def __wrt_vin(self):

        with open(self.fname+".vin",'w') as f:
            f.write(self.__wrt_vin_parameters())
            voro_dat=np.hstack((self.ind,self.points,self.w))
            np.savetxt(f,voro_dat,fmt='%d %.12f %.12f %.12f %.12f')

    def __del_vin(self):

        if os.path.exists(self.fname+".vin"):
            os.remove(self.fname+".vin")

    def __edge_cent(self,verts):

        n=len(verts)
        edges=np.empty([0,9])
        for i in range(n):
            edges=np.vstack((edges,np.hstack((0.5*(verts[i]+verts[(i+1)%n]),verts[i],verts[(i+1)%n]))))

        return edges

    def __euler(self,method='random'):
        
        if method=='random':
            eulers=Euler(self.n).angles

        return eulers

    def __wrt_euler(self):

        euler_line=''
        for eul in self.euler:
            euler_line=euler_line+str(eul)+' '

        return euler_line.rstrip()

    def __wrt_cuts(self):

        cuts_line=''
        for cut in self.cuts:
            cuts_line=cuts_line+str(cut)+' '

        return cuts_line.rstrip()

    def __read_euler(self,line):

        euler=np.empty([0,3],dtype=float)
        buff=line.rstrip()[1:-1].split('] [')
        for buffy in buff:
            if buff[0]!='':
                euler=np.vstack((euler,np.fromstring(buffy,dtype=float,sep=' ')))

        return euler

    def __read_cuts(self,line):

        cuts=np.empty([0,4],dtype=float)
        buff=line.rstrip()[1:-1].split('] [')
        for buffy in buff:
            if buff[0]!='':
                cuts=np.vstack((cuts,np.fromstring(buffy,dtype=float,sep=' ')))

        return cuts

    def __write_vol(self):
        
        with open(self.fname+".vor","a") as myfile:
            myfile.write(''.join(self.__vol_lines))

    def __write_hull(self):
        
        with open(self.fname+".vor","a") as myfile:
            myfile.write(''.join(self.__hull_lines))
            
    def __write_vor(self):
        
        lines=[]

        lines.append(str(self.n))
        lines.append(self.dim)

        lines.append(str(self.sizes))
        lines.append(str(self.periodicity))
        
        lines.append(str(self.dst)+' '+self.dst_type)
        lines.append(str(self.lloyd))
        
        lines.append(self.__wrt_cuts())
        
        lines.append(self.fname)
        lines.append(str(self.seed))
        lines.append(self.__wrt_euler())

        with open(self.fname+".vor","w") as f:
            f.write('\n'.join(lines))
            f.write('\n')

    def __read_voro_vor(self,fname):
        
        with open(fname+".vor","r") as myfile:
            lines=myfile.readlines()
        
        self.n=int(lines[0][:-1])
        self.dim=lines[1][:-1]

        self.sizes=np.fromstring(lines[2][1:-1],dtype=float,sep=' ')
        self.periodicity=np.fromstring(lines[3][1:-1],dtype=int,sep=' ')

        self.dst=float(lines[4].rstrip().split(' ')[0])
        self.dst_type=lines[4].rstrip().split(' ')[1]

        self.lloyd=int(lines[5][:-1])

        self.cuts=self.__read_cuts(lines[6])

        self.fname=lines[7][:-1]
        self.seed=int(lines[8][:-1])
        self.euler=self.__read_euler(lines[9])

    def __plotter(self,p,min_xyz,max_xyz):
        
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection

        fig = plt.figure()
        
        ax=fig.add_subplot(111,projection='3d', proj_type = ('ortho' if self.dim!='3d' else 'persp'))
        ax.add_collection3d(p)
        
        ax.set_xlim3d(1.1*min_xyz, 1.1*max_xyz)
        ax.set_ylim3d(1.1*min_xyz, 1.1*max_xyz)
        ax.set_zlim3d(1.1*min_xyz, 1.1*max_xyz)
        
        ax.xaxis.set(pane_color="white")
        ax.yaxis.set(pane_color="white")
        ax.zaxis.set(pane_color="white")
        
        ax.grid(False)
        
        ax.set_box_aspect([1,1,1])
        if self.dim!='3d':
            ax.view_init(360+90,90)

        plt.axis('off')
        
        plt.tight_layout()
        plt.show(block=False)

    def __poly_area(self,verts):

        n=len(verts)

        ref=verts.sum(axis=0)/n
        sorted_verts=sorted(verts,key=lambda x: -self.__angle(x,ref))

        psum=0
        nsum=0
        for i in range(n):
            psum+=verts[i][0]*verts[(i+1)%n][1]
            nsum+=verts[(i+1)%n][0]*verts[i][1]

        return abs(0.5*(psum-nsum))

    def __angle(self,p1,p2):

        k=(p2[1]-p1[1])/np.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2)

        if k>=0:
            return 2*np.pi-np.arcsin(k) if p2[0]>p1[0] else np.pi+np.arcsin(k)
        else:
            return np.arcsin(-k) if p2[0]>p1[0] else np.pi-np.arcsin(-k)
        
    def save(self):
        
        tmp=np.get_printoptions()['precision']
        np.set_printoptions(precision=16)
        
        self.__write_vor()
        self.__write_hull()
        self.__write_vol()
        
        np.set_printoptions(precision=tmp)

    def rerun(self,**kwargs):
        
        self.n=kwargs.get('n',self.n)
        self.dim=kwargs.get('dim',self.dim)
        self.sizes=kwargs.get('sizes',self.sizes)
        self.periodicity=kwargs.get('periodicity',self.periodicity)
        self.dst=kwargs.get('dst',self.dst)
        self.dst_type=kwargs.get('dst_type',self.dst_type)
        self.lloyd=kwargs.get('lloyd',self.lloyd)
        self.cuts=kwargs.get('cuts',self.cuts)
        self.fname=kwargs.get('fname',self.fname)
        self.mesh=[]

        self.__main()

    def create_mesh(self,mesh_input):

        self.mesh=create_mesh(self,mesh_input)

    def abaqus_export(self,abaqus_input):

        if self.mesh!=[]:
            abaqus_export(self,abaqus_input)
        else:
            print('Please create a mesh first!')

    def plot(self,hull=False):
        
        from matplotlib import cm
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        
        verts2plot=[]
        colarr=np.zeros(0)
        col=np.arange(len(self.conns)+1)/len(self.conns)
        np.random.shuffle(col)
        hull_stacked=np.empty([0,3])
        i=0
        for conn in self.conns:
            for facet in conn:
                if not hull:
                    verts2plot.append(np.array(self.verts[i])[facet])
                    colarr=np.append(colarr,col[i])
                else:
                    if self.dim=='3d':
                        facet_cent=np.sum(np.array(self.verts[i])[facet],axis=0)/facet.shape[0]
                        if (np.count_nonzero(np.sum(np.abs(self.hull_stacked-facet_cent),axis=1)<1e-5)>0):
                            verts2plot.append(np.array(self.verts[i])[facet])
                            colarr=np.append(colarr,col[i])
                    else:
                        edge_cent=self.__edge_cent(np.array(self.verts[i])[facet])
                        for ec in edge_cent:
                            if (np.count_nonzero(np.sum(np.abs(self.hull_stacked[:,0:2]-ec[0:2]),axis=1)<1e-5)>0):
                                verts2plot.append(np.vstack((ec[3:6],ec[6:9])))
                                colarr=np.append(colarr,col[i])
            i=i+1

        p=Poly3DCollection(verts2plot, edgecolors='k', linewidths=1, alpha=0.85,cmap=cm.rainbow)
        p.set_array(colarr)
      
        min_xyz=np.vstack(np.array(verts2plot,dtype="object").ravel()).min()
        max_xyz=np.vstack(np.array(verts2plot,dtype="object").ravel()).max()
        
        self.__plotter(p,min_xyz,max_xyz)

    def plot_vols(self):

        import matplotlib.pyplot as plt

        figure=plt.figure()
        ax=plt.axes()

        eq_r=np.power(self.vols*3.0/4.0/np.pi,1/3) if self.dim=='3d' else np.power(self.vols/self.z_size/np.pi,1/2)
        
        ax.hist(self.w,50,density=True,facecolor='gray',alpha=0.5,range=(np.min(eq_r),np.max(eq_r)))
        ax.hist(eq_r,50,density=True,facecolor='dimgray',alpha=0.8)

        ax.set_xlabel('eq radius')
        ax.set_ylabel('probability density')
        ax.set_xlim([np.min(eq_r)*0.8,np.max(eq_r)*1.2])

        plt.show(block=False)