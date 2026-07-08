#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 19:10:42 2023

@author: lara
"""

import json
from SideClasses import SideFunctions as SF
import numpy as np 

import matplotlib 
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from scipy.optimize import LinearConstraint
from scipy.optimize import minimize


import scipy 


class myArrow3D(FancyArrowPatch):
    
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)




class ReadLatticeParticle():
    
    
    def __init__(self, name , directory):
        
        file = open(directory + name +'.json', 'r')
        
        myDict = json.load(file)
        
        self.name = name
        self.directory=directory
        self.generators = myDict['Generators'] 
        self.edges = np.array(myDict['Edges'])
        self.n_edges = myDict['N_Edges']
        
        self.generators = np.array(self.generators)
        
        self.neighbors =np.array(myDict['Neighbors'])
        self.n_neighbors = myDict['N_Neighbors']
        
        
        self.neighbors_cartesian = np.zeros((self.n_neighbors,3))
        for i in range (self.n_neighbors):
            self.neighbors_cartesian[i] = self.lattice_to_cartesian(self.neighbors[i])
            norm = np.linalg.norm(self.neighbors_cartesian[i])
            self.neighbors_cartesian[i] = self.neighbors_cartesian[i]/norm
        
        self.n_orientations = myDict['N_Orientations']
        self.orientations = np.array(myDict['Orientations'])
        
        
        
        """Usually the only unique edge is edge=1, but in the case of the hexagon 
        we also need edge=4 to get all the structures"""
        self.face_type = myDict['Faces_type']
        self.unique_faces = np.unique(self.face_type)
        self.unique_edges = np.array([np.where(self.face_type==unique)[0][0] for unique in self.unique_faces])+1
        
                                
        
        """Total number of structures"""
        self.N_particles_type_max = 0
        self.all_N_structures = {}
        self.detail_all_N_structures = {}
        maxfound = False 
        i = 1
        while not maxfound:
            try :
                this_n_structures = myDict['N_Structures_'+str(i)]
                detail_this_n_structures = myDict['detail_N_Structures_'+str(i)]
            except :
                maxfound=True
            
            if not maxfound : 
                self.all_N_structures[i] = this_n_structures
                self.detail_all_N_structures[i] = detail_this_n_structures
                self.N_particles_type_max +=1
            i+=1
            
        self.read_structure_file()
        self.build_state_dict()
        self.compute_opposite_orientations()
        
        if name=='Triangular':
            self.map_old_to_new_empty_full()
                
    def read_structure_file(self):
        self.config_to_structure_dicts = {} # (o1, o2, edge) -> s
        self.structure_to_config_dicts = {} # s -> (o1, o2, edge=1), o1<o2
        self.structure_degeneracy_dicts = {} #how many configuratoin have the same s
        
        for ntype in range(1, self.N_particles_type_max):
            self.config_to_structure_dicts[ntype] = {}
            self.structure_to_config_dicts[ntype] = {}
            self.structure_degeneracy_dicts[ntype] = {}
        
            all_structures = SF.matrixToArray(self.directory+self.name+'_structures_'+str(ntype)+'particles.txt')
            for a_structure in all_structures :
                o1, o2, edge, s = a_structure
                
                self.config_to_structure_dicts[ntype][(o1, o2, edge)]=s
      
                if s in self.structure_degeneracy_dicts[ntype].keys():
                    self.structure_degeneracy_dicts[ntype][s]+=1  
                else :
                    self.structure_degeneracy_dicts[ntype][s] = 0
                    

                # if self.name=='Triangular' and s==71 and ntype==2:
                #     print('o1,o2,edge', o1, o2, edge, 's',s)
                        
                if (edge in self.unique_edges  or edge==0) :
                    
                    if not s in self.structure_to_config_dicts[ntype].keys():  
                        self.structure_to_config_dicts[ntype][s]=(o1, o2, edge)
                        
    def get_structure(self, n_particles_type,  o1,o2, edge, twoParticlesStructures=True):
        """Returns the indices as saved for python and not for c++, without the initial one particle structures"""
        N1 = self.detail_all_N_structures[n_particles_type][0]
        return(self.config_to_structure_dicts[n_particles_type][(o1, o2, edge)]-N1)
                    
    def get_particle_state(self, particle, orientation):
        if orientation==0:
            return(0)
        else :
            return(particle*self.n_orientations+orientation)
    
    def get_particle_id_orientation(self, particle_state):
        if particle_state ==0 :
            particle_id=0
            orientation=0
        else :
            particle_id = (particle_state-1)//self.n_orientations
            orientation = (particle_state-1)%self.n_orientations+1
        return(particle_id, orientation)
                    
    def build_state_dict(self):
        self.particle_orientation_to_state = {}
        self.state_to_particle_orientation = {}
        for particle in range (self.N_particles_type_max):
            for orientation in range (self.n_orientations+1):
                state = self.get_particle_state(particle, orientation)
                self.particle_orientation_to_state[(particle, orientation)] = state
                self.state_to_particle_orientation[state] = (particle, orientation)
        
        
    def lattice_to_cartesian(self, coordinate):
        """Transforms lattice coordinate to cartesian coordinate"""
        d =  len(self.generators)
        x,y,z = 0,0,0
        for i in range(d):
            dx, dy, dz = coordinate[i] * self.generators[i]
            x+=dx
            y+=dy
            z+=dz
        return(x,y,z)
    
    

    def lattice_resquare(self, coordinate, Lx, Ly):
        """In lattice coordinates"""
        x,y,z = coordinate
        y = y%Ly
        D=  y//2
        
        if x>Lx - D -1 :
            x = x-Lx
        if x<-D:
            x= x + Lx 
        return(np.array([x,y,z]))
    
    def find_opposite_edge_id(self, edge_id):
        """If edge_id is the id of the vector pointing towards the particle
        what is the id of the same vector pointing outwards the particle"""

        if edge_id ==0 :
            return(0)
        else :
            edge = self.neighbors[edge_id-1]
            opposite_edge = -edge
            searchZeros = np.abs(self.neighbors -opposite_edge)
            opposite_edge_id = np.where((searchZeros==0).all(axis=1))[0][0]
            opposite_edge_id +=1
            return(opposite_edge_id)


    def find_face_id(self, orientation, edge_id):
            """If particle is in given orientation, 
            what face is facing the edge number edge_id
            If the orientation is 0, then the face is -1"""
            if orientation==0 :
               face_id = 0
            else :
                face_id = np.where(self.orientations[orientation-1] == edge_id-1)[0][0]+1
            return(face_id)
    
    
    def get_faces_in_contact(self, structure, ntype):
        """Returns particle1, face1, particle2, face2"""        
        state1, state2, edge = self.structure_to_config_dicts[ntype][structure]  
        if edge!=0 :
            if state1!=0 :
                p1, o1 = self.state_to_particle_orientation[state1]
                face1 = self.find_face_id(o1, edge) 
            else :
                p1 = 0
                face1 = 0
            if state2 !=0:
                p2, o2 = self.state_to_particle_orientation[state2]               
                face2 = self.find_face_id(o2, self.find_opposite_edge_id(edge))
            else :
                p2 = 0
                face2 = 0
        else :
            p1, face1, p2, face2 = 0,0,0,0
        
        return(p1, face1, p2, face2)
    
    
    def show_structure(self, s=None, 
                       state1=None, state2=None, edge=None,
                       p1=None, p2=None,o1 = None, o2=None, 
                       ntype=1, axes=None, representation='arrowed', dimension=2, 
                       twoParticlesStructures=False):
        """Structure is either defined by s, by (state1, state2, edge), or by (p1, o1, p2, o2, edge)"""
        
        if s is not None and twoParticlesStructures:
            s = s+self.detail_all_N_structures[ntype][0]
        
        if o1 is not None :
            state1 = self.get_particle_state(p1, o1)
            state2 = self.get_particle_state(p2, o2)
            
            # print(state1, state2)
        
        if state1 is not None :
            s = self.config_to_structure_dicts[ntype][(state1, state2, edge)]
            
        
        state1, state2, edge = self.structure_to_config_dicts[ntype][s]
    
        p1, o1 = self.state_to_particle_orientation[state1]
        p2, o2 = self.state_to_particle_orientation[state2]
            
        
        
        print('p1,p2,o1,o2,edge', p1,p2,o1,o2,edge)
        size=0.95
        Particle = ParticleRepresentation(self)
        x1, y1, z1 = 0, 0, 0
        dx, dy, dz = self.neighbors_cartesian[edge-1]
        x2, y2, z2 = x1+dx, y1+dy, z1+dz
        # print(x1, x2, y1, y2)
        if axes is None :
            f,axes=plt.subplots(figsize=(3,3))
        print(axes, representation)
        Particle.plot_particle(x1, y1,z1,  size, p1, o1, dimension, axes=axes, representation=representation)
        Particle.plot_particle(x2, y2, z2, size, p2, o2, dimension, axes=axes, representation=representation)
        
        xmin, xmax =min(x1,x2), max(x1,x2)
        ymin, ymax =min(y1,y2), max(y1,y2)
        axes.set_xlim(xmin-size, xmax+size)
        axes.set_ylim(ymin-size, ymax+size)
        if twoParticlesStructures:
            stitle = s - self.detail_all_N_structures[ntype][0]
        else :
            stitle = s
        axes.set_title('s='+str(stitle))
        axes.set_aspect('equal', 'box')  
        axes.tick_params(axis='both',  bottom=False, left=False, labelbottom=False, labelleft=False)
        return(axes,s )
    
    
    
    def compute_opposite_orientations(self):
        """Creates a dict that for orientation o1 gives the reverse orientation """
        orientations = np.arange(self.n_orientations+1)
        opposites = - np.ones(self.n_orientations+1)
        for o1 in orientations :
            for o1op in range(o1+1) :
                s1 = self.config_to_structure_dicts[1][o1, 0, 1]
                s2 = self.config_to_structure_dicts[1][0, o1op, 1]
                if s1==s2 :
                    opposites[o1] = o1op
                    opposites[o1op] = o1
        
        self.opposite_orientations = opposites
    
    
    
    def get_matrix_structures(self, p1=0, p2=0, symmetric=True, n_type=1, emptyFullAlso=False):
        n = self.n_orientations + emptyFullAlso
        
        o1s = np.arange(1, n+1)
        if symmetric :
            o2s = self.opposite_orientations.copy()[1:]
        else :
            o2s = o1s
            
        if emptyFullAlso:
            o1s = np.concatenate((np.zeros(1), o1s))
            o2s = np.concatenate((np.zeros(1), o2s))
        
        structures_matrix = np.zeros((n,n), dtype=int)
        for i1 in range (n):
            for i2 in range (n):
                o1, o2 = o1s[i1], o2s[i2]
                state1, state2 = self.particle_orientation_to_state[p1, o1], self.particle_orientation_to_state[p2, o2]
                structures_matrix[i1, i2] = self.config_to_structure_dicts[n_type][state1, state2, 1]
        return(structures_matrix)
    
    
    def map_old_to_new_empty_full(self):
        matrice_structure = self.get_matrix_structures()
        
        faces_order_old = [3,2,1,6,5,4]
        # faces_order_old = [2,3,4,5,6,1]
        face_to_empty_structure_old = {}
        face_to_empty_structure_old[1]=6
        face_to_empty_structure_old[6]=3
        face_to_empty_structure_old[5]=2
        face_to_empty_structure_old[4]=1
        face_to_empty_structure_old[3]=4
        face_to_empty_structure_old[2]=5
        
        face_to_empty_structure_new = {}
        for o in range (1,7):
            s = self.config_to_structure_dicts[1][o,0,1]
            p,f,p,zero = self.get_faces_in_contact(s,1)
            face_to_empty_structure_new[f]= s
        
        self.old_to_new = {}
        self.new_to_old = {}
        
        for i in range(6):
            """which faces are in contact in the ith diagonal item of the ref matrix"""
            s_diag = matrice_structure[i][i]
            p1, face, p2, face = self.get_faces_in_contact(s_diag, 1)
            """To which empty-full strucutre this face correspond"""
            s_empty = face_to_empty_structure_new[face] -7
            
            """in the old matrix, the ith diagonal index correspond to old_face"""
            old_face = faces_order_old[i]
            """To which empty-full strucutre  old_face correspond"""
            old_s_empty = face_to_empty_structure_old[old_face]
            
            self.old_to_new[old_s_empty] = s_empty
            self.new_to_old[s_empty] = old_s_empty
    

       
            



    def full_structures_each_face(self, n_structures_type=1):
        s_faces = np.zeros((self.n_neighbors, self.n_neighbors), dtype=int)
        for s in range(self.all_N_structures[n_structures_type]):
            p1, f1, p2, f2 = self.get_faces_in_contact(s, n_structures_type)
            if f1*f2 !=0 :
                s_faces[f1-1][f2-1] = s
                s_faces[f2-1][f1-1] = s
        return(s_faces)
    
    
    
    def rotate_local_structure(self, local_structure, n):
        
        o, s1, s2, s3, s4, s5, s6 = local_structure
        structures = [s1, s2, s3, s4, s5, s6]
        
        o_new = (o-n-1)%6+1
        structures_new = structures[n:]+structures[:n]
        
        return([o_new]+structures_new)
    
    def local_structure_to_ref(self, local_structure):
        
        o, s1, s2, s3, s4, s5, s6 = local_structure
        n = (o-1 -1)%6+1
        return(self.rotate_local_structure(local_structure, n))
        
            
        
                    
class ParticleRepresentation():
    """
    This class is used to plot particles in 2D or 3D knowing the lattice they belong to
    You can have several choice for the representation, which is then an argument of the plotting function
    To plot a particle, you need to know its position, orientation, and type of particle
    """
    
    
    
    
    def __init__(self, lattice ):
        self.lattice = lattice
        self.n_faces = self.lattice.n_neighbors
        self.orientations = self.lattice.orientations
        
        self.contact_centers_polar = [2*np.pi*(k/self.n_faces) for k in range(self.n_faces)]    
        self.angular_increment =  self.contact_centers_polar[1]
        # the theta coordinate for the centers of the face of the polyhedra
        self.vertices_polar = np.array(self.contact_centers_polar) +self.angular_increment/2
        # the theta coordinate for the vertices of the face of the polyhedra
        
        self.cos_contacts = np.cos(self.contact_centers_polar)       
        self.sin_contacts = np.sin(self.contact_centers_polar)      
        self.cos_vertices = np.cos(self.vertices_polar)
        self.sin_vertices = np.sin(self.vertices_polar)
        
        
        
        
        
        
        """Cube vertices"""
        #front face 
        # 3   4
        # 1   2
        # back face :
        # 7   8
        # 5   6
        v1 = np.array([0,0,0])
        v2 = np.array([1,0,0])
        v3 = np.array([0,0,1])
        v4 = np.array([1,0,1])
        v5 = np.array([0,1,0])
        v6 = np.array([1,1,0])
        v7 = np.array([0,1,1])
        v8 = np.array([1,1,1])
        self.cubeVerticies = [v1, v2, v3, v4, v5, v6, v7, v8]
        
        
        # colormaps = [cm.viridis, cm.seismic]
        colormaps = [cm.viridis, cm.spring, cm.summer, cm.winter, cm.autumn,cm.cool]
        imax = [0.9, 0.8, 0.8, 0.8, 0.5, 0.8]
        self.myColors1 = [colormaps[i](np.linspace(0,imax[i], self.n_faces)) for i in range(len(colormaps))]
        
        self.myColors2 = ['gray', 'silver', 'royalblue', 'mediumorchid']
        self.myColors2 = ['silver', "goldenrod", 'mediumaquamarine', 'royalblue', 'mediumorchid']
        
        """FFL, LLK"""
        self.myColors3 = [[ 'dodgerblue','dodgerblue','gold','gold','gold','gold'], 
                          ['darkblue', 'darkblue','dodgerblue', 'dodgerblue','dodgerblue', 'dodgerblue']]
        self.myColors3 = [['darkblue','darkblue','dodgerblue','dodgerblue','gold','gold'], 
                          ['darkblue', 'darkblue','dodgerblue', 'dodgerblue','dodgerblue', 'dodgerblue']]
        
        
        
        self.myColorPatchy = ['royalblue', 'gold', 'salmon', 'forestgreen', 'rebeccapurple', 'darkorange']
        
        self.add_angle = self.n_faces==4
        
    def plot_2Dtriangular_patch(self,x0, y0, size,edge_id, color_id, axes):
        
        ratio_size = 0.5/self.cos_vertices[0]
        
        phi = self.contact_centers_polar[edge_id] -self.angular_increment/2
        
        path = [[x0,y0]]
        i1, i2 = edge_id, (edge_id+1)%self.n_faces
        path += [[x0+size*ratio_size*self.cos_vertices[i1], y0+size*ratio_size*self.sin_vertices[i1]]]
        path += [[x0+size*ratio_size*self.cos_vertices[i2], y0+size*ratio_size*self.sin_vertices[i2]]]
        path=np.array(path)
        
        h = mpatches.Polygon(path, color=self.myColors1[color_id[0]][color_id[1]], linewidth=0)#,ec='black')
        axes.add_artist(h)
        
    def plot_2Dtriangular_patch_triangular(self, x0, y0, size, edge_id, color_id, axes ):
        
        ratio_size = 0.5/self.cos_vertices[0] *1.5
        
        phi = self.contact_centers_polar[edge_id] -self.angular_increment/2
        
        path = [[x0,y0]]
        i1, i2 = (edge_id)%self.n_faces, (edge_id+2)%self.n_faces
        path += [[x0+size*ratio_size*self.cos_contacts[i1], y0+size*ratio_size*self.sin_contacts[i1]]]
        path += [[x0+size*ratio_size*self.cos_contacts[i2], y0+size*ratio_size*self.sin_contacts[i2]]]
        path=np.array(path)
        
        h = mpatches.Polygon(path, color=self.myColors3[color_id[0]][color_id[1]], linewidth=0)#,ec='black')
        axes.add_artist(h)
        
        
        
        
    def plot_2Dtriangle(self, x0, y0, size, particle, orientation, axes, order=1):
        if orientation > 0 :
            if orientation%2 ==0:
                i_edges = [0,2,4]
                
            else :
                i_edges = [1,3,5]
                
            for i_edge in i_edges :
                color_face = (self.orientations[1-orientation][i_edge]* order)%self.n_faces
                color_id = (particle, color_face)
                self.plot_2Dtriangular_patch_triangular(x0, y0, size,i_edge, color_id, axes)
        
        
        
        
        
    def plot_2Dcolored_particle(self, x0, y0, size, particle, orientation, axes, order=1):
        if orientation > 0 :
            for i_edge in range (self.n_faces):
                color_face = (self.orientations[1-orientation][i_edge]* order)%self.n_faces
                color_id = (particle, color_face)
                self.plot_2Dtriangular_patch(x0, y0, size,i_edge, color_id, axes)
                
                
    def plot_2Dpatchy_particle(self, x0, y0, size, particle, orientation, axes, patchID=[0,1,2,3,4,5], dotsize=100):
        
        
        if orientation >0 :
            
            radius = 0.5/self.cos_vertices[0] * size
            ratio_size = 3#0.5/self.cos_vertices[0] *1.5
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, orientation=phi, facecolor="silver", edgecolor='black')
            axes.add_artist(h)
            
            for i_edge in range (self.n_faces):
                color_id = patchID[self.orientations[1-orientation][i_edge]]
                x = x0+radius*self.cos_vertices[i_edge]
                y= y0+radius*self.sin_vertices[i_edge]
                axes.scatter(x,y,s=dotsize, color=self.myColorPatchy[color_id], zorder=1e4)
            
           
        
    
    def plot_2Darrowed_particle(self, x0,y0,size, particle, orientation, axes, edgecolor='none'):
        if orientation > 0 :
            radius = 0.5/self.cos_vertices[0] * size
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, 
                                        orientation=phi, facecolor=self.myColors2[particle],
                                        edgecolor=edgecolor)
            
            arrow_radius = 1*radius/2
            
            
            face_id = self.orientations[orientation-1][0] 
            xh,yh = self.cos_contacts[face_id],self.sin_contacts[face_id]
            xt, yt = -xh, -yh
            xa0, ya0 = x0+xt*arrow_radius , y0+yt*arrow_radius
            dxa, dya = (xh-xt)*arrow_radius, (yh-yt)*arrow_radius
            
            head_length = arrow_radius/2
            head_width= arrow_radius/3
            width=size/40
            
            arrow = mpatches.FancyArrow(xa0, ya0, dxa, dya, width=width, color='black', length_includes_head=True, head_width=head_width,head_length=head_length)
                
            axes.add_artist(h)
            axes.add_artist(arrow)
        
        if orientation==0:
            radius = 0.5/self.cos_vertices[0] * size
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, 
                                        orientation=phi, facecolor='white',
                                        edgecolor=edgecolor)
            axes.add_artist(h)
            
            
    def plot_2Dparticle_orientation(self, x0, y0, size, particle, orientation, axes, order):
        if orientation>0:
            radius = 0.5/self.cos_vertices[0] * size
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, 
                                        orientation=phi, facecolor=self.myColors1[particle][orientation-1])
                                        # ,edgecolor=edgecolor)
            axes.add_artist(h)
            
        
    def plot_2Dparticle(self, x0, y0, size, particle, orientation, axes, representation, order=1):
        if representation=="colored":
            self.plot_2Dcolored_particle(x0, y0, size, particle, orientation, axes, order=order)
        if representation=="arrowed":
            self.plot_2Darrowed_particle( x0,y0,size, particle, orientation, axes)
        if representation=="triangle":
            self.plot_2Dtriangle( x0,y0,size, particle, orientation, axes)
        if representation=="orientation":
            self.plot_2Dparticle_orientation(x0, y0, size, particle, orientation, axes, order=order)
            


    def plot_3Darrow(self, x0, y0, z0, size, particle, orientation, axes):
        if orientation >0 :  
            
            face1, face2, face3 = 1,2,3
            edges = []
            for face in [face1, face2, face3]:
                edges.append(self.orientations[orientation-1][face-1] + 1)
            edge1, edge2, edge3 = edges
            
            color = self.myColors2[particle]
            
            lw_0 = 1
            lw_ratio = 1
            for edge in edges[:2] :
                dx, dy, dz = self.lattice.neighbors_cartesian[edge-1]
                
                arrow_radius = size/lw_ratio/2
                x1,x2 = x0-dx*arrow_radius, x0+dx*arrow_radius
                y1,y2 = y0-dy*arrow_radius, y0+dy*arrow_radius
                z1,z2 = z0-dz*arrow_radius, z0+dz*arrow_radius
            
                lw = lw_0/lw_ratio
                lw_ratio+=1
                mutation_scale = 3*size
                # head_width = arrow_radius/3
                # head_length = arrow_radius/2
                arrow = myArrow3D([x1, x2], [y1,y2], [z1,z2], mutation_scale=mutation_scale, lw=lw,
                                  arrowstyle="simple", color=color)
                                  # , length_includes_head=True, 
                                  # head_width=head_width,head_length=head_length)
                axes.add_artist(arrow)


    def plot_Cube(self, x0, y0, z0, size, particle,  orientation , axes, alpha=0.95):
        """alpha is the transparency """
        if orientation > 0:
            
            """Get the colors of all the faces"""
            myEdges = self.lattice.orientations[orientation-1]
            faces_colors = []
            for edge in range(6):
                # print(myEdges, edge, myEdges-edge)
                face = np.where(np.array(myEdges)==edge)[0][0]
                faces_colors.append(self.myColors1[particle][face])
            
            """Get the coordinate of the vertices"""
            radius = size/2
            center= np.array([x0,y0,z0])
            vertices = self.cubeVerticies
            new_vertices = []
            for v in vertices :
                new_vertices.append( center-radius + 2*radius*v)
            v1, v2, v3, v4, v5, v6, v7, v8 = new_vertices
            
            face_1 = np.array([[v2,v6],[v4,v8]])
            face_2 = np.array([[v6,v5],[v8,v7]])
            face_3 = np.array([[v4,v8],[v3,v7]])
            
            face_4 = np.array([[v1,v5],[v3,v7]])
            face_5 = np.array([[v1,v2],[v3,v4]])
            face_6 = np.array([[v1,v2],[v5,v6]])
            
            faces = [face_1, face_2, face_3, face_4, face_5, face_6]
            
            """Plot each face"""
            for i in range(len(faces)):
                face = faces[i]
                X,Y,Z = face[:,:,0], face[:,:,1], face[:,:,2]
                color = faces_colors[i]
                color[3] = alpha
                axes.plot_surface(X, Y, Z ,color=color, edgecolor='none')
            
        
    def plot_3Dhexagon(self, x0, y0, z0, size, particle, orientation, axes, representation=None):
        """Plots a 3D hexagon in a 2d lattice, with deferent colors depending on if it's up or donw"""
        
        orientation_corresp = {}
        orientation_corresp[10]=1
        orientation_corresp[9]=2
        orientation_corresp[11]=3
        orientation_corresp[8]=4
        orientation_corresp[12]=5
        orientation_corresp[7]=6 
        Triangular = ReadLatticeParticle('Triangular', '/Users/lara/Documents/SelfAssembly2/Lattice/')
        newRepresentation = ParticleRepresentation(Triangular)
        if orientation in [1,2,3,4,5,6]:            
            newRepresentation.plot_2Dparticle(x0, y0, size, particle, orientation, axes, representation='arrowed')
        
        else :
            newRepresentation.plot_2Dparticle(x0, y0, size, particle+1, orientation_corresp[orientation], axes, representation='arrowed')
        
        
    

    def plot_3Dparticle(self, x0, y0, z0, size, particle, orientation, axes, representation):
        if representation=='arrowed':
            self.plot_3Darrow( x0, y0, z0, size, particle, orientation, axes)
        if representation=="colored":
             self.plot_Cube(x0, y0, z0, size, particle, orientation, axes)


    def plot_particle(self, x0, y0, z0, size, particle, orientation, dimension, axes, representation, order=1):
        if dimension==3 :
            self.plot_3Dparticle(x0, y0, z0, size, particle, orientation, axes, representation)
        if dimension==2:
            if self.lattice.name == 'Hexagonal' :
                """3D hexagons in 2 dimension"""
                self.plot_3Dhexagon(x0, y0, z0, size, particle, orientation, axes, representation=None)
            else :
                self.plot_2Dparticle(x0, y0, size, particle, orientation, axes, representation, order=order)
       

    def plot_one_contact( self, x0, y0, z0, i_edge, radius, structure, axes, LEL, lw=4, emax=None, color_map=None, zorder=None):
            """draw a line of a color corresponding to the energy between particle in x_real, y_real and it's neighbor following edge
            emax is given in absolute value, x0 and y0 are coordinates in real space (not lattice)
            Either LEL or color_map needs to be known to define the colors"""
            
            
            
            if color_map is None:
                e = LEL[structure]
                icol = (e-(-emax))/(2*emax) #number between 0 and 1 that gives the color
                color =  cm.seismic(icol)
            else:
                color = color_map[structure]
                
            
            cos1 = self.cos_vertices[i_edge-1]
            sin1 = self.sin_vertices[i_edge-1]
            cos2 = self.cos_vertices[i_edge]
            sin2 = self.sin_vertices[i_edge]
            
            x1, x2, y1, y2 = x0+cos1*radius, x0+cos2*radius, y0+sin1*radius, y0+sin2*radius
            axes.plot([x1,x2], [y1, y2], color=color, linewidth=lw,  solid_capstyle='round', zorder=zorder )
            # )




"""matrices used for reordering of the lel"""
permut = np.array([[0,1,0,0,0,0],
                   [0,0,1,0,0,0],
                   [0,0,0,1,0,0],
                   [0,0,0,0,1,0],
                   [0,0,0,0,0,1],
                   [1,0,0,0,0,0]])

permut_inv = np.linalg.inv(permut)

anti_diag = np.fliplr(np.eye(6))


class DesignLEL():
    """This class is used to design specific LELs
    The LEL array is built with the energy you choose 
    for a given relative orientation, or pair of particles"""
    
    
    def __init__(self, lattice, n_particles_type):
        
        self.lattice = lattice
        self.n_particles_type = n_particles_type
        self.n_structures = self.lattice.all_N_structures[n_particles_type]
        
        self.N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        self.compute_invariance() # see the function for explanation on invariances
        
        
    def LEL_sticky_isotropic(self,  stickiness):
        """Stickiness is a positive number"""
        
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            # try :
            o1, o2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            # except :
                # print(s)
            if (o1!=0 and o2!=0 and edge!=0):
                LEL[s] = -stickiness
        return(LEL)
    
    
    def LEL_predetermined(self, LELarg):
        """Takes a LEL that has the dimension of the 2 particles structures, 
        and returns the full LEL to write to C++"""
        LEL = np.zeros(self.n_structures)
        # N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        LEL[self.N1:] = LELarg
        return(LEL)
    
    def cut_LEL(self, LELarg):
        """Takes a LEL that has the dimension of the 1 and 2 particles structures
        that is used with C++ 
        and returns the LEL with only two particles structures"""
        # N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        return(LELarg[self.N1:])
    
    
    def LEL_two_levels(self, list_of_favored_structures):
        LEL = np.zeros(self.n_structures)
        LEL[list_of_favored_structures+self.N1] = 1
        
        return(LEL)
        
    
    def LEL_fixed_point(self, name, stickiness, repulsion):
        """liquid, gas,crystal,sponge,crystallite,double_fiber,micelle, cycle"""
        
        if name=='liquid':
            return(self.LEL_sticky_isotropic(stickiness))
        
        if name=='gas':
            return(np.zeros(self.n_structures))
        
        
        o_pairs = {}
        """Those lists are rundundant (1,1) and (4,4) are the same for example"""
        o_pairs['crystal'] = [(1,1),(2,2),(3,3),(4,4),(5,5),(6,6)] 
        o_pairs['sponge'] = [(2,5),(3,1),(4,6),(1,2),(5,4)]
        o_pairs['sponge_sparse'] = [(1,6),(5,2)]
        o_pairs['crystallite'] = [(1,5),(2,4),(4,2),(5,1)]
        o_pairs['double_fiber'] = [(2,5),(3,6),(1,1),(4,4)]
        o_pairs['fiber'] = [(1,1),(4,4)]
        o_pairs['fiber_dimer'] = [(5,2),(2,5)]
        o_pairs['micelle'] = [(2,2),(3,3),(5,5),(6,6),
                              # (5,1),(4,2),
                              # (4,3),(6,1),
                              # (5,4),(1,2),
                              # (3,1),(4,6)
                              (1,5),
                              (1,6),
                              (2,1),
                              (2,3)
                              ]
        o_pairs['cycle'] = [(1,2)]
        
        
        
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            if (state1!=0 and state2!=0 and edge!=0):
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                if (o1, o2) in o_pairs[name]:
                    LEL[s] = -stickiness
                else:
                    LEL[s]=repulsion
        return(LEL)
    
    
    def LEL_canonical(self, contactTopology, reference_face=1):
        """contactTopology is dimer, trimer, cycle, or fiber
        works only for hexagons so far"""
        f_pairs={}
        if reference_face == 1:
            f_pairs['dimer']=(1,1)
            f_pairs['trimer']=(1,2)
            f_pairs['cycle']=(1,3)
            f_pairs['fiber']=(1,4)
        else :
            f = reference_face
            f_pairs['dimer']=(f,f)
            f_pairs['trimer']=(f,f%6 + 1)
            f_pairs['cycle']=(f,(f+1)%6+1)
            f_pairs['fiber']=(f,(f+2)%6+1)
        
        
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            p1,f1,p2,f2 = self.lattice.get_faces_in_contact(s, self.n_particles_type)
            if p1==p2:
                # print((f1,f2), f_pairs[contactTopology])
                if (f1,f2)==f_pairs[contactTopology] or(f2,f1)==f_pairs[contactTopology]   :
                    
                    LEL[s]=1
        return(LEL)
        
    
    
    def LEL_homomeric_particles(self, attractions, repulsion, particlesEquivalent=True):
        """Similar particles are stacked whatever the orientation but repell other sorts of particle
        The attraction can be not the same depending the sort of particle : 
            - then attractions is a list of the list of length the number of particle type
            - otherwise attractions is a positive number
        Repulsion is a positive number"""
        
        if particlesEquivalent :
            attractions = np.ones(self.n_particles_type)*attractions
            
        
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            if (state1!=0 and state2!=0 and edge!=0):
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                if p1==p2:
                    LEL[s]=-attractions[p1]
                else:
                    LEL[s]=repulsion
        return(LEL)
    
   
    def areEquivalent(self, LEL1, LEL2):
        for n_perm in range(self.lattice.n_orientations):
            for d in range(2):
                LELcompare, n, d_chosen = self.reorder_LEL(LEL2, n_permutation=n_perm, direction=d)
                if np.all(LEL1==LELcompare):
                    return(True)
        return(False)
    
    def translate_LEL(self, LEL):
        """This functino only deals with full LEL (one and two particles structures)
        A two particle LEL only should first be processed with LEL_predetermined function"""
        
        newLEL = np.zeros(self.n_structures)
        # N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        n_total_faces = self.lattice.n_neighbors * self.n_particles_type
        J_faces = np.zeros(n_total_faces+1)
        for structure in range (self.N1, self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][structure]
            if state1 == 0 or state2 == 0:
                """this is the structure that we shift towards zero energ
                we also save the corresponding coupling to substract it"""
                p1, face1, p2, face2 = self.lattice.get_faces_in_contact(structure,self.n_particles_type)
                if face1!=0:
                    J_faces[p1*self.lattice.n_neighbors+ face1] = LEL[structure]
                elif face2!=0:
                    J_faces[p2*self.lattice.n_neighbor+face2] = LEL[structure]
                else :
                    J_faces[0] = LEL[structure]
                
                newLEL[structure] = 0
       
        for structure in range (self.N1, self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][structure]
            if state1 != 0 and state2 != 0:
                """this is the full-full structures that we keep non zero"""
                p1, face1, p2, face2 = self.lattice.get_faces_in_contact(structure,self.n_particles_type)
                
                J10 = J_faces[p1*self.lattice.n_neighbors+ face1] 
                J20 = J_faces[p2*self.lattice.n_neighbors+ face2] 
                J00 = J_faces[0]
                newLEL[structure] = LEL[structure] + J00 - J10 - J20
        
        return(newLEL)
    
    def get_energy_shift_for_translation(self, LELuntranslated, Nbonds, Nparticles):
        """This functino only deals with full LEL (one and two particles structures)
        A two particle LEL only should first be processed with LEL_predetermined function"""
        
        # N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        sum_J_full_empty = 0
        J00 = 0
        n_total_faces = self.lattice.n_neighbors * self.n_particles_type
        for structure in range (self.N1, self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][structure]
            if state1 == 0 or state2 == 0:
                """this is the structure that we shift towards zero energ
                we also save the corresponding coupling to substract it"""
                p1, face1, p2, face2 = self.lattice.get_faces_in_contact(structure,self.n_particles_type)
                
                if face1!=0:
                    sum_J_full_empty += LELuntranslated[structure] 
                elif face2!=0:
                    sum_J_full_empty += LELuntranslated[structure]
                else :
                    J00 = LELuntranslated[structure]
        energy_shift = J00 * (-Nbonds +n_total_faces*Nparticles) - sum_J_full_empty*Nparticles
        energy_shift = energy_shift / Nbonds
        return(energy_shift)
        
        
    
    def LEL_crystal_homomeric(self, attractions, repulsion, particlesEquivalent=True):
        """Similar particles form a crystal but repell other sorts of particle
        The attraction can be not the same depending the sort of particle : 
            - then attractions is a list of the list of length the number of particle type
            - otherwise attractions is a positive number
        Repulsion is a positive number"""
        if particlesEquivalent :
            attractions = np.ones(self.n_particles_type)*attractions
            
        
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            if (state1!=0 and state2!=0 and edge!=0):
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                if p1==p2 and o1==o2:
                    LEL[s]=-attractions[p1]
                else:
                    LEL[s]=repulsion
        return(LEL)
    
    def LEL_random(self, affinity, anisotropy):
        """all full-full structures are random number of mean=affinity, std=anisotropy"""
        n1, n2, n3 = self.lattice.detail_all_N_structures[self.n_particles_type]
        lel_fullfull = np.random.normal(affinity, anisotropy, n3)
        LEL = np.zeros(self.n_structures)
        LEL[n1+n2:]=lel_fullfull
        return(LEL)
    
    
    
    def LEL_triangle_one_particle(self, AA, AB, AC, BC, BB, CC, einf=100):
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            if (state1!=0 and state2!=0 and edge!=0):
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if (o1,o2)==(1,4) :
                    LEL[s] = AA
                elif (o1,o2)==(3,6) :
                    LEL[s] = BB
                elif (o1,o2)==(5,2) :
                    LEL[s] = CC
                elif (o1,o2)==(3,4) or (o1,o2)==(1,6) :
                    LEL[s] = AB
                elif (o1,o2)==(5,4) or (o1,o2)==(1,2) :
                    LEL[s] = AC
                elif (o1,o2)==(5,6) or (o1,o2)==(3,2) :
                    LEL[s] = BC
                else :
                    LEL[s] = einf
                    
        return(LEL)
    
    
    def LEL_two_levels_unique(self, n_favored, e_fav, e_unfav, index, local=False):
        """This function loos at pre computed non redundant two level LEL for hexagon"""
        if local :
             directory = "/Users/lara/Documents/LPTMS/"
        else :
             directory = "/home/lkoehler/"
        directory += 'SelfAssembly2/SavedFromSimu/230421EnumerationTwoLevels/'
        all_two_levels = SF.fileToPickle(directory+'all_two_levels')
        
        LELshort = all_two_levels[n_favored][index]
        
        LEL_full = LELshort[7:]
        i0 = np.where(LEL_full==0)[0]
        i1 = np.where(LEL_full==1)[0]
        LEL_full[i0]=e_unfav
        LEL_full[i1]=e_fav
        
        LELshort[7:] = LEL_full
        return(self.LEL_predetermined(LELshort))
            
        
    
    
    def LEL_triangle_two_particle(self, A1A1, A1B1, A1C1, B1C1, B1B1, C1C1, 
                                  A1A2, A1B2, A1C2, B1C2, B1B2, C1C2,
                                  A2C1, A2B1, B2C1,
                                  A2A2, A2B2, A2C2, B2C2, B2B2, C2C2,
                                  einf=100):
        LEL = np.zeros(self.n_structures)
        
        AAs = [A1A1,A1A2,A1A2,A2A2]
        BBs = [B1B1,B1B2,B1B2,B2B2]
        CCs = [C1C1,C1C2,C1C2,C2C2]
        ABs = [A1B1,A1B2,A2B1,A2B2]
        ACs = [A1C1,A1C2,A2C1,A2C2]
        BCs = [B1C1,B1C2,B2C1,B2C2]
        
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            if (state1!=0 and state2!=0 and edge!=0):
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if (o1,o2)==(1,4) :
                    if (p1,p2)==(0,0):
                        LEL[s] = AAs[0]
                    elif (p1,p2)==(0,1) :
                        LEL[s] = AAs[1]
                    elif (p1,p2)==(1,0) :
                        LEL[s] = AAs[2]
                    elif (p1,p2)==(1,1) :
                        LEL[s] = AAs[3]
                        
                elif (o1,o2)==(3,6) :
                    if (p1,p2)==(0,0):
                        LEL[s] = BBs[0]
                    elif (p1,p2)==(0,1) :
                        LEL[s] = BBs[1]
                    elif (p1,p2)==(1,0) :
                        LEL[s] = BBs[2]
                    elif  (p1,p2)==(1,1) :
                        LEL[s] = BBs[3]
                        
                elif (o1,o2)==(5,2) :
                    if (p1,p2)==(0,0):
                        LEL[s] = CCs[0]
                    elif (p1,p2)==(0,1) :
                        LEL[s] = CCs[1]
                    elif   (p1,p2)==(1,0) :
                        LEL[s] = CCs[2]
                    elif  (p1,p2)==(1,1) :
                        LEL[s] = CCs[3]
                        
                elif (o1,o2)==(3,4) or (o1,o2)==(1,6) :
                    if (p1,p2)==(0,0):
                        LEL[s] = ABs[0]
                    elif (p1,p2)==(0,1) :
                        LEL[s] = ABs[1]
                    elif (p1,p2)==(1,0) :
                        LEL[s] = ABs[2]
                    elif (p1,p2)==(1,1) :
                        LEL[s] = ABs[3]
                elif (o1,o2)==(5,4) or (o1,o2)==(1,2) :
                    if (p1,p2)==(0,0):
                        LEL[s] = ACs[0]
                    elif (p1,p2)==(0,1) :
                         LEL[s] = ACs[1]
                    elif (p1,p2)==(1,0) :
                        LEL[s] = ACs[2]
                    elif  (p1,p2)==(1,1) :
                        LEL[s] = ACs[3]
                elif (o1,o2)==(5,6) or (o1,o2)==(3,2) :
                    if (p1,p2)==(0,0):
                        LEL[s] = BCs[0]
                    elif (p1,p2)==(0,1) :
                        LEL[s] = BCs[1]
                    elif (p1,p2)==(1,0) :
                        LEL[s] = BCs[2]
                    elif (p1,p2)==(1,1) :
                        LEL[s] = BCs[3]
                else :
                    LEL[s] = einf
                    
        return(LEL)
    
                    
        
        
    def LEL_camembert2D(self, ecrystal, eline, sigma, sigma_unfav, einf):
        LEL = np.zeros(self.n_structures)
        
        
         
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
          
                      
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if o1==o2 :
                    """Crystal"""
                    LEL[s] = ecrystal
                
                elif (o1==1 and o2==6) or (o1==2 and o2==3):
                    """Defect line"""
                    LEL[s] = eline
                else :
                    """Forbidden contact"""
                    LEL[s] = einf
            
            elif state2==0 and state1!=0 :
                """Empty-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                if o1==4 or o1==5 :
                    LEL[s]= sigma
                else :
                    LEL[s] = sigma_unfav
                    
            elif state1==0 and state2!= 0 :
                 if o1==1 or o1==2 :
                    LEL[s]= sigma
                 else :
                    LEL[s] = sigma_unfav
        return(LEL)
    
    
    
    def LEL_defect_line(self, n_defect_line, ecrystal, eline, sigma, einf):
        
        if n_defect_line==6:
            o1_o2s = [(1,6), (2,3)]
        if n_defect_line==3:
            o1_o2s = [(2,6), (3,1)]
        if n_defect_line==1:
            o1_o2s = [(3,6), (2,5)]
            
        
        LEL = np.zeros(self.n_structures)
        #print(o1_o2s[0][0] ,o1_o2s[0][1], o1_o2s[1][0],o1_o2s[1][1])
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
          
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if o1==o2 :
                    """Crystal"""
                    LEL[s] = ecrystal
                
                elif (o1==o1_o2s[0][0] and o2==o1_o2s[0][1]) or (o1==o1_o2s[1][0] and o2==o1_o2s[1][1]):
                    """Defect line"""
                    LEL[s] = eline
                else :
                    """Forbidden contact"""
                    LEL[s] = einf
            
            elif state2==0 and state1!=0 :
                """Empty-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                
                LEL[s]= sigma
               
            elif state1==0 and state2!= 0 :
                LEL[s]= sigma
                 
        return(LEL)
    
    
    
    
    def LEL_square_camembert(self, ecrystal, eline, sigma, einf):
        
      
        o1_o2s = [(1,2), (1,4), (2,3), (4,3)]
        
        LEL = np.zeros(self.n_structures)
        #print(o1_o2s[0][0] ,o1_o2s[0][1], o1_o2s[1][0],o1_o2s[1][1])
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
          
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if o1==o2 :
                    """Crystal"""
                    LEL[s] = ecrystal
                
                elif (o1==o1_o2s[0][0] and o2==o1_o2s[0][1]) or (o1==o1_o2s[1][0] and o2==o1_o2s[1][1]):
                    """Defect line"""
                    LEL[s] = eline
                else :
                    """Forbidden contact"""
                    LEL[s] = einf
            
            elif state2==0 and state1!=0 :
                """Empty-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                
                LEL[s]= sigma
               
            elif state1==0 and state2!= 0 :
                LEL[s]= sigma
                 
        return(LEL)
    
        
    
    
    
    def get_x_y_camembert(self, LEL,n_line=6,fiber=False):
        
        if self.n_particles_type==1:
            if not fiber:
                LELex = self.LEL_defect_line(n_line, -1, -2, 0, 10)[self.N1:]
               
            else :
                LELex = self.LEL_fiber_zigzag_2D(-1, -2, 0, 0, 10)[self.N1:]
            s_line = np.where(LELex==-2)[0]
            s_crystal = np.where(LELex==-1)[0]
            # s_line = np.array( [12,15], dtype=int)
            # s_crystal =  np.array([7,14,20], dtype=int)
            
        if self.n_particles_type==2:
            s_line = np.array( [18,21], dtype=int)
            s_crystal = np.array( [13,20,26], dtype=int)
        s_sigma =  np.array(range(1,7), dtype=int)
        
        if  (not np.all(LEL[s_line]==LEL[s_line[0]]) or 
                not np.all(LEL[s_crystal]==LEL[s_crystal[0]]) or
                not np.all(LEL[s_sigma]==LEL[s_sigma[0]])
                
                ):
            print('this LEL does not correspond to camembert')
            return(False)
    
        
        eline = LEL[s_line[0]]
        ecrystal = LEL[s_crystal[0]]
        
        sigma = LEL[s_sigma[0]]
        x = eline-2*sigma
        y = ecrystal - 2*sigma
        return(x,y, 6*sigma)
    
    
    def get_colormap_camembert(self, color_crystal='mediumaquamarine',
                               color_line = 'mediumblue', 
                               color_surf = 'tomato', 
                               color_zero = 'white',
                               color_inf = 'darkred', n_line=6, fiber=False, square=False):
        
        if fiber:
            LEL = self.LEL_fiber_zigzag_2D(-1, -2, 0, 0, 10)[self.N1:]
            s_line = np.where(LEL==-2)[0]
            s_crystal = np.where(LEL==-1)[0]
           
            
        elif square:
            LEL = self.LEL_square_camembert(-1, -2, 0, 10)[self.N1:]
            s_line = np.where(LEL==-2)[0]
            s_crystal = np.where(LEL==-1)[0]
            
        
        elif self.n_particles_type==1:
            # if n_line==6:
            #     s_line = np.array( [12,15], dtype=int)
            #     s_crystal =  np.array([7,14,20], dtype=int)
         
            LEL = self.LEL_defect_line(n_line, -1, -2, 0, 10)[self.N1:]
            s_line = np.where(LEL==-2)[0]
            s_crystal = np.where(LEL==-1)[0]
            
        elif self.n_particles_type==2:
            s_line = np.array( [18,21], dtype=int)
            s_crystal = np.array( [13,20,26], dtype=int)
            
        s_sigma =  np.array(range(1,7), dtype=int)
        
        color_map = {}
        color_map[0] = color_zero
        for s in range(1,self.n_structures-self.N1):
            if s in s_line:
                color_map[s] = color_line
            elif s in s_crystal :
                color_map[s] = color_crystal
            elif s in s_sigma:
                color_map[s] = color_surf
            else :
                color_map[s] = color_inf
                
        return(color_map)        

        
        
    
    def LEL_camembert2D_and_passive(self, ecrystal, eline, sigma, sigma_unfav, einf, eneutral=100):
        """particles of type 0 are camembert, particles of type 2 interact with no one"""
        LEL = np.zeros(self.n_structures)
        
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
          
                      
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if p1==0 and p2 ==0:
                    
                    if o1==o2 :
                        """Crystal"""
                        LEL[s] = ecrystal
                    
                    elif (o1==1 and o2==6) or (o1==2 and o2==3):
                        """Defect line"""
                        LEL[s] = eline
                    else :
                        """Forbidden contact"""
                        LEL[s] = einf
                else :
                    LEL[s]=eneutral
            
            elif (state2==0 and state1!=0) :#or (state1==0 and state2!= 0) :
                """Empty-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                if p1==0 :
                    if o1==4 or o1==5 :
                        LEL[s]= sigma
                    else :
                        LEL[s] = sigma_unfav
                else :
                    LEL[s]=eneutral//2
                        
            elif state1==0 and state2!= 0 :
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                if p2==0:
                    if o2==1 or o2==2 :
                        LEL[s]= sigma
                    else :
                        LEL[s] = sigma_unfav
                else :
                    LEL[s]=eneutral
                   
                   
        return(LEL)
    
        
    
    
    
    
    
        
    def LEL_camembert3D(self, ecrystal, eline, sigma, sigma_unfav, einf, ez, chiral=True):
        LEL = np.zeros(self.n_structures)
        
        
         
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
          
                      
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if edge==1:
                    if (o1==o2) :
                        """Crystal"""
                        LEL[s] = ecrystal
                    
                    elif ((o1,o2)==(1,6)) or ((o1,o2)==(2,3)) :
                        """Defect line"""
                        LEL[s] = eline
                        
                    elif not chiral and ((o1, o2) in [(10,6),(1,9), (10,9),(7,3),(2,12),(7,12)]) :
                        """Defect line"""
                        LEL[s] = eline    
                        
                    elif not chiral and ((o1, o2) in [(1,10),(10,1),
                                                      (2,7),(7,2),
                                                      (3,12),(12,3)]) :
                        """Defect line"""
                        LEL[s] = ecrystal    
                        
                    else :
                        """Forbidden contact"""
                        LEL[s] = einf
                    
                else :
                    """Forbidden contact"""
                    LEL[s] = ez
            
                    
            
            elif state2==0 and state1!=0 :
                """Empty-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                if o1==4 or o1==5 :
                    LEL[s]= sigma
                else :
                    LEL[s] = sigma_unfav
                    
            elif state1==0 and state2!= 0 :
                 if o1==1 or o1==2 :
                    LEL[s]= sigma
                 else :
                    LEL[s] = sigma_unfav
        return(LEL)
    
    
    def LEL_change_population_two_particles(self,e1, e2 ):
        
        energies_particle = [e1,e2]
        
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            p1, o1 =  self.lattice.state_to_particle_orientation[state1]
            p2, o2 =  self.lattice.state_to_particle_orientation[state2]
            
            if (state1==0 and state2!=0) :
                LEL[s] = energies_particle[p2]
            elif (state2==0 and state1!=0) :
                LEL[s] = energies_particle[p1]
        return(LEL)
                
               
        
        
    
    def LEL_double_fiber(self, eline, edimer):
        
        o_pairs = {}
        o_pairs['line'] = [(1,1),(4,4)]
        o_pairs['dimer'] = [(2,5),(3,6)]
        
           
        LEL = np.zeros(self.n_structures)
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
            if (state1!=0 and state2!=0 and edge!=0):
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                if (o1, o2) in o_pairs['line']:
                    LEL[s] = eline
                elif (o1, o2) in o_pairs['dimer']:
                    LEL[s]=edimer
        return(LEL)
    
    
    
    
    def LEL_fiber_zigzag_2D(self, ecrystal, eline, sigma, sigma_unfav, einf):
        LEL = np.zeros(self.n_structures)
        
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
                      
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if o1==o2 :
                    """Crystal"""
                    LEL[s] = ecrystal
                
                elif ((o1,o2)==(1,4)) or ((o1,o2)==(2,5)) or ((o1,o2)==(4,1)) or ((o1,o2)==(3,6)):
                    """Defect line"""
                    LEL[s] = eline
                else :
                    """Forbidden contact"""
                    LEL[s] = einf
            
                    
            
            elif state2==0 and state1!=0 :
                """Empty-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                if o1==6 or o1==5 :
                    LEL[s]= sigma
                else :
                    LEL[s] = sigma_unfav
                    
            elif state1==0 and state2!= 0 :
                 if o1==3 or o1==2 :
                    LEL[s]= sigma
                 else :
                    LEL[s] = sigma_unfav
        return(LEL)
    
    def compute_invariance(self):
        """this function creates the invariances table
        invariances[i] is a vector such that the scalar product c.invariances[i] is always a constant
        c.invariances[0] should be one (conservation of the bonds)
        c.inverainces[i!=0] shoud be density/3/nparticles_type
           this is the conservation of the number of faces"""        
        
        
        invariances = np.zeros((self.lattice.n_neighbors*self.n_particles_type+1, self.n_structures-self.N1), dtype=int)
        
        for s in range(self.N1, self.n_structures):
             p1, f1, p2, f2 = self.lattice.get_faces_in_contact(s, self.n_particles_type)             
             """the faces in contact :  face 3 of particle 2 is labeled 2*n_faces + 3"""
             f1_tot = p1*self.lattice.n_neighbors + f1
             f2_tot = p2*self.lattice.n_neighbors + f2
             if not (f1_tot==0 and f1_tot==0):
                     invariances[f1_tot,s-self.N1]+=1
                     invariances[f2_tot,s-self.N1]+=1
        invariances[0,:]=np.ones(self.n_structures-self.N1)
        self.invariances = invariances
    
    
    
    
    def transform_in_matrix(self, LELshort, p1=0, p2=0, symmetric=True, emptyFullAlso=False):
        """Takes the short verision of LEL, without N1P structures"""
        LEL = self.LEL_predetermined(LELshort)
        structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, 
                                                               symmetric=symmetric, 
                                                               n_type=self.n_particles_type,
                                                               emptyFullAlso=emptyFullAlso)
        n, n = structures_matrix.shape
        matrix = np.zeros((n,n))
        for i in range(n):
            for j in range (n):
                matrix[i,j] = LEL[structures_matrix[i,j]]
                
        return(matrix)
    
    
    def transform_from_matrix(self, smallMatrix, p1=0, p2=0, symmetric=True):
        """Takes the nface*nface matrix and returns LEL"""
        
        LEL = np.zeros(self.n_structures)
        structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, symmetric=symmetric,n_type=self.n_particles_type)
        n, n = structures_matrix.shape
        for i in range(n):
            for j in range (n):
                LEL[structures_matrix[i,j]] = smallMatrix[i,j] 
        return(self.cut_LEL(LEL))
        
    
    
    
    def reorder_LEL(self, LELshort, n_permutation=None, direction=None, symmetric=True):
        """takes short LEL and return short LEL
        if n_permutation and direction are given, the reordering is predetermined"""
        matrix = self.transform_in_matrix(LELshort, symmetric=symmetric)
        empty_fulls = LELshort[1:self.lattice.n_neighbors+1]
        
        """Permutation of the matrix"""
        matrix = matrix.copy()
        empty_fulls.copy()
        
        if n_permutation is None :
            dim = matrix.shape[0]
            eig, eigv = np.linalg.eigh(matrix)
            eigv = eigv.T
            order = np.argsort(eig)
            # print('order',order)
            # print('eig', eig[order])
            eigv = eigv[order]
            eig = eig[order]
            
            target_best_i = 0
            
            """choosing the most negative lambda*u_i"""
            # mult = np.zeros((dim,dim))
            # for i in range(dim):
            #     mult[i] = eigv[i]*eig[i] #this is a vector
            # min_mult = np.argmin(mult)
            # absolute_max, best_i = np.unravel_index(min_mult, (dim,dim))
                
            # all_mult = np.unique(mult)
            # # print('min_mult and next', all_mult[0],all_mult[1])
            # abs_max_eigenvector = np.abs(eigv[absolute_max])
            # direction = np.argmin([abs_max_eigenvector[(best_i-1)%dim], abs_max_eigenvector[(best_i+1)%dim]]) #this is 0 or 1
           
            
            """choosing biggest projection on (1,0,0..) for the biggest eigenvector"""
            absolute_max = np.argmin(np.abs(eig))
            # eigenvector of the lowest eigenvalue :
            abs_max_eigenvector = np.abs(eigv[absolute_max])
            # print('eigv', eigv[absolute_max])
            # print('eigvAbs', abs_max_eigenvector)
            
            #We choose the permutation of the eigenvector 
            # with the biggest projection on [1,0,0,...] 
            best_i = np.argmin(abs_max_eigenvector) 
            # print('absolute_max', 'best_i',absolute_max,  best_i)
            
            # We choose the permutation of the eigenvector 
            # with the biggest projection on [0,1,0,...] 
            # direction is 0 or 1    
            direction = np.argmin([abs_max_eigenvector[(best_i-1)%dim], abs_max_eigenvector[(best_i+1)%dim]]) #this is 0 or 1
           
            
            """Choosing with the faces"""
            # face_stickiness = np.average(matrix, axis=1)
            # best_i = np.argmin(face_stickiness)
            
            # #direction is 0 or 1           
            # direction = np.argmin([face_stickiness[(best_i-1)%dim], face_stickiness[(best_i+1)%dim]]) 
            
            
            # print('best_i', best_i)
            # print('direction', direction)
            if direction==0 :
                n_permutation = np.abs(target_best_i - best_i)
            else :
                n_permutation = (dim-1) - np.abs(target_best_i - best_i)
                
    
        # we will do n_permutation cyclic permutation
        M = np.linalg.matrix_power(permut, n_permutation)
        Minv = np.linalg.inv(M)
        new_matrix = M.dot(matrix.dot(Minv))
        
        new_empty_fulls = M.dot(empty_fulls)
        #if needed, we also do the miror permutation
        if direction == 1 :
            new_matrix = anti_diag.dot(new_matrix.dot(anti_diag))
            new_empty_fulls = anti_diag.dot(empty_fulls)
        
        
        newLEL = np.zeros(len(LELshort))
        
        newLELfullfull =  self.transform_from_matrix(new_matrix, symmetric=symmetric)
        newLEL =newLELfullfull
        newLEL[0] = LELshort[0]
        newLEL[1:self.lattice.n_neighbors+1] = new_empty_fulls
        
        return(newLEL, n_permutation, direction)
    
    

    
    
    
    
    def remove_positive(self, short_reoredered_LEL, cVector):
        """No 1 particle structures, and LEL should already be translated with only full-full that are nonzero"""
        LEL = short_reoredered_LEL
        newLEL = LEL.copy()*0
        
        ind_plus = np.where(LEL>0)[0]
        ind_minus = np.where(LEL<0)[0]
        if len(ind_plus)>0:
            e_plus = np.sum(cVector[ind_plus]*LEL[ind_plus])
        else :
            e_plus = 0
        if len(ind_minus)>0:
            sum_c_minus = np.sum(cVector[ind_minus])
            newLEL[ind_minus] = LEL[ind_minus]+e_plus/sum_c_minus
        return(newLEL)
    
    
    def remove_positive2(self,short_reoredered_LEL, cVector):
        c = cVector.copy()
        LEL = short_reoredered_LEL.copy()
        order = np.argsort(LEL)
        s_max = order[-7]
        # s_max = np.argmax(LEL)
        c_sum = np.sum(c[7:])
        if LEL[s_max]>0:
            LEL[7:] = (LEL[7:]-LEL[s_max])/c_sum
        return(LEL)
            
    

    
    def convert_to_old_code(self, LELshort):
        m = self.transform_in_matrix(LELshort)
        m_indices_old= np.array([[18, 14, 10,  9, 12, 15],
           [14, 17, 13, 21,  8, 11],
           [10, 13, 16, 19, 20,  7],
           [ 9, 21, 19, 23, 24, 22],
           [12,  8, 20, 24, 26, 25],
           [15, 11,  7, 22, 25, 27]])
        
        
        LEL_old = np.zeros(28)
        LEL_old[0] = LELshort[0]
        
        for s in range (1,7):
            old_s = self.lattice.new_to_old[s]
            LEL_old[old_s] = LELshort[s]
        
        
        for i in range(6):
            for j in range(6):
                old_s = m_indices_old[i][j]
                coupling = m[i][j]
                LEL_old[old_s] = coupling
        return(LEL_old)
    
    def convert_from_old_code(self, LEL_old):
        m_indices_old= np.array([[18, 14, 10,  9, 12, 15],
           [14, 17, 13, 21,  8, 11],
           [10, 13, 16, 19, 20,  7],
           [ 9, 21, 19, 23, 24, 22],
           [12,  8, 20, 24, 26, 25],
           [15, 11,  7, 22, 25, 27]])
        m = np.zeros((6,6))
        
        LELshort = np.zeros(28)
        LELshort[0] = LEL_old[0]
        
        for old_s in range (1,7):
            s = self.lattice.old_to_new[old_s]
            LELshort[s] = LEL_old[old_s]
        
        
        for i in range(6):
            for j in range(6):
                m[i][j] = LEL_old[m_indices_old[i][j]]
        LELfullfull = self.transform_from_matrix(m)
        
        LELshort[7:] = LELfullfull[7:]
        return(LELshort)
        
        

class LELrepresentation():
    
    def __init__(self, lattice, n_particles_type):
        self.lattice = lattice
        self.n_particles_type = n_particles_type
        self.n_structures = self.lattice.all_N_structures[n_particles_type]
        
        self.designs = DesignLEL(self.lattice, self.n_particles_type)
        
        
    def get_colormatrixLEL(self, LEL, p1=0, p2=0,colormap=cm.seismic, 
                           colordict=None, maxi=None,
                           structures_in_grey=None,
                           grey_color='gray',
                           symmetric=True,
                           emptyFullAlso = False):
        """Takes the short verision of LEL, without N1P structures"""

        
        matrix = self.designs.transform_in_matrix(LEL, p1=p1, p2=p2, symmetric=symmetric, emptyFullAlso=emptyFullAlso)
        n = self.lattice.n_orientations + emptyFullAlso
        N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        
        if structures_in_grey is not None:
            if maxi is None and LEL is not None:
                maxi = max(-np.min(LEL), np.max(LEL))   
            structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, 
                                                                   symmetric=symmetric,
                                                                   n_type=self.n_particles_type,
                                                                   emptyFullAlso=emptyFullAlso)-N1
            colorMatrix = np.zeros((n,n,4))
            for i in range (n):
                for j in range (n):
                    s = structures_matrix[i,j]
                    if s in structures_in_grey:
                       colorMatrix[i,j]  =  matplotlib.colors.to_rgba(grey_color)
                    else :
                        if LEL is None :
                            colorMatrix[i,j]  =  matplotlib.colors.to_rgba('white')
                        else :
                            J = matrix[i][j]
                            level = 0.5 + J/maxi
                            colorMatrix[i,j] = colormap(level)
                            
        
        elif colordict is not None :
            structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, 
                                                                   symmetric=symmetric,
                                                                   n_type=self.n_particles_type,
                                                                   emptyFullAlso=emptyFullAlso)-N1
            colorMatrix = np.zeros((n,n,4))
            for i in range (n):
                for j in range (n):
                    colorMatrix[i,j] = matplotlib.colors.to_rgba(colordict[structures_matrix[i,j]])
        else :
            if maxi is None :
                maxi = max(-np.min(LEL), np.max(LEL))
            colorMatrix = np.zeros((n,n,4))
            for i in range (n):
                for j in range (n):
                    J = matrix[i][j]
                    level = 0.5 + J/maxi
                    colorMatrix[i,j] = colormap(level)
                
        return(colorMatrix)
    
    def get_colormatrixDensity(self, vector, p1=0, p2=0,
                               colormap=cm.seismic,maxi=None,
                               symmetric=True):
        """Takes the short verision of c, without N1P structures"""

        
        matrix = self.designs.transform_in_matrix(vector, p1=p1, p2=p2, symmetric=symmetric).copy()
        n = self.lattice.n_orientations
        N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        if maxi is None :
            maxi = np.max(vector[N1:])
        # print(maxi)
        colorMatrix = np.zeros((n,n,4))
        for i in range (n):
            for j in range (n):
                J = matrix[i][j]
                if i!=j:
                    J = J/2
                level = 0.5 - J/maxi
                colorMatrix[i,j] = colormap(level)

        return(colorMatrix)
    
    def plot_matrix(self, vector, vectorIsLEL=True, 
                    p1=0, p2=0,
                    axes=None, 
                    x0=0, y0=0, dx=1,dy=1,
                    fullMatrix = True, #if False, plot only the triangular inferior part
                    emptyFullAlso = False, #if true, plot also the structures empty-full
                    lw = 0.5,axes_off=True,
                    colormap=cm.seismic, 
                    colordict=None, #this is to colorcode the structures, not the energy level
                    structures_in_grey = None,
                    grey_color = 'gray',
                    arangeAxis = True,
                    maxi=None, 
                    margin=2, symmetric=True, 
                    plot_index_particle=False):
        """Takes the short verision of LEL, without N1P structures"""
        
        if vectorIsLEL :
            colorMatrix = self.get_colormatrixLEL(vector, p1=p1, p2=p2,
                                                  colormap=colormap, colordict=colordict, 
                                                  structures_in_grey=structures_in_grey,
                                                  grey_color=grey_color,
                                                  maxi=maxi, symmetric=symmetric, emptyFullAlso=emptyFullAlso)
        else  :
            colorMatrix = self.get_colormatrixDensity(vector, p1=p1, p2=p2,
                                                      colormap=colormap,maxi=maxi
                                                      , symmetric=symmetric)
                
        n,n, colordimension=colorMatrix.shape
    
        if axes is None :
            f,axes = plt.subplots()
    
        for i in range(n):
            if fullMatrix :
                jlist = range(0,n) 
            else :
                jlist = range(i,n)
            for j in jlist:
                color = colorMatrix[i,j]
                square = mpatches.RegularPolygon((x0+i*dx,y0-j*dy), 4, orientation=np.pi/4, 
                                                 radius=dx/np.sqrt(2),  color=color,ec='black', 
                                                 lw=lw)
                axes.add_artist(square)
    
        if arangeAxis :
            margin = margin
            axes.set_xlim(x0-margin, x0+(n-1)*dx+margin)
            axes.set_ylim(y0-(n-1)*dx - margin, y0+margin)
            axes.set_aspect('equal', 'box')  
            axes.set_xticks([])
            axes.set_yticks([])#,[])
            
        
        else:
            n1, n2 = -2, 7
            axes.set_xlim(x0+n1*dx, x0+n2*dx)
            axes.set_ylim(y0-n2*dx, y0-n1*dx)
            axes.set_aspect('equal', 'box')
        if axes_off :
            axes.axis('off')
            
        """plot matrix indices as particles"""
        if plot_index_particle:
            
            PR = ParticleRepresentation(self.lattice)
            
            dxP = dx
            rx = 0.7*dx
            ry=rx
            
            PR.plot_2Darrowed_particle(x0,y0+dxP, rx, 0, 0, axes, edgecolor='black')
            PR.plot_2Darrowed_particle(x0-dxP, y0, rx, 0, 0, axes, edgecolor='black')
            
            
            # plt.scatter(x0,y0+dxP, color='gray')
            if self.lattice.name=='Triangular':
                if emptyFullAlso:
                    o1s = [0,1,2,3,4,5,6]
                    o2s = [0,4,5,6,1,2,3]
                    my_range = range(1,7)
                else:
                    o1s = [1,2,3,4,5,6]
                    o2s = [4,5,6,1,2,3]
                    my_range = range(0,6)
            if self.lattice.name=='Square':
                if emptyFullAlso:
                    o1s = [0,1,2,3,4]
                    o2s = [0,3,4,1,2]
                    my_range = range(1,5)
                else:
                    o1s = [1,2,3,4]
                    o2s = [3,4,1,2]
                    my_range = range(0,4)
                
            for i in my_range:
                o = o2s[i]
                PR.plot_2Darrowed_particle(x0=x0+(i)*dx, y0=y0+dxP, size=rx, 
                                           orientation=o,particle=0, axes=axes, edgecolor='black')
            for j in my_range:
                o = o1s[j]
                PR.plot_2Darrowed_particle(x0=x0-dxP, y0=y0-(j)*dx, size=rx, 
                                           orientation=o, particle=0, axes=axes, edgecolor='black')
    
    
            
            #%%
class ClustersAndOrders():

    def __init__(self, name, index, directoryImage):
        """Takes a file data_image where position and orientations 
        of the particles were saved before
        Do analysis like measuring order parameters, clusters, etc"""
        
        self.name = name
        self.index = index
        self.directoryImage = directoryImage        
        data_image = SF.fileToPickle(directoryImage+name)[index]
        self.data_image = data_image
        self.NparticlesTot = data_image['NparticlesTot']
        self.particles_positions = data_image['particles_positions']
        self.particles_id = data_image['particles_id']
        self.particles_orientations = data_image['particles_orientations']
        self.dimension = data_image['dimension']
        
        self.limits = data_image['limits']
        
        self.lx = self.limits[0][1]-self.limits[0][0]
        self.ly = self.limits[1][1]-self.limits[1][0]
        self.lz = self.limits[2][1]-self.limits[2][0]
        
        #get lattice
        directoryLattice = '/Users/lara/Documents/SelfAssembly2/Lattice/'   
        self.lattice = ReadLatticeParticle(data_image['lattice'], directoryLattice)  
        
        
        self.distance_threshold = 1.1
        
        self.get_orientation_vectors()
        
        self.compute_distance_matrix()
        self.compute_adjacency_matrix()
        self.find_clusters()
        
        
    def get_orientation_vectors(self):
        generators = self.lattice.generators
        orientations = self.lattice.neighbors
        
        self.orientation_vectors = {}
        
        for o in range (len(orientations)):
            v = np.zeros(3)
            for d in range (len(generators)):
                v+= orientations[o][d]*generators[d]
            self.orientation_vectors[o+1]=v
                
        
        
    def measure_distance(self, i1, i2):
        """Measure distance between particles indexed i1 and i2"""
        x1,y1, z1 = self.particles_positions[0][i1],self.particles_positions[1][i1],self.particles_positions[2][i1]
        x2,y2, z2 = self.particles_positions[0][i2],self.particles_positions[1][i2],self.particles_positions[2][i2]
        
        dx = min( (x2-x1)%self.lx, (x1-x2)%self.lx )
        dy = min( (y2-y1)%self.ly, (y1-y2)%self.ly)
        if self.lz!=0:
            dz = min( (z2-z1)%self.lz, (z1-z2)%self.lz)
        else :
            dz=0
        d = np.sqrt(dx**2+ dy**2 + dz**2)
        
        return(d)
            
    def are_neighbors(self,i1, i2):
        """Test if particles indexed i1 and i2 are neighbors"""
        d = self.measure_distance(i1, i2)
        return(d<=self.distance_threshold)
    
    def compute_distance_matrix(self):
        self.distances  = np.zeros((self.NparticlesTot,self.NparticlesTot))
        for i in range (self.NparticlesTot):
            for j in range(self.NparticlesTot):
                self.distances [i,j] = self.measure_distance(i,j)
             
        
    def compute_adjacency_matrix(self):
        self.adjacency = np.zeros((self.NparticlesTot, self.NparticlesTot), dtype=bool)
        for i in range (self.NparticlesTot):
            for j in range(self.NparticlesTot):
                if self.distances[i,j]<= self.distance_threshold and self.distances[i,j]!=0:
                    self.adjacency[i,j] = True
       
    def get_neighbors(self, i):
        neighbors = np.where(self.adjacency[i])[0]
        return(neighbors)
       
    def find_clusters(self):
        
        never_considered = set()
        for i in range (self.NparticlesTot):
            never_considered.add(i)
        current_cluster = 0
        in_exploration = set()
        
        self.clusters = -np.ones(self.NparticlesTot, dtype=int)
        
        while not len(never_considered)==0:
            current_particle_ref = list(never_considered)[0]
            in_exploration.add(current_particle_ref)
            while not len(in_exploration)==0:
                current_particle = list(in_exploration)[0]
                
                self.clusters[current_particle] = current_cluster
                never_considered.remove(current_particle)
                in_exploration.remove(current_particle)
                
                for n in self.get_neighbors(current_particle):
                    if n in never_considered:
                        in_exploration.add(n)
            
            current_cluster+=1
        
        
        self.Nclusters = current_cluster+1
        self.Nparticles_per_cluster = np.zeros(self.Nclusters, dtype=int)
        for i in range (self.NparticlesTot):
            self.Nparticles_per_cluster[self.clusters[i]]+=1
        
        
        
        
        
        # # clusters = [[0]]
        # # to_explore = [queue.Queue()]
        
        # explored = np.zeros(self.NparticlesTot, dtype=bool)
        # reference = -np.ones(self.NparticlesTot, dtype=int)
        # reference[0] = 0
        
        # references_connexion = np.arange(self.NparticlesTot, dtype=int)
        # references_found = []
        # references_common = {}
        
        # self.clusters = -np.ones(self.NparticlesTot, dtype=int)
        # # clusters_found = 0
        
        # count = 0
        # while not np.all(explored) and count<self.NparticlesTot:
        #     count+=1
        #     next_to_explore = np.where(1-explored)[0][0] #first not explored index
        #     explored[next_to_explore] = True
                
        #     neighbors = self.get_neighbors(next_to_explore)
            
        #     if reference[next_to_explore] == -1:
        #         if np.all(reference[neighbors]==-1):
        #             # print('s1')
        #             # none of the neighbors were aldready assigned a reference
        #             # it means we found a new cluster
        #             this_ref = next_to_explore
                    
        #             reference[next_to_explore] = this_ref
        #             reference[neighbors] = this_ref
        #         else :
        #             # print('s2')
        #             identified_neighbor = neighbors[np.where(reference[neighbors]!=-1)[0][0]]
        #             # print('neighbors', next_to_explore, neighbors, 'identified',identified_neighbor)
        #             this_ref = reference[identified_neighbor]
                    
        #             reference[next_to_explore] = this_ref
        #             reference[neighbors] = this_ref
                    
        #     else :
        #         this_ref = reference[next_to_explore]
        #         # print('s3', next_to_explore,this_ref )
        #         # print('neighbors', next_to_explore, neighbors)
        #         for n in neighbors :
        #             other_ref = reference[n]
        #             other_refs = [other_ref]
                    
        #             if  other_ref!=-1 and other_ref !=this_ref:
        #                 deepest_ref = references_connexion[min(this_ref, other_ref)]
        #                 references_connexion[this_ref] = deepest_ref
        #                 references_connexion[other_ref] = deepest_ref
        #         reference[neighbors] = this_ref
        #     references_found.append(this_ref)
                    
        # print(references_connexion)
        # for ref in references_found:
        #     ref1 = ref
        #     ref2 = references_connexion[ref]
        #     identity = ref1==ref2
        #     while not identity:
                
        #         ref2 = ref1
        #         ref1 = references_connexion[ref1]
        #         identity = ref1==ref2
            
        #     references_common[ref]=ref1
                
                
                
        # print(references_common)
        # for r in references_common.keys():
        #     print(r, reference[references_common[r]])
            
            
        # for i in range (self.NparticlesTot):
        #     ref = reference[i]
        #     if ref in references_common.keys() :
        #         reference[i] = reference[references_common[ref]]
                
        # clusters_id = np.unique(reference)
     
                

        
        
        
                
        
    def plot(self, representation='colored', axes=None, dx=0, dy=0,
             add_cluster_colors=False):
        if not add_cluster_colors :
            plot_a_system(self.data_image,  representation=representation, axes=axes, dx=dx, dy=dy)
        else :
            colors_cluster = cm.hsv(np.linspace(0,1, self.Nclusters))
            f, axes = plt.subplots(1,2)
            plot_a_system(self.data_image,  representation=representation, axes=axes[0], dx=dx, dy=dy)
            for i in range(self.NparticlesTot):
                x,y = self.particles_positions[0][i], self.particles_positions[1][i]#, zs_plot[i] 
                axes[1].scatter(x,y, marker='o', color =colors_cluster[self.clusters[i]] )
                axes[1].text(x,y, str(i), fontsize=8)
           
            
           
    def measure_orientational(self):
        Ms = np.zeros((self.Nclusters,2))
        for i in range (self.NparticlesTot):
            c = self.clusters[i]
            o = self.particles_orientations[i]
            mx,my, mz = self.orientation_vectors[o]
            Ms[c][0]+=mx
            Ms[c][1]+=my
            
            
            
        return(Ms)
    
    
    
#%%

    
    
    
    
    
            
def plot_a_system(data_image,  representation='colored', axes=None, dx=0, dy=0, 
                  order=1, plot_contacts=False,color_map_contacts=None, LEL=None,
                  lw=2):
    
    
                
    NparticlesTot = data_image['NparticlesTot']
    xs_plot, ys_plot, zs_plot=  data_image['particles_positions'] 
    particle_ids = data_image['particles_id'] 
    orientations = data_image['particles_orientations'] 
            
    dimension = data_image['dimension']
    xlim, ylim, zlim =  data_image['limits'] 
    
    
    if dy!=0:
        ys_rel = ys_plot-ylim[0]
        new_ys = (ys_rel+dy)%(ylim[-1]-ylim[0])
        ys_plot = new_ys + ylim[0]
        
    if dx!=0:
        xs_rel = xs_plot-xlim[0]
        new_xs = (xs_rel+dx)%(xlim[-1]-xlim[0])
        xs_plot = new_xs + xlim[0]
       
    directoryLattice = '/Users/lkoehler/Documents/SelfAssembly2/Lattice/'    
    lattice = ReadLatticeParticle(data_image['lattice'] , directoryLattice)
    Particle = ParticleRepresentation(lattice)
    
    if axes is None :
        if dimension<3 :
            figure,axes = plt.subplots(figsize=(7,7))
        else :
            figure = plt.figure(figsize=(7,7))
            axes = plt.axes(projection='3d')
            
    if lattice.name=="Triangular":
        particle_size = 0.95
    elif lattice.name=="Square":
        particle_size=0.96
    for i in range(NparticlesTot):
        x,y,z = xs_plot[i], ys_plot[i], zs_plot[i] 
        
        Particle.plot_particle(x, y, z, particle_size, particle_ids[i],
                               orientations[i], dimension, axes, 
                               representation, order=order)
        
    if plot_contacts and 'contacts' in data_image.keys():
        contacts = data_image['contacts']
        for key in contacts.keys():
            x,y,z,i_edge = key
            structure = contacts[key]
            if structure <= lattice.n_orientations:
                zorder=NparticlesTot*3
            elif color_map_contacts[structure]=='mediumblue' or color_map_contacts[structure]=="mediumaquamarine":
                zorder = NparticlesTot
            else :
                zorder = NparticlesTot*2
            if lattice.name=="Triangular":
                radius = 0.5/np.cos(np.pi/6) #* particle_size
            elif lattice.name=="Square":
                radius = 0.5/np.cos(np.pi/4)
            Particle.plot_one_contact(x, y, z, i_edge, radius, structure, 
                                       axes, 
                                       LEL, lw=lw, emax=30, 
                                       color_map=color_map_contacts, 
                                       zorder=zorder)
        

        
           
            
    """Limits"""
    
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    
    if dimension>=3 :
        axes.set_zlim(zlim)
    else :
        axes.set_aspect('equal', 'box')  
    
    axes.tick_params(axis='both',  bottom=False, left=False, labelbottom=False, labelleft=False)
        

directoryLattice = "/Users/lkoehler/Documents/SelfAssembly2/Lattice/"
Triangular = ReadLatticeParticle('Triangular', directoryLattice)  
matrix= DesignLEL(Triangular,1).lattice.get_matrix_structures()-7

 
def c00(s,T, dense):
     if not dense:
         return(np.exp(s/T))
     else:
         return(0)
 
def c0a(s, a, T, dense):
     if not dense:
         return(np.exp((s+a)/T))
     else:
         return(0)
 
def cab(s,a,b, J, T):
     return( np.exp((s+a+b-J)/T ))
 
 
 
def J(i,j, LEL):
     return LEL[matrix[i,j]]
 
 
def c_sum(s,a_values, LEL,T, dense):
     n_tot=  c00(s,T, dense) 
     
     for i in range(6):
         n_tot+= c0a(s,a_values[i],T, dense)
         
     for i in range(6):
         for j in range(i,6):
             n_tot+= cab(s, a_values[i], a_values[j], J(i,j, LEL),T)
             
             
     return(n_tot)
     
     
def c_face(i, s, a_values, LEL, dense,T=1):
     
     na = c0a(s, a_values[i],T, dense) 
     
     for j in range(6):
         n = cab(s, a_values[i], a_values[j], J(i,j, LEL),T)
         if i!=j:
             na+=n
         else:
             na+=2*n
             
     return(na)
             
 
 
def func(sa_values, LEL, Nface, Nbonds,T, dense):
     
     res = np.zeros(7)
     s = sa_values[0]
     a_values = sa_values[1:7]
     
     res[0] = ( c_sum(s, a_values, LEL,T, dense) - 1 )
     for i in range(6):
         res[i+1] = ( c_face(i, s, a_values, LEL,T, dense) - Nface/Nbonds)
         
     return(res)
 
def sa_values_real(c, dense):
     if not dense:
         s = np.log(c[0] )
         a_values = np.log(c[1:7])-s 
         return np.concatenate((np.array([s]), a_values))
 
 
def c_b_func(sa_values, LEL, Nbonds,T, dense):
     
     c = np.zeros(28)
     
     s = sa_values[0]    
     a_values = sa_values[1:7]
     
     c[0] = c00(s,T, dense)
     for i in range(1,7):
         c[i] = c0a(s,a_values[i-1],T, dense)
         
     for i in range(6):
         for j in range(i,6):
             c[matrix[i,j]] = cab(s, a_values[i], a_values[j], J(i,j, LEL),T)
             
     return(c)  
 
    
def solve_c_boltzmann(designs, LEL, density, Nparticles, Nsites, sa_init,
                            dense = False, T=1, tol=1e-10):

   
    Nbonds = 3*Nsites
    Nface = Nparticles

    sa_values_found = scipy.optimize.fsolve(func, sa_init, 
                                                            args = (LEL, Nface, Nbonds, T, dense) ,
                                                            maxfev=2000,
                                                            full_output=False,
                                                            xtol=tol)
            
    sol_found = np.linalg.norm(func(sa_values_found, LEL, Nface, Nbonds, T, dense))<tol
    # 

    return(sa_values_found, sol_found)




def measure_c_boltzmann(designs, LEL, density, Nparticles, Nsites, c_init=None,
                            dense = False, T=1, tol=1e-10, N_steps=50, print_info=False):
    

    T = 1

    sa_init = sa_values_real(c_init, dense=False)
    sa_values_found, found = solve_c_boltzmann(designs, LEL, density, Nparticles, Nsites,sa_init=sa_init, tol=tol)
    n = 0
    dT = 0.5
    one_solution_found = found
    """if the solution was not found, we start at higher temperature and use the solution as 
    a strating point of the next optimization"""
    while not (found and T==1) and n < N_steps :
        
        if (not found):
            if one_solution_found:
                dT = dT/2
                T = T+dT
                
            else :
                T=2*T
        else:
            one_solution_found = found
            T = max(1, T-dT)
            sa_init = sa_values_found
        if print_info:
            print(n,found,T)
        sa_values_found, found = solve_c_boltzmann(designs, LEL, density, Nparticles, Nsites,sa_init=sa_init, T=T, dense=dense,tol=tol)
        n+=1
    c_b =   c_b_func(sa_values_found, LEL, 3*Nsites,T, dense=dense)
    # print("len",len(sa_values_found), len(c_b))
    return(c_b, (found and T==1))
    


    
    
    
def measure_c_min( designs, LEL, density, Nparticles, Nsites, 
                        c_init=None, dilute=True, 
                        maxiter=2000, disp=False, method="SLSQP",
                        add_surface_constraints = True
                        ):
        """Minimizes the enrgy under some constraints"""
        
        
        """dimension depends on if the system is dense or dilute"""
        if dilute or (not add_surface_constraints) :
            dim=28
        else :
            dim = 21
            LEL =  LEL[7:]
            if c_init is not  None :
                c_init = c_init[7:]
        
        
        
        """Initialisation"""
        if c_init is None :
            c_init = np.zeros(dim)
            c_init_is_zero = True
        else :
            c_init_is_zero = np.all(c_init==0)
        
        
        if dilute :
            """Compute the minium number of surface bonds of an aggregate 
            with nparticles in hexagonal lattice"""
            small_surfaces = np.array([0,6,10,12,14,16,18,18, 20, 22])        
            radius = lambda n: (-3 + np.sqrt(9+12*n)) / 6        
            if Nparticles>9:
                nsurf =  12*radius(Nparticles)+6
            else :
                nsurf =  small_surfaces[Nparticles]
                
            c_surf_min = nsurf/(3*Nsites)
            
            """c_surf_max"""
            c_surf_max = 2*density

            
            """sum(c_surf) is between c_surf_gas and c_surf_perfect_bulk"""      
            #     #the sum of the empty-full is constrained
            constraints_surface = np.concatenate((np.zeros(1), np.ones(6), np.zeros(21)))
            c_surf_perfect_bulk, c_surf_gas = c_surf_min,c_surf_max
        
        
        """Constraints"""
        allConstraints = []
        """sum(c) = 1 and sum(c_face_i)=density/3"""
        if dilute or (not add_surface_constraints):
            constraints_conservation = designs.invariances
        else :
            constraints_conservation = designs.invariances[:,7:]
        d=density/3
    
        if c_init_is_zero  :
            constraints_conservation_values_min = np.array([1,d,d,d,d,d,d])
        else :
            constraints_conservation_values_min = constraints_conservation.dot(c_init)
        constraints_conservation_values_max = constraints_conservation_values_min
       
        ConstraintsConservation = LinearConstraint(constraints_conservation, constraints_conservation_values_min, constraints_conservation_values_max)
        # ConstraintsLimits = LinearConstraint(constraints_limits,constraints_limits_values_min, constraints_limits_values_max )
    
        allConstraints.append(ConstraintsConservation)#,ConstraintsLimits]
        
        
        """all the cs are between 0 and 1"""
        # constraints_limits = np.eye(dim)
        # constraints_limits_values_min = np.zeros(dim)
        # constraints_limits_values_max = np.ones(dim)
        
        
        
        # """constraints on the empty full for dense system"""
        # if not dilute :
        #     constraints_limits_dense = np.zeros((7,28))
        #     for k in range (7):
        #         constraints_limits_dense[k,k] = 1
        #     constraints_limits_values_min_dense = np.zeros(7)#np.concatenate((np.zeros(7), np.zeros(21)))
        #     constraints_limits_values_max_dense = np.zeros(7)#np.concatenate((np.zeros(7), np.zeros(21))
           # if not dilute :
        #     ConstraintsLimitsDense = LinearConstraint(constraints_limits_dense,constraints_limits_values_min_dense, constraints_limits_values_max_dense )
        #     allConstraints.append(ConstraintsLimitsDense)
        
       
        
     
        if dilute and add_surface_constraints :
            ConstraintsLimitsSurf = LinearConstraint(constraints_surface,c_surf_perfect_bulk, c_surf_gas )
            allConstraints.append(ConstraintsLimitsSurf)
        
        
        """Objective function"""
        def func(c):
            return(np.sum(c*LEL))
        
        """Minimisation"""
        res = minimize(func, c_init,  method=method,
                       constraints=allConstraints,
                       options={'maxiter':maxiter, 'disp':disp}, bounds=[(0,1) for k in range(dim)]
                       )
        
        c_boltzmann = res['x']
        
        if not dilute and add_surface_constraints :
            c_boltzmann = np.concatenate((np.zeros(7),c_boltzmann))
        
        
        return(c_boltzmann)
    
    
def measure_frustration(c, LEL, designs, density, Nparticles, Nsites, alpha=0, dilute=True,
                        add_surface_constraints=True, T=0):
    """returns frustration (energy per link) and c_boltzmann)"""
    
    
    """Dwo initialisation from two different inital conditions and take the best"""
    if T==0:
        c_boltzmann1 = measure_c_min(designs, LEL, density, Nparticles, Nsites, c_init=c, dilute=dilute, add_surface_constraints=add_surface_constraints)
        c_boltzmann2 = measure_c_min(designs, LEL, density, Nparticles, Nsites, c_init=c*0, dilute=dilute,add_surface_constraints=add_surface_constraints)
        
        if c_boltzmann1.dot(LEL)<c_boltzmann2.dot(LEL):
            c_boltzmann=c_boltzmann1
        else: 
            c_boltzmann=c_boltzmann2
    
    else:
        c_boltzmann = measure_c_boltzmann(designs, LEL, density, Nparticles, Nsites, c_init=c,
                                          dense = not dilute)
 
        
    """Compute frustration"""
    Dc = c - c_boltzmann
    
    frustration = np.sum( Dc * LEL * c**alpha)
    # if frustration<0:
    #     c_boltzmann1 = measure_c_boltzmann(designs, LEL, density, Nparticles, Nsites, c_init=c, dilute=dilute, maxiter=10000)
    #     c_boltzmann2 = measure_c_boltzmann(designs, LEL, density, Nparticles, Nsites, c_init=None, dilute=dilute, maxiter=10000)
    
    
    return(frustration, c_boltzmann)
    
    
def measure_frustration_dense(c_dense, c_dilute, LEL, Nparticles, Nsites, Nsites_dense):
    """Returns e_dense-e, e_dense, and e (energy per particles) """
    e_dense = np.sum(c_dense*LEL) * 3 * Nsites_dense / Nparticles
    e = np.sum(c_dilute*LEL) * 3 * Nsites / Nparticles
    
   
    frustration = e_dense-e
    return(frustration, e_dense, e)

#%% Colormaps 

color_map_camembert = {1: 'tomato',
 2: 'tomato',
 3: 'tomato',
 4: 'tomato',
 5: 'tomato', 
 6: 'tomato',
 
 0: 'white',
 
 7: 'mediumaquamarine',
 14: 'mediumaquamarine',
 20: 'mediumaquamarine',
 
 12: 'mediumblue',
 15: 'mediumblue',

 8: 'darkred',
 9: 'darkred',
 10: 'darkred',
 11: 'darkred',
 13: 'darkred',
 16: 'darkred',
 17: 'darkred',
 18: 'darkred',
 19: 'darkred',
 21: 'darkred',
 22: 'darkred',
 23: 'darkred',
 24: 'darkred',
 25: 'darkred',
 26: 'darkred',
 27: 'darkred'}    

color_map_test= {}
for s in range (28):
    color_map_test[s] = cm.viridis(s/28)
    
    
#%% Load lattice files    
                
if __name__ == "__main__": 
    
    directory = '/Users/lkoehler/Documents/SelfAssembly2/Lattice/'   
    # Cubic = ReadLatticeParticle('Cubic', di∂rectory)  
    
    # FCC = ReadLatticeParticle('FaceCenteredCubic', directory)  
    # BCC = ReadLatticeParticle('BodyCentereCubic', directory)  
    Triangular = ReadLatticeParticle('Triangular', directory)  
    Square = ReadLatticeParticle('Square', directory)  
    # Hexagonal = ReadLatticeParticle('Hexagonal', directory)  


#%% Test particle representation
    

    """Test Square"""
    # myParticle = ParticleRepresentation(Square)
    # f, axes = plt.subplots()
    # # Psquare.plot_triangular_patch(0, 0, 1, 0, 0, axes1)
   
    # size = 0.9
    # for y in range(1):
    #     for x in range(1,myParticle.n_faces+1):
    #         plt.text(x, y-1, str(x) )
    #         myParticle.plot_2Dcolored_particle( x,y,size, y, x, axes)
    # # y=2
    # # for x in range(myParticle.n_faces):
    # #         myParticle.plot_arrowed_particle( 2*x,2*y,1, 0, x, axes)
    # y=3
    # for x in range(1,myParticle.n_faces+1):
    #         myParticle.plot_2Darrowed_particle( x,y,size, 1, x, axes)
            
    # axes.set_xlim(-2,myParticle.n_faces+2) 
    # axes.set_ylim(-2,6) 
    # axes.set_aspect('equal', 'box')  
    
    """Tests Triangular"""
    # myParticle = ParticleRepresentation(Triangular)
    # f, axes2 = plt.subplots()
    # # Psquare.plot_triangular_patch(0, 0, 1, 0, 0, axes1)
    # for y in range(2):
    #     for x in range(myParticle.n_faces):
    #         # plt.text(2*x, 2*y-0.5, str(x) )
    #         # myParticle.plot_2Dcolored_particle( x,y,1, y, x, axes2)
    #         myParticle.plot_particle(x, y, 1, 0.8, y, x, 2, axes2, representation='colored', order=1)
            
    # y=2
    # for x in range(myParticle.n_faces):
    #         myParticle.plot_2Darrowed_particle( x,y,0.8, 0, x, axes2)
    # axes2.set_xlim(-2,myParticle.n_faces+2) 
    # axes2.set_ylim(-2,4) 
    # axes2.set_aspect('equal', 'box')    
    
    
    """Test cube"""
    # myParticle = ParticleRepresentation(Cubic)
    # f = plt.figure()
    # axes = plt.axes(projection='3d')
    # x0,y0,z0 = 0,0,0
    # size=1
    
    # for x0 in range (1, 4,2):
    #     orientation  =x0//2
    #     myParticle.plot_3Dparticle(x0, y0, z0, size, 0, orientation, axes, 'arrowed')
    # axes.set_xlim([-1, x0+1])
    # # axes.set_ylim([-1, x0+1])
    # # axes.set_zlim([-1, x0+1])
    # # axes.set_aspect('equal', 'box')  
    
    # f = plt.figure()
    # axes = plt.axes(projection='3d')
    # for x0 in range (2, 25,2):
    #     orientation  = x0//2
    #     myParticle.plot_Cube(x0, y0, z0, size, 0, orientation, axes)
    # axes.set_xlim([-1, x0+1])
    # axes.set_ylim([-1, x0+1])
    # axes.set_zlim([-1, x0+1])
    # axes.set_xlabel('x')
    # axes.set_ylabel('y')
    # axes.set_zlabel('z')
    
    
    """Test triangle"""
    # myParticle = ParticleRepresentation(Triangular)
    
    # f,axes =  plt.subplots()
    # x0,y0,z0 = 0,0,0
    # size=1
    
    # for x0 in range (1, 7,1):
    #     orientation  =x0
    #     myParticle.plot_2Dtriangle(x0, y0, size, 0, orientation, axes)
    # axes.set_xlim([-1, x0+1])
    # axes.set_ylim([-1, x0+1])
    # axes.set_aspect('equal', 'box')  
    
    
    
    
    
    
    
#%% Test LEL design 
    # l = DesignLEL(Triangular, 2)
    # lel1 = l.LEL_sticky_isotropic(4)
    # lel2 = l.LEL_homomeric_particles(2, 3)
    
    
    # designs = DesignLEL(Triangular, 1)
    # lel2 = designs.LEL_fixed_point('sponge', 10, 2)
                
    # lel2 = designs.LEL_predetermined(np.ones(28))
    
    # l= DesignLEL(Hexagonal, 1)
    # ecrystal, eline, sigma, sigma_unfav, einf, ez =  -8, 0.5, 6, 6, 15, 2.1
    # lel1= l.LEL_camembert3D(ecrystal, eline, sigma, sigma_unfav, einf, ez)
    # lel2= l.LEL_camembert3D(ecrystal, eline, sigma, sigma_unfav, einf,ez, chiral=False)
                
#%% Test cube 
                
    
     
    # from mpl_toolkits import mplot3d

    # x = np.arange(0,2,1)
    # y = np.arange(0,2,1)
    # X,Y = np.meshgrid(x,y)
    # Z = 0*X
    
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    
    # ax.plot_surface(X, Y, Z ,color='blue', edgecolor='none')
    
    # plt.show()
        
 #%% Test ClustersAndOrders
 
# directoryImage = "/Users/lara/Documents/LPTMS/SelfAssembly2/SavedFromSimu/230613Renorm/ImagesData/"
# name = '230613Renorm_Af_-40_Ani_090_n_000'
# index = 0
# CO = ClustersAndOrders(name, index,directoryImage )
# CO.plot(add_cluster_colors=True)
# # print(CO.measure_orientational())

 
 