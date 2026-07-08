#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 17:15:01 2023

@author: lara
"""

import json
import numpy as np
import time
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize

from LatticeTools import ReadLatticeParticle
from LatticeTools import DesignLEL
from LatticeTools import ParticleRepresentation

from SideClasses import SideFunctions as SF
import subprocess


from matplotlib import cm
# import imageio

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 18})
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 

class ReadSystem():
    
    def __init__(self, name, directoryResults, processLocal=True):
        
        
        self.name = name
        
        
        file = open(directoryResults+name+'input_parameters.json', 'r')
        self.parameters =  json.load(file)
        if processLocal:
            self.directoryResults = self.parameters["directoryResultsRead"]
        else :
            self.directoryResults = self.parameters["directoryResults"]
    
        self.Nparticles = self.parameters['Nparticles']
        self.NparticleTypes = self.parameters['NparticleTypes']
        self.NparticlesTot = np.sum(self.Nparticles)
        self.Lx, self.Ly, self.Lz = self.parameters['Lx'],self.parameters['Ly'],self.parameters['Lz']
        
        if self.parameters['saveSnapshot']:
            self.sites_positions = SF.makeDictionaryPosition(self.directoryResults+"_sites_positions.txt")
        
        if processLocal:
            self.lattice = ReadLatticeParticle(self.parameters["lattice"], self.parameters["directoryInputsRead"])
        else :
            self.lattice = ReadLatticeParticle(self.parameters["lattice"], self.parameters["directoryInputs"])
        
        self.Particle = ParticleRepresentation(self.lattice)
        self.get_dimension()
        
        self.LEL = SF.fileToArray(self.directoryResults+'_LEL.txt')
        
        
    def get_dimension(self):
        sizes = [self.Lx, self.Ly, self.Lz]
        dimension = len(np.where(np.array(sizes)>1)[0])
        self.dimension = dimension
        
        
    def get_temperature_ramp(self, linear=True):
        
        Temperature_start = self.parameters['Temperature_start']
        Temperature_end = self.parameters['Temperature_end']
        NTemperatures =  self.parameters['Ntemperatures']
        
        if NTemperatures==1:
            return([Temperature_end])
        
        if linear :
            a = (Temperature_start-Temperature_end) / (NTemperatures-1) ;
            b = Temperature_start 
            Temperatures = -a*np.arange(0,NTemperatures)+ b 
        else :
            a = - (Temperature_start-Temperature_end) / np.log(NTemperatures) ;
            b = Temperature_start ;
            Temperatures = a*np.log((0,NTemperatures)+1)+ b 
        return(Temperatures)
        
        
        
        
    
    def plot_temperature_evolution(self, axes=None, return_array=False):
        if self.parameters["saveTemperatures"]:
            all_temperatures = np.zeros(0)
            
            """Annealing temperatures"""
            Temperatures = SF.fileToArray(self.directoryResults+'_Temperatures.txt')
            for t in Temperatures :
                all_temperatures = np.concatenate((all_temperatures, t*np.ones(self.parameters['Nsteps_per_T'])))
    
            """Quenching temperatures"""
            t0 = self.parameters["Temperature_quenching"]
            all_temperatures = np.concatenate((all_temperatures, t*np.ones(int(self.parameters['Nsteps_quenching']))))
        
            if axes is None :
                f, axes = plt.subplots()  
            else :
                f=None
            axes.plot(all_temperatures)
            axes.set_xlabel('Monte-Carlo steps')
            axes.set_ylabel('Temperature')
            if f is not None :
                f.tight_layout()
                
            if return_array :
                return(all_temperatures)
        else :
            print('Temperatures were not saved for this simulation')
    
    def plot_energy_evolution(self, axes=None, return_array=False):
        if self.parameters["saveEnergies"]:
            energies = SF.fileToArray(self.directoryResults+'_Energies.txt')/self.NparticlesTot
            
            if axes is None :
                f, axes = plt.subplots()
            else :
                f=None
            axes.plot(energies)
           
            
            axes.set_ylabel('Energy per particle')
            if f is not None :
                axes.set_xlabel('Monte-Carlo steps')
                f.tight_layout()
            
            if return_array:
                return(energies)
         
        else :
            print('Energies were not saved for this simulation')
        
        
    def get_MC_success(self):
        success = SF.fileToArray(self.directoryResults+'_MCsuccess.txt')
        return(success)
    
    def plot_energy_and_temperature(self):
        f, axes = plt.subplots(2,1, figsize=(12,8), sharex=True)
        self.plot_energy_evolution(axes[0])
        axes[0].tick_params(axis='x',  bottom=False, left=True, labelbottom=False, labelleft=True)
        
        self.plot_temperature_evolution(axes[1])
        plt.subplots_adjust(hspace=0)
        f.tight_layout()
        

            
    def plot_system_state(self, representation="colored",  Temperature=None, index=None,
                          figure=None, axes=None, origin=None, 
                          return_data = False, do_the_plot = True, order=1, 
                          return_origin=False ,
                          plot_contact=False, return_contact=False,
                          color_map_contacts = None):
        
        """Load the data about particles positions and orientations"""
        if Temperature is None :
            directory1 = self.directoryResults+"adressFulls_groundState.txt"
            directory2 = self.directoryResults+"orientationFulls_groundState.txt"
            
        else :
            T = str(Temperature)
            i = str(index)
            directory1 = self.directoryResults+"adressFulls_T" + T +"_step"+i+".txt"
            directory2 = self.directoryResults+"orientationFulls_T" + T +"_step"+i+".txt"
    
        addressOfParticles = SF.fileToArray(directory1)
        orientationOfParticles = SF.fileToArray(directory2)

        """Get the lattice positions and orientatins"""
        # lattice coordinate and orientation
        xs_l, ys_l, zs_l, states = np.zeros(self.NparticlesTot), np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot, dtype=int)
        for i in range (self.NparticlesTot):
             current_orientation=orientationOfParticles[i]
             current_address=addressOfParticles[i]
             
             latticePos = self.sites_positions[current_address]
             xs_l[i] = latticePos[0]
             ys_l[i] = latticePos[1]
             zs_l[i] = latticePos[2]
             states[i] = current_orientation
             
             
             
        """Add contacts"""
        if plot_contact or return_contact:
            contacts_temp = {}
            def add_contact(x_init, y_init, z_init, orientation):
                for p in [0,1]:
                    """p=0, the particle is at the origin of the edge vector
                    p=1, the particle it at the end of the edge vector"""
                    
                    for i_edge in range(self.lattice.n_edges) :
                        if p==1:
                            
                            x, y, z = np.array([x_init, y_init, z_init], dtype=int) - self.lattice.edges[i_edge]
                            
                            x,y,z = x%self.Lx, y%self.Ly, z%self.Lz
                        if p==0:
                            x,y,z  = x_init, y_init, z_init
                           
                        if (x, y,z, i_edge) not in contacts_temp :
                            contacts_temp[(x,y,z, i_edge)] = [0, 0]
                        
                        contacts_temp[(x,y,z,i_edge)][p] = orientation
                       
            
            for i in range(self.NparticlesTot):
                x,y,z = xs_l[i], ys_l[i], zs_l[i]
                orientation = states[i]
                add_contact(x,y,z, orientation)
                    
            contacts = {}
            for key in contacts_temp.keys():
                x,y,z, i_edge = key
                o1, o2 = contacts_temp[key]
                structure = self.lattice.get_structure(self.NparticleTypes, o1, o2, i_edge+1, twoParticlesStructures=True)
                
                contacts[x,y,z, i_edge]= structure
           
        """Transform the state id into particle and orientation"""
        particle_ids, orientations = np.zeros(self.NparticlesTot, dtype=int), np.zeros(self.NparticlesTot, dtype=int)
        for i in range(self.NparticlesTot):
            particle_ids[i], orientations[i] = self.lattice.get_particle_id_orientation(states[i])
     
        
        
             
        """Find the coordinates origin such that cluster of particles are in the middle"""
        if origin is None :
            origin = self.find_best_lattice_origin( xs_l.copy(), ys_l.copy(), zs_l.copy())
            
        x0, y0, z0 = origin
        ys_l = (ys_l-y0)%self.Ly
       
        if return_origin :
            return(origin)
        
             
        
        
        """Create the figure"""
        
        if do_the_plot and axes is None :
            if self.dimension<3 :
                figure,axes = plt.subplots(figsize=(7,7))
            else :
                figure = plt.figure(figsize=(7,7))
                axes = plt.axes(projection='3d')
       
        """Save the position in cartesian coordinates"""
        xs_plot, ys_plot, zs_plot = np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot)
        
             
                
             
        """Transform positions from lattice coordinate to cartesian coordinate"""
        # cartesian coordinate and orientation
        # xs_c, ys_c, zs_c = np.zeros(self.NparticlesTot), np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot)
        particle_size = 0.95
        for i in range(self.NparticlesTot):
            x,y,z = xs_l[i], ys_l[i], zs_l[i]
            
            if (self.lattice.name=='Triangular' or self.lattice.name=='Hexagonal'):
               x,y,z  = self.lattice.lattice_resquare([x,y,z],self.Lx, self.Ly )
        
            x,y,z = self.lattice.lattice_to_cartesian([x,y,z])
            x = (x-x0)%self.Lx
            
            xs_plot[i], ys_plot[i], zs_plot[i] = x,y,z
            
            if do_the_plot :
                self.Particle.plot_particle(x, y, z, particle_size, 
                                            particle_ids[i], orientations[i], 
                                            self.dimension, axes, 
                                            representation, order=order)
                
        """Also transform the contacts in cartesian coordinate"""
        
        if plot_contact or return_contact:
            contacts_cartesian = {}
            for key in contacts.keys():
                x,y,z,i_edge = key
                y =(y-y0)%self.Ly
                structure = contacts[key]
                if (self.lattice.name=='Triangular' or self.lattice.name=='Hexagonal'):
                   x,y,z  = self.lattice.lattice_resquare([x,y,z],self.Lx, self.Ly )
                x,y,z = self.lattice.lattice_to_cartesian([x,y,z])
                x = (x-x0)%self.Lx
                
                contacts_cartesian[x,y,z,i_edge]=structure
                
                if plot_contact:
                    # print(x,y,i_edge, structure)
                    self.Particle.plot_one_contact(x, y, z, i_edge, 0.5, structure, 
                                                   axes, 
                                                   self.LEL, lw=2, emax=30, 
                                                   color_map=color_map_contacts, 
                                                   zorder=None)
                    
            
            
           
            
        """Compute the limits"""        
        xLmin, xLmax = 0-particle_size/2, self.Lx-particle_size/2 
        yLmin, yLmax = 0-particle_size/2, self.Ly-particle_size/2        
        zLmin, zLmax = 0-particle_size/2, self.Ly-particle_size/2 
        
        xmin,y,z = self.lattice.lattice_to_cartesian((xLmin, 0, 0))
        xmax,y,z  = self.lattice.lattice_to_cartesian((xLmax, 0, 0))
        x,ymin,z = self.lattice.lattice_to_cartesian((0, yLmin, 0))
        x,ymax,z  = self.lattice.lattice_to_cartesian((0, yLmax, 0))
        x,y,zmin = self.lattice.lattice_to_cartesian((0, 0, zLmin))
        x,y,zmax  = self.lattice.lattice_to_cartesian((0, 0, zLmax))
        
        
        """Make the limits"""
        if do_the_plot :
            axes.set_xlim(xmin, xmax)
            axes.set_ylim(ymin, ymax )
            
            if self.dimension>=3 :
                axes.set_zlim(zmin, zmax)
            else :
                axes.set_aspect('equal', 'box')  
            
            axes.tick_params(axis='both',  bottom=False, left=False, labelbottom=False, labelleft=False)
        
        
        """Save to replot later without having to load everything"""
        if return_data :
            data_image = {}
            data_image['NparticlesTot'] = self.NparticlesTot
            data_image['particles_positions'] = [xs_plot, ys_plot, zs_plot]
            data_image['particles_id'] = particle_ids
            data_image['particles_orientations'] = orientations
            
            data_image['lattice'] = self.lattice.name
            data_image['dimension'] = self.dimension
            data_image['limits'] = [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
            
            if return_contact :
                data_image['contacts'] = contacts_cartesian
            
            return(data_image)
        
        
        
    def get_local_structures(self):
        
        """Load the data about particles positions and orientations"""
            
        directory1 = self.directoryResults+"adressFulls_groundState.txt"
        directory2 = self.directoryResults+"orientationFulls_groundState.txt"
            
       
        addressOfParticles = SF.fileToArray(directory1)
        orientationOfParticles = SF.fileToArray(directory2)

        """Get the lattice positions and orientatins"""
        # lattice coordinate and orientation
        xs_l, ys_l, zs_l, states = np.zeros(self.NparticlesTot), np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot),np.zeros(self.NparticlesTot, dtype=int)
        for i in range (self.NparticlesTot):
             current_orientation=orientationOfParticles[i]
             current_address=addressOfParticles[i]
             
             latticePos = self.sites_positions[current_address]
             xs_l[i] = latticePos[0]
             ys_l[i] = latticePos[1]
             zs_l[i] = latticePos[2]
             states[i] = current_orientation
             
             
             
        """Add contacts"""
        contacts_temp = {}
        def add_contact(x_init, y_init, z_init, orientation):
            for p in [0,1]:
                """p=0, the particle is at the origin of the edge vector
                p=1, the particle it at the end of the edge vector"""
                
                for i_edge in range(self.lattice.n_edges) :
                    if p==1:
                        
                        x, y, z = np.array([x_init, y_init, z_init], dtype=int) - self.lattice.edges[i_edge]
                        
                        x,y,z = x%self.Lx, y%self.Ly, z%self.Lz
                    if p==0:
                        x,y,z  = x_init, y_init, z_init
                       
                    if (x, y,z, i_edge) not in contacts_temp :
                        contacts_temp[(x,y,z, i_edge)] = [0, 0]
                    
                    contacts_temp[(x,y,z,i_edge)][p] = orientation
                    #contacts_temp[(x_init,y_init,z_init,i_edge+3)][int(not p)] = orientation
                   
        
        for i in range(self.NparticlesTot):
            x,y,z = xs_l[i], ys_l[i], zs_l[i]
            orientation = states[i]
            add_contact(x,y,z, orientation)
                
        contacts = {}
        
        for i in range(self.NparticlesTot):
            x,y,z = xs_l[i], ys_l[i], zs_l[i]
            orientation = states[i]
            contacts[x,y,-1] = orientation
            
        for key in contacts_temp.keys():
            x,y,z, i_edge = key
            o1, o2 = contacts_temp[key]
            structure = self.lattice.get_structure(self.NparticleTypes, o1, o2, i_edge+1, twoParticlesStructures=True)
            
            #contacts[x,y,z, i_edge]= structure
            contacts[x,y, i_edge]= structure
            
            if o2!=0:
                xprime, yprime, zprime = np.array([x, y, z], dtype=int) + self.lattice.edges[i_edge]
                xprime,yprime,zprime = xprime%self.Lx, yprime%self.Ly, zprime%self.Lz
                #contacts[xprime, yprime, zprime, i_edge+3] = structure         
                contacts[xprime, yprime, i_edge+3] = structure       
                
                
        
        
        local_structures = {}
        for i in range(self.NparticlesTot): 
            x,y,z = xs_l[i], ys_l[i], zs_l[i]
            ls = [contacts[x,y,-1]]+[contacts[x,y,i_edge] for i_edge in range(6)]
            local_structures[x,y] =  self.lattice.local_structure_to_ref(ls)
        
       
        
        return(local_structures)
    
    
    
          
       
        
       
        
        
    def find_best_lattice_origin(self, xs, ys, zs):
        """From a list of poistions of particles (x, y), gets the translation x0 and y0 
        that needs to be done so that the cluster will be represented in the middle"""
        
        
        """Look for y with no particle to put the zero there"""
        yempty = list(np.arange(0, self.Ly))
        yoccupied = np.unique(ys)
        for myY in yoccupied :
           try: yempty.remove(myY)
           except: donothing=0
           
        if len(yempty)==0:
            """We will look for a line where there is no cluster going up"""
            notFound = True
            i=0
            while notFound and i<self.Ly:
                myY = ys[i]
                
                inds = np.where(ys==myY)[0]
                indsUP = np.where(ys==(myY+1)%self.Ly)[0]
                
                xoccupied = xs[inds]
                xups = xs[indsUP] #x coordinate of particles on the line above
              
                forbidUP = np.unique(np.concatenate((xoccupied, (xoccupied-1)%self.Lx))) 
                #coordinates that should not be occupied on the line above 
                
                intUP = list(set(list(forbidUP)) & set(list(xups))) 
                notFound = ( len(intUP) != 0 )
                if notFound : i+=1
            
            y0 = ys[i]
            
            
        else:
            """We will look for the y in the middle of an empty zone to be the new zero"""
           
            yempty = np.sort(yempty)
            y0 = self.find_middle_bigger_zone(yempty, self.Ly)
       
              
        ys = (ys-y0)%self.Ly
        
        
        """The y has been dealt with, let's look for the x if necessary"""
        """Let's make a square"""
        
        
        for i in range (self.NparticlesTot):
            xs[i] , ys[i], zs[i] = self.lattice.lattice_resquare([xs[i],ys[i],zs[i]],self.Lx, self.Ly )
        

        """We modify also the xs"""
        emptyxs = list(np.arange(self.Lx+1) -1 - self.Ly//4)
        fullxs = np.unique(np.concatenate((xs[np.where(ys==self.Ly//2)[0]] , xs[np.where(ys==self.Ly//2+1)[0]])))
        # print('full x ', fullxs)
        for xoccupied in fullxs :
            try: emptyxs.remove(xoccupied)
            except: donothing=0
   
        
        emptyxs = np.array(emptyxs)+self.Ly//4+1
        x0 = self.find_middle_bigger_zone(emptyxs, self.Lx+1) 

        return(x0,y0, 0)
        
    def find_middle_bigger_zone(self, emptylist, L):
        """This is a private function used to find a good middle of the coordinates to represent the system nicely"""
        if len(emptylist)==0:
            return(0)
        elif len (emptylist)==1:
            return(emptylist[0])
        else:
            zonestart, zonewidth = [emptylist[0]],[]
            zonesize=0
            if (0 in emptylist) and (L-1 in emptylist) :
                lastIsSplitted = True
            else :
                lastIsSplitted = False
            for i in range (1,len(emptylist)):
                # print(i, emptylist[i], zonestart, zonewidth,zonesize)
                if emptylist[i%len(emptylist)]!=(emptylist[(i-1)%len(emptylist)]+1)%L:
                    zonestart.append(emptylist[i%len(emptylist)])
                    zonewidth.append(zonesize)
                    zonesize=1
                else :
                    zonesize+=1
            zonewidth.append(zonesize)
            if lastIsSplitted :
                zonestart = np.array(zonestart)[1:]
                zonewidth[-1]+=zonewidth[0]
                zonewidth = np.array(zonewidth)[1:]
            
            if len(zonestart)!=0:
                y0 = (zonestart[np.argmax(zonewidth)] + np.max(zonewidth)//2 )%L
            else :
                y0=0
            return(y0)
              
    
    def plot_one_structure_s(self, s=None, 
                       state1=None, state2=None, edge=None,
                       p1=None, p2=None,o1 = None, o2=None, 
                       representation='arrowed' ):
        
        """Structure is either defined by s, by (state1, state2, edge), or by (p1, o1, p2, o2, edge)"""
        
        
        axes, s = self.lattice.show_structure(s=s, 
                       state1=state1, state2=state2, edge=edge,
                       p1=p1, p2=p2,o1 = o1, o2=o2, 
                       ntype=self.NparticleTypes,  representation=representation)
        axes.set_title('s='+str(s)+ ', e='+str(self.LEL[s]))
            
    
    def save_movie(self, directory_save=None, representation="colored",
                   showTemperature=True, plot_contact=False,color_map_contacts=None,
                   save_video_only = True, index_temperature_start=None):
        
        if directory_save is None :
         
            directory_save = '/Users/lkoehler/Documents/SelfAssembly2/Movies/'+(self.name.split('/')[0])
            subprocess.run(['mkdir', directory_save])
    
        
        try :
            if not save_video_only :
                open(directory_save + '/00000.png')
            else  :
                open(directory_save + '/Video.gif')
            file_already_exist = True
        except :
            file_already_exist = False
        
        if not file_already_exist :
            temperatures = range(self.parameters['Ntemperatures'])
            if index_temperature_start is not None:
                temperatures = temperatures[index_temperature_start:]
            if self.parameters["Nmes_per_T"] ==1 :
                indices = [self.parameters["Nsteps_per_T"]-1]
            else :
                freq_mes = self.parameters["Nsteps_per_T"] // self.parameters["Nmes_per_T"]
                i0 = self.parameters["Nsteps_per_T"] -1 - (self.parameters["Nmes_per_T"]-1)*freq_mes
                indices = np.arange(i0, self.parameters["Nsteps_per_T"]+1,freq_mes)
                # indices = np.arange(freq_mes, self.parameters["Nsteps_per_T"]+1,freq_mes)
                # indices[-1]-=1
            image_index = 0
            
            # x0,y0 = self.plot_ground_state(onlyGetPositionsCoordinates=True, getInfoOrigin=True)
            
            if self.parameters["saveTemperatures"] :
                Temperatures =  SF.fileToArray(self.directoryResults+'_Temperatures.txt')
                Temperatures = np.concatenate((Temperatures,[self.parameters["Temperature_quenching"]]))
            else :
                Temperatures= self.get_temperature_ramp()
                
            origin = self.plot_system_state(return_origin=True,do_the_plot=False)
            
            frames = []
            
            for t in temperatures :
                for i in indices : 
                    self.plot_system_state(representation=representation, 
                                           Temperature=t, index=i,
                                           origin=origin,plot_contact=plot_contact,
                                           color_map_contacts=color_map_contacts)
                    if showTemperature :#and self.parameters["saveTemperatures"]:
                        if self.dimension==2 :
                            plt.text(self.Lx-1, 0.1,'T=%.3f'%Temperatures[t], fontsize=12, ha='right')
                        else:
                            plt.title('T=%.2f'%temperatures[t], fontsize=12, ha='right')
                    
                    if not save_video_only:
                        plt.savefig(directory_save+'/'+str(image_index).zfill(5), dpi=300, bbox_inches='tight')
                        plt.close()
                    else :
                        plt.savefig(directory_save+'/temp.png', dpi=300, bbox_inches='tight')
                        plt.close()
                        image = imageio.imread(directory_save+'/temp.png')
                        frames.append(image)
                        
                    image_index+=1
        else :
            print('Files already exist in ',directory_save,'we wont overload them')
                
    
        imageio.mimsave(directory_save+'/Video.gif', # output gif
                        frames,          # array of input frames
                        duration = 10)   
    
ROUNDLEL=12
class ReadSeveralSystems():

    def __init__(self, name, local, processLocal=True, Nsimu=None, n_init=0, PKS=False):
        print("local",local)
        if local :
            directoryResults = '/Users/lkoehler/Documents/SelfAssembly2/ClusterMove/Results/'
        else :
            directoryResults = "/Users/lkoehler/Documents/LPTMSscratch/SelfAssembly2/Results/"
        if not processLocal :
            directoryResults = "/Scratch/lkoehler/SelfAssembly2/Results/"
            
        if PKS:
            directoryResults = '/Users/lkoehler/Documents/PKSdata/SelfAssembly2/Results/'
            if not processLocal : 
                directoryResults = '/data/biophys/lkoehler/SelfAssembly2/Results/'
        
        self.name = name
        self.processLocal = processLocal
        self.n_init = n_init
        
        print('Load ', directoryResults+self.name)
        
        if Nsimu is None :
            allSimuFound = False
            index = 0
            while not allSimuFound and index<1e6:
                try :   
                    LEL =SF.fileToArray(directoryResults+name+'/'+str(self.n_init+index)+'_LELchosen.txt')    
                    index+=1
                except :
                    allSimuFound = True
            self.Nsimu = index
        else :
            self.Nsimu = Nsimu
            
            
                
        print(directoryResults)
        file = open(directoryResults+name+'/'+str(n_init)+'input_parameters.json', 'r')
                
        
        """We assume parameters of the first file will stay the same"""
        self.parameters =  json.load(file)
        self.directoryResults = self.parameters["directoryResultsRead"]
        if self.directoryResults[-1]!='/':
            self.directoryResults = directoryResults+self.directoryResults.split('/')[-2]
    
        print(self.Nsimu, 'Nsimu in this folder')
        self.Nparticles = self.parameters['Nparticles']
        self.NparticleTypes = self.parameters['NparticleTypes']
        self.NparticlesTot = np.sum(self.Nparticles)
        self.Lx, self.Ly, self.Lz = self.parameters['Lx'],self.parameters['Ly'],self.parameters['Lz']
        
        self.Nsites = self.Lx*self.Ly*self.Lz
        self.density = self.NparticlesTot/self.Nsites
        # print(self.directoryResults)
        # if self.parameters['saveSnapshot']:
        #     self.sites_positions = SF.makeDictionaryPosition(self.directoryResults+"_sites_positions.txt")
        
        if processLocal :
            self.lattice = ReadLatticeParticle(self.parameters["lattice"], self.parameters["directoryInputsRead"])
        else :
            self.lattice = ReadLatticeParticle(self.parameters["lattice"], self.parameters["directoryInputs"])
        
        self.Particle = ParticleRepresentation(self.lattice)

        self.Naverage = int(self.parameters["Naverage"])
        self.Nsteps = int(self.parameters["Nsteps_per_T"]*self.parameters["Ntemperatures"])
        
        self.N1particleStructures = self.lattice.detail_all_N_structures[self.NparticleTypes][0]
        
        self.load_all_LELs()
        self.cluster_statistics_loaded = False #this willl become True if you call the fucntion load_cluster_statistics()
        
    def load_all_LELs(self) :
        t0 = time.time()
        self.LELs = []
        p = 0.01
        for index in range(self.Nsimu):
            if index/self.Nsimu > p :
                """Print information about the time it takes to load"""
                t = time.time()
                dt = str(round((t-t0)/60, 2))
                print(str(int(p*100))+'% loaded in '+dt+' min')
                p+=0.2
            LEL = SF.fileToArray(self.directoryResults+'/'+str(index+self.n_init)+'_LELchosen.txt') 
            #print(self.directoryResults+'/'+str(index+self.n_init)+'_LELchosen.txt')
            self.LELs.append(LEL[self.N1particleStructures:])
        self.LELs = np.array(self.LELs)
        
    def load_temperature_protocol(self, index):
        if self.parameters["saveTemperatures"]:
            all_temperatures = np.zeros(0)
            
            """Annealing temperatures"""
            Temperatures = SF.fileToArray(self.directoryResults+'/'+str(index)+'_Temperatures.txt')
            for t in Temperatures :
                all_temperatures = np.concatenate((all_temperatures, t*np.ones(self.parameters['Nsteps_per_T'])))
    
            """Quenching temperatures"""
            t0 = self.parameters["Temperature_quenching"]
            all_temperatures = np.concatenate((all_temperatures, t*np.ones(self.parameters['Nsteps_quenching'])))
        
            return(all_temperatures)
        else :
            print('Temperatures were not saved for this simulation')
            return(None)
        
        
    def load_temperatures(self, index):
        if self.parameters["saveTemperatures"]:
            Temperatures = SF.fileToArray(self.directoryResults+'/'+str(index)+'_Temperatures.txt')
            return(Temperatures)
        
    def load_statistics(self, toload, index):
        """toload is "c", "d", "Cab" 
        This gives the txt file without transformation
        c and d will be of dimension N1+Nstructures, ans not just Nstructures
        To remove the N1 first useless values in the c, use self.measure() """
        
        if self.parameters["saveEnergies"]:
            vector = SF.fileToArray(self.directoryResults+'/'+str(index)+'_'+toload+'.txt')
            return(vector)


    def load_energies(self, index):
        if self.parameters["saveEnergies"]:
            energies = SF.fileToArray(self.directoryResults+'/'+str(index)+'_Energies.txt')
            energies = energies/self.NparticlesTot
            return(energies)
        else :
            print('Energies were not saved for this simulation')
            return(None)
        
            
    def load_success_rate(self, index):
        if self.parameters["saveMCsuccess"]:
            success = SF.fileToArray(self.directoryResults+'/'+str(index)+'_MCsuccess.txt')
            #energies = energies/self.NparticlesTot
            return(success)
        else :
            print('MC success rate were not saved for this simulation')
            return(None)
        
        
    def find_indices(self, LEL=None, precision=1e-4):
        """Find the indices corresponding to the required LEL"""

        LEL = np.round(LEL,ROUNDLEL )
        
        allLEL = np.array(self.LELs)
        
        if len(allLEL)>0:
            
            searchZeros = np.abs(allLEL -LEL)
            indices = np.where((searchZeros<precision).all(axis=1))[0]
            
            if len(indices)==0:
                print('No simulation found for LEL = ', LEL)
            
            return(indices)    
        else:
            return([])
        
    def get_local_structures(self, index):
        
        s = ReadSystem(self.name+'/'+str(index+self.n_init), self.directoryResults[:-len(self.name)],
                       processLocal = self.processLocal)
        return(s.get_local_structures())
            
        
        
        
    def load_cluster_statistics(self):
        
        if not self.cluster_statistics_loaded:
     
            self.Nclusters = []
            self.AllclusterSizes = []
            self.AllclusterVacancies= []
            self.AllclusterVacanciesSize = []
            self.AllclusterInnerSurface = []
            self.AllclusterOuterSurface = []
            self.AllclusterAdresses = []
               
                
            for index in range(self.Nsimu) :
                self.Nclusters.append([])
                self.AllclusterSizes.append([])
                self.AllclusterVacancies.append([])
                self.AllclusterVacanciesSize.append([])
                self.AllclusterInnerSurface.append([])
                self.AllclusterOuterSurface.append([])
                self.AllclusterAdresses.append([])
                
                for imes in range(self.parameters['Ncluster_measurements']):
                    m = SF.matrixToArray(self.directoryResults+'/'+str(self.n_init+index)+"_Nfull_Nvac_SizeVac_InSurf_OutSurf_"+str(imes)+".txt")
                    Nfulls = list(m[:,0])
                    Nvacancies = list(m[:,1])
                    SizeVacancies = list(m[:,2])
                    NinnerSurface =list( m[:,3])
                    NouterSurface = list(m[:,4])
                    oneAdressPerCluster = [hex(adress) for adress in m[:,5]]
                    Ncluster = len(Nfulls)                
                    self.Nclusters[index]+=[Ncluster]
                    self.AllclusterSizes[index]+=Nfulls
                    self.AllclusterVacancies[index]+=Nvacancies
                    self.AllclusterVacanciesSize[index]+=SizeVacancies
                    self.AllclusterInnerSurface[index]+=NinnerSurface
                    self.AllclusterOuterSurface[index]+=NouterSurface
                    self.AllclusterAdresses[index]+=oneAdressPerCluster
            self.cluster_statistics_loaded
            
    def measure_size_distribution(self, LEL):
        """Returns the distribution array where distribution[k] is the number of observed clusters of size k+1"""
        indices = self.find_indices(LEL)
        sizes = []
        for i in indices :
            sizes+=list(self.AllclusterSizes[i])
        sizes = np.array(sizes)
        distribution = np.zeros(self.Nparticles)
        for s in sizes :
            distribution[s-1]+=1
        return(distribution)
            
    def measure(self, toMeasure, indices):
        """toMeasure = Cab, c, d"""
        
        elements = []
        for index in indices:
            if toMeasure != 'Cab':
                element = SF.fileToArray(self.directoryResults+'/'+str(self.n_init+index)+'_'+toMeasure+'.txt')
                element = element[self.N1particleStructures:]
                elements.append(element)
                
            else :
                element = SF.matrixToArray(self.directoryResults+'/'+str(self.n_init+index)+'_'+toMeasure+'.txt')
                element = element[self.N1particleStructures:,self.N1particleStructures:]
                # C++ does not normalize the matrix
                element = element/(self.Naverage*(3*self.Nsites))**2
                
                #symmetrize the matrix 
                CT = (element.copy()).T
                np.fill_diagonal(CT, np.zeros(len(CT)))
                element = element+CT                  
                
                elements.append(element)
        
        elements=np.array(elements, dtype=float)
        
        return(np.mean(elements, axis=0), np.std(elements, axis=0))
    
    

    
    
    def measure_c_boltzmann(self, designs, LEL, c_init):
            
        if c_init is None :
            c_init = np.zeros(28)
        
        
        """sum(c) = 1 and sum(c_face_i)=density/3"""
        constraints_conservation = designs.invariances
        d=self.density/3
        constraints_conservation_values_min = np.array([1,d,d,d,d,d,d])
        constraints_conservation_values_max = constraints_conservation_values_min
        
        """all the cs are between 0 and 1"""
        constraints_limits = np.eye(28)
        constraints_limits_values_min = np.zeros(28)
        constraints_limits_values_max = np.ones(28)
        
        """sum(c_surf) is between c_surf_gas and c_surf_perfect_bulk"""
        constraints_surface = np.concatenate((np.zeros(1), np.ones(6), np.zeros(21)))
        c_surf_perfect_bulk, c_surf_gas = self.c_surf_min(self.NparticlesTot),self.c_surf_max()
        
        
        
        
        # """stack constraints"""
        # A = np.concatenate((constraints_limits, constraints_conservation))
        # lb = np.concatenate((constraints_limits_values_min, constraints_conservation_values_min))
        # ub = np.concatenate((constraints_limits_values_max, constraints_conservation_values_max))
        # print(A, lb, ub)
        # scipy.optimize
        ConstraintsConservation = LinearConstraint(constraints_conservation, constraints_conservation_values_min, constraints_conservation_values_max)
        ConstraintsLimits = LinearConstraint(constraints_limits,constraints_limits_values_min, constraints_limits_values_max )
        ConstraintsLimitsSurf = LinearConstraint(constraints_surface,c_surf_perfect_bulk, c_surf_gas )
        
        
        def func(c):
            return(np.sum(c*LEL))
        
        res = minimize(func, c_init,  
                       constraints=[ConstraintsConservation,ConstraintsLimits, ConstraintsLimitsSurf] )
        
        c_boltzmann = res['x']
        
        return(c_boltzmann)#['x']
    
  
    def measure_cluster_information(self, toMeasure, indices, ponderate=True):
        """Returns average and std of some observable of the clusters, over all clusters and all realizations 
        with the simulation at index indices
        toMeasure = fullVolume, totalVolume, porosity, areaToVolumeRatio, sphericity, holeSize, holesPerParticle
        if ponderate, we ponderage the average by the number of particles per cluster """
        
        sizes = []
        vacancies = []
        outer_surface = []
        nvacancies = []
        sizeMax = []
        
        for i in indices :
            sizes+=list(self.AllclusterSizes[i])
            vacancies+=list(self.AllclusterVacanciesSize[i] )  
            nvacancies+=list(self.AllclusterVacancies[i] )  
            outer_surface += list(self.AllclusterOuterSurface[i])
            sizeMax += [np.max(self.AllclusterSizes[i])]
        
        sizes = np.array(sizes)
        vacancies = np.array(vacancies)
        nvacancies = np.array(nvacancies)
        outer_surface = np.array(outer_surface)
        sizeMax = np.array(sizeMax)
        
        if toMeasure == 'fullVolume':
            elements = sizes
        elif toMeasure == 'totalVolume' :
            elements = sizes + vacancies
        elif toMeasure == 'sizeMax' :
            elements = sizeMax
            
        elif toMeasure == 'porosity':
            elements = vacancies/ sizes
        elif toMeasure ==  'holeSize' :
            ids_holes = np.where(nvacancies>0)[0]
            if len(ids_holes>0):
                elements = vacancies[ids_holes] / nvacancies[ids_holes]
                sizes = sizes[ids_holes]
            else :
                elements = 0
        elif toMeasure == 'holesPerParticle' :
            elements = nvacancies / sizes
        elif toMeasure == 'savRatio':
            elements = outer_surface / (3*(sizes + vacancies))
            #this is a ratio between number of bonds, and there is 3 bonds per particle
        elif toMeasure == 'sphericity':
            volume = sizes + vacancies
            radius = lambda v: (-3 + np.sqrt(9+12*v)) / 6 
            small_surfaces = np.array([0,6,10,12,14,16,18,18, 20, 22])
            Nmax_out = np.zeros(len(volume))
            for i in range(len(volume)):
                v=volume[i]
                if v>9:
                    Nmax_out[i] = 12*radius(v)+6
                else :
                    Nmax_out[i] = small_surfaces[v]
            elements = 1/(outer_surface/Nmax_out)
        if toMeasure == 'sizeMax' or (not ponderate) :
            average, std = np.mean(elements), np.std(elements)
        
        else :
            average = np.sum(elements*sizes)/np.sum(sizes)
            std = np.sqrt( np.sum(elements**2*sizes)/np.sum(sizes) - average**2 )
            
        return(average, std)
    
    
    
    def plot_state(self,index, representation='colored', Temperature=None, 
                   index_system=None, figure=None, axes=None, return_data=False, 
                   do_the_plot=True, order=1,
                   plot_contact =False, return_contact=False,
                   color_map_contacts=None, 
                   origin=None):
        s = ReadSystem(self.name+'/'+str(index+self.n_init), self.directoryResults[:-len(self.name)],
                       processLocal = self.processLocal)
        data_image = s.plot_system_state(representation=representation, Temperature=Temperature, 
                            index=index_system, figure=figure, axes=axes, 
                            return_data=return_data, do_the_plot=do_the_plot, order=order,
                            plot_contact=plot_contact, return_contact=return_contact,
                            color_map_contacts=color_map_contacts, origin=origin)
        if return_data:
            return(data_image)
        
        
    def save_images_files(self, directoryImage, save_contacts=False):
        
        try :
            all_images = SF.fileToPickle(directoryImage+self.name)
        except :
            all_images = {} 
        p=0.01
        for i in range (self.Nsimu):            
            if (i/self.Nsimu)>p:
                print(str(round(i/self.Nsimu*100,1))+'%')
                p+=0.2
            if i+self.n_init not in all_images.keys():
                try :
                    data = self.plot_state(i, return_data=True, do_the_plot=False, return_contact=save_contacts)
                    plt.close()
                    all_images[i+self.n_init] = data
            
                    SF.pickleToFile(directoryImage+self.name, all_images)
                except: 
                    donothing=0
        
        

    
if __name__ == "__main__": 
    directory = '/Users/lara/Documents/SelfAssembly2/Results/'
    name = '230210TriangularHomomeric3particles_2/'
    name = '230210SquareHomomeric4particles/'
    name = '230210TriangularHomomeric4particles/'
    name = '230210SquareCrystal1p/'
    name = '230210TriangularCrystal4p/'
    name = '230210TriangularRandom4p/'
    name = '230210CubicCrystal2p/'
    name = '230210HexagonalCrystal3p_2/'
    
    # name = '230212Crystal3p/'
    # name = '230212CrystalAndMutation4/'
    # name = '230212Crystal/'
    
    # name = '230212CubicCrystal/'
    # name = '230212CubicCrystalMutation/'
    name = '230212CubicCrystal3p/'
    
    # name = '230213TestMartin/'
    
    # name = '230212TriangularRandom3p_af2_8/'
    name = '230217TestTriangle_2p_2/'
    name = '230217Triangle_ABC_AAC/'
    name = '230217Crystal_L010_d005/2'
    
    name = '230224TestCamembert/0'
    name = '230224TestCamemberthexagonal/0'
    name = '230224TestCamemberthexagonal2/0'
    
    name = '230228CamembertCalib_T0_10_n_10/0'
    name = '230302TestFiber2/0'
    name = '230310Camembert3D_chiral/0'
    
    
    
    directory = '/Users/lara/Documents/LPTMSscratch/SelfAssembly2/Results/'
    # s = ReadSystem(name, directory)
    # s.plot_system_state()
    
    name = '230310Camembert3D_achiral/1'
    # s2 = ReadSystem(name, directory)
    # s2.plot_system_state()
    # s.plot_energy_and_temperature()
    # s.plot_system_state('colored')
    
    
    # s.plot_system_state('arrowed')
    
    # s.plot_one_structure_s(edge=1, p1=0, p2=1, o1=1, o2=1)
    
    # s.plot_one_structure_s(edge=2, p1=0, p2=1, o1=2, o2=2)
    
    
    #%% Figure camembert 3D 
    # name='230302TestFiber'
    
    
    # # R = ReadSeveralSystems('230310Camembert3D_achiral', local=False)
    # f, axes = plt.subplots()
    # i=4
    # R.plot_state(i, axes=axes, figure=f)  
    # axes.text(0.03,0.9,'b)',transform = axes.transAxes, color='black', fontsize=16)
    # d = '/Users/lara/Documents/Thèse/Camembert/Figures/3Dcamembert_230310Camembert3D_achiral_'+str(i)+'.pdf'
    # plt.savefig(d, bbox_inches='tight')
    
    
    
    
    
    
    #%% debug
    
    # name='230302TestFiber2'
    
    # R = ReadSeveralSystems('230315Crystallite', local=False)
    # R1 = ReadSeveralSystems(name, local=True)
    # d1 = R1.NparticlesTot/R1.Nsites
    # t1 = np.sum(R1.load_statistics('c', 0)[7:])
    # R1.load_cluster_statistics()
    #%% Verifiy translation 
    
    if False :
        R = ReadSeveralSystems('230316VerifyTranslation', local=False)
        designs = DesignLEL(R.lattice, R.NparticleTypes)
        Nsystems = 5*2
        
        f,axes = plt.subplots(1,2)
        
        colors= ['darkorange', 'forestgreen', 'royalblue', 'deeppink']
        
        for i in range (4):
            try :
                lel = R.LELs[i*Nsystems]
                lelT = designs.cut_LEL(designs.translate_LEL(designs.LEL_predetermined(lel)))
                inds = R.find_indices(lel)
                indsT = R.find_indices(lelT, 1e-3)
                
                c,dc = R.measure('c', inds)
                cT,dcT = R.measure('c', indsT)
                
                axes[0].plot(lel, color=colors[i], marker='+')
                axes[0].plot(lelT, color=colors[i], linestyle='--', marker='x')
                axes[1].errorbar(np.arange(1,28), c[1:], dc[1:]/np.sqrt(len(inds)), color=colors[i], marker='+')
                axes[1].errorbar(np.arange(1,28),cT[1:],dcT[1:]/np.sqrt(len(indsT)), color=colors[i], linestyle='--',marker='x')
            
                E = np.sum(c*lel)
                ET = np.sum(cT*lelT)
                # print("E = "+str(E))
                # print("E' = "+str(ET))
                
                print("E'-E = "+str(ET-E))
                shift = designs.get_energy_shift_for_translation(designs.LEL_predetermined(lel), 3*R.Nsites, R.NparticlesTot)
                print("E'-E predicted = ",shift)
            except:
                donothing=0
    
        axes[0].set_xlabel('Structures')
        axes[0].set_ylabel('LEL')
        axes[1].set_ylabel('c')
        plt.tight_layout()

 
    #%%
    
color_map_camembert = {1: 'tomato',
 2: 'tomato',
 3: 'tomato',
 4: 'tomato',
 5: 'tomato', 
 6: 'tomato',
 
 0: 'white',
 
 7: 'skyblue',
 14: 'skyblue',
 20: 'skyblue',
 
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
    