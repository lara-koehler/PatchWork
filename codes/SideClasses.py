#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:19:08 2020

@author: lara
"""



import numpy as np 
import scipy as sp 
import pickle
import matplotlib
import colorsys

from time import sleep
# import pandas as pd

# import psutil



import ast


class SideFunctions : 
    def __init__(self):
        self.nothing = 0
        
        


        
    def listToString(aList):
        res=""
        for i in range (len(aList)):
            res+=str(aList[i])
        return res
    
    def stringToList(aString):
        res=[]
        for s in aString :
            res.append(int(s))
        return(res)
    
    def stringToList2(aStringSeparatedWithSpaces):
        res=[]
        line = aStringSeparatedWithSpaces.split()
        for i in range(len(line)):
            res.append(float(line[i]))
        return(np.array(res))
    
    
    
    """!!!"""
    """ READING """
    """!!!"""
    
    def fileToMatrixOfArrays(filename):
        file = open(filename, 'r')
        line = file.readline()
        line=line.split('] [')
        res = []
        value=line[-1]
        i=0
    
        while (value !="#END" and i<1e10) :
            res.append([])
            #print('line', line)
            for j in range (len(line)):
                vector=line[j]
                # print('vector', vector)
                if vector[0]=='[':
                    vector=vector[1:]
                if vector[-3:][0]==']':
                    vector = vector[:-3]
                #print(vector, vector[-3:], vector[-3:][0]=='] ')
                res[i].append(SideFunctions.stringToList2(vector))
                
                
            line=file.readline()
            line=line.split('] [')
            value=line[-1]
            res[i]= np.array(res[i], dtype=object)
            i+=1 
        file.close()
        return (res)
    
    

    def tryeval(val):
      try:
        val = ast.literal_eval(val)
      except ValueError:
        pass
      return val

    
    def matrixToArray(filename):
        
        
        file = open(filename, 'r')
        line = file.readline()
        line=line.split()
        res = []
        value=line[-1]
        i=0
        isFloat = True
        isInt = True
    
        while (value !="#END" and i<1e10) :
            res.append([]);
            
            for j in range (len(line)):
                # try :
                res[i].append(SideFunctions.tryeval(line[j]))
                # except :
                #     isInt = False
                #     try :
                #         res[i].append(float(line[j]))
                #     except :
                #         isFloat = False
            
                
            line=file.readline()
            line=line.split()
            value=line[-1]
            
            res[i]=np.array(res[i])
            
            i+=1 
        #print("isInt or isFloat: "+str(isInt or isFloat))
        if isFloat or isInt:
            l1 = len(res)
            if l1>0 :
                l2 = len(res[0])
                res2 = np.zeros((l1,l2), type(res[0][0]))
                for i in range (l1):
                    for j in range(l2):
                        res2[i][j]=res[i][j]
            else:
                res2=res
        else:
            res2=res
        file.close()
        return (res2)
    
    def fileToPickle(filename):
        if filename[-3:]!='pkl':
            filename+='.pkl'
            

        file = open(filename, "rb")
        myObject = pickle.load(file)

        file.close()
        return(myObject)
    
    def fileToDict(filename):
        res = {}
        b=0

        f=open(filename, 'r')
        line = f.readline()
        while line!= "#END":
            line=line.split()
            key = line[0]
            attribute = line[1]
            line = f.readline()
            isInt, isFloat = True, True
            try :
                if key !='sequence':
                    attribute = int(attribute)
                else :
                    isInt = False
            except:
                isInt = False
            if (not isInt):
                try :
                    if key !='sequence':
                        attribute = float(attribute)
                    else :
                        isFloat=False
                except :
                    isFloat = False
                    
            res[key]=attribute
            # print(key)
        return(res)
        
        
    def fileToArray(filename):
        file = open(filename, 'r')
        line = file.readline()
        line=line.split()
        res = []
        value=line[0]
        i=1
        isInt, isFloat = True, True
        while (value !="#END" and i<1e10) :
            isInt, isFloat = True, True
            try :
                value = int(value)
            except :
                isInt = False
            if not isInt :
                try : 
                    value = float(value)
                except :
                    isFloat = False
            
            res.append(value)
            value=line[i]
            
            i+=1
        if isInt or isFloat:
            res=np.array(res)
        file.close()
        return (res)
    
    def fileToArrayFast(filename):
        a = pd.read_csv(filename, header=None, sep=' ')
        b = np.array(a.values[0][:-1])
        return(b)

    def makeDictionaryPosition(filename):
        file = open(filename, 'r')
        l=file.readline()
        positionDictionary = {}
        count=0
        while (l!="#END" and count<1e10):
            count+=1
            current_site = l.split()
            adress = current_site[0]
            #print(adress)
            position=[int(current_site[1]),int(current_site[2]), int(current_site[3])]
            positionDictionary[adress]=position;
            l=file.readline()
        #print(count)
        file.close()
        return(positionDictionary)
    
    def readClusterPosition(filename):
        file = open(filename, 'r')
        l=file.readline()
        clusters=[]
        count=0
        while (l!="#END" and count<1e10):
            current_cluster = l.split()
            clusters.append([])
            for i in range (len(current_cluster)):
                clusters[count].append(current_cluster[i])
            l=file.readline()
            count+=1
        #print(count)
        file.close()
        return(clusters)
        
    
    
    
    
    """ !!! """
    """ WRITING """
    """ !!! """
    
    
    # def wait_if_has_handle(fpath):
    #     """check if a file is opened by another process"""
    #     wait = False
    #     for proc in psutil.process_iter():
    #         try:
    #             for item in proc.open_files():
    #                 if fpath == item.path:
    #                     return True
    #         except Exception:
    #             pass
        
    #     if wait :
    #         sleep(0.1)
    
        
    
    def matrixToFile(filename, matrix):
        SideFunctions.wait_if_has_handle(filename)
            
        file = open(filename, 'w')
        string =""
        for i in range (len(matrix)):
            line = matrix[i]
            for j in range (len(line)):
                string+=str(matrix[i][j])+" "
            string+="\n"
        string+="#END"
        file.write(string)
        file.close()    
        
    def threeDArrayToFile(filename, matrix):
        SideFunctions.wait_if_has_handle(filename)
           
        file = open(filename, 'w')
        string =""
        for i in range (len(matrix)):
            line = matrix[i]
            for j in range (len(line)):
                string+='['
                theArray = matrix[i][j]
                for k in range (len(theArray)):
                    string+=str(matrix[i][j][k])+" "
                string+=']'
            string+="\n"
        string+="#END"
        file.write(string)
        file.close()
        
        
        
    def pickleToFile(filename, obj):
        # SideFunctions.wait_if_has_handle(filename)
        
        filename+='.pkl'
        file = open(filename, "wb")
        pickle.dump(obj, file)
        file.close()
        
        

    
    def arrayToFile(filename, array):
        SideFunctions.wait_if_has_handle(filename)
            
        file = open(filename, 'w')
        string =""
        for i in range (len(array)):
            string+=str(array[i])+" "
        string+="#END"
        file.write(string)
        file.close()


    def dictToFile(dictionnary, filename):
        SideFunctions.wait_if_has_handle(filename)
           
        f=open(filename, 'w')
        for key in dictionnary :
            f.write(key+" "+str(dictionnary[key])+"\n")
        f.write("#END")
        f.close()



        
    def normalizeArray(array):
        maxi = np.max(array)
        mini = np.min(array)
        res = array-mini
        if (maxi-mini)!=0:
            res= res/(maxi-mini)
        return(res)
        
        
        
    def propagate_error(result, values,errors):
        resultError = 0
        for i in range (len(values)):
            resultError += (errors[i]/values[i])**2
            print(errors[i]/values[i])
        print(result*np.sqrt(resultError))
        print("\n")
        return(result*np.sqrt(resultError))
        
        
        
    def converged_measure(aList, ratio=20):
        if len(aList)<3:
            return(False)
        lastDelta = np.linalg.norm(aList[-1]-aList[-2])
        norm = np.mean(np.linalg.norm(np.array(aList), axis=0))
        totalDelta = np.linalg.norm(aList[-1]-aList[0])
        
        #print(lastDelta/totalDelta, 1./ratio)
        return(lastDelta/totalDelta < 1./ratio)
    
    
    def moving_average(a, n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
        
    
    
    def arrayInArrayList (array, arrayList):
        n = len(arrayList)
        if n==0:
            return(False)
        else:
            for i in range(n):
                if (array==arrayList[i]).all():
                    return(True)
            return(False)
        
    def search_change_sign(x,y):
        n = len(x)
        previous_sign = np.sign(y[0])
        xchanges = []
        for i in range(1,n):
            if np.sign(y[i]) != previous_sign:
                x1, x2 = x[i-1], x[i]
                y1, y2 = y[i-1], y[i]
                if y2!=0:
                    r = np.abs(y1)/np.abs(y2)
                    # (x0-x1)/(x2-x0) = r (thales theorem)
                    x0 = (r*x2+x1)/(1+r)
                else :
                    x0 = x2
                # print('y1, y2, x1, x2, x0', y1, y2, x1, x2, x0)
                xchanges.append(x0)
            previous_sign = np.sign(y[i])
        return(xchanges)
    
    
    def search_change_sign_fluctu2(x,y):
        n = len(x)
        previous_sign = np.sign(y[0])
        xchanges = []
        for i in range(1,n):
            if np.sign(y[i]) != previous_sign:
                x1, x2 = x[i-1], x[i]
                y1, y2 = y[i-1], y[i]
                if y2!=0:
                    r = np.abs(y1)/np.abs(y2)
                    # (x0-x1)/(x2-x0) = r (thales theorem)
                    x0 = (r*x2+x1)/(1+r)
                else :
                    x0 = x2
                # print('y1, y2, x1, x2, x0', y1, y2, x1, x2, x0)
                xchanges.append(x0)
            previous_sign = np.sign(y[i])
            
        n_candidates = len(xchanges)
        if n_candidates >1 :
            distances = np.zeros(n_candidates)
            for i in range (n_candidates):
                if i==0 :
                    distances[i]= np.abs(xchanges[0]-x[0])+np.abs(xchanges[1]-xchanges[0])
                elif i!=n_candidates -1 :
                    distances[i] = np.abs(xchanges[i]-xchanges[i-1])+np.abs(xchanges[i]-xchanges[i+1])
                else :
                    distances[i]= np.abs(xchanges[i]-x[-1])+np.abs(xchanges[i]-xchanges[i-1])
            print('distances', distances)
            xchange = xchanges[np.argmax(distances)]
        else :
            xchange = xchanges[0]
            
        return(xchange)
    
    
    
    def search_change_sign_fluctu(x, y):
        nmax = 15
        single_zero_found = False
        ntest = 0
        nav = 1
        xchange = 0
        
        while (not single_zero_found) and (ntest < nmax) :
            # print(nav, ntest)
            x, y = np.array(x), np.array(y )
            if nav ==1 :
                averaged_y, averaged_x = y, x
            else :
                averaged_x, averaged_y = SideFunctions.average_signal(x, y, nav)
            xchanges = SideFunctions.search_change_sign(averaged_x, averaged_y)
            
            if len(xchanges) !=1:
                nav+=1
                ntest+=1
            else :
                xchange = xchanges[0]
                single_zero_found = True
        return(xchange)
    
    def average_signal(x, y, nav):
        averaged_y = y[:-nav+1]
        for k in range (1, nav):
            iend = -nav+1+k
            if iend !=0:
                averaged_y += y[k:iend]
            else :
                averaged_y += y[k:]
        averaged_x = x[nav//2-1:-nav//2]
        return(averaged_x, averaged_y)
    
    def scale_lightness(colorName, addedLightness):
        """AddedLightness is a number betweeen 0 and 1 """
        scale_l = 1+addedLightness
        rgb = matplotlib.colors.ColorConverter.to_rgb(colorName)
        # convert rgb to hls
        h, l, s = colorsys.rgb_to_hls(*rgb)
        # manipulate h, l, s values and return as rgb
        return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)
    
    
    def devide_line(x, y, N):
        newx, newy = [], []
        for i in range (len(x)-1):
            x1,y1 = x[i], y[i]
            x2,y2 = x[i+1], y[i+1]
            newx+= list(np.linspace(x1, x2, int(round( N/len(x))), endpoint=False))
            newy+= list(np.linspace(y1, y2, int(round( N/len(x))), endpoint=False))
        newx+=[x2]
        newy+=[y2]
        return(np.array(newx), np.array(newy))
    
    
    
    def find_middle_bigger_zone(emptylist, L):
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
        
        
            
    
class myLinearAlgebra():
    
    def __init__(self):
        self.nothing =  0
        
    def complement(subspace, precision=1e-7):
        """all that is below precisoin*maximalvalue is considered to be zero"""
        dimension = len(subspace[0])
        size = len(subspace)
        space = np.concatenate((subspace, np.zeros(((dimension-size), dimension))))
        complement = sp.linalg.null_space(space, rcond=precision).T
        return(complement)
        
        
    def isIncluded(vector, subspace, precision= 1e-7):
        """all that is below precisoin*maximalvalue is considered to be zero"""
        """a vector is indluded in a subspace if it has no components in its complementary"""
        theComplement = myLinearAlgebra.complement(subspace, precision)
        allSpace = np.concatenate((subspace, theComplement))
        components = np.linalg.pinv((allSpace).T).dot(vector.T)
        componentsOutside = components[len(subspace):]
         
        return ((np.abs(componentsOutside/max(components))<precision).all())
    
    
    
    def projection(vector, subspace, precision=1e-7):
        
        theComplement = myLinearAlgebra.complement(subspace, precision)
        allSpace = np.concatenate((subspace, theComplement))
        components = np.linalg.pinv((allSpace).T).dot(vector.T)
        projection = components[:len(subspace)]

        return(projection)
    
    
    def orthogonalVector(vector):
        """Will genereate an orthogonal vector to the given one, of same dimension"""
        N=len(vector)
        res = np.ones(N)
        if not (vector==0).all():
            nonzeroIndex = np.where(vector!=0)[0][0]
            
            xi = vector[nonzeroIndex]
            res[nonzeroIndex] = -(np.sum(vector)-xi)/xi
        return(res)
            
    
    def areColinear(vector1, vector2, precision=1e-7):
        precision0 = max(np.max(vector1), np.max(vector2))*precision
        n = len(vector2)
        k = 0
        
        for i in range(n):
            if abs(vector2[i])<precision0:
                if abs(vector1[i])>=precision0:
                    """v1 is non zero where v2 is zero"""
                    return(False)
            else:
                current_k = 1.*vector1[i]/vector2[i]
                
                if k==0:
                    k=current_k
                elif abs(k-current_k)>precision:   
                    """the proportionality coefficient is not the same than before"""
                    return(False)
        return(True)
        
        
        
    
        
        



            
        
