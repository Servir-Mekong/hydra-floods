from __future__ import print_function,division
import itertools
import numpy as np
from scipy import ndimage
from sklearn.cluster import KMeans
import multiprocessing as mp


class starfm(object):
    
    def __init__(self,img,fslv,cslv,windowSize=33,nClusters=3):
        self.img = img
        self.fslv = fslv
        self.cslv = cslv
        
        if windowSize%2 == 0:
            print("Selected window size of {0} is an even number, \
                  forcing window size to be odd number by subtracting 1")
            windowSize = windowSize-1
        self.w = windowSize
        
        if type(nClusters) != int:
            print("Number of clusters specified is not an interger value, \
                  forcing vlaue to int type")
            nClusters = int(nClusters)
            
        self.m = nClusters
        self.sigma = fslv.std()
        self.center= int(np.around(self.w/2.))-1
        
    def __str__(self):
        return "Object for applying the STARFM algorithm"
        
    def clusterWindow(self,windowArray):
        shape = windowArray.shape
        # reformat data array to 1-d for clutering
        stack = windowArray.ravel()
        # initialize k-means algorithm
        kmeans = KMeans(init='k-means++', n_clusters=2, n_init=10)
        # apply kmeans algorithm
        kmeans.fit(stack)
        clusters = kmeans.cluster_centers_.squeeze()[:, 0]
        # get labeled data based on kmeans centroids
        c_idx = kmeans.labels_
        cluster = c_idx.reshape(shape)
        val = cluster[self.center,self.center]
        
        # return reshaped 2-D array of clustered pixels
        return np.ma.masked_where(cluster!=val,windowArray)
        
    def likePixels(self,windowArray):
        threshold = self.sigma * (2/float(self.m))
        center = windowArray[self.center,self.center]
        diff = np.abs(windowArray-center)
        like = diff<=threshold
        
        # return masked array of like pixels
        return np.ma.masked_where(like==False,windowArray)
    
    def calcWeights(self,imgWindow,slaveWindow):        
        noMask = np.where(slaveWindow.mask==False)
        weights = np.zeros_like(imgWindow)
        centerIdx = self.center
        iters = noMask[0].size
        for i in range(iters):
            
            ri = stats.pearsonr([imgWindow[self.center,self.center],imgWindow[noMask[0][i],noMask[1][i]]],
                                [slaveWindow[self.center,self.center],slaveWindow[noMask[0][i],noMask[1][i]]])[0]
            di = 1 + (np.sqrt((centerIdx-noMask[0][i])**2 + (centerIdx-noMask[1][i])**2) / (self.w/2.))
            
            Di = (1-ri)*di
            
            weights[noMask[0][i],noMask[1][i]] = Di
            
        outWeights = (1/weights)/weights.sum()
        
        return np.ma.masked_where(slaveWindow.mask==True,outWeights)
        
    
    def calcConversion(self,imgWindow,slaveWindow):
        c = imgWindow[np.where(slaveWindow.mask==False)].ravel()
        f = slaveWindow[np.where(slaveWindow.mask==False)].ravel()
        coefs = stats.linregress(c,f)
        
        return coefs[0]
    
    def _mapper(self,args):
        i,j,s,f,c = list(args)
        clusterMask = self.clusterWindow(slvWindow)
        likeMask = self.likePixels(clusterMask)
        distanceWeights = self.calcWeights(imgWindow,likeMask)
        V = self.calcConversion(imgWindow,likeMask)

        conversion = distanceWeights*V*imgWindow

        downscaledImg[i,j] = self.slv[i,j] + conversion.sum()
        
        return 
        
    
    def downscale(self):
        yRatio = self.slv.shape[0]/float(self.img.shape[0])
        xRatio = self.slv.shape[1]/float(self.img.shape[1])
        
        imgRescale = ndimage.zoom(self.img,[yRatio,xRatio],order=0)
        downscaledImg = np.zeros_like(imgRescale)
        
        args = []
        
        for i in range(imgRescale.shape[0]):
            ymin = int(np.around(i-(self.w/2.)))
            ymax = int(np.around(i+(self.w/2.)))
            if ymin<0:
                ymin=0
            if ymax>imgRescale.shape[0]:
                ymax = imgRescale.shape[0]
            for j in range(imgRescale.shape[1]):
                xmin = int(np.around(j-(self.w/2.)))
                xmax = int(np.around(j+(self.w/2.)))
                if xmin < 0:
                    xmin=0
                if xmax > imgRescale.shape[1]:
                    xmax = imgRescale.shape[1]
                if imgWindow.mask[i,j] != True:
                    fslvWindow = self.fslv[ymin:ymax,xmin:xmax]
                    cslvWindow = self.cslv[ymin:ymax,xmin:xmax]
                    imgWindow = imgRescale[ymin:ymax,xmin:xmax]

                    arg ={'i':i, 'j':j, 'Image':imgWindow, 'HiresSlave':fslvWindow,
                          'CoarseSlave':cslvWindow}
                    args.append(arg)
                
        processors = mp.cpu_count()
        pool = mp.Pool(processors)
        results = pool.map_async(self._mapper,args)
        results.close()
        out = result.join()
        print(out)
        
        return downscaledImg
    


def strum():
    print("Algorithm not yet implemented...check back soon")
    
    return


def ustarfm():
    print("Algorithm not yet implemented...check back soon")
    
    return

def bathtub(dem,fraction,probablistic=False,demStdDev=4.2):
    
    
    
    return