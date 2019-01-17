from __future__ import print_function,division
import itertools
import numpy as np
import xarray as xr
import multiprocessing as mp

import rastersmith as rs


class STARFM(object):
    def __init__(self,gr,windowSize=33,A=0.5):
        self.gr = gr

        if windowSize%2 == 0:
            print("Selected window size of {0} is an even number, \
                  forcing window size to be odd number by subtracting 1")
            windowSize = windowSize-1
        self.w = windowSize

        self.center= int(np.around((self.w-1)/2.))
        self.A = A
        self.buff = int((self.w - 1) / 2)

    def calcWeights(self,imgWindow,slaveWindow,logWeight=False):
        weights = imgWindow.copy()
        weights[:,:,:,:,:] = 1

        yPos = np.array(np.where(weights[:,:,0,:,0])[0]).reshape(weights.shape)
        xPos = np.array(np.where(weights[:,:,0,:,0])[1]).reshape(weights.shape)

        dijk = np.sqrt((self.center-xPos)**2 + (self.center-yPos)**2)
        Dijk = 1 + (dijk/self.A)

        Sijk = np.abs(imgWindow-slaveWindow)

        Tijk = self.tDiff
        if Tijk < 1:
            Tijk = 1

        if logWeight:
            Cijk = np.log(Sijk) * np.log(Tijk) * np.log(Dijk)
        else:
            Cijk = Sijk * Tijk * Dijk

        Wijk = (1/Cijk) / np.sum(1/Cijk)

        return Wijk

    def getDownscaled(self,coarse,fine,slave,pos):
        self.yEdges = [pos[0]-self.buff,pos[0]+self.buff]
        self.xEdges = [pos[1]-self.buff,pos[1]+self.buff]
        # print(pos)

        if ((pos[0] < self.center) | (pos[1] < self.center)) |\
           ((self.gr.dims[0] - pos[0] < self.center) | (self.gr.dims[1] - pos[1] < self.center)):

           result = np.nan

        else:
            # M_t0
            cWindow = coarse.isel(dict(lat=slice(self.yEdges[0],self.yEdges[1]),
                                       lon=slice(self.xEdges[0],self.xEdges[1]))
                                 )
            # L_t0
            fWindow = fine.isel(dict(lat=slice(self.yEdges[0],self.yEdges[1]),
                                     lon=slice(self.xEdges[0],self.xEdges[1]))
                               )
            # M_tk
            sWindow = slave.isel(dict(lat=slice(self.yEdges[0],self.yEdges[1]),
                                      lon=slice(self.xEdges[0],self.xEdges[1]))
                                )

            # Wijk calculation
            weights = self.calcWeights(fWindow,sWindow)#.where(clusters==0)

            # predict the
            outWindow = (fWindow.isel(time=0)/np.sum(weights)) + cWindow.isel(time=0) - sWindow.isel(time=0)

            cntrPx = outWindow.isel(dict(lat=self.center,lon=self.center)).values

            try:
                result = cntrPx.astype(np.float)
            except TypeError:
                result = np.nan

        return result


    def apply(self,coarse,fine,slave,parallel=False):
        self.tDiff = np.abs((coarse.attrs['date']-slave.attrs['date']).days)

        outRast = fine.copy()

        iters = [[i,j] for i in range(self.gr.dims[0]) for j in range(self.gr.dims[1])]

        if parallel:
            ncores = mp.cpu_count() - 1
            p = mp.Pool(ncores)
            args = ((coarse,fine,slave,x) for x in iters)
            pred = p.map(self._applyParallel,args)

        else:
            pred = map(lambda x: self.getDownscaled(coarse,fine,slave,x),iters)

        out = np.array(list(pred)).reshape(outRast.shape)
        outRast[:,:,:,:,:] = out

        outRast.coords['time'] = coarse.coords['time']

        return outRast

    def _applyParallel(self,args):
        coarse,fine,slave,x = args
        return self.getDownscaled(coarse,fine,slave,x)


class STRUM(object):
    def __init__():
        print("Algorithm not yet implemented...check back soon")

        return


class USTARFM(object):
    def __init__():
        print("Algorithm not yet implemented...check back soon")

        return

class Bathtub(object):
    def __init__(self,gr,probablistic=False,demStdDev=4.2,nIter=100):
        self.gr = gr

        return

    def selGridByCell(self,i,j,fWater,hand):
        cell = fWater.isel(dict(lat=i,lon=j))
        x = cell.coords['lon'].values
        y = cell.coords['lat'].values
        res = fWater.attrs['resolution']
        return hand.sel(dict(lat=slice(y+res,y-res),lon=slice(x-res,x+res)))

    def grids2raster(self,fWater,raster,grids,yRange,xRange):
        remap = raster.copy()
        remap[:,:,:,:,:] = -32768
        # print(remap.shape)
        cnt = 0
        # print(grids[0].shape)
        for i in yRange:
            for j in xRange:
                cell = fWater.isel(dict(lat=i,lon=j))
                x = cell.coords['lon'].values
                y = cell.coords['lat'].values
                res = fWater.attrs['resolution']
                remap.sel(dict(lat=slice(y+res,y-res),lon=slice(x-res,x+res)))[:,:,:,:,:] = grids[cnt].values
                cnt += 1

        remap.coords['band'] = ['water','mask']

        return remap

    def fillGrid(self,grid,fraction,seed=None):

        water = grid.copy()
        if np.isnan(fraction):
            final = water.copy()
            final[:,:,:,:,:] = 0

        else:
            water.values = np.zeros(water.shape)

            elv = range(15)

            results = list(map(lambda x: xr.where(grid==x,x+1,water),elv))
            result = xr.concat(results,dim='z')

            depth = result.sum(dim='z')
            depth = xr.where(depth==0,30,depth)

            simFracs = list(map(lambda x: depth.where(depth<=x).count().values / float(depth.size),elv))
            diff = np.array(list(map(lambda x: simFracs[x] - fraction, elv)))
            absDiff = np.abs(diff)

            maxDepth = absDiff.argmin()
            resid = diff[maxDepth-1]

            simWater = depth <= maxDepth

            dFlat = depth.isel(dict(time=0,band=0))

            if resid > 0:
                rIdx = np.where(dFlat == maxDepth)
                rand = np.zeros(dFlat.shape) + rIdx[0].size

                rand[rIdx[0],rIdx[1]] = np.random.random(rIdx[0].size) * rIdx[0].size

                thresh = rIdx[0].size * (fraction / simFracs[maxDepth])

                residWater = rand<thresh

            else:
                residWater = np.zeros(dFlat.shape)

            residWater = xr.DataArray(residWater,coords={'lat':simWater.coords['lat'].values,
                                                         'lon':simWater.coords['lon'].values},
                                                 dims=['lat','lon'])

            final = simWater.astype(np.bool) | residWater.astype(np.bool)
            final = final.expand_dims('z').transpose('lat','lon','z','band','time')

        return final


    def apply(self,fWater,hand,seedRaster=None,parallel=False):

        xDim = range(len(fWater.coords['lon'].values))
        yDim = range(len(fWater.coords['lat'].values))

        idxs = [[i,j] for i in yDim for j in xDim]

        fracs = list(map(lambda i: fWater.isel(dict(lat=i[0],lon=i[1])).values.flatten()[0], idxs))
        grids = list(map(lambda i: self.selGridByCell(i[0],i[1],fWater,hand),idxs))

        args = list(zip(*[grids,fracs]))

        if parallel==True:
            ncores = mp.cpu_count() - 1
            p = mp.Pool(ncores)
            waterGrids = p.map(self._applyParallel,args)

        else:
            waterGrids = list(map(lambda x: self.fillGrid(x[0],x[1]),args))

        waterMap = self.grids2raster(fWater,hand,waterGrids,yDim,xDim)
        waterMap = waterMap.isel(band=slice(0,1))
        waterMap.coords['bands'] = 'water'

        return waterMap

    def _applyParallel(self,args):
        grid,frac = args
        return self.fillGrid(grid,frac)
