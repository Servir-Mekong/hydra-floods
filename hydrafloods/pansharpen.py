from __future__ import print_function,division
import numpy as np
import pywt
from scipy import ndimage,signal
from matplotlib import colors
from sklearn.decomposition import PCA
from . import io,utils

def _interp(r,g,b,i,method='cubic'):
    code = {'cubic':3,'bilinear':1,'nearest':0}
    R = ndimage.zoom(r,2,order=code[method])
    G = ndimage.zoom(g,2,order=code[method])
    B = ndimage.zoom(b,2,order=code[method])
    I = ndimage.zoom(i,2,order=code[method])
    
    return R,G,B,I

def _rescale(data,maxVal,minVal,stretch=False,inverse=False):
    if stretch:
        if inverse:
            out = ((maxVal-minVal) * (maxVal-data)) / (maxVal-minVal)
        else:
            out = ((maxVal-minVal) * (data-minVal)) / (maxVal-minVal)
    else:
        if inverse:
            out = (maxVal-data) / (maxVal-minVal)
        else:
            out = (data-minVal) / (maxVal-minVal)
            
    return out


def gramSchmidt():
    
    return


def hsv(R,G,B,I,P,weight=0.3,S=None):
    '''Function to apply hue-saturation-intensity pansharpening algorithm
    
    Args:
        R (ndarray): 2-d array representing the red channel of image.
        G (ndarray): 2-d array representing the green channel of image.
        B (ndarray): 2-d array representing the blue channel of image.
        I (ndarray): 2-d array representing the infrared channel of image.
        P (ndarray): 2-d array representing the panchromatic band of image.
        weight (float): Value representing the weight of IR channel in the pan image. 

    Returns:
        panImg (ndarrayd): 3-d array of the pan sharpened RGB image
    '''
    if not (S is None):
        P = utils.hist_match(P,S)
    
    R,G,B,I = _interp(R,G,B,I)
    rgb = np.dstack([R,G,B])
    hsv = colors.rgb_to_hsv(rgb)
    
    data_max, data_min = hsv[:,:,2].max(),hsv[:,:,2].min()
    intensity = _rescale((P - I*weight),data_max,data_min)

    outHsv = np.dstack([hsv[:,:,0],hsv[:,:,1],intensity])

    panImg = colors.hsv_to_rgb(outHsv)
    
    return panImg


def wavelet(R,G,B,I,P):
    
    interped = _interp(R,G,B,I)
    
    multiSpec = []
    panCoeffs = pywt.dwt2(P, 'haar')
    
    for i in range(len(interped)):
        wave = pywt.dwt2(interped[i], 'haar')
        coeffs = (wave[0],panCoeffs[1])
        multiSpec.append(pywt.idwt2(coeffs, 'haar'))
    
    panImg = np.dstack(multiSpec)
    
    return _rescale(panImg,panImg.max(),panImg.min())
    

def brovey(R,G,B,I,P,weights=None,S=None):
    '''Function to apply brovey pansharpening algorithm
    
    Args:
        R (ndarray): 2-d array representing the red channel of image.
        G (ndarray): 2-d array representing the green channel of image.
        B (ndarray): 2-d array representing the blue channel of image.
        I (ndarray): 2-d array representing the infrared channel of image.
        P (ndarray): 2-d array representing the panchromatic band of image.
        weights (list, ndarray): List representing the weights of channels in 
                                 the pan image in order or RGBI. 

    Returns:
        panImg (ndarray): 3-d array of the pan sharpened RGB image
    '''
    if not (S is None):
        P = utils.hist_match(P,S)
        
    R,G,B,I = _interp(R,G,B,I)
    panImg = np.zeros([I.shape[0],I.shape[1],4])
    
    if weights:
        if type(weights) == str:
            try:
                DNF = (P - (weights[3] * I)) / ((weights[0] * R) + (weights[1] * G) + (weights[2] * B))

                panImg[:,:,2] = B * DNF
                panImg[:,:,1] = G * DNF
                panImg[:,:,0] = R * DNF
                panImg[:,:,3] = I * DNF
            except IndexError:
                raise IndexError('weights list has wrong number of elements')
        else:
            raise TypeError('weights keyword requires a list')
    
    else:
        panImg[:,:,0] = _rescale(R / ((B + G + R) * P),R.max(),R.min())
        panImg[:,:,1] = _rescale(G / ((B + G + R) * P),G.max(),G.min())
        panImg[:,:,2] = _rescale(B / ((B + G + R) * P),B.max(),B.min())
        panImg[:,:,3] = _rescale(I / ((B + G + R) * P),I.max(),I.min())
    
    return panImg


def pca(R,G,B,I,P,n_components=4,S=None):
    '''Function to apply principal component analysis pansharpening algorithm
    
    Args:
        R (ndarray): 2-d array representing the red channel of image.
        G (ndarray): 2-d array representing the green channel of image.
        B (ndarray): 2-d array representing the blue channel of image.
        I (ndarray): 2-d array representing the infrared channel of image.
        P (ndarray): 2-d array representing the panchromatic band of image.
        n_components (int): Integer value representing the number of components to solve 

    Returns:
        panImg (ndarray): 3-d array of the pan sharpened RGB image
    ''' 
    if not (S is None):
        P = utils.hist_match(P,S)
   
    R,G,B,I = _interp(R,G,B,I)
    stack = np.vstack([R.ravel(),G.ravel(),B.ravel(),I.ravel()]).T
    
    mu = np.mean(stack, axis=1)
    
    pca_ = PCA(n_components=n_components)
    pca_.fit(stack)
    
    transform = pca_.transform(stack)
    data_min,data_max = transform[:,0].min(),transform[:,0].max()
    
    transform[:,0] = _rescale(P.ravel(),data_max,data_min,stretch=True)

    Xhat = np.dot(transform[:,:n_components], pca_.components_[:n_components,:])
    Xhat += mu[:,np.newaxis]
    
    panImg = Xhat.reshape([P.shape[0],P.shape[1],Xhat.shape[1]])
    
    return _rescale(panImg,panImg.max(),panImg.min())

    
if __name__ == "__main__":
    infile = '/home/ubuntu/hydra/data/VNP09GA_NRT.A2018053.h27v07.001.h5'
    
    viirs = io.readViirs(infile)
    
    try:
        pan = hsv(viirs)
    except:
        raise RuntimeError('Opps...something went wrong with module')