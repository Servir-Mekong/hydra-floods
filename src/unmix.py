from __future__ import print_function,division
import numpy as np
import pysptools.util as util
import pysptools.eea as eea
import pysptools.abundance_maps as amp

def get_endmembers(data, n_members):
    print('Endmembers extraction with NFINDR')
    nfindr = eea.NFINDR()
    U = nfindr.extract(data, n_members, maxit=10, normalize=True, ATGP_init=True)
    # return an array of endmembers
    return U


def abundance_map(data, U):
    print('Abundance maps generation with NNLS')
    nnls = amp.NNLS()
    amaps = nnls.map(data, U, normalize=True)
    # return an array of abundance maps
    return amaps


def linear(img,endmembers,scale=True):
    
    if type(endmembers) != np.ndarray:
        endmembers = np.array(endmembers)
        
    if img.shape[2] != endmembers.shape[1]:
        raise IndexError('The array dimensions of channels and endmembers do not match')
        
    inverse = np.linalg.pinv(endmembers).T
    
    unmixed = np.matmul(inverse,img)
    
    return


def non_linear():
    
    return