import numpy as np

def cc(image1, image2, nodata):
    try:
    ccoeff = np.corrcoef(image1, image2)
    return ccoeff[0][1]
    except:
    return nodata

def rmse(tru, est):
    err = np.sqrt((est - tru)**2) 
    return err.mean()


def ssim():
    
    return


def uiqi():
    
    return


def ergas():
    
    return
