"""
This file takes the height data of Hawaii and perform gradient calcualtion
and plotting the result in a contour graph.
Note that this file takes the code given in lab 3 solution and we did
grouping and modification to gurantee correctness.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import os
import urllib.request
import zipfile

def getdata():
    """
    This function retrieve the island data from the given website,
    then we sort the raw height data into the a 2d array
    """
    file="N19W156.hgt"
    url='https://dds.cr.usgs.gov/srtm/version2_1/SRTM3/Islands/{}.zip'.format(file)
    dest="{}.zip".format(file)
    if not os.path.exists("{}.zip".format(file)):
        urllib.request.urlretrieve(url, dest)
        
    data=zipfile.ZipFile(dest,'r').read(file)
    size=1201
    hgt=np.frombuffer(data, dtype='>h').reshape(size,size)[::-1]
    w=np.where(hgt<=-32760)
    fhgt=np.array(hgt, dtype=float)
    window=5
    for loc in zip(*w):
        around=fhgt[loc[0]-window : loc[0]+window, loc[1]-window : loc[1]+window]
        avg=np.sum(around*(around>0))/np.sum(around>0)
        fhgt[loc]=avg
    return fhgt, file
    
def grad1D(data, dx):
    """
    this function compute the gradient along one direction
    """
    delta=np.zeros_like(data)
    delta[1:-1]=(data[2:]-data[:-2])/(2*dx)
    delta[0]=(data[1]-data[0])/dx
    delta[-1]=(data[-1]-data[-2])/dx
    return delta

def gradient(data, dx=420, dy=420):
    """
    this function compute the gradient in both x and y direction.
    """
    dwdx=grad1D(data.T, dx).T
    dwdy=grad1D(data, dy)
    return dwdx, dwdy

def find_grid(name, size=1201):
    """
    this function takes the 2d array and convert the x and y into 
    longitude and latitude for graphing correctness.
    """
    loc=dict()
    for l in["N","E","S","W"]:
        if l in name:
            loc[l]=name.index(l)
    if len(loc)>2:
        print("Oh dear, weird location")
        return None
    else:
        first=loc.get("N", loc.get("S",0))
        second=loc.get("W", loc.get("E",0))
        latstart, lonstart=(
                float(name[first+1: second]),
                float(name[second+1: name.index(".")]),
        )
        if "W" in loc:
            lonstart =- lonstart
        if "S" in loc:
            latstart =- latstart
            
        lats=np.arange(latstart, latstart+1,1/(size))
        lons=np.arange(lonstart, lonstart+1,1/(size))
        return lons, lats
    
def make_figure(ax, xval, yval, data, sampleX, sampleY, title,cmap,vmin=None,vmax=None):
    """
    this function produces the contour graphs given the gradients.
    """
    d=data[sampleY, sampleX]
    vmin=vmin or d.min()
    vmax=vmax or d.max()
    ax.pcolormesh(xval[sampleX], yval[sampleY], d, cmap=cmap, vmin=vmin,vmax=vmax)
    ax.set_title(title)
    ax.set_ylabel("Latitude")
    ax.set_xlabel("Longitude")

