"""
This file takes an import from lab3 which calculated gradients in both x and y
direction and create a plot according to the value of gradient.
Then we will performa fourier transfrom on the raw height data and try to 
comput the gradient via fourier transfrom and compare the results from both methods.
"""
import matplotlib.pyplot as plt
import numpy as np
from time import time
from lab03_mapping import getdata, gradient, make_figure, find_grid

fhgt, file = getdata()

def lab3gradients(fhgt, file): 
    """
    calling functions from lab3 which does the following:
        1. compute gradients in x and y direction
        2. translating the position into longitude and latitude
        3. plot both gradients
    """
    dwdx, dwdy = gradient(fhgt)
    xloc, yloc = find_grid(file)
    fig, axs=plt.subplots(2,1, figsize=(5,12))
    axs=axs.ravel()
    sample=slice(None,None,1)
    plots=dict(dwdx=[dwdx,'gray',-.03,.03, axs[0]],
               dwdy=[dwdy,'gray',-.03,.03,axs[1]],)
    for k in["dwdx","dwdy"]:
        d=dict(zip(["data","cmap","vmin","vmax","ax"],plots[k]))
        make_figure(d["ax"],xloc, yloc, d["data"],sample,sample,k,d["cmap"],d["vmin"],d["vmax"])
    
def fft_dwdx(data):
    """
    This function performs fourier transform on each row and multiply the results
    with the derivative coeffient which result in gradients
    """
    result = []
    for i in range(len(data)):
        ft_data = np.fft.fft(data[i])
        k_coeff = np.fft.fftfreq(len(data[i]))
        diff_coeff = 1j * 2 * np.pi * k_coeff
        ft_diff = ft_data * diff_coeff
        grad = np.fft.ifft(ft_diff)
        result.append(grad.real)
    return np.asarray(result)

def fft_dwdy(data):
    """
    This function uses the fft_dwdx to perform calcualtion along the coloumns
    by transposing the data and transposing back the gradients
    """
    result = fft_dwdx(np.transpose(data))
    return result.T

def partCgradients(fhgt, file):
    """
    calling functions from lab3 and fft_dwdx, fft_dwdy which does the following:
        1. compute gradients in x and y direction with fft_dwdx and fft_dwdy
        2. translating the position into longitude and latitude
        3. plot both gradients
    """
    dwdx, dwdy = fft_dwdx(fhgt), fft_dwdy(fhgt)
    xloc, yloc = find_grid(file)
    fig, axs=plt.subplots(2,1, figsize=(5,12))
    axs=axs.ravel()
    sample=slice(None,None,1)
    plots=dict(dwdx=[dwdx,'gray',-10,10, axs[0]],
               dwdy=[dwdy,'gray',-10,10, axs[1]],)
    for k in["dwdx","dwdy"]:
        d=dict(zip(["data","cmap","vmin","vmax","ax"],plots[k]))
        make_figure(d["ax"],xloc, yloc, d["data"],sample,sample,k,d["cmap"],d["vmin"],d["vmax"])

#call to get the answers for part c
lab3gradients(fhgt, file)        
partCgradients(fhgt, file)

def padded_fft_dwdx(data, padSize):
    """
    This function first pad the data with given size, then it performs fourier 
    transform on each row and multiply the results with the derivative 
    coeffient which result in gradients
    """
    padded = np.pad(data, (padSize, padSize), 'constant')
    paddedResult = fft_dwdx(padded)
    return paddedResult[padSize:padSize+len(data),padSize:padSize+len(data)]

def padded_fft_dwdy(data, padSize):
    """
    This function uses the padded_fft_dwdx to perform calcualtion along the 
    coloumns by transposing the data and transposing back the gradients
    """
    padded = np.pad(data, (padSize, padSize), 'constant')
    paddedResult = fft_dwdx(np.transpose(padded))
    return paddedResult.T[padSize:padSize+len(data),padSize:padSize+len(data)]

def partD(fhgt):
    """
    calling functions padded_fft_dwdx, padded_fft_dwdy which does the following:
        1. compute gradients in x and y direction for every pad size
        2. compute the time taken
        3. plot the time vs pad size
    """
    padSizeRange = np.arange(1,100,1)
    timeResults = []
    totalSize = []
    N = len(fhgt)
    for i in padSizeRange:
        start = time()
        padded_fft_dwdx(fhgt, i)
        padded_fft_dwdy(fhgt, i)
        end=time()
        diff=end-start
        timeResults.append(diff)
        totalSize.append((N+2*i))
    fig, Graph = plt.subplots()
    Graph.plot(totalSize, timeResults)
    Graph.set(xlabel='Padding Size', ylabel='Time',
           title='Frequency and amplitude')
    Graph.grid()
    fig.savefig("q3c.png")
    plt.show()
    
#call to get the answer for part d
partD(fhgt)
    