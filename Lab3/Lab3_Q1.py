import struct
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve

#Read the binary file
rawData = open('N19W156.hgt', 'rb')

def binaryReader(rawValues):
    #initialize list to store results into 2D array
    results = []
    for i in range(1201):
        #initialize list to store a row
        temp = []
        for j in range(1201):
            buf = rawValues.read(2)
            value = struct.unpack('>h', buf)[0]
            if value < -50:
                 #Invalid points
                 value = -50
            temp.append(value)
        #store every row to form a 2D array
        results.append(temp)
    return np.asarray(results)

def adjustArrayValue(array):
    #filtering 3x3 kernal to find avg from neighbours
    kernel = 1/8 * np.asarray([[1,1,1],[1,0,1],[1,1,1]])
    #filter the whole array with the kernal
    c = convolve(array, kernel, mode='constant')
    for y in range(len(array)):
        for x in range(len(array[0])):
            if array[y][x] == -50:
                #replace the invalid data points with avg of neighbours
                array[y][x] = c[y][x]
    return array
                
def findXGradients(data, delta):
    #initialize list to store results into 2D array
    results = []
    for row in range(len(data)):
        #initialize list to store a row
        temp = []
        #Forward Difference for first element
        first = (data[row][1]-data[row][0])/delta
        if first > 0.035 or first < -0.035:
                first = 0.035 if first>0.035 else -0.035
        temp.append(first)
        for i in range(1,len(data[0])-1):
            #central Difference for every element other than first and last in same row
            dx = (data[row][i+1]-data[row][i-1])/(2.0*delta)
            #limit the max and min to achieve higher contrast
            if dx > 0.035 or dx< -0.035:
                dx = 0.035 if dx>0.035 else -0.035
            temp.append(dx)
        #Backward Difference for last element
        last = (data[row][-1]-data[row][-2])/delta
        if last > 0.035 or last< -0.035:
                last = 0.035 if last>0.035 else -0.035
        temp.append(last)
        #store every row to form a 2D array
        results.append(temp)
    return np.asarray(results)

def findYGradients(data, delta):
    #turn columns into rows
    transposedData = np.transpose(data)
    #turn rows back to columns
    return np.transpose(findXGradients(transposedData,delta))

def calculateI(x, y, theta):
    return -((np.cos(theta)*x)+(np.sin(theta)*y))/(np.sqrt(x**2 + y**2)+1)

#set up constants and calculations
h = 420 #increment (m)
wArray = binaryReader(rawData)
wArray = adjustArrayValue(wArray)
dfdx = findXGradients(wArray,h)
dfdy = findYGradients(wArray,h)
I = calculateI(dfdx,dfdy, np.pi)

#plotting
fig1 = plt.figure(figsize=(6, 3.2))
ax = fig1.add_subplot(111)
ax.set_title('contourmap of height(w)\nfrom sealevel (m)')
ax.set_xlabel("x coordinate (420m)")
ax.set_ylabel("y coordinate (420m)")
plt.imshow(wArray)
ax.set_aspect('equal')
cax = fig1.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
fig1.savefig("q1_w.png")

fig2 = plt.figure(figsize=(6, 3.2))
ax = fig2.add_subplot(111)
ax.set_title('contourmap of intensity(I)\nof illumination (lx)')
ax.set_xlabel("x coordinate (420m)")
ax.set_ylabel("y coordinate (420m)")
plt.imshow(I)
ax.set_aspect('equal')
cax = fig2.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
fig2.savefig("q1_I.png")
plt.show()
