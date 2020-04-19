"""
This file uses fourier transform from numpy module on data filemsl_rems.csv
which recorded surface air pressure measurements from the Mars over 2
Mars years. We will explore the how to extract information through
Fourier transfrom
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

#importing data from csv
marsData = np.genfromtxt('msl_rems.csv', delimiter=',',skip_header=1)
Time,L_S,A,B,C,D,E = marsData.T


def cosine_wave(time, amplitude, period, phase):
    """Generate a cosine wave with known amplitude, period, and phase"""
    return amplitude*np.cos((time/period + phase) *2*np.pi)

def fourierAnalysis(data,dt):
    """
    Perform Fourier Transform and return the fft coeffiecients,
    frequency, amplitude, and phase
    """
    fftData = np.fft.fft(data)
    fftData = fftData[1::]
    freqs = np.fft.fftfreq(len(data),dt)
    freqs = freqs[1::]
    amplitude = 2 * np.abs(fftData)/len(data)
    phase = np.arctan2(fftData.imag, fftData.real)
    
    ln = len(data)//2
    amplitude = amplitude[:ln]
    phase = phase[:ln]
    freqs = freqs[:ln]
    
    return fftData, freqs, amplitude, phase

def printAnswers(peaksIndices, amplitude, freqs, phase):
    """
    This function takes the indices of peaks and prints out the 
    amplitude, period, and phase of the peaks.
    """
    for i in peaksIndices:
        print("Peak {}:".format(np.where(peaksIndices == i)[0]+1))
        print("amplitudes:", amplitude[i])
        print("period:", 1/freqs[i], "sol")
        print("phase:", phase[i], "rad")
        print()
    
def filterbyMinPeriod(peakIndices, freqs, amplitude, minimum, num):
    """
    This function takes the indices of peaks and filter out for peaks with 
    periods over a minimum threshold. Then, it will print out the 
    amplitude, period, and phase of the peaks.
    """
    temp = []
    for i in peakIndices:
        if 1/freqs[i] >= minimum:
            temp.append(i)
    result = sorted(temp, key=lambda x: amplitude[x], reverse=True)
    return result[:num]

def filterbyMaxPeriod(peakIndices, freqs, amplitude, maximum, num):
    """
    This function takes the indices of peaks and filter out for peaks with 
    periods under a maximum threshold. Then, it will print out the 
    amplitude, period, and phase of the peaks.
    """
    temp = []
    for i in peakIndices:
        if 1/freqs[i] <= maximum:
            temp.append(i)
    result = sorted(temp, key=lambda x: amplitude[x], reverse=True)
    return result[:num]


## PART A
#set arbituary time from 0s to 100s for testing
dt = 0.1
timeRange = np.arange(0,100,dt)
#construction of our test case with amplitude of 3, period of 20, and phase of 0.25
y = cosine_wave(timeRange, 3, 20, 0.25)
#fourier transformed data
fy, freqs, amplitude, phase = fourierAnalysis(y, dt)
#calculating amplitude, period, phase from fft data
print("fft amplitude:", np.max(amplitude), "original amplitude:",3)
print("fft period:", 1/freqs[np.argmax(amplitude)], "original period:",20)
print("fft phase:", phase[np.argmax(amplitude)]/(2*np.pi), "original phase:",np.pi/2)
#plotting the results
fig, graph = plt.subplots()
graph.plot(timeRange, y)
graph.set(xlabel='Time', ylabel='Amplitude',
      title='Frequency and amplitude')
graph.grid()
plt.show()

fig, graph = plt.subplots()
graph.plot(1/freqs, amplitude)
graph.set(xlabel='Frequency', ylabel='Amplitude',
      title='Frequency and amplitude')
graph.grid()
fig.savefig("q2a.png")
plt.show()


## PART B
#fourier transform data set A
fy, freqs, amplitude, phase = fourierAnalysis(A, 1/24)
#find the peaks in amplitude and printing respective informations:
peaksIndices = find_peaks(amplitude)[0][::-1]
printAnswers(peaksIndices, amplitude, freqs, phase)

#plotting amplitude and phase over period in log scale
fig, graph = plt.subplots()
graph.plot(1/freqs, amplitude)
graph.set(xlabel='Period in log scale (sols)', ylabel='Pressure Amplitude (Pa?)',
       title='Period and amplitude',xscale="log")
graph.grid()
fig.savefig("q2b_amplitude.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(1/freqs, phase)
graph.set(xlabel='Period in log scale (sols)', ylabel='Phase (rad)',
       title='Period and phase',xscale="log")
graph.grid()
fig.savefig("q2b_phase.png")
plt.show()


## PART C
#fourier transform data set B
fy, freqs, amplitude, phase = fourierAnalysis(B, 1/24)
#find the peaks in amplitude and printing respective informations:
peaksIndices = find_peaks(amplitude)[0][::-1]
printAnswers(peaksIndices, amplitude, freqs, phase)

#plotting amplitude and phase over period in log scale
fig, graph = plt.subplots()
graph.plot(1/freqs, amplitude)
graph.set(xlabel='Period in log scale (sols)', ylabel='Pressure Amplitude',
       title='Period and amplitude', xscale="log")
graph.grid()
fig.savefig("q2c_amplitude.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(1/freqs, phase)
graph.set(xlabel='Period in log scale (sols)', ylabel='Phase (rad)',
       title='Period and phase',xscale="log")
graph.grid()
fig.savefig("q2c_phase.png")
plt.show()


## PART D
#fourier transform data set C
fy, freqs, amplitude_C, phase = fourierAnalysis(C, 1/24)
#find the peaks in amplitude and printing respective informations:
peaksIndices_C = find_peaks(amplitude_C)[0][::-1]
printAnswers(peaksIndices_C, amplitude_C, freqs, phase)

#plotting amplitude and phase over period in log scale
fig, graph = plt.subplots()
graph.plot(1/freqs, amplitude_C)
graph.set(xlabel='Period in log scale (sols)', ylabel='Pressure Amplitude',
       title='Period and amplitude', xscale="log")
graph.grid()
fig.savefig("q2d_amplitude.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(1/freqs, phase)
graph.set(xlabel='Period in log scale (sols)', ylabel='Phase (rad)',
       title='Period and phase',xscale="log")
graph.grid()
fig.savefig("q2d_phase.png")
plt.show()


## PART E
#fourier transform data set D
fy, freqs, amplitude_D, phase = fourierAnalysis(D, 1/24)
#find the peaks in amplitude and printing respective informations:
#using information from part d to kinda cheat
mask = []
for peak in peaksIndices_C:
    mask.append(np.where(np.logical_and(amplitude_D <= amplitude_C[peak]*1.1,
                               amplitude_D >= amplitude_C[peak]*0.9))[0])
mask = [item for sublist in mask for item in sublist]
peaksIndices = find_peaks(amplitude_D, distance=16047)[0][::-1]
printAnswers(mask, amplitude_D, freqs, phase)

#plotting amplitude and phase over period in log scale
fig, graph = plt.subplots()
graph.plot(1/freqs, amplitude_D)
graph.set(xlabel='Period in log scale (sols)', ylabel='Pressure Amplitude (Pa?)',
       title='Period and amplitude', xscale="log")
graph.grid()
fig.savefig("q2e_amplitude.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(1/freqs, phase)
graph.set(xlabel='Period in log scale (sols)', ylabel='Phase (rad)',
       title='Period and phase',xscale="log")
graph.grid()
fig.savefig("q2e_phase.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(Time, D)
graph.set(xlabel='time (s)', ylabel='Pressure (Pa)',
       title='Time series of D')
graph.grid()
fig.savefig("q2e_timeseries.png")
plt.show()

## PART F
#fourier transform data set E
fy, freqs, amplitude, phase = fourierAnalysis(E, 1/24)
#find the peaks in amplitude and printing respective informations:
peaksIndices = find_peaks(amplitude)[0][::-1]
#finding peaks over 100 sols and under 1 sol
filtered100 = filterbyMinPeriod(peaksIndices, freqs, amplitude, 100, 4)
filtered1 = filterbyMaxPeriod(peaksIndices, freqs, amplitude, 1.1, 4)
print("largest  waves  with periods longer than 100 sols:")
printAnswers(filtered100, amplitude, freqs, phase)
print("largest  waves  with periods shorted than 1 sols:")
printAnswers(filtered1, amplitude, freqs, phase)

#plotting amplitude and phase over period in log scale
fig, graph = plt.subplots()
graph.plot(1/freqs, amplitude)
graph.set(xlabel='Period in log scale (sols)', ylabel='Pressure Amplitude (Pa?)',
       title='Period and amplitude', xscale="log")
graph.grid()
fig.savefig("q2f_amplitude.png")
plt.show()

fig, graph = plt.subplots()
graph.plot(1/freqs, phase)
graph.set(xlabel='Period in log scale (sols)', ylabel='Phase (rad)',
       title='Period and phase',xscale="log")
graph.grid()
fig.savefig("q2f_phase.png")
plt.show()
