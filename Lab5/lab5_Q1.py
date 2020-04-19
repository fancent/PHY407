"""
This file takes the closing value for each business day from late 2014 until 
2019 of the S&P 500stock index. By performing fourier transform, we hope to
extract information that will help us to understand the composition of the
stock market prices.
"""
import numpy as np
import matplotlib.pyplot as plt
## PART A
#import data
data = np.genfromtxt("sp500c.csv", delimiter=',', skip_header=1)
day = [i[0] for i in data]
close = [i[1] for i in data]

#plotting raw data
fig, graph = plt.subplots()
graph.plot(day, close)
graph.set(xlabel='day from late 2014 (day)', ylabel='closing value ($)',
       title='Raw data day vs closing value')
graph.grid()
fig.savefig("Q1a.jpg")
plt.show()

## PART B
#performing fourier transform and extracting the correct fruqeuncy for rfft
close_fft = np.fft.rfft(close)
N = len(close_fft)
freqs = np.fft.fftfreq(2*N, 1.0)
freqs = freqs[:N]

#plotting the result
fig, graph = plt.subplots()
graph.plot(freqs[1:], close_fft[1:])
graph.set(xlabel='Frequency (1/day)', ylabel='amplitude',
       title='fft frequency vs amplitude')
graph.grid()
fig.savefig("Q1b.jpg")
plt.show()

## PART C
#Test for inverse
iclose_fft = np.fft.irfft(close_fft, len(close))
#Taking the first 10% of the fft data
first_10p = np.zeros(N, float).astype("complex128")
first_10p[0 : int(N*0.1)] = (close_fft[0 : int(N*0.1)])
close_first10 = np.fft.irfft(first_10p,len(close))

#plotting the inverse of 10% of fourier transform data
fig, graph = plt.subplots()
graph.plot(day, close, label="oringinal data")
graph.plot(day, close_first10, label="inversed fft data")
graph.set(xlabel='days from late 2014 (day)', ylabel='closing value ($)',
       title='Original data vs inverse of 10% transformed data')
graph.grid()
graph.legend()
fig.savefig("Q1c.jpg")
plt.show()

#Part D
#period = 130 days, so freq = 1/130 = 0.0077
#picking the frequency up till 1/130   
cleanedFreq = np.where(freqs<=1/130)[0]
sixMonths = np.zeros(N, float).astype("complex128")
sixMonths[cleanedFreq[0] : cleanedFreq[-1]] = close_fft[cleanedFreq[0] : cleanedFreq[-1]]
close_sixMonths = np.fft.irfft(sixMonths,len(close))

#plotting the inverse of frequency up till 1/130 
fig, graph = plt.subplots()
graph.plot(day, close, label="oringinal data")
graph.plot(day, close_sixMonths, label="inversed fft data")
graph.set(xlabel='days from late 2014 (day)', ylabel='closing value ($)',
       title='Original data vs inverse of 6 months transformed data')
graph.grid()
graph.legend()
fig.savefig("Q1d.jpg")
plt.show()

#Part E
#similar to above where we restrict our frequency to 1/5
cleanedFreq_1week = np.where(freqs<=1/5)[0]
oneweek = np.zeros(N, float).astype("complex128")
oneweek[cleanedFreq_1week[0] : cleanedFreq_1week[-1]] = close_fft[cleanedFreq_1week[0] : cleanedFreq_1week[-1]]
close_oneWeek = np.fft.irfft(oneweek,len(close))

#plotting the inverse of frequency up till 1/5
fig, graph = plt.subplots()
graph.plot(day, close, label="oringinal data")
graph.plot(day, close_oneWeek, label="inversed fft data")
graph.set(xlabel='days from late 2014 (day)', ylabel='closing value ($)',
       title='Original data vs inverse of 1 week transformed data')
graph.grid()
graph.legend()
fig.savefig("Q1e.jpg")
plt.show()
#plotting the histogram of filtered data
n_bins = 20
histogram, axs = plt.subplots()
axs.hist(close_oneWeek, bins=n_bins)
axs.set(xlabel='Closing value ($)', ylabel='Frequency',
       title='Histogram of the filtered data')
histogram.savefig("Q1e_histogram.png")
plt.show()