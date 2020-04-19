import numpy as np
import random

#import the data. Note that the txt file has to be in the same folder
data = np.loadtxt('cdata.txt')

#function for calculating standard deviation by using eq (5) from lab 2
def std5(val):
    mean = (1./len(val))*np.sum(val) # 1st pass to find mean
    sumOfDifference = 0
    for i in val: # 2nd pass to find difference
        sumOfDifference += (i - mean) ** 2
    ans5 = np.sqrt((1./(len(val)-1)) * sumOfDifference)
    return ans5

#function for calculating standard deviation by using eq (6) from lab 2
def std6(val):
    sumForMean = 0 #this is for calculating the mean 
    sumOfSquare = 0 #this is for calculating the sum of sqaured element
    for i in val: #only pass through data set
        sumForMean += i
        sumOfSquare += (i ** 2)
        if sumOfSquare < sumForMean: #check to prevent sqrt negative number
            print('you are wrong, but the code will still run')
            sumOfSquare += sumForMean #will make it larger than sum for mean but std will be wrong
    ans6 = np.sqrt((1./(len(val)-1)) * (sumOfSquare - len(val)*((sumForMean/len(val)) ** 2)))
    return ans6

#using numpy.std method as our reference answer
correctAnswer = np.std(data, ddof=1)

#using (x-y)/y to calculate relative error from lab 2
print('The relative error of method 5: ', np.abs(std5(data) - correctAnswer)/correctAnswer)
print('The relative error of method 6: ', np.abs(std6(data) - correctAnswer)/correctAnswer)

#constants for generating set 1
mean1, sigma1, n1 = (0., 1., 2000)
#constants for generating set 2
mean2, sigma2, n2 = (1.e7, 1., 2000)

#set 1 data
data1= np.random.normal(mean1, sigma1, n1)
#set 2 data
data2= np.random.normal(mean2, sigma2, n2)

#using numpy.std method to calculate respective answers
correctAnswer1 = np.std(data1, ddof=1)
correctAnswer2 = np.std(data2, ddof=1)

#using same method as above to calculate relative error for set 1
print('The relative error of method 5 of set 1: ', np.abs(std5(data1) - correctAnswer1)/correctAnswer1)
print('The relative error of method 6 of set 1: ', np.abs(std6(data1) - correctAnswer1)/correctAnswer1)

#using same method as above to calculate relative error for set 2
print('The relative error of method 5 of set 2: ', np.abs(std5(data2) - correctAnswer2)/correctAnswer2)
print('The relative error of method 6 of set 2: ', np.abs(std6(data2) - correctAnswer2)/correctAnswer2)

def std6Improved(val):
    sumForMean = 0 #same setup as above
    sumOfSquare = 0
    adjust = random.choice(val) #randomly pick an element to shift data set close to 0 and 1
    for i in val:
        sumForMean += i - adjust #shifting each data point
        sumOfSquare += ((i - adjust) ** 2) #shifting each data point
        if sumOfSquare < 0:
            print('you are wrong, but the code will still run')
            np.abs(sumOfSquare)
    #no need to account the adjustment into std calculation because shifting data set does not effect value of std
    ans6 = np.sqrt((1./(len(val)-1))*(sumOfSquare - (sumForMean ** 2)/len(val)))
    return ans6

#using same method as above to calculate relative error for improved eq 6
print('The relative error of method 6 with improved method: ', np.abs(std6Improved(data) - correctAnswer)/correctAnswer)
