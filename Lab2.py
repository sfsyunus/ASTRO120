import numpy as np
import scipy.misc
import os
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib import rcParams
from scipy.signal import argrelextrema
from cycler import cycler

rcParams['axes.labelsize'] = 18
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15
rcParams['figure.titlesize'] = 15
rcParams['savefig.dpi'] = 600

pathCal = "Data/Calibration/"
pathMulti = "Data/Multi/"
pathPlots = "Data/Plots/" 


def loadData(path, file):
	return np.transpose(np.genfromtxt(path + file + ".txt", dtype='float', skip_header=17, skip_footer=1))


def getBias():
	#returns avg value 134.2559668
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + "Bias_multi/"
	dataMean = np.array([])
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float'))
		mean = np.mean(data[1])
		dataMean = np.append(dataMean, mean)
	avg = np.mean(dataMean)
	arr_avg = np.zeros(2048)
	arr_avg += avg
	arr_avg[0] = 0
	return arr_avg

def getBias2():
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + "Bias_multi/"
	dataMean = np.zeros(2048)
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float'))
		dataMean += data[1]
	Bias = dataMean/100
	return Bias

def getOffset():
	#returns avg value 132.29318359
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + "Darks_multi/"
	dataMean = np.array([])
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))
		mean = np.mean(data[1])
		dataMean = np.append(dataMean, mean)
	avg = np.mean(dataMean)
	arr_avg = np.zeros(2048)
	arr_avg += avg
	arr_avg[0] = 0
	return arr_avg

def getOffset2():
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + "Darks_multi/"
	dataMean = np.zeros(2048)
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float',  skip_header=17, skip_footer=1))
		dataMean += data[1]
	darks = dataMean/100
	return darks

def plotOffsets():
	darks = getOffset2()
	darks = darks[1:]
	bias = getBias2()
	bias = bias[1:]
	avgDarks = getOffset()
	avgDarks = avgDarks[1:]
	avgBias = getBias()
	avgBias = avgBias[1:]
	plt.plot(darks, color = 'b', label = 'darks')
	plt.plot(bias, color = 'g', label = 'bias')
	plt.plot(avgDarks, color = 'r', label = 'mean darks')
	plt.plot(avgBias, color = 'k', label = 'mean bias')
	plt.legend()
	plt.xlim(0,2048)
	plt.xlabel("Pixels")
	plt.ylabel("Intensity")
	plt.title("Comparing bias and dark current")


def getMeanNeon():
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + "Neon_multi/"
	dataMean = np.zeros(2048)
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))
		dataMean += data[1]
	NeonMean = dataMean/100
	Neon = NeonMean - getOffset2()
	return Neon


def getMeanFluor():
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + "Fluor_multi/"
	dataMean = np.zeros(2048)
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))
		dataMean += data[1]
	FluorMean = dataMean/100
	Fluor = FluorMean - getOffset2()
	return Fluor

def loadNew(path, file):
	raw_data = loadData(path,file)
	arr_avg = getOffset2()
	data = raw_data
	data[1] = raw_data[1] - arr_avg
	return data

def DrawSpectrum(path, file):
	data = loadNew(path, file)
	x = data[0]
	y = data[1]
	plt.plot(x,y, color = "red", label = "Manufacturer's calibration")
	#plt.xlim(0,2048)
	plt.ylim(-50,2700)
	plt.title("Solar Spectrum")
	plt.xlabel("Pixels")
	plt.ylabel("Intensity [ADU]")
	plt.legend()
	#plt.savefig(pathPlots + "")

'''
def getPeaks(path, file):
#finds Neon peaks
	data = loadNew(path, file)
	x = data[0]
	y = data[1]
	peaks = np.array([])
	peak_values = np.array([])
	peak_index = np.array([])
	for idx, element in enumerate(y):
		if element > 200:
			if element > 5 + y[idx-1] and element > 5 + y[idx+1]:
				peak_values = np.append(peak_values, element)	
				peak_index = np.append(peak_index, idx)
	return data, peak_index, peak_values
'''

#for fluorescent lights use getPeaks: order = 50 and getCentroid: radius = 30
#for neon lamp use getPeaks: order = Null and getCentroid: radius = 5


def getPeaks(path, file):
#finds fluorescent light peaks
	data = loadNew(path,file)
	y = data[1]
	peak_values = y[argrelextrema(y, np.greater, order = 50)[0]]
	peak_index = argrelextrema(y, np.greater, order = 50)[0]
	return data, peak_index, peak_values


def getCentroids(path, file):
	data, peak_index, peak_values = getPeaks(path, file)
	counts = data[1]
	centroids = np.array([])
	for i in peak_index:
		for idx, element in enumerate(counts):
			if i == idx:
				x = np.arange(idx - 5, idx + 5)
				y = counts[idx - 5:idx + 5]
				centroid = (np.sum(x*y))/(np.sum(y))
				centroids = np.append(centroids, centroid)
	return centroids


def centroidError(path, file):
	data, peak_index, peak_values = getPeaks(path, file)
	counts = data[1]
	centroidErrors = np.array([])
	for i in peak_index:
		for idx, element in enumerate(counts):
			if i == idx:
				x = np.arange(idx - 5, idx + 5)
				y = counts[idx - 5:idx + 5]
				var = np.var(x)
				sumI = np.sum(y)
				centroidError = var/sumI
				centroidErrors = np.append(centroidErrors, centroidError)
	totError = np.mean(centroidErrors)
	return centroidErrors, totError

'''
#Neon Centroid Poisition Errors	
0.00035454,  0.00328467,  0.00230041,  0.00627453,  0.0057188 ,
0.00211554,  0.00140891,  0.00096095,  0.00221666,  0.00258267,
0.00142721,  0.00152799,  0.00089949,  0.0007456 ,  0.00104103,
0.0020007 ,  0.00161239,  0.0008856 ,  0.001424  
Mean error = 0.0020411421452219852
#################################
#Hg Centroid Position Errors
0.0439742 ,  0.0857321 ,  0.00558742,  0.00076541,  0.0018887 ,
0.00036238,  0.00144338,  0.00027584,  0.0024551 ,  0.00855587,
0.01233682
Mean error =  0.014852474765414484
'''

def NeonPixelWavelength(path, file):
	neon_centroids = getCentroids(path, file)
	neon_centroids = np.delete(neon_centroids, (np.r_[1,3:7,15,18]))
	neon_pixel = np.array([1358.6, 1417.8, 1547.1, 1560.6, 1596.4, 1629.5, 1675.4, 1708.6, 1721.7, 1794.0, 1858.8, 1915.1])
	#known pixels from USB 2000 handout for book keeping
	neon_wave = np.array([585.2487, 594.48342, 614.30626, 616.35939, 621.72812, 626.6495, 633.44278, 638.29917, 640.2248, 650.65281, 659.89529, 667.82762])
	#well known Ne I wavelengths
	##########################################
	plt.scatter(neon_wave, neon_centroids)
	plt.title("Neon Lamp known wavelengths vs pixel centroids")
	plt.xlabel("Wavelength [nm]")
	plt.ylabel("Pixel Number")
	return neon_centroids, neon_centroids.size

def FluorPixelWavelength(path, file):
	fluor_centroids = getCentroids(path, file)
	fluor_centroids = np.delete(fluor_centroids, (np.r_[1,4,6:11]))
	fluor_pixel = np.array([101.7, 310.2, 479.1, 1115.2])
	#known pixels from USB 2000 handout for book keeping
	fluor_wave = np.array([365.0153, 404.6563, 435.8328, 546.0735])
	#well known Hg I wavelengths
	##########################################
	plt.scatter(fluor_wave, fluor_centroids)
	plt.title("Hg I known wavelengths vs pixel centroids")
	plt.xlabel("Wavelength [nm]")
	plt.ylabel("Pixel Number")
	return fluor_centroids

def plotCentroid(path, file):
	centroids = getCentroids(path, file)
	DrawSpectrum(path, file)
	for idx,element in enumerate(centroids):
		if idx == 0:
			plt.axvline(element, linestyle = '-.', color = "purple", label = "Centroids")
		else:
			plt.axvline(element, linestyle = '-.', color = "purple")
	plt.legend(loc = 1)
	#plt.savefig(pathPlots + file)



NeWavelengths = np.array([585.2487, 594.48342, 614.30626, 616.35939, 621.72812, 626.6495, 633.44278, 638.29917, 640.2248, 650.65281, 659.89529, 667.82762])
NeCentroids = np.array([1360.33424625, 1419.35164458, 1548.52808431, 1561.91607869, 1598.08561939, 1631.11361148, 1677.00470621, 1710.11949035, 1723.06204597, 1795.37574513, 1860.3151117, 1916.89366745])
HgWavelengths = np.array([365.0153, 404.6563, 435.8328, 546.0735])
HgCentroids = np.array([104.26733898, 313.93123276, 482.31453689, 1110.05392963])
centroids = np.array([104.26733898, 313.93123276, 482.31453689, 1110.05392963, 1360.33424625, 1419.35164458, 1548.52808431, 1561.91607869, 1598.08561939, 1631.11361148, 1677.00470621, 1710.11949035, 1723.06204597, 1795.37574513, 1860.3151117, 1916.89366745])
wavelengths = np.array([365.0153, 404.6563, 435.8328, 546.0735, 585.2487, 594.48342, 614.30626, 616.35939, 621.72812, 626.6495, 633.44278, 638.29917, 640.2248, 650.65281, 659.89529, 667.82762])

def fit(x,y):
	#returns the same coefficients
	covm = np.cov(x, y)
	m = covm[0,1]/covm[0,0]
	c = np.mean(y) - m * np.mean(x)
	N = x.size
	y1 = m*x + c
	plt.plot(x,y1)
	res = y - y1
	var = np.var(res)
	mErr = N*var/(N*np.sum(x**2) - np.sum(x)**2 )
	cErr = var*np.sum(x**2)/(N*np.sum(x**2) - np.sum(x)**2)
	'''
	5.9811567411272906,
 	-2113.5157532705007,
	0.0033236364872263357,
 	1151.0024193760516
 	'''
	return m, c, mErr, cErr


def quadFit(x,y):
	#x is centroids y is wavelengths
	# 5.98115674, -2113.51575327
	#1.66943579e-01,   3.53700691e+02
	#x = np.arange(0,2048)
	coeffsLin = np.polyfit(x,y,1)
	m = coeffsLin[0]
	c = coeffsLin[1]
	y1 = m*x + c
    ################################################
    # 3.39499984e-03,   2.45091086e+00,  -1.23712295e+03
    #-1.53324772e-05,   1.98209060e-01,   3.44190745e+02
	coeffsQuad = np.polyfit(x,y,2)
	a = coeffsQuad[0]
	b = coeffsQuad[1]
	c = coeffsQuad[2]
	y2 = a*x*x + b*x + c
	################################################
	'''
	plt.plot(x,y1, '-.', linewidth = "2", color = "red", label = "Linear Fit")
	plt.plot(x,y2, '--', linewidth = "2", color = "purple", label = "Quadratic Fit")
	plt.title("Linear and Quadratic Fit for Pixel to Wavelenghth Calibration")
	plt.xlabel("Wavelength [nm]")
	plt.ylabel("Pixel Number")
	plt.xlim(350,700)
	plt.ylim(0,2048)
	plt.scatter(NeWavelengths,NeCentroids, s = 30, color = "blue", label = "Ne I Centroids")
	plt.scatter(HgWavelengths,HgCentroids, s = 30, color = "green", label = "Hg I Centroids")
	plt.legend(loc = 2)
	'''
	return coeffsLin, coeffsQuad


def plotResiduals():
	coeffsLin, coeffsQuad = quadFit(wavelengths,centroids)
	#x is centroids value, y is wavelengths
	m = coeffsLin[0]
	c = coeffsLin[1]
	y1 = m*wavelengths + c
	y1Ne = NeCentroids - (m*NeWavelengths + c)
	y1Hg = HgCentroids - (m*HgWavelengths + c)
	res1 = centroids - y1
	sigma1 = np.std(res1)
	#################
	a = coeffsQuad[0]
	b = coeffsQuad[1]
	d = coeffsQuad[2]
	y2 = a*wavelengths**2 + b*wavelengths + d
	y2Ne = NeCentroids -  (a*NeWavelengths**2 + b*NeWavelengths + d)
	y2Hg = HgCentroids - (a*HgWavelengths**2 + b*HgWavelengths + d)
	res2 = centroids - y2
	sigma2 = np.std(res2)
	#################
	#manufacturer's coeff
	e = -1.4913E-5
	f = 0.19743726
	g = 344.311846
	y3 = e*wavelengths**2 + f*wavelengths + g
	res3 = centroids - y3
	sigma3 = np.std(res3)
	#################
	fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1, sharex=True, sharey=False, squeeze=True)
	ax1.scatter(NeWavelengths, NeCentroids, color = "blue", label = "Ne I")
	ax1.scatter(HgWavelengths, HgCentroids, color = "green", label = "Hg I")
	ax1.plot(wavelengths, y2, color = "red", label = "Quadratic Fit")
	ax1.set_ylabel("Pixel [pixels]")
	ax1.legend(loc = 4)
	##
	ax2.scatter(NeWavelengths, y1Ne, color = "blue", label = "Ne I")
	ax2.scatter(HgWavelengths, y1Hg, color = "green", label = "Hg I")
	ax2.set_ylabel("Pixel Error [pixels]")
	ax2.legend(loc = 4)
	##
	ax3.scatter(NeWavelengths, y2Ne, color = "blue", label = "Ne I")
	ax3.scatter(HgWavelengths, y2Hg, color = "green", label = "Hg I")
	ax3.set_ylabel("Pixel Error [pixels]")
	ax3.legend(loc = 4)
	##
	plt.xlim(350, 700)
	plt.xlabel("$\lambda$ [nm]")
	#sigma1, sigma2 = 3.5322767535816961, 0.2450458984350104, 0.26729315568560091
	return sigma1, sigma2, sigma3

def pixToWave(pixel):
	#coeffsLin, coeffsQuad = quadFit(centroids,wavelengths)
	coeffsLin, coeffsQuad = quadFit(wavelengths, centroids)
	a2 = coeffsQuad[0]
	a1 = coeffsQuad[1]
	a0 = coeffsQuad[2]
	#wave = a2*pixels**2 + a1*pixels + a0
	wave = (-a1 + np.sqrt(a1**2 - 4*a2*(a0-pixel)))/(2*a2)
	return wave

def plotVar():
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/Data/Various/Colors/"
	labels = np.array(['Incandescent','Red reflected', 'Blue reflected'])
	for i, filename in enumerate(os.listdir(newPath)):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))
		x = pixToWave(data[0])
		y = data[1]
		#plt.rc('axes', prop_cycle=(cycler('label', ['Incandescent','Red reflected', 'Blue reflected'])))
		plt.rc('axes', prop_cycle=(cycler('color', ['k','r','b'])))
		plt.plot(x, y, label = labels[i])
		plt.xlim(pixToWave(0), pixToWave(2048))
		plt.xlabel("$\lambda$ [nm]")
		plt.ylabel("Intensity")
	plt.legend(loc = 2)

def timeSeries(path, pixN):
	#Plot 1 looks at Ne, Pixel 1360 at 585nm
	#Plot 2 looks at Hg, Pixel 1110 at 542nm
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + path + '/'
	fileN = np.arange(0,100)
	pixelVals = np.array([])
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))
		y = data[1]
		pixelVal = y[pixN]
		pixelVals = np.append(pixelVals, pixelVal)
	plt.plot(fileN, pixelVals)
	plt.xlabel("Sample Number")
	plt.ylabel("Intensity [ADU]")
	plt.title("Pixel 1541 response for Hg I Spectrum")
	return np.std(pixelVals)

def powerSpec(path):
	nyqFreq = 0.5/0.1
	c = 3E8
	pixels = np.arange(0,2048)
	inten = np.zeros(2048)
	wave = np.zeros(2048)
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + path + '/'
	for filename in os.listdir(newPath):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))
		wave += data[0]
		inten += data[1]
	meaninten = inten/100
	power = meaninten**2/2048
	#freq = c*100/wave
	freq = 	np.linspace(0,5,2048)
	plt.plot(freq/nyqFreq, power, color = 'b')
	plt.xlabel(r"$\nu$/$\nu_N$")
	plt.ylabel("Average power [ADU$^2$/sample]")
	plt.title("Average Power Series of Fluorescent light Spectrum")
	return freq.size, power.size


def noise(path):
	newPath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab2/" + pathMulti + path + '/'
	ADU = np.zeros(2048)
	pixArr = np.zeros([100,2048])
	varADU = np.array([])
	stdADU = np.array([])
	for i,filename in enumerate(os.listdir(newPath)):
		data = np.transpose(np.genfromtxt(newPath + filename, dtype='float', skip_header=17, skip_footer=1))	
		ADU += data[1]
		for j, pixVal in enumerate(data[1]):
			pixArr[i,j] += pixVal
	for k in np.arange(0,2048):
		var = np.var(pixArr[:,k])
		varADU = np.append(varADU, var)
		std = np.std(pixArr[:,k])
		stdADU = np.append(stdADU, std)
	ADU = ADU/100
	ADU0 = getBias2()
	#############
	x = ADU - ADU0
	y = varADU
	#y2 = stdADU
	plt.scatter(x,y, s = 10, color = "r")
	coeffsLin, coeffsQuad = quadFit(x,y)
	m = coeffsLin[0]
	c = coeffsLin[1]
	y1 = m*x + c
	a2 = coeffsQuad[0]
	a1 = coeffsQuad[1]
	a0 = coeffsQuad[2]
	y2 = a2*x**2 + a1*x + a0
	plt.semilogy(x,y1, '-.', color = 'b', lw = 2, label = "Linear fit")
	#plt.plot(x,y2, ':', color = 'g', lw = 2, label = "Second order polynomial fit")
	#plt.xlim(0,4000)
	#plt.ylim(0,700)
	plt.xlabel("$ADU - ADU_0$ [ADU]")
	plt.ylabel("Variance [ADU$^2$]")
	plt.title("Variance vs Mean plot for Fluorescent light Spectrum")
	plt.legend(loc = 2)
	return x,y, coeffsLin, coeffsQuad

'''
For Fluorescent linear fit
m = 0.090328135535942281,
c = 6.3608292582303747,
mErr = 4.1281609207483945e-06,
cErr = 1.6020592348379281
For Neon linear fit
m = 1.2535067053969415,
c =  31.822063990942091,
mErr =  0.00016786999221332321,
cErr = 5.5057831523320031
'''