import numpy as np
from matplotlib import pyplot as plt
import astropy.io.fits as pf
import os
import scipy.misc
import scipy.stats as stats
import matplotlib as mpl
from matplotlib import rcParams
from scipy.signal import argrelextrema
from cycler import cycler
import scipy.ndimage as ndi
import urllib as url
import string as str
from scipy import signal
import pdb
from numpy.linalg import inv, det, norm
import sys

rcParams['axes.labelsize'] = 18
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15
rcParams['figure.titlesize'] = 20
rcParams['savefig.dpi'] = 600


#Echelle orders 
neonEch = {'1':[40,90], '2':[120,170], '3':[190,230], '4':[260,300], '5':[320,370], '6':[380,420], '7':[440,480], '8':[490,530], '9':[540,580]}
ledEch = {'1':[310,370]}
halEch = {'1':[40, 80], '2':[110, 160], '3':[180, 230], '4':[250, 300], '5':[310, 360],
        '6':[380, 420], '7':[430, 480], '8':[485, 530], '9':[535, 580], '10':[580, 630], 
        '11':[630, 670], '12':[675, 715], '13':[718, 755], '14':[757, 793], '15':[795, 830]}
sunEch = {'1':[52, 96], '2':[126, 170], '3':[196, 240], '4':[261, 304], '5':[322, 365],
        '6':[380, 420], '7':[435, 480], '8':[485, 530], '9':[535, 575], '10':[580, 621], 
        '11':[625, 665], '12':[666, 707]}




def getDark1():
	#average over 1s exposure dark files
	path = 'data/dark_1s/'

	files = os.listdir(path)
	#create an empty array of image size
	file = pf.open(path + files[0])
        data = file[0].data
        x_size = data.shape[1]
        y_size = data.shape[0]
        dark = np.zeros((y_size, x_size))
	#loop over all darks and get sum of each pixel
	for i in range(len(files)):
            file = pf.open(path + files[i])
            data = file[0].data
            data = data[::-1,::-1]
            dark += data
	#get average over all pixels
	dark /= len(files)
	'''
	im = plt.imshow(dark, cmap = "gray_r", aspect = "equal", vmin=np.median(dark)-10, vmax=np.median(dark)+10, interpolation = "bilinear",  origin = "lower")
	clb =plt.colorbar(im)
	clb.ax.set_ylabel('Counts [ADU]')
	plt.xlabel("x [pixels]")
	plt.ylabel("y [pixels]")
	plt.title("Dark Current")
	'''
	return dark, files

def getDark10():
	#average over 1s exposure dark files
	path = 'data/dark_10s/'

	files = os.listdir(path)
	#create an empty array of image size
	file = pf.open(path + files[0])
        data = file[0].data
        x_size = data.shape[1]
        y_size = data.shape[0]
        dark = np.zeros((y_size, x_size))
	#loop over all darks and get sum of each pixel
	for i in range(len(files)):
            file = pf.open(path + files[i])
            data = file[0].data
            data = data[::-1,::-1]
            dark += data
	#get average over all pixels
	dark /= len(files)
	'''
	im = plt.imshow(dark, cmap = "gray_r", aspect = "equal",vmin=np.median(dark)-20, vmax=np.median(dark)+20, interpolation = "bilinear",  origin = "lower")
	clb =plt.colorbar(im)
	clb.ax.set_ylabel('Counts [ADU]')
	plt.xlabel("x [pixels]")
	plt.ylabel("y [pixels]")
	plt.title("Dark Current")
	'''
	return dark

def getNeon():
    path = 'data/neon_1s/'
    darks = getDark1()

    files = os.listdir(path)
   
    file = pf.open(path + files[0])
    temp = file[0].data
    y_size = temp.shape[0]
    x_size = temp.shape[1]
    arr = np.zeros((y_size, x_size))
    
    for i in range(len(files)):
        file = pf.open(path + files[i])
        data = file[0].data
        data = data[::-1, ::-1]
        arr += data
	
    avgNeon = arr/len(files)
    avgNeon = avgNeon - darks
    
    im = plt.imshow(avgNeon, cmap = "gray_r", aspect = "equal", vmin = 7*np.mean(avgNeon), vmax=15*np.mean(avgNeon), origin = 'lower')
    clb =plt.colorbar(im)
    clb.ax.set_ylabel('Counts [ADU]')
    plt.xlabel("x [pixels]")
    plt.ylabel("y [pixels]")
    plt.title("Neon lamp")
	
    return avgNeon


def getHal():
    path = 'data/halogen_10s/'
    darks = getDark10()

    files = os.listdir(path)
   
    file = pf.open(path + files[0])
    temp = file[0].data
    y_size = temp.shape[0]
    x_size = temp.shape[1]
    arr = np.zeros((y_size, x_size))
    
    for i in range(len(files)):
        file = pf.open(path + files[i])
        data = file[0].data
        data = data[::-1, ::-1]
        arr += data
	
    avgHal = arr/len(files)
    avgHal = avgHal - darks
    '''
    im = plt.imshow(avgHal, cmap = "gray_r", aspect = "equal", vmin = 0, vmax = np.mean(avgHal), origin='lower')
    clb =plt.colorbar(im)
    clb.ax.set_ylabel('Counts [ADU]')
    plt.xlabel("x [pixels]")
    plt.ylabel("y [pixels]")
    plt.title("Halogen lamp")
	'''
    return avgHal

def getLED():
    path = 'data/635nm_10s/'
    darks = getDark10()

    files = os.listdir(path)
   
    file = pf.open(path + files[0])
    temp = file[0].data
    y_size = temp.shape[0]
    x_size = temp.shape[1]
    arr = np.zeros((y_size, x_size))
    
    for i in range(len(files)):
        file = pf.open(path + files[i])
        data = file[0].data
        data = data[::-1, ::-1]
        arr += data
	
    avgLED = arr/len(files)
    avgLED = avgLED - darks
    '''
    im = plt.imshow(avgLED, cmap = "gray_r", aspect = "equal", vmin = np.mean(avgLED)/50, vmax = 50*np.mean(avgLED), origin='lower')
    clb =plt.colorbar(im)
    clb.ax.set_ylabel('Counts [ADU]')
    plt.xlabel("x [pixels]")
    plt.ylabel("y [pixels]")
    plt.title("653nm LED")
	'''
    return avgLED



def plot2():
    path = 'data/neon_1s/'
    darks = getDark1()

    files = os.listdir(path)
   
    file = pf.open(path + files[0])
    temp = file[0].data
    y_size = temp.shape[0]
    x_size = temp.shape[1]
    arr = np.zeros((y_size, x_size))
    
    for i in range(len(files)):
        file = pf.open(path + files[i])
        data = file[0].data
        data = data[::-1]
        arr += data
	
    ccdim = arr/len(files)
    ccdim = ccdim - darks
    

    dimensions = ccdim.shape
    midpoint_x = int(np.floor(dimensions[1]/2.0))
    
    #pick a slice from the center of the CCD image
    vertical_slice = [ccdim[y,midpoint_x] for y in np.arange(dimensions[0])]

    plot = plt.figure()
    plt.plot(np.arange(dimensions[0]),vertical_slice,'g',color="red")
    plt.xlabel('Pixel Number', size=15)
    plt.xlim([0,dimensions[0]])
    plt.ylim([0,np.max(vertical_slice)*1.1])
    plt.ylabel('Count [ADU]', size=15)
    
    return ccdim, dimensions, midpoint_x


def getSpectrum(img, echelle, orderstr):
#make ccd image 1d
    yRange = echelle[orderstr] # sets y range of echelle order
    imgLim = img[yRange[0]:yRange[1],:] # selects array values for specified y range
    xRange = np.arange(len(imgLim[0,:]))
    intensity = np.array([])
    for i in xRange:
        intensity = np.append(intensity, np.mean(imgLim[:,i])) # getting average of intensity values
   	
    return intensity, xRange #returns the intensity data and x values (pixels)

def plotSpectrum(img, echelle, orderstr):
#make ccd image 1d
    yRange = echelle[orderstr] # sets y range of echelle order
    imgLim = img[yRange[0]:yRange[1],:] # selects array values for specified y range
    xRange = np.arange(len(imgLim[0,:]))
    intensity = np.array([])
    for i in xRange:
        intensity = np.append(intensity, np.mean(imgLim[:,i])) # getting average of intensity values
   	
    plt.plot(xRange, intensity, label = "Order %s" %orderstr)
    plt.legend(loc = "best", ncol = 2)
    plt.xlim(0,1048)
    plt.ylim(-100, 12500)
    plt.xlabel('x Pixel Number')
    plt.ylabel('Counts [ADU]')
    
    return intensity, xRange #returns the intensity data and x values (pixels)


def plotNeonSpectrum():
	orders = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
	flux = np.array([])
	xRange = np.array([])
	for order in orders:
		intensity, x = plotSpectrum(getNeon(), neonEch, order)
		flux = np.append(flux, intensity)
		xRange = np.append(xRange, x)
	return flux, xRange

def plotLEDSpectrum():
	orders = ['1']
	flux = np.array([])
	xRange = np.array([])
	for order in orders:
		intensity, x = plotSpectrum(getLED(), ledEch, order)
		flux = np.append(flux, intensity)
		xRange = np.append(xRange, x)
	return flux, xRange


def getPeaks(img, echelle, orderstr):
#finds Neon peaks
	flux, x = getSpectrum(img, echelle, orderstr)
	peaks = np.array([])
	peak_values = np.array([])
	peak_pix = np.array([])
	for idx, element in enumerate(flux):
		if element > 100 and idx < 1040:
			if element > 5 + flux[idx-1] and element > 5 + flux[idx+1]:
				if idx not in [15, 477, 152]:
					peak_values = np.append(peak_values, element)	
					peak_pix = np.append(peak_pix, idx)
	return peak_pix, peak_values



def getCentroids(img, echelle, orderstr):
	#find centroids given peaks from the spectrum

	flux, x = getSpectrum(img, echelle, orderstr)
	peak_pix, peak_values = getPeaks(img, echelle, orderstr)

	centroids = np.array([])
	centroidErrors = np.array([])

	for i in peak_pix:
		for idx, element in enumerate(flux):
			if i == idx:
				x = np.arange(idx - 10, idx + 10)
				y = flux[idx - 10:idx + 10]
				centroid = (np.sum(x*y))/(np.sum(y))
				centroids = np.append(centroids, centroid)
				var = np.var(x)
				sumI = np.sum(y)
				centroidError = var/sumI
				centroidErrors = np.append(centroidErrors, centroidError)
			
	totError = np.mean(centroidErrors)

	return centroids, centroidErrors, totError


def NeonCent():
	orders = ['1', '2', '3', '4', '5', '6', '7', '8', '9']

	centroids = np.array([])
	centroidErrors = np.array([])

	for order in orders:
		x,y = plotSpectrum(getNeon(), neonEch, order)
		cent, centErr, totErr = getCentroids(getNeon(), neonEch, order)
		centroids = np.append(centroids, cent)
		centroidErrors = np.append(centroidErrors, centErr)
		
	totError = np.mean(centroidErrors)
	return centroids, centroidErrors, totError



def plotCent(img, echelle, orderstr):
	centroids, centroidErrors, totError = getCentroids(img, echelle, orderstr)
	I, x = plotSpectrum(img, echelle, orderstr)

	for idx,element in enumerate(centroids):
		if idx == 0:
			plt.axvline(element, linestyle = '-.', color = "purple", label = "Centroids")
		else:
			plt.axvline(element, linestyle = '-.', color = "purple")
	plt.legend(loc = 1)
	#plt.savefig(pathPlots + file)


#### all parameters for fit
neonWav = np.array([585.24878, 588.1895, 
					588.1895, 594.48342, 597.46276, 602.99969, 
					607.43377, 609.61631, 614.30626, 616.35939, 621.72812, 
					621.72812, 626.6495, 630.4789, 633.44278, 638.29917, 640.2246, 
					638.29917, 640.2246, 650.65281, 653.28822, 659.8953, 
					659.89529, 667.82762, 671.7043, 
					692.94673, 
					703.24131, 717.39381, 
					724.51666, 747.24386])

neonCent = np.array([714.3921528, 869.8691231, 
					83.64854374, 386.8675134, 538.6042546, 814.6641011,
					220.3124687, 323.0120201, 547.7210432, 648.0203816, 916.8037346,
					84.47574946, 308.0277121, 485.715391, 625.8109609, 861.0368293, 956.4791089, 
					13.14578582, 96.98266278, 564.15251918, 686.3188149, 1001.753767, 
					113.05757, 456.3820864, 629.1587041, 
					647.2545527, 
					157.2681132, 745.3153014,
					82.79411475, 865.2195724])

#### fit parameters by echelle order
neonWavs = {'9':[585.24878, 588.1895], '8':[588.1895, 594.48342, 597.46276, 602.99969], 
			'7':[607.43377, 609.61631, 614.30626, 616.35939, 621.72812], 
			'6':[621.72812, 626.6495, 630.4789, 633.44278, 638.29917, 640.2246], 
			'5':[638.29917, 640.2246, 650.65281, 653.28822, 659.8953], 
			'4':[659.89529, 667.82762, 671.7043], '3':[692.94673], 
			'2':[703.24131, 717.39381], '1':[724.51666, 747.24386]}

neonCents = {'9':[714.3921528, 869.8691231], '8':[83.64854374, 386.8675134, 538.6042546, 814.6641011], 
			'7':[220.3124687, 323.0120201, 547.7210432, 648.0203816, 916.8037346], 
			'6':[84.47574946, 308.0277121, 485.715391, 625.8109609, 861.0368293, 956.4791089], 
			'5':[13.14578582, 96.98266278, 564.15251918, 686.3188149, 1001.753767], 
			'4':[113.05757, 456.3820864, 629.1587041], '3':[647.2545527], 
			'2':[157.2681132, 745.3153014], '1':[82.79411475, 865.2195724]}



def fit(orderstr):
	
	x = np.array(neonWavs[orderstr])
	y = np.array(neonCents[orderstr])
	
	#x = wavelengths
	#y = centroids

	covm = np.cov(x, y)
	m = covm[0,1]/covm[0,0]
	c = np.mean(y) - m * np.mean(x)

	N = x.size
	y1 = m*x + c
	res1 = y - y1
	var = np.var(res1)
	mErr = N*var/(N*np.sum(x**2) - np.sum(x)**2 )
	cErr = var*np.sum(x**2)/(N*np.sum(x**2) - np.sum(x)**2)
	

	coeffsQuad, cov = np.polyfit(x,y,2, cov=True)
	a2 = coeffsQuad[0]
	a1 = coeffsQuad[1]
	a0 = coeffsQuad[2]
	y2 = a2*x*x + a1*x + a0
	res2 = y - y2
	quadErr = np.sqrt(np.diag(cov))
	'''
	x_extra = np.append(x,  x[-1:])
	y_extra = np.append(y, y[-1:])
	weights = [1.0, 1.0, 1.0, 1.0, 1.0, sys.float_info.epsilon]

	fit_extra, cov_extra = np.polyfit(x_extra, y_extra, 2, w=weights, cov=True)
	quadErr_extra = np.sqrt(np.diag(cov_extra))


	fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1, sharex=True, sharey=False, squeeze=True)
	
	ax1.scatter(x, y, color = "blue", label = "Measured centroids")
	ax1.plot(x, y1, linestyle='--', marker='x', linewidth = "2", color = "red", label = "Linear Fit")
	ax1.set_ylabel("Pixel [pixels]")
	ax1.legend(loc = 4)
	ax1.set_ylim(0,1024)
	##
	ax2.scatter(x, y, color = "blue", label = "Measured centroids")
	ax2.plot(x, y2, linestyle='--',marker='v', linewidth = "2", color = "green", label = "Quadratic Fit")
	ax2.set_ylabel("Pixel [pixels]")
	ax2.legend(loc = 4)
	ax2.set_ylim(0,1024)
	##
	plt.xlabel("$\lambda$ [nm]")
	plt.xlim(637,661)
	'''
	return  x, y, y1, y2, m, c, mErr, cErr, coeffsQuad, quadErr #, fit_extra, quadErr_extra


def plotResiduals(orderstr):
	x, y, y1, y2, m, c, mErr, cErr, coeffsQuad, quadErr = fit(orderstr)

	res1 = y - y1
	sigma1 = np.std(res1)
	#################
	res2 = y - y2
	sigma2 = np.std(res2)
	#################
	fig, (ax1, ax2, ax3) = plt.subplots(nrows = 3, ncols = 1, sharex=True, sharey=False, squeeze=True)
	
	ax1.scatter(x, y, color = "blue", label = "Ne I")
	ax1.plot(x, y1, color = "red", label = "Linear Fit")
	ax1.set_ylabel("Pixel [pixels]")
	ax1.legend(loc = 4)
	##
	ax2.scatter(x, res1, color = "red", label = "Linear Residuals")
	ax2.set_ylabel("Pixel Error [pixels]")
	ax2.legend(loc = 4)
	##
	ax3.scatter(x, res2, color = "green", label = "Quadratic Residuals")
	ax3.set_ylabel("Pixel Error [pixels]")
	ax3.legend(loc = 4)
	##
	#plt.xlim(350, 700)
	plt.xlabel("$\lambda$ [nm]")
	


def pixToWave(orderstr, pixels):
	x, y, y1, y2, m, c, mErr, cErr, coeffsQuad, quadErr = fit(orderstr)

	a2 = coeffsQuad[0]
	a1 = coeffsQuad[1]
	a0 = coeffsQuad[2]
	
	wave1 = (pixels-c)/m
	wave2 = (-a1 + np.sqrt(a1**2 - 4*a2*(a0-pixels)))/(2*a2)
	
	return wave1


def getSun(nthImg):
    path = 'data/solar/Nov-28-2017/transit1/'
    #Transit 1: start at 17 end at 53
    #Transit 2: start at 12 end at 29 but seems too short
    #Transit 3: start at 25 end at 60

    files = os.listdir(path)
   
    file = pf.open(path + files[nthImg])
    data = file[0].data
    data = data[::-1, ::-1]
   
    
    #create an empty array of image size
    temp = pf.open(path + files[1])
    datemp = temp[0].data
    x_size = datemp.shape[1]
    y_size = datemp.shape[0]
    darkarr = np.zeros((y_size, x_size))
    #loop over all darks and get sum of each pixel
    for i in range(1,16):
        file = pf.open(path + files[i])
        dark = file[0].data
        darkarr += dark
	#get average over all pixels
	darkarr /= len(files)

	sun_data = data - darkarr

    return file, sun_data

def plotSunCCD(nthImg):
    file, sun_data = getSun(nthImg)
    im = plt.imshow(sun_data, cmap = "gray", aspect = "equal", origin = 'lower')
    clb =plt.colorbar(im)
    clb.ax.set_ylabel('Counts [ADU]')
    plt.xlabel("x [pixels]")
    plt.ylabel("y [pixels]")
    plt.title("Solar Spectrum")


def plotSunSpec(nthImg, orderstr):
	file, sun_data = getSun(nthImg)
	sunFlux, sunWavs = getSpectrum(sun_data, sunEch, orderstr)

	sunWavs = pixToWave(orderstr, sunWavs)
	sunFlux /= np.max(sunFlux)
	
	sunFlux = sunFlux[8:]
	sunWavs = sunWavs[8:]
	
	plt.plot(sunWavs, sunFlux, label = "Sun Spectrum")
	plt.axvline(656.25, color='k', linestyle='-.', label = 'H-alpha')
	plt.legend(loc = "best")
	#plt.xlim(0,1048)
	#plt.ylim(-100, 12500)
	plt.xlabel('$\lambda$ [nm]')
	plt.ylabel('Relative Flux')
	
	return sunFlux, sunWavs


def plotLimbs(nLimb, nCent, orderstr):
	file1, sun_limb = getSun(nLimb)
	file2, sun_cent = getSun(nCent)
	
	sunFluxL, sunXL = getSpectrum(sun_limb, sunEch, orderstr)
	sunFluxC, sunXC = getSpectrum(sun_cent, sunEch, orderstr)

	sunWavsL = pixToWave(orderstr, sunXL)
	sunFluxL /= np.max(sunFluxL)
	sunWavsC = pixToWave(orderstr, sunXC)
	sunFluxC /= np.max(sunFluxC)

	
	sunFluxL = sunFluxL[8:]
	sunWavsL = sunWavsL[8:]
	sunFluxC = sunFluxC[8:]
	sunWavsC = sunWavsC[8:]

	plt.plot(sunWavsC, sunFluxC, color = 'red', label = "Sun Center")
	plt.plot(sunWavsL, sunFluxL, color = 'blue', label = "Sun Limb")
	plt.legend(loc = "best", ncol = 2)
	plt.xlabel('$\lambda$ [nm]')
	plt.ylabel('Relative Flux')
	
	return sunWavsL, sunFluxL, sunWavsC, sunFluxC



def plotIntensity():
    path = 'data/solar/Nov-28-2017/transit1/'
	
    files = os.listdir(path)
   	#get array size
    file = pf.open(path + files[1])
    temp = file[0].data
    y_size = temp.shape[0]
    x_size = temp.shape[1]
    flux = np.array([])
    time = np.array([])

    #get darks
    darkarr = np.zeros((y_size, x_size))
    #loop over all darks and get sum of each pixel
    for i in range(1,16):
        file = pf.open(path + files[i])
        dark = file[0].data
        darkarr += dark
	#get average over all pixels
	darkarr /= len(files)
    
    #loop over all files and get sum of flux and times
    for i in range(1, len(files)):
        file = pf.open(path + files[i])
        time = np.append(time, file[0].header['jd'])
        data = file[0].data
        data = data[::-1,::-1]
        fSun = np.sum(data - darkarr)
        flux = np.append(flux, fSun)
    flux /= flux.max()
    
    #get time in seconds
    time = (time - time[0])*86400

    #bounds - left:16, right: 52
    plt.axvline(time[16], color='k', linestyle='--')
    plt.axvline(time[52], color='k', linestyle='--')
    plt.scatter(time, flux)
    plt.xlim(-5,290)
    plt.ylim(-0.03,1.05)
    plt.xlabel("Time [s]")
    plt.ylabel("Relative Flux")
    # Duration of transit: 137.00002133846283 seconds
    return flux, time



def getSunPeaks(nthImg,orderstr):
#finds fluorescent light peaks
	file, data = getSun(nthImg)
	y, x = getSpectrum(data, sunEch, orderstr)
	peak_values = y[argrelextrema(y, np.less, order = 25)[0]]
	peak_index = argrelextrema(y, np.less, order = 25)[0]
	return peak_index, peak_values, x, y


def getSunCentroids(nthImg, orderstr):
	#find centroids given peaks from the spectrum

	peak_pix, peak_values, x_pos, flux = getSunPeaks(nthImg, orderstr)

	centroids = np.array([])
	centroidErrors = np.array([])

	for i in peak_pix:
		for idx, element in enumerate(flux):
			if i == idx and i < 1037 and i > 10:
				x = np.arange(idx - 10, idx + 10)
				y = flux[idx - 10:idx + 10]
				centroid = (np.sum(x*y))/(np.sum(y))
				centroids = np.append(centroids, centroid)
				var = np.var(x)
				sumI = np.sum(y)
				centroidError = var/sumI
				centroidErrors = np.append(centroidErrors, centroidError)
			
	totError = np.mean(centroidErrors)

	return centroids#, centroidErrors, totError


def plotSunCent(nthImg, orderstr):
	centroids = getSunCentroids(nthImg, orderstr)
	flux, wavs = plotSunSpec(nthImg, orderstr)
	centroids = pixToWave(orderstr, centroids)

	for idx,element in enumerate(centroids):
		if idx == 0:
			plt.axvline(element, linestyle = '-.', color = "purple", label = "Centroids")
		else:
			plt.axvline(element, linestyle = '-.', color = "purple")
	plt.legend(loc = 1)


def window(data):
	return np.hanning(len(data))*data


def quadFit(x, y):

	coeffsQuad, cov = np.polyfit(x,y,2, cov=True)
	a2 = coeffsQuad[0]
	a1 = coeffsQuad[1]
	a0 = coeffsQuad[2]
	y2 = a2*x*x + a1*x + a0
	res2 = y - y2
	quadErr = np.sqrt(np.diag(cov))
	
	return  x, y2, quadErr


def process(nthImg, orderstr):
	#get the spectrum
	file, sun_data = getSun(nthImg)
	sunFlux, sunX = getSpectrum(sun_data, sunEch, orderstr)

	sunWavs = pixToWave(orderstr, sunX)
	
	sunFlux = sunFlux[8:]
	sunWavs = sunWavs[8:]

	x, y, Err = quadFit(sunWavs, sunFlux)

	sunFlux = sunFlux - y
	sunFlux /= np.max(sunFlux)
	sunFlux = window(sunFlux)
	return sunFlux, sunWavs

def plotLimbsN(nLimb, nCent, orderstr):
	
	sunFluxL, sunWavsL = process(nLimb, orderstr)
	sunFluxC, sunWavsC = process(nCent, orderstr)

	plt.plot(sunWavsC, sunFluxC, color = 'red', label = "Sun Center")
	plt.plot(sunWavsL, sunFluxL, color = 'blue', label = "Sun Limb")
	plt.legend(loc = "best", ncol = 2)
	plt.xlim(620,642.5)
	plt.xlabel('$\lambda$ [nm]')
	plt.ylabel('Relative Flux')
	
	return sunWavsL, sunFluxL, sunWavsC, sunFluxC

def random():
	sunWavsL1, sunFluxL1, sunWavsC1, sunFluxC1 = plotLimbs(18, 35, '6')
	sunWavsL2, sunFluxL2, sunWavsC2, sunFluxC2 = plotLimbsN(18, 35, '6')


	fig, (ax1, ax2) = plt.subplots(nrows = 2, ncols = 1, sharex=True, sharey=False, squeeze=True)

	ax1.plot(sunWavsC1, sunFluxC1, color = 'red', label = "Sun Center")
	ax1.plot(sunWavsL1, sunFluxL1, color = 'blue', label = "Sun Limb")
	ax1.set_ylabel('Relative Flux')
	ax1.legend(loc = 'best')
	##
	ax2.plot(sunWavsC2, sunFluxC2, color = 'red', label = "Sun Center")
	ax2.plot(sunWavsL2, sunFluxL2, color = 'blue', label = "Sun Limb")
	ax2.set_ylabel('Relative Flux')
	ax2.legend(loc = 'best')
	##
	plt.xlim(620,642.5)
	plt.xlabel("$\lambda$ [nm]")


def plotProcessed(nthImg, orderstr):
	sunFlux, sunWavs = process(nthImg, orderstr)
	centroids = getSunCentroids(nthImg, orderstr)
	centroids = pixToWave(orderstr, centroids)

	for idx,element in enumerate(centroids):
		if idx == 0:
			plt.axvline(element, linestyle = '-.', color = "purple", label = "Centroids")
		else:
			plt.axvline(element, linestyle = '-.', color = "purple")
	plt.legend(loc = 1)

	plt.plot(sunWavs, sunFlux, label = "Sun Spectrum")
	plt.legend(loc = "best", ncol = 2)
	#plt.xlim(0,1048)
	#plt.ylim(-100, 12500)
	plt.xlabel('$\lambda$ [nm]')
	plt.ylabel('Relative Flux')


def crosscorr(nlimb, ncent, orderstr):
	centy, centx = process(ncent, orderstr)
	limby, limbx = process(nlimb, orderstr)
	#use peak at 548 and in range 512
	#centy = centy[::-1]
	#limby = limby[::-1]
	cross_corr = np.correlate(centy[548-512:548+512], limby[548-512:548+512], 'same')
	# Create the pixel axis for computing the lag
	n_pixels = len(cross_corr)
	if n_pixels % 2: # Number of pixels is odd
		shift_axis = np.arange(-n_pixels/2+1, n_pixels/2+1)
	else: # Number of pixels is even
		shift_axis = np.arange(-n_pixels/2, n_pixels/2)

	midPoint = n_pixels/2.
	x = np.arange(-10, 10)
	y = cross_corr[midPoint-10:midPoint+10]
	centroid = (np.sum(x*y))/(np.sum(y))
	'''
	plt.plot(shift_axis, cross_corr, 'gx-', label = 'Cross Correlation')
	plt.xlabel('Lag [$pixels$]')
	plt.ylabel('Cross-correlation [$ADU^2$]')
	plt.axvline(centroid, linestyle = '-.', color = "k", label = 'Centroid: -0.311 pixels' )	
	plt.legend(loc='best')
	'''
	return cross_corr, shift_axis, centroid



def pixShift():
	path = 'data/solar/Nov-28-2017/transit1/'

	files = os.listdir(path)

	pixels = np.array([])
	times = np.array([])

	for i in np.arange(17,54):
		file = pf.open(path + files[i])
		times = np.append(times, file[0].header['jd'])
		corr, shift, cent = crosscorr(i, 35, '6')
		pixels = np.append(pixels, cent)

	temp = pf.open(path+files[1])
	time0 = temp[0].header['jd']
	times = (times - time0)*86400
	pixels = pixels - pixels[len(pixels)/2]
	
	coeffs, cov = np.polyfit(times, pixels, 1, cov = True)
	offsetFit = coeffs[0]*times + coeffs[1]
	Err = np.sqrt(np.diag(cov))

	plt.plot(times, pixels, 'x')
	plt.plot(times, offsetFit)
	plt.xlabel('Time [s]')
	plt.ylabel('Pixel Shift [pixels]')
	# m = 0.0010861299873041041
	# mErr = 7.82625590e-05
	# c = -0.14828319964233436
	# cErr = 1.07434099e-02

	return pixels, times, Err


dpdt = 0.0010861299873041041
disp = 0.02118744
time = 137.00002133846283
c = 299792
l0 = 631.685
xi = np.deg2rad(1.26)
eta = np.deg2rad(73.94)
dec = np.deg2rad(-21.2533)
diam = 0.009431390819379791 #radians
'''
dpdt = 0.00109
disp = 0.02188348
time = 133
c = 299792
l0 = 649.6020659880282
xi = np.deg2rad(3.36)
eta = np.deg2rad(226.33)
dec = np.deg2rad(-21.2533)
diam = 0.009431390819379791 #radians
'''
def vrot():
	vrot1 = (dpdt*disp*time*c)/(l0*np.cos(xi)*np.cos(eta))
	vrot2 = vrot1/(15*np.cos(dec))
	return vrot1

def Rsun():
	v=vrot()
	rsun = (1.9*26.24*24.*60.*60.)/(2.*np.pi)
	return rsun

def AU():
	r = Rsun()
	d = 2*r/np.sin(diam)
	return d

