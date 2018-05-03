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
import pdb
from numpy.linalg import inv, det, norm

rcParams['axes.labelsize'] = 18
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15
rcParams['figure.titlesize'] = 20
rcParams['savefig.dpi'] = 600

#path on home computer
datapath = "C:/Users/Sameen_Yunus/Desktop/Fall17/Astro120/Lab3/data/"
#path on UGAstro
#datapath = "/home/global/optical_2017"

#10/10/2017
path1 = "10-10-2017/0032.fts"
darks1 = "10-10-2017/darks/" #22 files

#10/18/2017
path2 = "10-18-2017/0133.fts"
darks2 = "10-18-2017/darks/" #11 files

#10/25/2017
path3 = "10-25-2017/0052.fts"
darks3 = "10-25-2017/darks/" #9 files


def loadFits(file):
	#needs full argument of path and filename
	hdulist = pf.open(file)
	data = np.array(hdulist[0].data)
	header = hdulist[0].header
	return data, header


def getDark():
	#average over files
	pathdarks = datapath+darks1
	files = os.listdir(pathdarks)
	#create an empty array of image size
	size = loadFits(pathdarks+files[0])[0].shape
	dark = np.zeros(size)
	#loop over all darks and get sum of each pixel
	for file in files:
		dark += loadFits(pathdarks+file)[0]
	#get average over all pixels
	dark /= len(files)
	'''
	im = plt.imshow(dark, cmap = "gray_r", aspect = "equal", vmin = 600, vmax= 750,interpolation = "bilinear",  origin = "lower")
	clb =plt.colorbar(im)
	clb.ax.set_ylabel('Counts [ADU]')
	plt.xlabel("x Pixel [pixel units]")
	plt.ylabel("y Pixel [pixels units]")
	plt.title("Dark Current")
	'''
	return dark


def plotField():
	path = datapath+path1
	data, header = loadFits(path)
	data = data - getDark()
	data = np.rot90(data)
	data = np.flipud(data)
	'''
	plt.imshow(data, origin = "lower", cmap = "gray_r")
	clb = plt.colorbar()
	clb.ax.set_ylabel('Counts [ADU]')
	plt.clim(0,100)
	plt.xlabel("x Pixel [pixel units]")
	plt.ylabel("y Pixel [pixels units]")
	plt.title("Dark subtracted 25 phocaea field")
	'''
	return data, header


def getStars(image, max_stars, star_radius, gaussian_sigma):
    #
    # James R. Graham 11/7/2011 University of Toronto
    #
    # image - image
    # max_stars - number of stars to find
    # star_radius - radius of exclusion zone
    # gaussian_sigma - sigma in the Gaussian smoothing function 
    #
    # Step 1) Find the stars in an image starting with the brightest peak.
    # Step 2) Remove that star and the search for the next brightest peak.

    filtered_image = ndi.filters.gaussian_filter(image, gaussian_sigma) # filter the image
    
    starxy = np.zeros([max_stars,2])

    x = np.arange(np.float(image.shape[1]))  # row/column order
    y = np.arange(np.float(image.shape[0]))
    xx, yy = np.meshgrid(x, y)
    
    for star in np.arange(max_stars):
        coordinate = np.unravel_index(np.argmax(filtered_image), image.shape)  # only pick out one value

        starxy[star,0] = coordinate[1]   # x - row/column order
        starxy[star,1] = coordinate[0]   # y 

        r2 = (xx-coordinate[1])**2. + (yy-coordinate[0])**2.
       
        filtered_image[np.where(r2 <= star_radius**2.)] = -999.0 # zero out that part of the array 

    return starxy, filtered_image
	
def plotStars():
	ccd_data, header = loadFits(datapath+path1)
	darkSubData = ccd_data - getDark()
	image = darkSubData - np.median(darkSubData)
	image = np.rot90(image)
	image = np.flipud(image)
	star_radius = 18. # radius of region to blank out when removing star from image
	gaussian_sigma = 4. # smoothing rms to use when search for stars
	max_stars = 15 # maximum number of stars to locate

	# find the brightest stars and return smoothed image

	starxy, filtered_image = getStars(image, max_stars, star_radius, gaussian_sigma)
	'''
	plt.figure()
	img = plt.imshow(image, origin = "lower", cmap="gray_r") #Set colormap

	plt.plot(starxy[:,0], starxy[:,1], 'ro', markersize=10, markerfacecolor='none')

	clb = plt.colorbar()
	clb.ax.set_ylabel('Counts [ADU]')
	plt.clim(0,50)

	plt.xlabel("x Pixel [pixel units]")
	plt.ylabel("y Pixel [pixels units]")
	plt.xlim(0,1336)
	plt.ylim(0,2004)
	'''
	return starxy, image, filtered_image


def getCentroids(image, starxy, star_radius, annulus_width):
    # James R. Graham 11/7/2011 University of Toronto
    #
    # Measure star centroids
    #
    # starrad  radius of sky aperture 
    # nsky     number of pixels insky annulus relative to star 

    x = np.arange(0, image.shape[1])  # row/column order
    y = np.arange(0, image.shape[0])

    sky_radius = np.sqrt(annulus_width+1)*star_radius

    xx, yy = np.meshgrid(x, y)

    x_centroids = np.array([]) # x-centroid
    y_centroids = np.array([]) # y-centroid
    starflux = np.array([]) # star counts in aperture
    rms_x = np.array([]) # rms width in x
    rms_y = np.array([]) # rms width in y

    i = 1
    for star in starxy: 

        r2 = (xx-star[0])**2. + (yy-star[1])**2.

        wstar = np.where(  r2 <= star_radius**2.)
        wsky  = np.where( (r2 >  star_radius**2.) & (r2 < sky_radius**2.) )

        # measure the centroid 

        medsky = np.median(image[wsky])

        # print 'Star %d'%i,star,' Median sky = ',medsky

        # compute the moments

        si   = np.sum((image[wstar] - medsky))						#sum over I
        six  = np.sum((image[wstar] - medsky)*xx[wstar])/si 		#x centroid
        siy  = np.sum((image[wstar] - medsky)*yy[wstar])/si 		#y centroid
        six2 = np.sum((image[wstar] - medsky)*xx[wstar]**2.)/si  	#x squared centroid
        siy2 = np.sum((image[wstar] - medsky)*yy[wstar]**2.)/si 	#y squared centroid

        rms_x     = np.append(rms_x, np.sqrt(six2 - six**2. )) 		#like the variance in x
        rms_y     = np.append(rms_y, np.sqrt(siy2 - siy**2. )) 		#like the variance in y	

        x_centroids = np.append(x_centroids, six)
        y_centroids = np.append(y_centroids, siy)
        starflux = np.append(starflux,si) 
        i += 1

    return x_centroids, y_centroids, rms_x, rms_y, starflux, medsky


def plotCentroids():
	star_radius = 20. # radius of circular region to measure centroid
 	annulus_width = 10. # size of sky region Area(sky) = nsky x Area(star)
 	starxy, image, filtered_image = plotStars()
	x_centroids, y_centroids, rms_x, rms_y, starflux, medsky= getCentroids(image, starxy, star_radius, annulus_width)

	#filtered_image = ndi.filters.gaussian_filter(image, gaussian_sigma) # filter the image
	'''
	plt.figure()
	img = plt.imshow(image, origin = "lower", cmap="gray_r") #Set colormap

	plt.plot(x_centroids, y_centroids, 'rx', markersize = 8, mew=2)

	clb = plt.colorbar()
	clb.ax.set_ylabel('Counts [ADU]')
	plt.clim(0,50)

	plt.xlabel("x Pixel [pixel units]")
	plt.ylabel("y Pixel [pixels units]")
	plt.xlim(0,1336)
	plt.ylim(0,2004)
	'''
	return x_centroids, y_centroids, rms_x, rms_y, starflux, medsky



def usno(radeg,decdeg,fovam,epoch):

    # James R. Graham 2013/10/13 UC Berkeley
    
    # get USNO B-1 stars centered at
    # radeg and decdeg (J2000.0) in degrees
    # centered in a square field of view (arc min) 
    #
    # Corrects for proper motion to current epoch
  
    a1 = 'http://webviz.u-strasbg.fr/viz-bin/asu-tsv/?-source=USNO-B1'
    a2 = '&-c.ra={:4.6f}&-c.dec={:4.6f}&-c.bm={:4.7f}/{:4.7f}&-out.max=unlimited'.format(radeg,decdeg,fovam,fovam)
    

    a0 = a1+a2

    print 'Calling Vizier',a0

    f = url.urlopen(a0)

	# Read from the object, storing the page's contents in 's'.
    s = f.read() 
    f.close()

    sl = s.splitlines()
    sl = sl[45:-1]       # get rid of header  - updated Oct 2013

    name = np.array([])  # star name
    raD  = np.array([])  # RA in degrees
    decD  = np.array([])  # DEC in degrees
    rmag = np.array([])  # rmage

    for k in sl:
        kw = k.split('\t')

        ded0 = float(kw[2])

        pmrad = float(kw[6])/3600e3/np.cos(np.deg2rad(ded0))   # convert from mas/yr to deg/year
        pmded = float(kw[7])/3600e3   

        name  = np.append(name,kw[0])
        raD   = np.append(raD,float(kw[1]) + pmrad*(epoch-2000.0)     ) 
        decD   = np.append(decD,float(kw[2]) + pmded*(epoch-2000.0)     )
 
        if kw[12] != '     ':             # case when no mag is reported
            rmag = np.append(rmag,float(kw[12])) 
        else:
            rmag = np.append(rmag,np.nan) 

    return name,raD,decD,rmag


def plotUSNO():
	
	###10-10-2017
	ras = '20:49:06.93'	 #20:49:06.93 +09:22:39.8
	des = '+09:22:39.8'

	###10-18-2017
	#ras = '21:00:21.40'  #21:00:21.40 +07:18:19.9
	#des = '+07:18:19.9'

	###10-25-2017
	#ras = '21:10:44.45'  #21:10:44.45 +05:50:29.5                 
	#des = '+05:50:29.5'	 
	

	#convert to degrees
	radeg = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
	dsgn = np.sign(float(des[0:3]))
	dedeg = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.	
	fovam = 22.0 # size of square search field in arc min
	name,rad,ded,rmag = usno(radeg,dedeg,fovam,2017.10)

	ra = np.deg2rad(rad)
	dec = np.deg2rad(ded)
	ra0 = np.deg2rad(radeg)
	dec0 = np.deg2rad(dedeg)
	#convert to cartesian
	X = -(np.cos(dec)*np.sin(ra-ra0))/(np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0)+np.sin(dec)*np.sin(dec0))
	Y = -(np.sin(dec0)*np.cos(dec)*np.cos(ra-ra0) - np.cos(dec0)*np.sin(dec))/(np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0)+np.sin(dec)*np.sin(dec0))

	#X = X[::-1]
	#Y = Y[::-1]
	#convert to pixels
	f = 6300 #in mm
	p = 0.018 #in mm
	x0 = 668 + 38
	y0 = 1002 - 121
	x = f*(X/p)+x0
	y = f*(Y/p)+y0

	w = np.where(rmag < 12.5)[0] # select only bright stars r < 12 mag.
	
 	'''
	#plot USNO in degrees
	plt.plot(rad[w],ded[w],'g.')
	plt.locator_params(axis='x',nbins=4) # Nice tick marks
	plt.locator_params(axis='y',nbins=4)
	plt.tick_params('x',pad=10)
	plt.xlabel('RA [Deg]')
	plt.ylabel('Dec [Deg]')
	plt.ticklabel_format(useOffset=False)
	plt.axis('scaled') # Make the aspect ratio correct
	ax = plt.gca()
	ax.set_xlim(ax.get_xlim()[::-1]) # reverse the x-axis direction
	x0,x1 = ax.get_xlim()
	y0,y1 = ax.get_ylim()
	ax.set_aspect(abs(x1-x0)/abs(y1-y0))
	'''
	'''
	#plot USNO in cartesian
	plt.plot(X[w],Y[w], 'r.')
	plt.locator_params(axis='x',nbins=4) # Nice tick marks
	plt.locator_params(axis='y',nbins=4)
	plt.tick_params('x',pad=10)
	plt.xlabel('X')
	plt.ylabel('Y')
	plt.ticklabel_format(useOffset=False)
	plt.axis('scaled') # Make the aspect ratio correct
	ax = plt.gca()
	'''
	'''
	plt.plot(x[w],y[w], 'ro')
	plt.locator_params(axis='x',nbins=4) # Nice tick marks
	plt.locator_params(axis='y',nbins=4)
	plt.tick_params('x',pad=10)
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	plt.ticklabel_format(useOffset=False)
	plt.xlim(0,1336)
	plt.ylim(0,2004)
	plt.axis('scaled') # Make the aspect ratio correct
	'''
	return X[w], Y[w], x[w], y[w]
	
	

def plotUSNOcent():
	x_cent, y_cent, rms_x, rms_y, starflux = plotCentroids()
	X_USNO, Y_USNO, x_USNO, y_USNO = plotUSNO()
	plt.plot(x_cent, y_cent, 'rx', mew=2, label = 'CCD')
	plt.plot(x_USNO, y_USNO, 'b+', mew=2, label = 'USNO')
	#plt.locator_params(axis='x',nbins=4) # Nice tick marks
	#plt.locator_params(axis='y',nbins=4)
	#plt.tick_params('x',pad=10)
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	#plt.ticklabel_format(useOffset=False)
	#plt.axis('scaled') # Make the aspect ratio correct
	plt.legend(loc=1, numpoints=1)
	plt.xlim(0,1336)
	plt.ylim(0,2004)
	return x_cent, y_cent, x_USNO, y_USNO

'''
#10-10-2017
x_cent = np.array([440.1549305, 1192.812012, 177.6366562, 982.5542522, 819.1538439, 1163.422387, 685.0084847, 892.635618])
y_cent = np.array([290.2608806, 1367.694057, 1250.174549, 1680.239492, 1711.615728, 908.5897848, 193.4551821, 785.6629279])
x_USNO = np.array([434.8567728, 1190.899814, 186.7396168, 986.9460518, 825.83967, 1155.996295, 674.2936448, 887.1398716])
y_USNO = np.array([291.0542848, 1345.822837, 1242.361642, 1657.108089, 1689.759048, 893.4802231, 192.6672854, 775.2251887])
X_USNO = np.array([-0.00215183779, 0.00000828518262, -0.00286074395, -0.000574439852, -0.0010347438, -0.0000914391581, -0.00146773244, -0.000859600367])
Y_USNO = np.array([-0.00291127, 0.00010235, -0.00019325, 0.00099174, 0.00108503, -0.00119006, -0.00319238, -0.00152793])

#10-18-2017
x_cent = np.array([41.76101999, 944.6712395, 1323.842533, 952.1251489, 1069.197712, 75.46685817, 184.9256435, 703.8967158, 321.2540469])
y_cent = np.array([257.4653784, 1265.531207, 968.785584, 1870.479906, 1585.180593, 56.61733393, 729.5358071, 1431.610438, 1543.922738])
x_USNO = np.array([1071.193596, 944.0541046, 710.1575105, 332.066081, 187.5743626, 958.8523097, 71.00893862, 40.22862684, 1314.728384])
y_USNO = np.array([1558.59207, 1244.166521, 1412.945311, 1527.156176, 723.6234258, 1842.045431, 60.55878998, 259.537365, 946.4033153])
X_USNO = np.array([-0.00162792, 0.00095444, 0.00201351, 0.00099672, 0.0013177, -0.00153997, -0.00120693, 0.00028616, -0.0007941])
Y_USNO = np.array([-0.00212989, 0.00068333, -0.00016742, 0.00239156, 0.00158169, -0.0026984, -0.00080393, 0.00116556, 0.00149187])
'''
#10-25-2017
x_cent = np.array([1262.873865, 450.0135574, 354.8243826, 834.7322436, 812.2458684, 944.2793293, 497.9318032, 39.0939353])
y_cent = np.array([652.3405811, 226.3386064, 1588.119371, 1591.093729, 164.9852918, 105.6467709, 1574.073768, 1912.22766819])
x_USNO = np.array([1246.264979, 438.9642436, 361.9333422, 835.1990386, 796.0157889, 925.3072942, 502.8703393, 53.71593487])
y_USNO = np.array([642.7665454, 232.5224186, 1578.368407, 1574.962937, 168.3386572, 107.8043005, 1562.663722, 1902.337174])
X_USNO = np.array([0.00154361, -0.00076296, -0.00098305, 0.00036914, 0.00025719, 0.00062659, -0.00058037, -0.00186367])
Y_USNO = np.array([-0.00068067, -0.00185279, 0.00199248, 0.00198275, -0.00203618, -0.00220913, 0.00194761, 0.00291811])


def leastSq():
	#x/y is dependent variable: pixels
	#X/Y is independent variable: coords
	matx = x_cent
	maty = y_cent
	matX = X_USNO
	matY = Y_USNO
	f = 6300
	p = 0.018
	# c = ( a11 a12 x0)
	# d = ( a21 a22 y0)
	ones = np.ones(matX.size)
	B = np.column_stack((f*matX/p, f*matY/p, ones))
	# c = (B(transpose)B)(inverse)B(transpose)matx
	
	invTransB = np.dot(inv(np.dot(B.T, B)),B.T)
	c = np.dot(invTransB, matx)
	d = np.dot(invTransB, maty)

	return c, d

def trans():
	c, d = leastSq()
	f = 6300
	p = 0.018
	a11 = c[0]
	a12 = c[1]
	a21 = d[0]
	a22 = d[1]
	x0 = c[2]
	y0 = d[2]

	#transformation matrix using plate constants defined above
	T = np.array(([f*a11/p, f*a12/p, x0], [f*a21/p, f*a22/p, y0], [0,0,1]))
	fP = np.sqrt(det(T))
	return fP, T

def resid():
	fP, T = trans()
	X = np.column_stack((X_USNO, Y_USNO, np.ones(X_USNO.size)))
	x = np.column_stack((x_cent, y_cent, np.ones(x_cent.size)))
	xtrans = np.dot(T, X.T)
	xfit = xtrans.T
	xUSNOfit = xfit[:,0]
	yUSNOfit = xfit[:,1]
	xRes = x_cent - xUSNOfit
	yRes = y_cent - yUSNOfit

	xRMS = np.sqrt((xRes**2)/xRes.size)
	yRMS = np.sqrt((yRes**2)/yRes.size)

	xRMS1 = np.sqrt(np.sum(xRes**2)/xRes.size)
	yRMS1 = np.sqrt(np.sum(yRes**2)/yRes.size)
	'''
	plt.plot(x_cent, xRes, 'rv', label = 'x residuals', markersize=10, markeredgecolor='r')
	plt.plot(y_cent, yRes, 'b^', label = 'y residuals', markersize=10, markeredgecolor='b')
	plt.legend(loc=1, numpoints=1)
	plt.xlim(0,2004)
	plt.ylim(-1,1)
	plt.xlabel('x or y [pixel]')
	plt.ylabel('Residuals [pixel]')
	'''
	return xUSNOfit, yUSNOfit, xRMS, yRMS

def plotUSNOfit():
	xfit, yfit, xRMS, yRMS = resid()
	x_cent, y_cent, rms_x, rms_y, starflux = plotCentroids()
	
	plt.plot(x_cent, y_cent, 'rx', mew=2, label = 'CCD')
	plt.plot(xfit, yfit, 'b+', mew=2, label = 'USNO')
	#plt.locator_params(axis='x',nbins=4) # Nice tick marks
	#plt.locator_params(axis='y',nbins=4)
	#plt.tick_params('x',pad=10)
	plt.xlabel('x [pixel]')
	plt.ylabel('y [pixel]')
	#plt.ticklabel_format(useOffset=False)
	#plt.axis('scaled') # Make the aspect ratio correct
	plt.legend(loc=1, numpoints=1)
	plt.xlim(0,1336)
	plt.ylim(0,2004)
	return x_cent, y_cent
	
#pixel coords of 25 Phocaea
Ax1 = 637.7872011	
Ex1 = 4.78677207
Ay1 = 840.5330361	#1010
Ey1 = 4.98936914
Ax2 = 613.9706193	
Ex2 = 4.75305943
Ay2 = 1014.453312	#1018
Ey2 = 4.50672353
Ax3 = 565.2971165	
Ex3 = 3.13883316
Ay3 = 677.1151069	#1025
Ey3 = 2.86866458

#degree coords of 25 Phocaea
Adx1 = 312.278874999999	
Edx1 = 5.452346762
Ady1 = 9.37772222222222	#1010
Edy1 = 0.161212818
Adx2 = 315.088390676933	
Edx2 = 5.501385047
Ady2 = 7.30477951429136	#1018
Edy2 = 0.125044919
Adx3 = 317.7095104122325	
Edx3 = 5.546693279
Ady3 = 5.80804693302507	#1025
Edy3 = 0.099488791

#radian coords of 25 Phocaea
Arx1 = np.deg2rad(312.278874999999)
Erx1 = np.deg2rad(5.452346762)
Ary1 = np.deg2rad(9.37772222222222)	#1010
Ery1 = np.deg2rad(0.161212818)
Arx2 = np.deg2rad(315.088390676933)	
Erx2 = np.deg2rad(5.501385047)
Ary2 = np.deg2rad(7.30477951429136)	#1018
Ery2 = np.deg2rad(0.125044919)
Arx3 = np.deg2rad(317.7095104122325)		
Erx3 = np.deg2rad(5.546693279)
Ary3 = np.deg2rad(5.80804693302507)	#1025
Ery3 = np.deg2rad(0.099488791)


stdxRMS = 0.028879822238979301
stdyRMS = 0.065964808192441729
def PixtoCel():
	#convert from pixels to standard coordinates to celestial coordinates using the linear least squares fit
	fP, T = trans()
	#x = TX where x is pixels and X is standard
	xMat = np.array([Ax3, Ay3, 1])
	XMat = np.dot(inv(T),xMat)
	
	xEMat = np.array([stdxRMS, stdyRMS, 1])
	XEMat = np.dot(inv(T),xEMat)

	X = XMat[0]
	Y = XMat[1]
	XE = XEMat[0]
	YE = XEMat[1]


	#now convert from standard coords to RA and Dec
	###10-10-2017
	#ras = '20:49:06.93'	 #20:49:06.93 +09:22:39.8
	#des = '+09:22:39.8'

	###10-18-2017
	#ras = '21:00:21.40'  #21:00:21.40 +07:18:19.9
	#des = '+07:18:19.9'

	###10-25-2017
	ras = '21:10:44.45'  #21:10:44.45 +05:50:29.5                 
	des = '+05:50:29.5'	 
	

	#convert to degrees
	rad0 = 15*(float(ras[0:2]) + float(ras[3:5])/60. + float(ras[6:])/3600.)
	dsgn = np.sign(float(des[0:3]))
	decd0 = float(des[0:3]) + dsgn*float(des[4:6])/60. + dsgn*float(des[7:])/3600.	

	ra0 = np.deg2rad(rad0)
	dec0 = np.deg2rad(decd0)

	#ra and dec in radians
	raR = ra0 + np.arctan(-X/(np.cos(dec0)-Y*np.sin(dec0)))
	decR = np.arcsin((np.sin(dec0)+Y*np.cos(dec0))/np.sqrt(1+X**2+Y**2))
	raER = ra0 + np.arctan(-XE/(np.cos(dec0)-YE*np.sin(dec0)))
	decER = np.arcsin((np.sin(dec0)+YE*np.cos(dec0))/np.sqrt(1+XE**2+YE**2))


	#ra and dec in degrees
	raD = np.rad2deg(raR)
	decD = np.rad2deg(decR)
	raED = np.rad2deg(raER)
	decED = np.rad2deg(decER)

	#find error in asteroid position

	return rad0, decd0, raD, decD, raER, decER

def deg2HMS(ra='', dec='', round=False):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    if round:
      decS = int((abs((dec-deg)*60)-decM)*60)
    else:
      decS = (abs((dec-deg)*60)-decM)*60
    DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, decS)
  
  if ra:
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    if round:
      raS = int(((((ra/15)-raH)*60)-raM)*60)
    else:
      raS = ((((ra/15)-raH)*60)-raM)*60
    RA = '{0}{1} {2} {3}'.format(rs, raH, raM, raS)
  
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC



#Use Ax1, Ax2, Ax3 and Ay1, Ay2, Ay3 to compute parallax
def parallax():

	#first epoch
	X1 =	0.952806104140452
	Y1 =	0.285128234518564
 	Z1 =	0.123459343396735
 	R1 = np.array([X1, Y1, Z1])
 	#second epoch
 	X2 =	0.899585999316349
	Y2 =	0.40193079893623
	Z2 =	0.174090123545291
	R2 = np.array([X2, Y2, Z2])
	#third epoch
	X3 =	0.838830658413973
	Y3 =	0.49796692555458
	Z3 =	0.215723949818143
	R3 = np.array([X3, Y3, Z3])

	#define unit vectors s1, s2 and s3
	x1 = np.cos(Arx1)*np.cos(Ary1)
	y1 = np.sin(Arx1)*np.cos(Ary1)
	z1 = np.sin(Ary1)
	s1 = np.array([x1, y1, z1])
	xe1 = np.cos(Erx1)*np.cos(Ery1)
	ye1 = np.sin(Erx1)*np.cos(Ery1)
	ze1 = np.sin(Ery1)
	se1 = np.array([xe1, ye1, ze1])

	x2 = np.cos(Arx2)*np.cos(Ary2)
	y2 = np.sin(Arx2)*np.cos(Ary2)
	z2 = np.sin(Ary2)
	s2 = np.array([x2, y2, z2])
	xe2 = np.cos(Erx2)*np.cos(Ery2)
	ye2 = np.sin(Erx2)*np.cos(Ery2)
	ze2 = np.sin(Ery2)
	se2 = np.array([xe2, ye2, ze2])

	x3 = np.cos(Arx3)*np.cos(Ary3)
	y3 = np.sin(Arx3)*np.cos(Ary3)
	z3 = np.sin(Ary3)
	s3 = np.array([x3, y3, z3])
	xe3 = np.cos(Erx3)*np.cos(Ery3)
	ye3 = np.sin(Erx3)*np.cos(Ery3)
	ze3 = np.sin(Ery3)
	se3 = np.array([xe3, ye3, ze3])


	#define epoch times 
	
	T1 = '04:08:36.27' #on 10/10
	T2 = '04:17:01.90' #on 10/18 so +7
	T3 = '03:45:09.32' #on 10/25 so +15

	t1 = float(T1[0:2])/24. + float(T1[3:5])/1440. + float(T1[6:])/86400.
	t2 = 7 + float(T2[0:2])/24. + float(T2[3:5])/1440. + float(T2[6:])/86400.
	t3 = 15 + float(T3[0:2])/24. + float(T3[3:5])/1440. + float(T3[6:])/86400.

	tau1 = t2 - t1
	tau3 = t3 - t2
	#s dot is velocity
	v2 = tau3*(s2-s1)/(tau1*(tau1+tau3)) + tau1*(s3-s2)/(tau3*(tau1+tau3))
	ve2 = tau3*(se2-se1)/(tau1*(tau1+tau3)) + tau1*(se3-se2)/(tau3*(tau1+tau3))
	# s double dot is acceleration
	a2 = 2*(s3-s2)/(tau3*(tau1+tau3)) - 2*(s2-s1)/(tau1*(tau1+tau3))
	ae2 = 2*(se3-se2)/(tau3*(tau1+tau3)) - 2*(se2-se1)/(tau1*(tau1+tau3))
	
	return s1, s2, s3, v2, a2, R2

def iterate():
	s, v, a, R = parallax()
	k = 0.017202098950

	r = 1.878477410275053

	rho = (k**2)*(1/norm(R)**3 - 1/r**3)*(np.dot(v, np.cross(R, s)))/(np.dot(v, np.cross(a, s)))
	
	r = np.sqrt(rho**2 + norm(R)**2 + 2*rho*np.dot(R, s))
	return rho, r

