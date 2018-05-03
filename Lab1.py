import numpy as np
import scipy.misc
import matplotlib.pyplot as plt 

#variables
file_name = 'pmt_17083_170831_1956_40_max.csv'
t = np.loadtxt(file_name, delimiter = ',', dtype = 'int32', usecols=[1])
dt = t[1:] - t[0:-1]
dt_gated = dt[dt > 3000]

'''
#Raw event time plot
plot = plt.figure()
plt.plot(t, color="red")
plt.xlabel("Event Number", size = 15)
plt.ylabel("Time [clock ticks]", size = 15)
plt.title("Raw event time in clock ticks vs Event number", size = 15)
plt.tight_layout()
plt.show()



#Time interval plot
plot = plt.figure()
plt.plot(dt, ",", color="red")
plt.xlabel('Event number', size = 15)
plt.ylabel('Interval [clock ticks]', size = 15)
plt.title("Time Interval vs Event Number", size = 15)
plt.show()



figa, (plt1, plt2) = plt.subplots(2,1)

#Mean interval per chunk of nstep
nstep = 1000
marr = np.array([]) 			# define an empty array
charr = np.arange(0, dt.size, nstep) # define chunks array
i = np.arange(dt.size)
for j in i[0::nstep]: 			# iterate in strides of nstep
	m = np.mean(dt[j:j+nstep]) 	# compute the mean
	marr = np.append(marr,m) 	# append the new value
plt1.plot(charr, marr,'o')
plt1.set_xlabel('Start Index', size = 15)
plt1.set_ylabel('Mean Interval [clock ticks]', size = 15)
plt1.set_title("Mean Interval for 1000 events", size = 15)

nstep = 100
marr = np.array([]) 			# define an empty array
charr = np.arange(0, dt.size, nstep) # define chunks array
i = np.arange(dt.size)
for j in i[0::nstep]: 			# iterate in strides of nstep
	m = np.mean(dt[j:j+nstep]) 	# compute the mean
	marr = np.append(marr,m) 	# append the new value
plt2.plot(charr, marr,'o')
plt2.set_xlabel('Start Index', size = 15)
plt2.set_ylabel('Mean Interval [clock ticks]', size = 15)
plt2.set_title("Mean Interval for 100 events", size = 15)
plt.tight_layout()
plt.show()



#Mean intervals against number of intervals averaged

nstep = 100
marr = np.array([]) 			# define an empty array
charr = np.arange(0, dt.size, nstep) # define chunks array
i = np.arange(dt.size)
for j in i[0::nstep]: 			# iterate in strides of nstep
	m = np.mean(dt[0:j+nstep]) 	# compute the mean
	marr = np.append(marr,m) 	# append the new value
plt.plot(charr, marr,'o')
plt.xlabel('Number of Intervals', size = 15)
plt.ylabel('Mean Interval [clock ticks]', size = 15)
plt.title("Mean Intervals against Number of Intervals Averaged in steps of 100", size = 15)
plt.show()



#plot standard deviation against 1/sqrtN

sdarr = np.array([]) 					# define an empty array for mean
Narr = np.arange(10,1000,10)		 	# define N array	
sqrtN = np.array([])					# define empty sqrt array


i = np.arange(dt.size)
for j in range(1,100):
        nstep = j*10
        marr = np.array([])				# mean array has to be initiated each time
        for k in i[0::nstep]:
                m = np.mean(dt[k:k+nstep])
                marr = np.append(marr,m)
        mu = np.sum(marr)/np.float(marr.size)
        sample_std = np.sqrt(np.sum((marr - mu)**2.)/(np.float(marr.size)-1.))
        sdarr = np.append(sdarr, sample_std)
        n = 1/np.sqrt(j+nstep)
        sqrtN = np.append(sqrtN, n)

SDOM = np.std(dt)
#plt.plot(Narr, sdarr,'.')
plt.plot(sqrtN, sdarr, '.', color = 'green')
plt.plot(sqrtN, SDOM*sqrtN, color = 'red')
plt.xlabel('$1/\sqrt{N}$', size = 15)
plt.ylabel('Standard Deviation of the Mean [ticks]', size = 15)
plt.xlim(0.03,0.31)
plt.title("Standard Deviation versus square root N", size = 15)
plt.tight_layout()
plt.show()




#Create Histograms
N = 500

#Gate the afterpulse events
dt_gated = dt[dt>3000]
# define the lower and upper bin edges and bin width
bw = (dt_gated.max()-dt_gated.min())/(N-1.)
binl = dt_gated.min() + bw * np.arange(N)
# define the array to hold the occurrence count
bincount = np.array([])
# loop through the bins
for bin in binl:
	count = np.where((dt_gated >= bin) & (dt_gated < bin+bw))[0].size
	bincount = np.append(bincount,count)
#compute bin centers for plotting
binc = binl + 0.5*bw
plt.figure()
plt.plot(binc,bincount,drawstyle='steps-mid')
plt.xlabel('Interval [ticks]', size = 15)
plt.xlim(0, 8E6)
plt.ylabel('Frequency', size = 15)
plt.title('Frequency/Interval for 50 bins, $\\tau = 1.04\\times10^6$')


#poisson distribution
def poisson( x, tau ):
	return (1/tau)*np.exp(-x/tau)
tau = np.mean(dt_gated)
N = dt_gated.size*poisson(binc,tau)*bw

plt.plot(binc, N)
#plt.semilogy(binc, N, nonposy = "clip")
plt.show()





#Changing brightness of LED
#load files
marr = np.array([])
sdarr = np.array([])

files = [1956,1957,1959,2001,2003,2007]
for file in files:
	t_new = np.loadtxt('varied_intensity/pmt_17083_170831_%d.csv' %(file), delimiter = ',', dtype = 'int32', usecols=[1])
	dt_new = t_new[1:] - t_new[0:-1]
	dt_newgate = dt_new[dt_new>3000]

	mean = np.mean(dt_newgate)
	marr = np.append(marr, mean)
	std = np.std(dt_newgate)
	sdarr = np.append(sdarr, std)

plt.plot( marr, sdarr, 'o')
plt.plot([0, 1.6*10**7], [0, 1.6*10**7])
plt.xlabel('Interval Sample Mean [ticks]', size = 15)
plt.ylabel('Interval Standard Deviation [ticks]', size = 15)
plt.title('$\mu = \sigma$ Plot')


plt.show()
plt.ticklabel_format(style='sci',scilimits=(0,0))
plt.tight_layout()
'''







'''
fig1, (plota, plotb, plotc) = plt.subplots(3,1)


# Reconstruct the time series (without the 32-bit jumps)
t1 = np.cumsum(dt_gated)
plota.plot(t1)
plota.set_xlabel('Time [ticks]', size = 15)
plota.set_ylabel('Event Number', size = 15)
plota.set_title('Sum of event times versus Photon arrivals')
#plt.ticklabel_format(style='sci',scilimits=(0,0))

def poisson_dist (events, x):
	mean = np.mean(events)
	N = events.size
	return np.exp(-mean)*mean**(x)*N/scipy.misc.factorial(x)

t1_hist = plt.hist(t1, bins = 3000)[0]

plotb.hist(t1, bins=3000, histtype='step', color="red")
plotb.set_xlim(0,2.5E9)
plotb.set_xlabel('Time [ticks]', size = 15)
plotb.set_ylabel('Counts per Bin', size = 15)
plotb.set_xlim(0, 2.5E9)
plotb.set_title('Number of Events per Bin')


x = np.linspace(0, t1_hist.max(), 3000)
y = poisson_dist(t1_hist, x)
plotc.plot(x,y, color = 'blue')
plotc.hist(t1_hist, bins=t1_hist.max(), histtype='step', align='left', color="red")
plotc.set_xlim([-0.5,t1_hist.max()-0.5])
plotc.set_xlabel('Photons per Bin', size = 15)
plotc.set_ylabel('Frequency', size = 15)
plotc.set_title('Frequency of occurence of photons arriving per bin', size = 15)
plt.show()
plt.tight_layout()
'''



fig1, (plota, plotb, plotc) = plt.subplots(3,1)


# Reconstruct the time series (without the 32-bit jumps)
t1 = np.cumsum(dt_gated)
plota.plot(t1)
plota.set_xlabel('Time [ticks]', size = 15)
plota.set_ylabel('Event Number', size = 15)
plota.set_title('Sum of event times versus Photon arrivals')
#plt.ticklabel_format(style='sci',scilimits=(0,0))

# define the lower and upper bin edges and bin width

t1 = np.cumsum(dt_gated)
t1_gated = t1[t1<2.5E9]
N = t1_gated.size
bw = (t1_gated.max()-t1_gated.min())/(N-1.)
binl = t1_gated.min() + bw * np.arange(N)
# define the array to hold the occurrence count
bincount = np.array([])
# loop through the bins
for bin in binl:
	count = np.where((t1_gated >= bin) & (t1_gated < bin+bw))[0].size
	bincount = np.append(bincount,count)


plotb.plot(t1_gated,bincount,drawstyle='steps-mid')
plotb.set_xlabel('Time [ticks]', size = 15)
plotb.set_ylabel('Counts per Bin', size = 15)
plotb.set_xlim(0, 2.5E9)
plotb.set_ylim(0, 6)
plotb.set_title('Number of Events per Bin')


#plot histogram of events per time bin
def poisson_dist (events, x):
	mean = np.mean(events)
	N = events.size
	return np.exp(-mean)*mean**(x)*N/scipy.misc.factorial(x)



freqarr = np.zeros(18)
counts = np.arange(0,18)
for bins in bincount:
	for n in counts:
		if bins == n:
			freqarr[n] += 1

yarr = np.array([])
for x in range(0,18): 
  y = poisson_dist(bincount,x)
  yarr = np.append(yarr, y)

x = np.arange(0,counts.size)
y1 = yarr

plotc.plot(x,y1, color = 'red', linewidth = 2)
plotc.plot(counts, freqarr, drawstyle = 'steps-mid', color = 'blue', linewidth = 2)
plotc.set_xlabel('Photons per Bin', size = 15)
plotc.set_ylabel('Frequency', size = 15)
plotc.set_xlim(0,6)
plotc.set_ylim(-50,1000)
plotc.set_title('Frequency of occurence of photons arriving per bin', size = 15)
plt.tight_layout()
plt.show()



#plt.ticklabel_format(useOffset = False)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))


