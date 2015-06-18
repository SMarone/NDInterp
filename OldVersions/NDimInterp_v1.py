from numpy import linalg as LA
import numpy as np
from  scipy import spatial
from matplotlib import cm
import matplotlib.pyplot as plt
import time
import math
from mpl_toolkits.mplot3d import Axes3D


## Notes:
## x and y currently have similar range, can be changed later
## n dimensional will be next step
print '3 Dimensional KD-Tree with Linear Interpolation'
print ' '

'''
# This was updated with FindNormal function since it is much slower as dimensions increase
def FindNorm(vects):
	dimension = len(vects[1,:])
	concvects = np.concatenate((vects, vects), 1)
	a = np.ones(dimension)
	b = np.ones(dimension)
	for d in range(dimension):
		for n in range(dimension-1):
			a[d] *= concvects[n,d+n+1]
			b[d] *= concvects[n,d-n-1]
	norm = (a - b)/LA.norm(a - b)
	return norm
'''

def FindNormal(vects, r, da, db):
	# It is known that the given array will be of size [dims-1, dims]
	dims = len(vects[1,:])
	# Dims is the same number as neighbors
	proda = np.prod(vects[r,da], axis=1)
	prodb = np.prod(vects[r,db], axis=1)
	norm = (proda - prodb)/LA.norm(proda - prodb)
        return norm


def CreateCrateFunct():
	maxi = 512
	mini = -512
	stp = 6
	# Setup Function Data
	x = range(mini,maxi,stp)
	y = range(mini,maxi,stp)
	# Egg Crate Function
	z = -x * np.sin(np.sqrt(abs(x - y - 47))) - (y+47) * np.sin(np.sqrt(abs(x/2 + y + 47)))
	return x, y, z, maxi, mini


def FunctChk(x, y, z):
	# Crate Visual Check
	x, y = np.meshgrid(x, y)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()


def CreateExp(mini, maxi, exppts):
	# Create Random Experimental Values
	xe = np.random.rand(exppts) * (maxi-mini) + mini
	ye = np.random.rand(exppts) * (maxi-mini) + mini
	xem, yem = np.meshgrid(xe, ye)
	ze = -xe * np.sin(np.sqrt(abs(xe - ye - 47))) - (ye+47) * np.sin(np.sqrt(abs(xe/2 + ye + 47)))
	return xe, ye, ze


def ExpChk():
	# Create Surface plot to Check Unstructured Experimental Data
	fige = plt.figure()
	axe = fige.gca(projection='3d')
	#surf = ax.plot_surface(xem, yem, ze, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	surfe = axe.tricontour(xe.flatten() , ye.flatten() , ze.flatten(), 100, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	fige.colorbar(surfe, shrink=0.5, aspect=5)
	plt.show()


def ExpintoKD(xe, ye, ze, exppts):
	# Form experimental data into KDTree
	numleafs = 8
	#Data Structure for KD is np.array[totalpoints,dimensions]
	#xe = np.append(xe,[3.0])
	#ye = np.append(ye,[4.0])

	# Organize experimental Data
	ExpData = np.transpose(np.array([xe,ye]))
	#print ExpData
	#print np.transpose(ze)

	# Make into Tree
	leafsz = math.ceil(exppts/numleafs)
	KData = spatial.KDTree(ExpData,leafsize=leafsz)

	# Form into dictionary as well
	expdict = {}
	for e in range(exppts):
		expdict["%s" % KData.data[e,:]] = ze[e]

	# Check
	#print ExpData[0:11,0:2]
	#print xe[0:11]
	#print ye[0:11]
	#print leafsz
	#print exppts//leafsz
	#chk = np.array([[0.0,0.0],[512.0,512.0]])
	return expdict, KData


def WriteExp(xe, ye, ze, filename):
	# Write out data in NumPy array of x, y, z
	fe = open(filename, 'w')
	expdata = np.array[xe, ye, ze]
	f.write(expdata)
	f.close()
	return


def ReadExp(filename):
	# Read in NumPy array of x, y, z data
	fe = open(filename, 'r')
	expdata = f.read()
	f.close
	return expdata

def ManInput():
	# Input Variables
	numin = int(input ('Enter the number of points to interpolate: '))
	xin = []
	yin = []
	
	for n in range(numin):
		i = n+1
		print 'Variable %s/%s: ' % (i, numin)
		print ' '
		xin = np.append(xin, float(input('Enter the x value: ')))
		yin = np.append(yin, float(input('Enter the y value: ')))

	return xin, yin


def RandTrains(mini, maxi, trnpts):
	xt = np.random.rand(trnpts) * (maxi-mini) + mini
	yt = np.random.rand(trnpts) * (maxi-mini) + mini
	trndata = np.transpose(np.array([xt,yt]))
	return trndata


def LHTrains(mini, maxi, trnpts):
	# STEPHEN - Improve this later to work with any problem
	mi = int(math.ceil(mini / 100.0)) * 100
	ma = int(math.ceil(maxi / 100.0)) * 100
	step = (ma-mi-100)/(trnpts-1)
	xt = np.array(range(mi,ma,step))
	yt = np.array(range(mi,ma,step))
	# Latin Hypercube it up
	np.random.shuffle(xt)
	np.random.shuffle(yt)
	trndata = np.transpose(np.array([xt,yt]))
	#print xt
	#print yt
	#print TrnData
	return trndata


def TwoNeighbors(expdict, KData, trndata):
	exppts = len(KData.data[:,0])
	trnpts = len(trndata[:,0])
	trndict = {}
	trnz = []
	# KData query takes (location,#ofneighbors) 
	# Determine closest points to ing data
	ndist, nloc = KData.query(trndata,2)
	# Loop through each point to decide z value
	for n in range(len(nloc[:,0])):
		# Grab point information from KData
		pta = KData.data[nloc[n,0],:]
		za = expdict['%s' % pta]
		ptb = KData.data[nloc[n,1],:]
		zb = expdict['%s' % ptb]
		dx = (trndata[n,0] - pta[0])/(ptb[0] - pta[0])
		dy = (trndata[n,1] - pta[1])/(ptb[1] - pta[1])
		trnz.append(((zb-za)/2) * (dy + dx)+za)
		trndict['%s' % trndata[n,:]] =  trnz[n] 
	return trnz, trndict


def ThreeNeighbors(expdict, KData, trndata):
	exppts = len(KData.data[:,0])
	trnpts = len(trndata[:,0])
	trndict = {}
	trnz = []
	# KData query takes (location,#ofneighbors) 
	# Determine closest points to ing data
	ndist, nloc = KData.query(trndata,3)
	# Loop through each point to decide z value
	for n in range(len(nloc[:,0])):
		# Get distance weights
		ndist[n,:] = ndist[n,:]/sum(ndist[n,:])
		# Grab point information from KData
		pta = KData.data[nloc[n,0],:]
		za = expdict['%s' % pta]
		ptb = KData.data[nloc[n,1],:]
		zb = expdict['%s' % ptb]
		ptc = KData.data[nloc[n,2],:]
		zc = expdict['%s' % ptc]
		# Determine the change in z wrt x and y on average
		# Order will be similar to [x or y, a b or c]

		#dzd_ab = (zb-za)*np.array([1/(ptb[0] - pta[0]), 
		#	                   1/(ptb[1] - pta[1])])
		#dzd_bc = (zc-zb)*np.array([1/(ptc[0] - ptb[0]), 
		#	                   1/(ptc[1] - ptb[1])])
		#dzd_ca = (za-zc)*np.array([1/(pta[0] - ptc[0]), 
		#	                   1/(pta[1] - ptc[1])])
		# Get percentage effect of z wrt x and y
		#deff = abs(dzd_ab * ndist[n,0] + 
		#	    dzd_bc * ndist[n,1] + 
		#	    dzd_ca * ndist[n,2])
		#deff = deff / (deff[0] + deff[1])
		# find x and y distance from each point
		deff = [0.5, 0.5]
		d_ab = np.array([(trndata[n,0] - pta[0])/(ptb[0] - pta[0]), 
			         (trndata[n,1] - pta[1])/(ptb[1] - pta[1])])
		d_bc = np.array([(trndata[n,0] - ptb[0])/(ptc[0] - ptb[0]), 
			         (trndata[n,1] - ptb[1])/(ptc[1] - ptb[1])])
		d_ca = np.array([(trndata[n,0] - ptc[0])/(pta[0] - ptc[0]), 
			         (trndata[n,1] - ptc[1])/(pta[1] - ptc[1])])
		# Find 3 possible z value with x and y effects on z factored in
		zd_ab = ((zb-za)*(deff[0]*d_ab[0] + deff[1]*d_ab[1]) + za)
		zd_bc = ((zc-zb)*(deff[0]*d_bc[0] + deff[1]*d_bc[1]) + zb)
		zd_ca = ((za-zc)*(deff[0]*d_ca[0] + deff[1]*d_ca[1]) + zc)
		# Use average z location
		#trnz.append((zd_ab+zd_bc+zd_ca)/3)
		trnz.append(zd_ab*ndist[n,0]+zd_bc*ndist[n,1]+zd_ca*ndist[n,2])
		trndict['%s' % trndata[n,:]] = trnz[n]
	return trnz, trndict


def N_Neighbors(expdata, trndata, numleafs):
	# Expdata has n dimensions, trndata has n-1 dimensions because last
	# dimension is the dependant one we are solving for.  

	# Need same number of neighbors as problem dimensions

	trnpts = len(trndata[:,0])
	exppts = len(expdata[:,0])
	neighs = len(expdata[0,:])
	if (numleafs > exppts):
		# Can not have more leaves than experimental data points. Also 
		# need to ensure there will be enough neighbors in each leaf
		numleafs = exppts // (neighs*2)
		print 'Number of leafs has been changed to %s' % numleafs
		if (numleafs > exppts):
			print 'The problem has too few experimental points for'
			print 'its number of dimensions.'
			raise SystemExit

	print numleafs
	print neighs

	# Make into Tree
	leafsz = math.ceil(exppts/numleafs)
	# Can only use first n-1 dimensions of expdata since that is all
	# which is inside the trndata
	KData = spatial.KDTree(expdata[:,:-1],leafsize=leafsz)

	# KData query takes (data, #ofneighbors) 
	# Determine closest points to ing data
	ndist, nloc = KData.query(trndata ,neighs)
	# Loop through each point to decide z value
	for t in range(trnpts):
		## STEPHEN - you can probably remove the npts middle man here
		# Grab point information from KData
		#npts = np.empty((0,neighs), float)
		#vect = np.empty((0,neighs), float)
		#npts = np.append(npts, [expdata[nloc[t,0],:]], 0)
		nvect = np.empty((0,neighs), float)
		for n in range(neighs-1):
			# Creates array [neighbor, dimension]
			#npts = np.append(npts, [expdata[nloc[t,n+1],:]], 0)
			#vect = np.append(vect, [npts[n+1,:] - npts[n,:]], 0)
			nvect = np.append(nvect, [expdata[nloc[t,n+1],:] - 
				                  expdata[nloc[t,n],:]], 0)
		if debug:
			print ' '
			print npts
			print expdata[nloc[t,1],:]
			print npts[1,1]
			print vect
			print ' '

		normal = FindNorm(nvect)
		# Set the trndata to a temporary array so it has the same
		# dimensions as the experimental data with z as 0 for now
		temp = np.concatenate((trndata[t,:], np.array([0])), 1)
		print temp
		#pta = KData.data[nloc[n,0],:]
		#pta = np.r_[pta, expdict['%s' % pta]]
		#trnz.append(zd_ab*ndist[n,0]+zd_bc*ndist[n,1]+zd_ca*ndist[n,2])
		#trndict['%s' % trndata[n,:]] = trnz[n]
	#return trnz, trndict


# Run Code
print 'Running Code'
numleafs = 8
debug = False
count = 0
mini = -512
maxi = 512
exppts = 10
trnpts = 5
xe, ye, ze = CreateExp(mini, maxi, exppts)
expdata = np.transpose(np.array([xe, ye, ze]))
trndata = LHTrains(mini, maxi, trnpts)
t1 = time.time()
#trnz, trndict = TwoNeighbors(expdict, KData, trndata)
t2 = time.time()
#trnzh, trndicth = ThreeNeighbors(expdict, KData, trndata)
t3 = time.time()
N_Neighbors(expdata, trndata, numleafs)
t4 = time.time()

'''
if debug:
	print '** KData **'
	print KData.data
	ndist, nloc = KData.query(trndata,2)
	print '** ndist **'
	print ndist
	print '** nloc **'
	print nloc
	print '** expdict **'
	print expdict
	print '** trndata **'
	print trndata
	#print '** trnz Two Neighbors **'
	#print trnz
	#print '** trndict Two Neighbors **'
	#print trndict
	#print '** trnz Three Neighbors **'
	#print trnzh
	#print '** trndict Three Neighbors **'
	#print trndicth

for n in range(trnpts):
	xck = trndata[n,0]
	yck = trndata[n,1]
	zck = -xck * np.sin(np.sqrt(abs(xck - yck - 47))) - (yck+47) * np.sin(np.sqrt(abs(xck/2 + yck + 47)))
	twon = abs(trnz[n]-zck)
	threen = abs(trnzh[n]-zck)
	if (twon > threen):
		count += 1
	if debug:
		print "** Two Neighbor vs Three Neighbor Values"
		print twon, threen , (twon > threen)

print ' '
print 'Count', count, 'out of', trnpts, 'points are better using 3 neighbors.'
print ' '
print "## Two Neighbors Time ##"
print t2-t1
print "## Three Neighbors Time ##"
print t3-t2
'''

## Ignore below this, just junk I used as I got more comfortable with python

# Manual input test
#xin, yin = ManInput()
#print xin
#print yin

# return arrays of [neighbordistances], [neighborlocations] respectively
#print KData.query(chk)
#print KData.query(chk[0],3)

# Me just getting comfy with these arrays
#print 'Experimental Data Array'
#print expdata
#print 'X, Y, and Z for Experimental Data'
#print xe, ye, ze
#print xe[3], ye[5], ze[8]
#print expdata[3,0], expdata[5, 1], expdata[8,2]
#expdict, KData = ExpintoKD(xe, ye, ze, exppts)
#trndata = RandTrains(mini, maxi, trnpts)

