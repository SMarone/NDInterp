from numpy import linalg as LA
import numpy as np
from  scipy import spatial
from matplotlib import cm
import matplotlib.pyplot as plt
import time
import math
from mpl_toolkits.mplot3d import Axes3D
import copy as cp

## Notes:
print 
print '^---- N Dimensional Linear Interpolation ----^'
print 

# Setting up the original problem ==============================================

def CreateCrateFunct():
	maxi = 512
	mini = -512
	stp =6
	print 'Setting up Crate Function'
	print '-Maximum x and y:', maxi
	print '-Minimum x and y:', mini
	print
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	# Setup Function Data
	x = np.arange(mini,maxi,stp)
	y = np.arange(mini,maxi,stp)
	x, y = np.meshgrid(x, y)
	# Egg Crate Function
	z = -x * np.sin(np.sqrt(abs(x - y - 47))) - (y+47) * np.sin(np.sqrt(abs(x/2 + y + 47)))
	# Crate Visual Check
	surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	fig.suptitle('Actual Function Representation', fontsize=14, fontweight='bold')
	ax.set_xlabel('X - Independent Variable')
	ax.set_ylabel('Y - Independent Variable')
	ax.set_zlabel('Z - Dependent Variable')
	#plt.show()
	return x, y, z, maxi, mini


def CreatePWFunct():
	maxi = 500
	mini = -500
	stp = 5
	print 'Setting up 3 Plane Piecewise Function'
	print '-Maximum x and y:', maxi
	print '-Minimum x and y:', mini
	print
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	# Setup Function Data
	x = np.arange(mini,maxi,stp)
	y = np.arange(mini,maxi,stp)
	x, y = np.meshgrid(x, y)
	x = x.flatten()
	y = y.flatten()
	z = np.empty(len(x))
	# Make Piecewise Function with 3 Planes
	count = 0
	for i in range(len(x)):
		if (x[i] > 0):
			if (y[i] > 0):
				z[count] = x[i] + 8*y[i]
				count +=1
			else:
				z[count] = 5*x[i] + y[i]
				count +=1
		else:
			z[count] = 5*x[i] + 6*y[i]
			count +=1
	fdata = np.transpose(np.array([x, y, z]))
	# PW Visual Check
	surf = ax.tricontour(x.flatten() , y.flatten() , z.flatten(), 100, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	fig.colorbar(surf, shrink=0.5, aspect=5)
	fig.suptitle('Actual Function Representation', fontsize=14, fontweight='bold')
	ax.set_xlabel('X - Independent Variable')
	ax.set_ylabel('Y - Independent Variable')
	ax.set_zlabel('Z - Dependent Variable')
	#plt.show()
	return fdata, maxi, mini

# Generate and Check Training Data =============================================

def CreateTrainCrate(mini, maxi, numtrn):
	## Create Random Training Values
	xe = np.random.rand(numtrn) * (maxi-mini) + mini
	ye = np.random.rand(numtrn) * (maxi-mini) + mini
	## Alter training data for more useable format
	ze = -xe * np.sin(np.sqrt(abs(xe - ye - 47))) - (ye+47) * np.sin(np.sqrt(abs(xe/2 + ye + 47)))
	traindata = np.transpose(np.array([xe, ye, ze]))
	return traindata


def CreateTrainPW(mini, maxi, numtrn):
	## Create Random Training Values
	xe = np.random.rand(numtrn) * (maxi-mini) + mini
	ye = np.random.rand(numtrn) * (maxi-mini) + mini
	# Ensure There are points on edge of PW problem
	xe[0] = mini/2
	xe[1] = mini
	xe[2] = 0
	xe[3] = maxi/2
	xe[4] = maxi
	xe[5:10] = 0
	ye[0:5] = 0
	ye[5] = mini/2
	ye[6] = mini
	ye[7] = maxi/2
	ye[8] = maxi
	## Alter training data for more useable format
	ze = np.empty(len(xe))
	# Make Piecewise Function with 3 Planes
	count = 0
	for i in range(len(xe)):
		if (xe[i] > 0):
			if (ye[i] > 0):
				ze[count] = xe[i] + 8*ye[i]
				count +=1
			else:
				ze[count] = 5*xe[i] + ye[i]
				count +=1
		else:
			ze[count] = 5*xe[i] + 6*ye[i]
			count +=1
	traindata = np.transpose(np.array([xe, ye, ze]))
	return traindata

def DataPlot(xe, ye, ze, dname):
	## Create Surface plot to Check Unstructured Training Data
	fige = plt.figure()
	axe = fige.gca(projection='3d')
	## Nice Surface plot does not work because of unstructured data
	# surf = ax.plot_surface(xem, yem, ze, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	surfe = axe.tricontour(xe.flatten() , ye.flatten() , ze.flatten(), 100, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	fige.colorbar(surfe, shrink=0.5, aspect=5)

	fige.suptitle('%s Representation' % dname, fontsize=14, fontweight='bold')
	axe.set_xlabel('X - Independent Variable')
	axe.set_ylabel('Y - Independent Variable')
	axe.set_zlabel('Z - Dependent Variable')
	#plt.show()

# Generate and Check Check Data =============================================

def RandPreds(mini, maxi, numprd):
	## Creates random locations for the predicted points
	xt = np.random.rand(numprd) * float((maxi-mini)) + float(mini)
	yt = np.random.rand(numprd) * float((maxi-mini)) + float(mini)
	prddata = np.transpose(np.array([xt,yt]))
	return prddata


def LHPreds(traindata, numprd, res):
	## This takes in training data and sets up predict data that is dispersed
	# in the form of a latin hypercube
	## Dim of Check data is -1 of Train data since it does not include
	# the dependent value you are solving for
	dim = len(traindata[0,:])
	## The res variable is an integer that proclaims resolution of 
	# problem step, default 100 for crate function
	mini = traindata.min(axis=0)
	mini = mini[:-1]
	maxi = traindata.max(axis=0)
	maxi = maxi[:-1]
	## Adjust data to be integer forms of minimum and maximum
	mi = np.ceil(mini / float(res)) * int(res)
	ma = np.floor(maxi / float(res)) * int(res)
	## Floor and ceiling to avoid going out of bounds
	p = np.zeros((numprd,dim-1), dtype="float")
	for i in range(dim-1):
		## Latin Hypercube it up
		p[:,i] = np.linspace(mi[i],ma[i],numprd)
		np.random.shuffle(p[:,i])
	return p

def PlotPred2D(p, mini, maxi):
	xp = p[:,0]
	yp = p[:,1]
	fig = plt.figure()
	plt.scatter(xp, yp)
	fig.suptitle('2D Predicted Point Locations', fontsize=14, fontweight='bold')
	plt.xlabel('X Value')
	plt.xlim((mini,maxi))
	plt.ylabel('Y Value')
	plt.ylim((mini,maxi))

# Interpolation Main Body Functions ============================================

def PrepareTrain(traindata):
	## Extra Inputs for Find Normal
	## The dimension of the problem is the same as the length of any row of data
	dims = len(traindata[0,:])
	## Number of row vectors needed for interpolation always dimensions - 1
	r = np.arange(dims-1, dtype="int")
	## Diagonal counters are found here in order to speed up finding each normal
	da = np.zeros((dims, dims-1), dtype="int")
	db = np.zeros((dims, dims-1), dtype="int")
	for i in xrange(dims):
	    da[i] = (r+1+i) % dims
	    db[i] = (i-r-1) % dims
	return dims, r, da, db


def FindNormal(vects, point, r, da, db):
	## It is known that the given array will be of size [dims-1, dims]
	## First two vectors are values which must be multiplied together for
	# an exterior product
	proda = np.prod(vects[r,da], axis=1)
	prodb = np.prod(vects[r,db], axis=1)
	norm = (proda - prodb)
	## The pc is the constant of the n dimensional plane function.  Uses 
	# point which is from the closes neighbor
	pc = np.dot(point,norm)
	return norm, pc


def N_Neighbors(traindata, prddata, numleaves):
	## Traindata has n dimensions, prddata has n-1 dimensions because last
	# dimension is the dependant one we are solving for.  
	## Get the values for finding each normal
	dims, r, da, db = PrepareTrain(traindata)
	maxi = traindata.max(axis=0)
	mini = traindata.min(axis=0)
	numprd = len(prddata[:,0])
	prddict = {}
	prdz = np.zeros((numprd), dtype="float")
	numtrn = len(traindata[:,0])
	if (numleaves > numtrn):
		## Can not have more leaves than training data points. Also 
		# need to ensure there will be enough neighbors in each leaf
		numleaves = numtrn // (dims*2)
		print 'Number of leaves has been changed to %s' % numleaves
		if (numleaves > numtrn):
			print 'The problem has too few training points for'
			print 'its number of dimensions.'
			raise SystemExit

	## Make into Tree
	leavesz = math.ceil(numtrn/numleaves)
	## Can only use first n-1 dimensions of traindata since that is all
	# which is inside the prddata
	KData = spatial.KDTree(traindata[:,:-1],leafsize=leavesz)

	## KData query takes (data, #ofneighbors) to determine closest 
	# training points to predicted data
	ndist, nloc = KData.query(prddata ,dims)
	print 'Nearest Neighbor (NN) KDTree Results'
	print '-Number of Neighbors per Predicted Point:', dims
	print '-Farthest Neighbor Distance:', np.max(ndist)
	print
	for t in range(numprd):  ## Go through each predict point 
		if (np.all(abs(prddata[t,:-1] - traindata[nloc[t,0],:-1])) > .001):
			nvect = np.empty((0,dims), float)
			for n in range(dims-1): ## Go through each neighbor
				## Creates array [neighbor, dimension] from KData results
				'''
				if (np.all(abs(traindata[nloc[t,n+1],:] - 
				               traindata[nloc[t,  n],:])) < 0.00000001):
					print 'Imagine this. . .'
					print [traindata[nloc[t,n+1],:] - traindata[nloc[t,  n],:]]
					print 'WTH'
				'''
				nvect = np.append(nvect, [traindata[nloc[t,n+1],:] - 
					                      traindata[nloc[t,  n],:]], 0)

			point = traindata[nloc[t,0],:]
			normal, pc = FindNormal(nvect, point, r, da, db)
			## Set the prddata to a temporary array so it has the same
			# dimensions as the training data with z as 0 for now
			prdtemp = np.concatenate((prddata[t,:], np.array([0])), 1)
			if (normal[-1] == 0):
				print 'Error, dependent variable has infinite answers.'
				raise SystemExit
			prdz[t] = (np.dot(prdtemp,normal)-pc)/-normal[-1]

		else:
			# Closest neighbor overlaps predicted point
			print 'That just happened'
			print (prddata[t,:-1] - traindata[nloc[t,0],:-1])
			prdz[t] = traindata[nloc[t,0],-1]
		prddict['%s' % prddata[t,:]] = prdz[t]
	return prdz, prddict


def N_WeightNeighbors(traindata, prddata, numleaves, N, L):
	## Traindata has n dimensions, prddata has n-1 dimensions because last
	# dimension is the dependant one we are solving for.  
	## Get the values for finding each normal
	numprd = len(prddata[:,0])
	prddict = {}
	prdz = np.zeros((numprd), dtype="float")
	numtrn = len(traindata[:,0])
	
	if (numleaves > numtrn):
		## Can not have more leaves than training data points. Also 
		# need to ensure there will be enough neighbors in each leaf
		numleaves = numtrn // (len(traindata[0,:])*2)
		print 'Number of leaves has been changed to %s' % numleaves
		if (numleaves > numtrn):
			print 'The problem has too few training points for'
			print 'its number of dimensions.'
			raise SystemExit

	## Make into Tree
	leavesz = math.ceil(numtrn/numleaves)
	KData = spatial.KDTree(traindata[:,:-1],leafsize=leavesz)

	## KData query takes (data, #ofneighbors) to determine closest 
	# training points to predicted data
	ndist, nloc = KData.query(prddata ,N)
	print 'Weighted Nearest Neighbor (WNN) KDTree Results'
	print '-Number of Neighbors per Predicted Point:', N
	print '-Farthest Neighbor Distance:', np.max(ndist)
	print
	for t in range(numprd):  ## Go through each predict point 
		if (np.all(abs(prddata[t,:-1] - traindata[nloc[t,0],:-1])) > .001):
			sumdist = np.sum(1/(ndist[t,:]**L))
			wt = np.sum(traindata[nloc[t,:],-1]/(ndist[t,:]**L))
			prdz[t] = wt/sumdist

		else:
			# Closest neighbor overlaps predicted point
			print 'That just happened'
			print (prddata[t,:-1] - traindata[nloc[t,0],:-1])
			prdz[t] = traindata[nloc[t,0],-1]
		prddict['%s' % prddata[t,:]] = prdz[t]
	return prdz, prddict


def HermFunct(y, mu, tension, bias):
	## Should  have 4 neighbor values (y) and percent distance along line 
	# from y1 and y2 to get the wanted value (mu).  Bias and tension 
	# are optional values.

	mu2 = mu * mu
	mu3 = mu2 * mu

	m0  = (y[1]-y[0])*(1+bias)*(1-tension)/2
	m0 += (y[2]-y[1])*(1-bias)*(1-tension)/2
	m1  = (y[2]-y[1])*(1+bias)*(1-tension)/2
	m1 += (y[3]-y[2])*(1-bias)*(1-tension)/2
	a0 =  2*mu3 - 3*mu2 + 1
	a1 =    mu3 - 2*mu2 + mu
	a2 =    mu3 -   mu2
	a3 = -2*mu3 + 3*mu2

	return (a0*y[1] + a1*m0 + a2*m1 + a3*y[2])


def N_HermNeighbors(traindata, prddata, numleaves, N, tension, bias):
	## Traindata has n dimensions, prddata has n-1 dimensions because last
	# dimension is the dependant one we are solving for.  
	## N neighbors will be found but only 4 will be used per dimension
	dims = len(traindata[0,:])
	numprd = len(prddata[:,0])
	prddict = {}
	prdz = np.zeros((numprd), dtype="float")
	numtrn = len(traindata[:,0])
	if (numleaves > numtrn):
		## Can not have more leaves than training data points. Also 
		# need to ensure there will be enough neighbors in each leaf
		numleaves = numtrn // (dims*2)
		print 'Number of leaves has been changed to %s' % numleaves
		if (numleaves > numtrn):
			print 'The problem has too few training points for'
			print 'its number of dimensions.'
			raise SystemExit

	## Make into Tree
	leavesz = math.ceil(numtrn/numleaves)
	KData = spatial.KDTree(traindata[:,:-1],leafsize=leavesz)

	## KData query takes (data, #ofneighbors) to determine closest 
	# training points to predicted data
	ndist, nloc = KData.query(prddata ,N)
	print 'Hermite Neighbors (HN) Interpolation KDTree Results'
	print '-Number of Neighbors per Predicted Point:', N
	print '-Farthest Neighbor Distance:', np.max(ndist)
	print
	for t in range(numprd):  ## Go through each predict point 
		if (np.all(abs(prddata[t,:-1] - traindata[nloc[t,0],:-1])) > .001):
			## This sorts the rows by the values in the column specified. 
			## The neighvals dim is [neigh#, dimension] 
			neighvals = traindata[nloc[t,:],:]
			orgneighs = neighvals[neighvals[:,0].argsort(),:]
			
			i = 1
			while True:
				i += 1
				if (i == N-2):
					i -= 1
					break
				elif (prddata[t,0] < orgneighs[i,0]):
					break

			diff = ((prddata[t,0] - orgneighs[i-1,0]) / 
			        (orgneighs[i+2,0]-orgneighs[i-1,0]))
			prdz[t] = HermFunct(orgneighs[(i-2):(i+2),-1], diff, tension, bias)
			'''
			diff = ((prddata[t,0] - orgneighs[0,0]) / 
			        (orgneighs[3,0]-orgneighs[0,0]))
			prdz[t] = HermFunct(orgneighs[:4,-1], diff, tension, bias)
			'''
			for D in range(dims-2):
				## Dim started above to initiate prdz, then continues to finish
				dim = D+1
				orgneighs = neighvals[neighvals[:,dim].argsort(),:]
				
				i = 1
				while True:
					i += 1
					if (i == N-2):
						i -= 1
						break
					elif (prddata[t,dim] < orgneighs[i,dim]):
						break

				## Place of prddata along mid neighbors
				diff = ((prddata[t,dim] - orgneighs[i-1,dim]) / 
				        (orgneighs[i+2,dim]-orgneighs[i-1,dim]))
				prdz[t] += HermFunct(orgneighs[(i-2):(i+2),-1], diff, tension, bias)
				'''
				diff = ((prddata[t,dim] - orgneighs[0,dim]) / 
			    	    (orgneighs[3,dim]-orgneighs[0,dim]))
				prdz[t] += HermFunct(orgneighs[:4,-1], diff, tension, bias)
				'''
				prdz[t] = prdz[t]/2 ## Average with previous
		else:
			# Closest neighbor overlaps predicted point
			print 'That just happened'
			print (prddata[t,:-1] - traindata[nloc[t,0],:-1])
			prdz[t] = traindata[nloc[t,0],-1]
		prddict['%s' % prddata[t,:]] = prdz[t]
	return prdz, prddict


def N_CosNeighbors(traindata, prddata, numleaves, N, tension, bias):
	## Traindata has n dimensions, prddata has n-1 dimensions because last
	# dimension is the dependant one we are solving for.  
	## N neighbors will be found but only 4 will be used per dimension
	dims = len(traindata[0,:])
	numprd = len(prddata[:,0])
	prddict = {}
	prdz = np.zeros((numprd), dtype="float")
	numtrn = len(traindata[:,0])
	if (numleaves > numtrn):
		## Can not have more leaves than training data points. Also 
		# need to ensure there will be enough neighbors in each leaf
		numleaves = numtrn // (dims*2)
		print 'Number of leaves has been changed to %s' % numleaves
		if (numleaves > numtrn):
			print 'The problem has too few training points for'
			print 'its number of dimensions.'
			raise SystemExit

	## Make into Tree
	leavesz = math.ceil(numtrn/numleaves)
	KData = spatial.KDTree(traindata[:,:-1],leafsize=leavesz)

	## KData query takes (data, #ofneighbors) to determine closest 
	# training points to predicted data
	ndist, nloc = KData.query(prddata ,N)
	print 'Cosine Neighbors (CN) Interpolation KDTree Results'
	print '-Number of Neighbors per Predicted Point:', N
	print '-Farthest Neighbor Distance:', np.max(ndist)
	print
	for t in range(numprd):  ## Go through each predict point 
		if (np.all(abs(prddata[t,:-1] - traindata[nloc[t,0],:-1])) > .001):
			## This sorts the rows by the values in the column specified. 
			## The neighvals dim is [neigh#, dimension] 
			neighvals = traindata[nloc[t,:],:]
			orgneighs = neighvals[neighvals[:,0].argsort(),:]
			
			i = 0
			while True:
				i += 1
				if (i == N-1):
					i -= 1
					break
				elif (prddata[t,0] < orgneighs[i,0]):
					break

			diff = ((prddata[t,0] - orgneighs[i-1,0]) / 
			        (orgneighs[i,0]-orgneighs[i-1,0]))
			mu2 = (1-np.cos(diff*np.pi))/2
			prdz[t] = orgneighs[i-1,-1]*(1-mu2) + orgneighs[i,-1]*mu2
			'''
			diff = ((prddata[t,0] - orgneighs[0,0]) / 
			        (orgneighs[1,0]-orgneighs[0,0]))
			prdz[t] = y[0]*(1-diff) + y[1]*diff
			'''
			for D in range(dims-2):
				## Dim started above to initiate prdz, then continues to finish
				dim = D+1
				orgneighs = neighvals[neighvals[:,dim].argsort(),:]
				
				i = 0
				while True:
					i += 1
					if (i == N-1):
						i -= 1
						break
					elif (prddata[t,dim] < orgneighs[i,dim]):
						break

				## Place of prddata along mid neighbors
				diff = ((prddata[t,dim] - orgneighs[i-1,dim]) / 
				        (orgneighs[i,dim]-orgneighs[i-1,dim]))
				mu2 = (1-np.cos(diff*np.pi))/2
				prdz[t] += orgneighs[i-1,-1]*(1-mu2) + orgneighs[i,-1]*mu2
				'''
				diff = ((prddata[t,dim] - orgneighs[0,dim]) / 
			    	    (orgneighs[1,dim]-orgneighs[0,dim]))
				prdz[t] += y[0]*(1-diff) + y[1]*diff
				'''
				prdz[t] = prdz[t]/2 ## Average with previous
		else:
			# Closest neighbor overlaps predicted point
			print 'That just happened'
			print (prddata[t,:-1] - traindata[nloc[t,0],:-1])
			prdz[t] = traindata[nloc[t,0],-1]
		prddict['%s' % prddata[t,:]] = prdz[t]
	return prdz, prddict

# Error Check ==================================================================

def PlotErrorCrate(prddata, prdz, title):
	# Plot data against values from crate problem
	xc = prddata[:,0] 
	yc = prddata[:,1]
	zc = -xc * np.sin(np.sqrt(abs(xc - yc - 47))) - (yc+47) * np.sin(np.sqrt(abs(xc/2 + yc + 47)))
	point = 1+np.arange(len(zc), dtype='float')
	error = abs(((zc-prdz))/zc) * 100
	avg = np.average(error)
	'''
	for i in range(len(error)):
		if (error[i] > 99):
			print '**High Percent Error found of', error[i], '**'
			print '-Location:', xc[i], yc[i]
			print '-Calculated Value:', zc[i]
			print '-Found Value:', prdz[i] 
			print
	'''
	print '%s Information' % title 
	print '-Max Percent Error:', np.max(error)
	print '-Average Percent Error:', avg
	print 
	fig = plt.figure()
	plt.plot(point, error, 'k')
	fig.suptitle('%s Per Point' % title, fontsize=14, fontweight='bold')
	plt.xlabel('Point Identifier')
	plt.ylabel('Percent Error')
	plt.ylim((0,(avg*2)))


def PlotErrorPW(prddata, prdz, title):
	# Plot data against values from crate problem
	xc = prddata[:,0] 
	yc = prddata[:,1]
	zc = np.empty(len(xc))
	# Make Piecewise Function with 3 Planes
	count = 0
	for i in range(len(xc)):
		if (xc[i] > 0):
			if (yc[i] > 0):
				zc[count] = xc[i] + 8*yc[i]
				count +=1
			else:
				zc[count] = 5*xc[i] + yc[i]
				count +=1
		else:
			zc[count] = 5*xc[i] + 6*yc[i]
			count +=1
	point = 1+np.arange(len(zc), dtype='float')
	error = abs((zc-prdz)/zc) * 100
	avg = np.average(error)
	'''
	for i in range(len(error)):
		if (error[i] > 99):
			print '**High Percent Error found of', error[i], '**'
			print '-Location:', xc[i], yc[i]
			print '-Calculated Value:', zc[i]
			print '-Found Value:', prdz[i] 
			print
	'''
	print '%s Information' % title 
	print '-Max Percent Error:', np.max(error)
	print '-Average Percent Error:', avg
	print 

	fig = plt.figure()
	plt.plot(point, error, 'k')
	fig.suptitle('%s Per Point' % title, fontsize=14, fontweight='bold')
	plt.xlabel('Point Identifier')
	plt.ylabel('Percent Error')
	plt.ylim((0,avg))

# Run Code =====================================================================

## Problem function type, choices: Crate=Egg Crate; PW=Piecewise
PT = 'PW'
## Predicted data selection type, choices: LH=Latin Hypercube; Rand=random
PDT = 'LH'   

if (PT == 'Crate'):
	x, y, z, maxi, mini = CreateCrateFunct()
elif (PT == 'PW'):
	fdata, maxi, mini = CreatePWFunct()
else:
	print 'Incorrect Problem Type'
	raise SystemExit

## Misc. Code Inputs
WN = 5 ## # of neighbors to solve for, affects WN.  Default 5
HN = 8 ## # of neighbors to solve for, affects HN.  Min of 5, default 8
CN = 5 ## # of neighbors to solve for, affects CN.  Min of 3, default 5
L = 5 ## Value which alters weighted neighbors function dependance to distance
      # can go from 1 to inf. Default 5
res = 10 ## Resolution of rounding for the Latin Hypercube DOE. Default 10
numleaves = 100  ## # of leaves for KD Tree, must be < numtrn. Default 100
numprd = 1000 ## # of predicted points, increase to raise plot resolution. Default 1000
numtrn = 500000 ## # of training points to be create, increase to lower error. Default 500000
tension = 0 ## For Herm Neighs, 1 is high, 0 default, -1 is low
bias = 0 ## For Herm Neighs, 0 default, + or - closer to first seg or sec resp.

print 'Interpolation Input Values'
print '-KDTree Leaves:', numleaves
print '-Number of Training Points:', numtrn
print '-Number of Predicted Points:', numprd
print
print '^---- Running Code ----^'
print


if (PT == 'Crate'):
	traindata = CreateTrainCrate(mini, maxi, numtrn)
elif (PT == 'PW'):
	traindata = CreateTrainPW(mini, maxi, numtrn)
DataPlot(traindata[:,0], traindata[:,1], traindata[:,2], 'Training Data')

if (PDT == 'Rand'):
	prddata = RandPreds(mini, maxi, numprd)
elif (PDT == 'LH'):
	prddata = LHPreds(traindata, numprd, res)
else:
	print'Incorrect Predicted Data Type'
	raise SystemExit
PlotPred2D(prddata, mini, maxi)

wpd = cp.deepcopy(prddata)
cpd = cp.deepcopy(prddata)
hpd = cp.deepcopy(prddata)

t1 = time.time()
zch, zchd = N_Neighbors(traindata, prddata, numleaves)
t2 = time.time()
zchW, zchWd = N_WeightNeighbors(traindata, wpd, numleaves, WN, L)
t3 = time.time()
zchH, zchHd = N_HermNeighbors(traindata, hpd, numleaves, HN, tension, bias)
t4 = time.time()
zchC, zchCd = N_CosNeighbors(traindata, cpd, numleaves, CN, tension, bias)
t5 = time.time()

print 'Results Found and Plotting Now'
print '-NN Interpolation Run Time:', (t2-t1)
print '-WNN Interpolation Run Time:', (t3-t2)
print '-CN Interpolation Run Time:', (t5-t4)
print '-HN Interpolation Run Time:', (t4-t3)
print
DataPlot(prddata[:,0], prddata[:,1], zch, 'NN Predicted Data')
DataPlot(wpd[:,0], wpd[:,1], zchW, 'WNN Predicted Data')
DataPlot(hpd[:,0], hpd[:,1], zchH, 'HN Predicted Data')
DataPlot(cpd[:,0], cpd[:,1], zchC, 'CN Predicted Data')

if (PT == 'Crate'):
	PlotErrorCrate(prddata, zch, 'NN Error') 
	PlotErrorCrate(prddata, zchW, 'WNN Error') 
	PlotErrorCrate(prddata, zchH, 'HN Error') 
	PlotErrorCrate(prddata, zchC, 'CN Error') 
elif (PT == 'PW'):
	PlotErrorPW(prddata, zch, 'NN Error') 
	PlotErrorPW(prddata, zchW, 'WNN Error') 
	PlotErrorPW(prddata, zchH, 'HN Error') 
	PlotErrorPW(prddata, zchC, 'CN Error') 

'''
with open("LV2_Error.txt", "a") as efile:
    efile.write("\nRun Times\n")
    efile.write("\n-LN Interpolator:")
    efile.write(str(t2-t1))
    efile.write("\n-WN Interpolator:")
    efile.write(str(t3-t2))
    efile.write("\n-CN Interpolator:")
    efile.write(str(t5-t4))
    efile.write("\n-HN Interpolator:")
    efile.write(str(t4-t3))
'''
#plt.show()
