import numpy as np
from  scipy import spatial
from matplotlib import cm
import matplotlib.pyplot as plt
import time
import math
from mpl_toolkits.mplot3d import Axes3D

## Notes:
print 
print '^---- N Dimensional Interpolation ----^'
print 

## Written by Stephen Marone at NASA GRC
## Nomenclature:
# LN-LinearNearest Neighbor,   WN-Weighted Nearest Neighbor, 
# CN-Cosine Nearest Neighbor,  HN-Hermite Nearest Neighbor
# prd - Predicted Points: Points which you want the values at
# trn - Training Points:  Points with known values 
# Ind - Independent: Dimensions which should be known
# Dep - Dependant:   Dimensions which are a function of Ind

# Setting up the original problem ==============================================

class N_Data(object):
	## Main data type with N dimension.  Only needed if you want to use the 
	# example functions, check for error, or plot.
	def __init__(self, funct='None', dtype='None', dims=3):
		## The funct is Crate or PW for Piecewise.
		## The dtype is cart for cartesian or rand for random
		self.funct = funct
		self.dtype = dtype
		self.dims = dims

	def AssignInd(self, points):
		self.points = points

	def AssignDep(self, val):
		self.values = val

	def CreateInd(self, mini, maxi, numpts):
		## Setup independent data for each function via distribution
		if (self.dtype == 'rand'):
			self.points = np.random.rand(numpts,self.dims-1) * \
										(maxi-mini) + mini
			if (self.funct == 'PW'):
				## This ensures points exist on discontinuities
				self.points[0,0] = mini/2
				self.points[1,0] = mini
				self.points[2,0] = 0
				self.points[3,0] = maxi/2
				self.points[4,0] = maxi
				self.points[5:10,0] = 0
				self.points[0:5,1] = 0
				self.points[5,1] = mini/2
				self.points[6,1] = mini
				self.points[7,1] = maxi/2
				self.points[8,1] = maxi

		elif (self.dtype == 'LH'):
			points = np.zeros((numpts,self.dims-1), dtype="float")
			for d in range(self.dims-1):
				## Latin Hypercube it up
				points[:,d] = np.linspace(mini,maxi,numpts)
				np.random.shuffle(points[:,d])
			self.points = points
			
		elif (self.dtype == 'cart'):
			if (self.dims == 3):
				stp = (maxi-mini)/(np.sqrt(numpts)-1)
				x = np.arange(mini,(maxi+stp),stp)
				y = np.arange(mini,(maxi+stp),stp)
				x, y = np.meshgrid(x, y)
				self.points = np.transpose(np.array([x.flatten(), y.flatten()]))
			else:
				print 'Can not currently do cartesian in any dimension \
						except for 3'
				raise SystemExit
			'''
			## Could not figure out why np would not do a meshgrid in N dims
			## Defaults for crate: mini=-512,maxi=512,numpts=81
			stp = (maxi-mini)/(np.sqrt(numpts)-1)
			xm = np.meshgrid(*[np.arange(mini,(maxi+stp),stp) for \
					d in range(dims)], sparse=False, indexing='ij')
			'''
		else:
			print 'Distribution type not found.'
			raise SystemExit

	def CreateDep(self):
		## This not only creates z but the data subtype as well
		
		## Check if points have been created
		try:
			self.points
		except NameError:
			print 'Independent variables must be created first'
			print
			return

		if (self.funct == 'Crate'):
			x = self.points[:,0]
			y = self.points[:,1]
			## Egg Crate Function
			z = -x * np.sin(np.sqrt(abs(x - y - 47))) - \
			    (y+47) * np.sin(np.sqrt(abs(x/2 + y + 47)))
			self.values = z
		elif (self.funct == 'PW'):
			x = self.points[:,0]
			y = self.points[:,1]
			z = np.empty(len(x))
			## Piecewise Function with 3 Planes
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
			self.values = z
		else:
			print 'Function type not found.'
			raise SystemExit

	def PlotResults(self, title):
		## Plots the value results of a 3 dimensional problem

		## Ensure points and values have been created
		try:
			self.points
			self.values
		except NameError:
			print 'Points or values have not been set.'
			return

		x = self.points[:,0]
		y = self.points[:,1]
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		surf = ax.tricontour(x.flatten(), y.flatten(), self.values, 
							 100, rstride=1, cstride=1, 
							 cmap=cm.coolwarm, linewidth=0, 
							 antialiased=False)
		fig.colorbar(surf, shrink=0.5, aspect=5)
		fig.suptitle('%s' % title, 
					 fontsize=14, fontweight='bold')
		ax.set_xlabel('X - Independent Variable')
		ax.set_ylabel('Y - Independent Variable')
		ax.set_zlabel('Z - Dependent Variable')

	def PlotPoints(self, title):
		## Dim-1 plot of point locations for visualization

		## Ensure points have been created
		try:
			self.points
		except NameError:
			print 'Points have not been set.'
			return

		fig = plt.figure()
		if (self.dims == 3):
			plt.scatter(self.points[:,0], self.points[:,1])
			fig.suptitle('%s' % title, 
							fontsize=14, fontweight='bold')
			plt.xlabel('X Value')
			plt.ylabel('Y Value')
		elif (self.dims == 4):
			ax = fig.gca(projection='3d')
			ax.scatter(self.points[:,0], self.points[:,1], 
						self.points[:,2], c=c, marker=m)
			fig.suptitle('%s' % title, 
							fontsize=14, fontweight='bold')
			ax.set_xlabel('X Value')
			ax.set_ylabel('Y Value')
			ax.set_zlabel('Z Value')

	def FindError(self, title, plot=True, check=False):
		## Find and plot error of the found values
		
		## Ensure points and values have been created
		try:
			self.points
			self.values
		except NameError:
			print 'Points or values have not been set.'
			return

		xc = self.points[:,0]
		yc = self.points[:,1]
		loc = 1+np.arange(len(self.values), dtype='float')
		
		if (self.funct == 'Crate'):
			zc = -xc * np.sin(np.sqrt(abs(xc - yc - 47))) - \
					(yc+47) * np.sin(np.sqrt(abs(xc/2 + yc + 47)))
		elif (self.funct == 'PW'):
			# Plot data against values from crate problem
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
		else:
			print 'No function available for comparison.'
			print
			return

		error = abs((zc-self.values)/(zc+0.000000000001)) * 100
		avg = np.average(error)
		'''
		if check:
			for i in range(len(error)):
				if (error[i] > 100):
					print '**High Percent Error found of', error[i], '**'
					print '-Location:', xc[i], yc[i]
					print '-Calculated Value:', zc[i]
					print '-Found Value:', self.values[i] 
					print
		'''
		print '%s Information' % title
		print '-Max Percent Error:', np.max(error)
		print '-Average Percent Error:', avg
		print

		if plot:
			fig = plt.figure()
			plt.plot(loc, error, 'k')
			fig.suptitle('%s Per Point' % title, fontsize=14, fontweight='bold')
			plt.xlabel('Point Identifier')
			plt.ylabel('Percent Error')
			plt.ylim((0,((avg+0.00000000001)*2)))

		self.error = error
		self.actualvalues = zc

class InterpBase(object):
	## The shared aspects of all interpolators
	def __init__(self, TrainPts, TrainVals, NumLeaves=8):
		## TrainPts and TrainVals are the known points and their
		# respective values which will be interpolated against.
		# The IType is the type of interpolation that the user wants.
		self.tp = TrainPts
		self.tv = TrainVals
		self.dims = len(TrainPts[0,:]) + 1
		self.ntpts = len(TrainPts[:,0])
		self.nl = NumLeaves

	def FindNeighs(self, PrdPts, N=5): 
		## Uses a KD Tree to set the distances, locations, and values 
		# of the closest neighbors to each point in the 
		nppts = len(PrdPts[:,0])
		prdz = np.zeros(nppts, dtype="float")
		numtrn = self.ntpts
		numleaves = self.nl

		if (numleaves > (numtrn/N)):
			## Can not have more leaves than training data points. Also 
			# need to ensure there will be enough neighbors in each leaf
			numleaves = numtrn // (N*1.5)
			print 'Number of leaves has been changed to %s' % numleaves
			if (numleaves > numtrn):
				print 'The problem has too few training points for'
				print 'its number of dimensions.'
				raise SystemExit

		## Make into training data into a Tree
		leavesz = math.ceil(numtrn/numleaves)
		KData = spatial.KDTree(self.tp,leafsize=leavesz)

		## KData query takes (data, #ofneighbors) to determine closest 
		# training points to predicted data
		ndist, nloc = KData.query(PrdPts ,N)
		return nppts, ndist, nloc

class LNInterp(InterpBase):
	def __class__(self, PrdPts):
		## This method uses linear interpolation by defining a plane with
		# a set number of nearest neighbors to the predicted
		dims = self.dims
		
		## Find them neigbors
		nppts, ndist, nloc = self.FindNeighs(PrdPts, N=self.dims)

		print 'Weighted Nearest Neighbors (WNN) KDTree Results'
		print '-Number of Neighbors per Predicted Point:', N
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		## Extra Inputs for Finding the normal are found below
		prdz = np.zeros((nppts), dtype="float")

		## Number of row vectors needed always dimensions - 1
		r = np.arange(dims-1, dtype="int")
		
		## Diagonal counters are found here to speed up finding each normal
		da = np.zeros((dims, dims-1), dtype="int")
		db = np.zeros((dims, dims-1), dtype="int")
		for i in xrange(dims):
		    da[i] = (r+1+i) % dims
		    db[i] = (i-r-1) % dims
		
		for t in range(self.nppts):  ## Go through each predict point 
			if (np.all(abs(PrdPts - self.tp[nloc[t,0],:])) > 0.000001):
				nvect = np.empty((0,dims), float)
				## Planar vectors need both dep and ind dimensions
				trnd = np.concatenate((self.tp[nloc[t,:],:] ,
										self.tv[nloc[t,:],:]),
				 							axis=1)
				for n in range(dims-1): ## Go through each neighbor
					## Creates array[neighbor, dimension] from NN results
					nvect = np.append(nvect, [trnd[(n+1),:] - trnd[n,:]], 0)
	
				## It is known that found array will be size [dims-1, dims]
				# and can be used below to in an exterior product.
				normal = np.prod(nvect[r,da], axis=1) - \
						 np.prod(nvect[r,db], axis=1)

				## The pc is the constant of the n dimensional plane. 
				## It uses the point which is from the closest neighbor
				pc = np.dot(traindata[nloc[t,0],:],normal)

				## Set the prddata to a temporary array so it has the same
				# dimensions as the training data with z as 0 for now
				prdtemp = np.concatenate((prddata[t,:], np.array([0])), 1)
				
				if (normal[-1] == 0):
					print 'Error, dependent variable has infinite answers.'
					raise SystemExit
				prdz[t] = (np.dot(prdtemp,normal)-pc)/-normal[-1]
	
			else:
				## Closest neighbor overlaps predicted point
				#print 'That just happened'
				#print (prddata[t,:-1] - traindata[nloc[t,0],:-1])
				prdz[t] = self.tv[nloc[t,0]]
		
			#prddict['%s' % prddata[t,:]] = prdz[t] #from when I wanted a dict
		return prdz
	
class WNInterp(InterpBase):
	## Weighted Neighbor Interpolation
	def __class__(self, PrdPts, N=8, DistEff=5):

		nppts, ndist, nloc = self.FindNeighs(PrdPts, N=N)
		
		prdz = np.zeros((nppts), dtype="float")

		print 'Weighted Nearest Neighbors (WNN) KDTree Results'
		print '-Number of Neighbors per Predicted Point:', N
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		for t in range(nppts):  ## Go through each predict point 
			if (np.all(abs(PrdPts[t,:] - self.tp[nloc[t,0],:])) > 0.000001):
				sumdist = np.sum(1/(ndist[t,:]**DistEff))
				wt = np.sum(self.tv[nloc[t,:]]/(ndist[t,:]**DistEff))
				prdz[t] = wt/sumdist
	
			else:
				## Closest neighbor overlaps predicted point
				#print 'That just happened'
				#print (prddata[t,:-1] - traindata[nloc[t,0],:-1])
				prdz[t] = self.tv[nloc[t,0]]
		return prdz

class CNInterp(InterpBase):
	## Cosine Neighbor Interpolation
	def __class__(self, PrdPts, N=5, tension=0, bias=0):

		nppts, ndist, nloc = self.FindNeighs(PrdPts, N=N)

		prdz = np.zeros((nppts), dtype="float")
		
		print 'Cosine Nearest Neighbors (CN) Interpolation KDTree Results'
		print '-Number of Neighbors per Predicted Point:', N
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		for t in range(nppts):  ## Go through each predict point 
			if (np.all(abs(PrdPts[t,:] - self.tp[nloc[t,0],:])) > 0.000001):
				## This sorts the rows by the values in the column specified. 
				## The neighvals dim is [neigh#, dimension] 

				neighvals = np.concatenate((self.tp[nloc[t,:],:] ,
											self.tv[nloc[t,:],:]),
				 							axis=1)
				orgneighs = neighvals[neighvals[:,0].argsort(),:]
				
				i = 0
				while True:
					i += 1
					if (i == N-1):
						i -= 1
						break
					elif (PrdPts[t,0] < orgneighs[i,0]):
						break
	
				diff = ((PrdPts[t,0] - orgneighs[i-1,0]) / 
				        (orgneighs[i,0]-orgneighs[i-1,0]))
				mu2 = (1-np.cos(diff*np.pi))/2
				prdz[t] = orgneighs[i-1,-1]*(1-mu2) + orgneighs[i,-1]*mu2
				'''
				diff = ((PrdPts[t,0] - orgneighs[0,0]) / 
				        (orgneighs[1,0]-orgneighs[0,0]))
				prdz[t] = y[0]*(1-diff) + y[1]*diff
				'''
				for D in range(self.dims-2):
					## Dim started above to initiate prdz, then continues to finish
					dim = D+1
					orgneighs = neighvals[neighvals[:,dim].argsort(),:]
					
					i = 0
					while True:
						i += 1
						if (i == N-1):
							i -= 1
							break
						elif (PrdPts[t,dim] < orgneighs[i,dim]):
							break
	
					## Place of prddata along mid neighbors
					diff = ((PrdPts[t,dim] - orgneighs[i-1,dim]) / 
					        (orgneighs[i,dim]-orgneighs[i-1,dim]))
					mu2 = (1-np.cos(diff*np.pi))/2
					prdz[t] += orgneighs[i-1,-1]*(1-mu2) + orgneighs[i,-1]*mu2
					'''
					diff = ((PrdPts[t,dim] - orgneighs[0,dim]) / 
				    	    (orgneighs[1,dim]-orgneighs[0,dim]))
					prdz[t] += y[0]*(1-diff) + y[1]*diff
					'''
					prdz[t] = prdz[t]/2 ## Average with previous
			else:
				# Closest neighbor overlaps predicted point
				#print 'That just happened'
				#print (PrdPts[t,:-1] - self.tp[nloc[t,0],:-1])
				prdz[t] = self.tv[nloc[t,0]]
		return prdz

class HNInterp(InterpBase):
	## Hermite Neighbor Interpolation
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
	
	def __class__(self, PrdPts, N=8, tension=0, bias=0):
		## Traindata has n dimensions, prddata has n-1 dimensions because last
		# dimension is the dependant one we are solving for.  
		## N neighbors will be found but only 4 will be used per dimension

		nppts, ndist, nloc = self.FindNeighs(PrdPts, N=N)

		prdz = np.zeros((nppts), dtype="float")
		
		print 'Hermite Nearest Neighbors (HN) Interpolation KDTree Results'
		print '-Number of Neighbors per Predicted Point:', N
		print '-Farthest Neighbor Distance:', np.max(ndist)
		print
		for t in range(numprd):  ## Go through each predict point 
			if (np.all(abs(prddata[t,:-1] - traindata[nloc[t,0],:-1])) > .001):
				## This sorts the rows by the values in the column specified. 
				## The neighvals dim is [neigh#, dimension] 

				neighvals = np.concatenate((self.tp[nloc[t,:],:] ,
											self.tv[nloc[t,:],:]),
				 							axis=1)
				orgneighs = neighvals[neighvals[:,0].argsort(),:]
				
				i = 1
				while True:
					i += 1
					if (i == N-2):
						i -= 1
						break
					elif (PrdPts[t,0] < orgneighs[i,0]):
						break
	
				diff = ((PrdPts[t,0] - orgneighs[i-1,0]) / 
				        (orgneighs[i+2,0]-orgneighs[i-1,0]))
				prdz[t] = self.HermFunct(orgneighs[(i-2):(i+2),-1],
										 diff, tension, bias)
				'''
				diff = ((PrdPts[t,0] - orgneighs[0,0]) / 
				        (orgneighs[3,0]-orgneighs[0,0]))
				prdz[t] = self.HermFunct(orgneighs[:4,-1],
										 diff, tension, bias)
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
						elif (PrdPts[t,dim] < orgneighs[i,dim]):
							break
	
					## Place of prddata along mid neighbors
					diff = ((PrdPts[t,dim] - orgneighs[i-1,dim]) / 
					        (orgneighs[i+2,dim]-orgneighs[i-1,dim]))
					prdz[t] += self.HermFunct(orgneighs[(i-2):(i+2),-1],
											  diff, tension, bias)
					'''
					diff = ((PrdPts[t,dim] - orgneighs[0,dim]) / 
				    	    (orgneighs[3,dim]-orgneighs[0,dim]))
					prdz[t] += self.HermFunct(orgneighs[:4,-1],
											  diff, tension, bias)
					'''
					prdz[t] = prdz[t]/2 ## Average with previous
			else:
				# Closest neighbor overlaps predicted point
				print 'That just happened'
				print (PrdPts[t,:-1] - self.tp[nloc[t,0],:-1])
				prdz[t] = self.tv[nloc[t,0]]
		return prdz

check = N_Data('Crate', 'cart')
check.CreateInd(-512,512,16641)
check.CreateDep()
check.PlotResults('Checked Data')
#print check.funct
#print check.dtype
#print len(check.points)
#print check.points
check.FindError('Checked Data')

train = N_Data('PW', 'rand')
train.CreateInd(-500,500,50000)
train.CreateDep()
train.PlotResults('Training Data')
#print len(train.points)
#print train.points
#print train.dims
train.FindError('Training Data', False, True)

pred = N_Data('Crate', 'LH', 3)
pred.CreateInd(-512,512,81)
pred.AssignDep(check.values[:81])
pred.PlotResults('Predicted Data')
pred.PlotPoints('Point Locations')
#print pred.points.shape
#print len(pred.points)
#print pred.points
#print pred.dims
pred.FindError('Bad Stuff', True, True)

bs = N_Data('PW', 'LH', 7)
#print bs.dims

plt.show()

## Side note, PrdPts are predicted points technically although it's not really 
# that simple. They are better named pride points.  Everyone needs pride points,
# but their value is unknown usually and can only be determined when related to
# those around you.
